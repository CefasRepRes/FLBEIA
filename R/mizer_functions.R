#############
##
############

project_dyn <- function (object, effort, t_max = 100, dt = 0.1, t_save = 1, 
          t_start = 0, initial_n, initial_n_pp, append = TRUE, progress_bar = TRUE, inter_mat,kappa_dyn,
          ...) 
{
  validObject(object)
  if (is(object, "MizerSim")) {
    params <- object@params
    no_t <- dim(object@n)[[1]]
    initial_n <- params@initial_n
    initial_n[] <- object@n[no_t, , ]
    initial_n_pp <- params@initial_n_pp
    initial_n_pp[] <- object@n_pp[no_t, ]
    initial_n_other <- object@n_other[no_t, ]
    t_start <- as.numeric(dimnames(object@n)[[1]][[no_t]])
  }
  else if (is(object, "MizerParams")) {
    params <- object
    if (missing(initial_n)) 
      initial_n <- params@initial_n
    if (missing(initial_n_pp)) 
      initial_n_pp <- params@initial_n_pp
    initial_n_other <- params@initial_n_other
  }
  else {
    stop("The `object` argument must be either a MizerParams or a MizerSim object.")
  }
  params@initial_n[] <- initial_n
  params@initial_n_pp[] <- initial_n_pp
  params@initial_n_other <- initial_n_other
  no_sp <- length(params@w_min_idx)
  assert_that(is.array(initial_n), is.numeric(initial_n), are_equal(dim(initial_n), 
                                                                    c(no_sp, length(params@w))))
  assert_that(is.numeric(initial_n_pp), length(initial_n_pp) == 
                length(params@w_full))
  assert_that(is.null(initial_n_other) || is.list(initial_n_other))
  other_names <- names(params@other_dynamics)
  if (length(other_names) > 0) {
    if (is.null(names(initial_n_other))) {
      stop("The initial_n_other needs to be a named list")
    }
    if (!setequal(names(initial_n_other), other_names)) {
      stop("The names of the entries in initial_n_other do not match ", 
           "the names of the other components of the model.")
    }
  }
  if (missing(effort)) 
    effort <- params@initial_effort
  if (is.null(dim(effort))) {
    no_gears <- dim(params@catchability)[1]
    if ((length(effort) > 1) & (length(effort) != no_gears)) {
      stop("Effort vector must be the same length as the number of fishing gears\n")
    }
    gear_names <- dimnames(params@catchability)[[1]]
    effort_gear_names <- names(effort)
    if (length(effort) == 1 && is.null(effort_gear_names)) {
      effort_gear_names <- gear_names
    }
    if (!all(gear_names %in% effort_gear_names)) {
      stop("Gear names in the MizerParams object (", 
           paste(gear_names, collapse = ", "), ") do not match those in the effort vector.")
    }
    time_dimnames <- c(t_start, t_start + t_max)
    effort <- t(array(effort, dim = c(no_gears, 2), dimnames = list(gear = effort_gear_names, 
                                                                    time = time_dimnames)))
  }
  no_gears <- dim(params@catchability)[1]
  if (dim(effort)[2] != no_gears) {
    stop("The number of gears in the effort array (length of the second dimension = ", 
         dim(effort)[2], ") does not equal the number of gears in the MizerParams object (", 
         no_gears, ").")
  }
  gear_names <- dimnames(params@catchability)[[1]]
  if (!all(gear_names %in% dimnames(effort)[[2]])) {
    stop("Gear names in the MizerParams object (", 
         paste(gear_names, collapse = ", "), ") do not match those in the effort array.")
  }
  effort <- effort[, gear_names, drop = FALSE]
  if (is.null(dimnames(effort)[[1]])) {
    stop("The time dimname of the effort argument must be numeric.")
  }
  time_effort <- as.numeric(dimnames(effort)[[1]])
  if (any(is.na(time_effort))) {
    stop("The time dimname of the effort argument must be numeric.")
  }
  if (is.unsorted(time_effort)) {
    stop("The time dimname of the effort argument should be increasing.")
  }
  t_start <- time_effort[[1]]
  t_end <- time_effort[[length(time_effort)]]+(1) ### run the last year!!
  t_max <- t_end - t_start
  if (t_max < t_save) {
    t_save <- t_max
  }
  if ((t_save < dt) || !isTRUE(all.equal((t_save - round(t_save/dt) * 
                                          dt), 0))) 
    stop("t_save must be a positive multiple of dt")
  skip <- round(t_save/dt)
  t_dimnames <- seq(t_start, t_end, by = t_save)
 # browser()
  sim <- MizerSim(params, t_dimnames = t_dimnames)
  resource_dynamics_fn <- get(params@resource_dynamics)
  other_dynamics_fns <- lapply(params@other_dynamics, get)
  rates_fns <- lapply(params@rates_funcs, get)
  if (progress_bar == TRUE) {
    pb <- progress::progress_bar$new(format = "[:bar] :percent ETA: :eta", 
                                     total = length(t_dimnames), width = 60)
    pb$tick(0)
  }
  if (is(progress_bar, "Progress")) {
    progress_bar$set(message = "Running simulation", 
                     value = 0)
    proginc <- 1/length(t_dimnames)
  }
  t <- t_start
  i_save_time <- 2
  i_effort <- 1
  current_effort <- effort[1, ]
  t_next_effort <- time_effort[[2]] - 1e-08
  n_list <- list(n = initial_n, n_pp = initial_n_pp, n_other = initial_n_other)
  sim@n[1, , ] <- initial_n
  sim@n_pp[1, ] <- initial_n_pp
  sim@n_other[1, ] <- initial_n_other
  sim@effort[1, ] <- current_effort
  r_w <- (params@w_full < params@resource_params$w_pp_cutoff)
  params@interaction[,] <- inter_mat[,,1]## to add the interaction matrix
  params@cc_pp[r_w] <- kappa_dyn[1]*(params@w_full[r_w])^(-params@resource_params$lambda) ## resource dynamics
  for (i in 2:length(t_dimnames)) {
    n_list <- project_simple(params, n = n_list$n, n_pp = n_list$n_pp, 
                             n_other = n_list$n_other, t = t, dt = dt, steps = skip, 
                             effort = current_effort, resource_dynamics_fn = resource_dynamics_fn, 
                             other_dynamics_fns = other_dynamics_fns, rates_fns = rates_fns)
    t <- t + t_save
    if (is(progress_bar, "Progress")) {
      progress_bar$inc(amount = proginc)
    }
    if (progress_bar == TRUE) {
      pb$tick()
    }
    sim@n[i, , ] <- n_list$n
    sim@n_pp[i, ] <- n_list$n_pp
    sim@n_other[i, ] <- n_list$n_other
    sim@effort[i, ] <- current_effort
    if (t >= t_next_effort) {
      i_effort <- i_effort + 1
      current_effort <- effort[i_effort, ]
      params@interaction[,] <- inter_mat[,,i_effort] ## to add the interaction matrix
      params@cc_pp[r_w] <- kappa_dyn[i_effort]*(params@w_full[r_w])^(-params@resource_params$lambda) ## resource dynamics
      t_next_effort <- ifelse(i_effort < length(time_effort), 
                              time_effort[[i_effort + 1]] - 1e-08, t_max + 
                                10000)
    }
  }
  if (is(object, "MizerSim") && append) {
    no_t_old <- dim(object@n)[1]
    no_t <- length(t_dimnames)
    new_t_dimnames <- c(as.numeric(dimnames(object@n)[[1]]), 
                        t_dimnames[2:no_t])
    new_sim <- MizerSim(params, t_dimnames = new_t_dimnames)
    old_indices <- 1:no_t_old
    new_indices <- seq(from = no_t_old + 1, length.out = no_t - 
                         1)
    new_sim@n[old_indices, , ] <- object@n
    new_sim@n[new_indices, , ] <- sim@n[2:no_t, , ]
    new_sim@n_pp[old_indices, ] <- object@n_pp
    new_sim@n_pp[new_indices, ] <- sim@n_pp[2:no_t, ]
    new_sim@n_other[old_indices, ] <- object@n_other
    new_sim@n_other[new_indices, ] <- sim@n_other[2:no_t, 
                                                  ]
    new_sim@effort[old_indices, ] <- object@effort
    new_sim@effort[new_indices, ] <- sim@effort[2:no_t, ]
    return(new_sim)
  }
  return(sim)
}
#############
##
############

getCatch <- function(sim,t_save=0.1){
  yield <- getYield(sim)
  t(simplify2array(by(yield,floor(as.numeric(rownames(yield))-t_save),colMeans)))
}
#############
##
############


getCatchGear_gov <- function(sim,t_save=0.1){
  yield <- getYieldGear_gov(sim)
  t(simplify2array(by(yield,floor(as.numeric(rownames(yield))-t_save),colMeans)))
}
#############
##
############


getYieldGear_gov<-function (sim) 
{
  biomass <- sweep(sim@n, 3, sim@params@w * sim@params@dw, 
                   "*")
  f_gear <- getFMortGear_gov(sim)
  yield <- apply(f_gear  * biomass, c(1, 2), sum)
  return(yield)
}
#############
##
############


getFMortGear_gov <- function(sim){
  effort <- sim@effort[,71]
  sel <- sim@params@selectivity[71,,]
  cat <- sim@params@catchability[71,]
  
  ret <- array(0,dim=c(length(effort),length(cat),ncol(sel)))
  for(i in 1:length(effort)){
    ret[i,,] <- effort[i] * sel * cat
  }
  return(ret)
}

########
##
########


update_waa <- function(waa,grow,spec,dt,params){
  g_fn <- stats::approxfun(c(params@w, params@species_params$w_inf[spec]),
                           c(grow[spec,], 0))
  myodefun <- function(t, state, parameters) {
    return(list(g_fn(state)))
  }
  deSolve::ode(y = waa[spec,], times = c(0,dt), func = myodefun)[2,-1]
}

##############
##
##############

project_dyn_waa <- function (object, effort, t_max = 100, dt = 0.1, t_save = 1, 
                         t_start = 0, initial_n, initial_n_pp, append = TRUE, progress_bar = TRUE, inter_mat,kappa_dyn,waa_initial,nages=20, ...) 
{
  validObject(object)
  if (is(object, "MizerSim")) {
    params <- object@params
    no_t <- dim(object@n)[[1]]
    initial_n <- params@initial_n
    initial_n[] <- object@n[no_t, , ]
    initial_n_pp <- params@initial_n_pp
    initial_n_pp[] <- object@n_pp[no_t, ]
    initial_n_other <- object@n_other[no_t, ]
    t_start <- as.numeric(dimnames(object@n)[[1]][[no_t]])
  }
  else if (is(object, "MizerParams")) {
    params <- object
    if (missing(initial_n)) 
      initial_n <- params@initial_n
    if (missing(initial_n_pp)) 
      initial_n_pp <- params@initial_n_pp
    initial_n_other <- params@initial_n_other
  }
  else {
    stop("The `object` argument must be either a MizerParams or a MizerSim object.")
  }
  params@initial_n[] <- initial_n
  params@initial_n_pp[] <- initial_n_pp
  params@initial_n_other <- initial_n_other
  ##
  if(missing(waa_initial)){
    waa_initial <- matrix(1:nages,nrow(params@species_params),nages,byrow=T)
  }
  ##
  no_sp <- length(params@w_min_idx)
  assert_that(is.array(initial_n), is.numeric(initial_n), are_equal(dim(initial_n), 
                                                                    c(no_sp, length(params@w))))
  assert_that(is.numeric(initial_n_pp), length(initial_n_pp) == 
                length(params@w_full))
  assert_that(is.null(initial_n_other) || is.list(initial_n_other))
  other_names <- names(params@other_dynamics)
  if (length(other_names) > 0) {
    if (is.null(names(initial_n_other))) {
      stop("The initial_n_other needs to be a named list")
    }
    if (!setequal(names(initial_n_other), other_names)) {
      stop("The names of the entries in initial_n_other do not match ", 
           "the names of the other components of the model.")
    }
  }
  if (missing(effort)) 
    effort <- params@initial_effort
  if (is.null(dim(effort))) {
    no_gears <- dim(params@catchability)[1]
    if ((length(effort) > 1) & (length(effort) != no_gears)) {
      stop("Effort vector must be the same length as the number of fishing gears\n")
    }
    gear_names <- dimnames(params@catchability)[[1]]
    effort_gear_names <- names(effort)
    if (length(effort) == 1 && is.null(effort_gear_names)) {
      effort_gear_names <- gear_names
    }
    if (!all(gear_names %in% effort_gear_names)) {
      stop("Gear names in the MizerParams object (", 
           paste(gear_names, collapse = ", "), ") do not match those in the effort vector.")
    }
    time_dimnames <- c(t_start, t_start + t_max)
    effort <- t(array(effort, dim = c(no_gears, 2), dimnames = list(gear = effort_gear_names, 
                                                                    time = time_dimnames)))
  }
  no_gears <- dim(params@catchability)[1]
  if (dim(effort)[2] != no_gears) {
    stop("The number of gears in the effort array (length of the second dimension = ", 
         dim(effort)[2], ") does not equal the number of gears in the MizerParams object (", 
         no_gears, ").")
  }
  gear_names <- dimnames(params@catchability)[[1]]
  if (!all(gear_names %in% dimnames(effort)[[2]])) {
    stop("Gear names in the MizerParams object (", 
         paste(gear_names, collapse = ", "), ") do not match those in the effort array.")
  }
  effort <- effort[, gear_names, drop = FALSE]
  if (is.null(dimnames(effort)[[1]])) {
    stop("The time dimname of the effort argument must be numeric.")
  }
  time_effort <- as.numeric(dimnames(effort)[[1]])
  if (any(is.na(time_effort))) {
    stop("The time dimname of the effort argument must be numeric.")
  }
  if (is.unsorted(time_effort)) {
    stop("The time dimname of the effort argument should be increasing.")
  }
  t_start <- time_effort[[1]]
  t_end <- time_effort[[length(time_effort)]] + 1
  t_max <- t_end - t_start
  if (t_max < t_save) {
    t_save <- t_max
  }
  if ((t_save < dt) || !isTRUE(all.equal((t_save - round(t_save/dt) * 
                                          dt), 0))) 
    stop("t_save must be a positive multiple of dt")
  skip <- round(t_save/dt)
  t_dimnames <- seq(t_start, t_end, by = t_save)
  # browser()
  sim <- MizerSim(params, t_dimnames = t_dimnames)
  waa_save <- array(0,dim = c(length(t_dimnames),nrow(waa_initial),ncol(waa_initial)))
  mort_save <- sim@n
  mort_save[] <- 0
  #### add names
  resource_dynamics_fn <- get(params@resource_dynamics)
  other_dynamics_fns <- lapply(params@other_dynamics, get)
  rates_fns <- lapply(params@rates_funcs, get)
  if (progress_bar == TRUE) {
    pb <- progress::progress_bar$new(format = "[:bar] :percent ETA: :eta", 
                                     total = length(t_dimnames), width = 60)
    pb$tick(0)
  }
  if (is(progress_bar, "Progress")) {
    progress_bar$set(message = "Running simulation", 
                     value = 0)
    proginc <- 1/length(t_dimnames)
  }
  t <- t_start
  i_save_time <- 2
  i_effort <- 1
  current_effort <- effort[1, ]
  #t_next_effort <- time_effort[[2]] - 1e-08
  #
  t_next_effort <- ifelse(i_effort < length(time_effort), 
                          time_effort[[i_effort + 1]] - 1e-08, t_max + 
                            10000)
  n_list <- list(n = initial_n, n_pp = initial_n_pp, n_other = initial_n_other)
  sim@n[1, , ] <- initial_n
  sim@n_pp[1, ] <- initial_n_pp
  sim@n_other[1, ] <- initial_n_other
  sim@effort[1, ] <- current_effort
  waa_save[1,,] <- waa_initial
  waa <- waa_initial
  mort_save[1,,] <- 0
  r_w <- (params@w_full < params@resource_params$w_pp_cutoff)
  params@interaction[,] <- inter_mat[,,1]## to add the interaction matrix
  params@cc_pp[r_w] <- kappa_dyn[1]*(params@w_full[r_w])^(-params@resource_params$lambda) ## resource dynamics
  for (i in 2:length(t_dimnames)) {
    ##
    
    ## 
    n_list <- project_simple_waa(params, n = n_list$n, n_pp = n_list$n_pp, 
                             n_other = n_list$n_other, t = t, dt = dt, steps = skip, 
                             effort = current_effort, resource_dynamics_fn = resource_dynamics_fn, 
                             other_dynamics_fns = other_dynamics_fns, rates_fns = rates_fns,waa=waa)
    t <- t + t_save
    if (is(progress_bar, "Progress")) {
      progress_bar$inc(amount = proginc)
    }
    if (progress_bar == TRUE) {
      pb$tick()
    }
    sim@n[i, , ] <- n_list$n
    sim@n_pp[i, ] <- n_list$n_pp
    sim@n_other[i, ] <- n_list$n_other
    sim@effort[i, ] <- current_effort
    waa_save[i,,] <- n_list$waa
    mort_save[i,,] <- n_list$rates$mort
    waa <- n_list$waa
    #print(waa[1,1])
    if (t >= t_next_effort) {
      i_effort <- i_effort + 1
      params@interaction[,] <- inter_mat[,,i_effort] ## to add the interaction matrix
      params@cc_pp[r_w] <- kappa_dyn[i_effort]*(params@w_full[r_w])^(-params@resource_params$lambda) ## resource dynamics
      current_effort <- effort[i_effort, ]
      t_next_effort <- ifelse(i_effort < length(time_effort), 
                              time_effort[[i_effort + 1]] - 1e-08, t_max + 
                                10000)
      waa <- cbind(params@species_params$w_min,waa[,-(ncol(waa))])
    }
  }
  if (is(object, "MizerSim") && append) {
    no_t_old <- dim(object@n)[1]
    no_t <- length(t_dimnames)
    new_t_dimnames <- c(as.numeric(dimnames(object@n)[[1]]), 
                        t_dimnames[2:no_t])
    new_sim <- MizerSim(params, t_dimnames = new_t_dimnames)
    old_indices <- 1:no_t_old
    new_indices <- seq(from = no_t_old + 1, length.out = no_t - 
                         1)
    new_sim@n[old_indices, , ] <- object@n
    new_sim@n[new_indices, , ] <- sim@n[2:no_t, , ]
    new_sim@n_pp[old_indices, ] <- object@n_pp
    new_sim@n_pp[new_indices, ] <- sim@n_pp[2:no_t, ]
    new_sim@n_other[old_indices, ] <- object@n_other
    new_sim@n_other[new_indices, ] <- sim@n_other[2:no_t, 
    ]
    new_sim@effort[old_indices, ] <- object@effort
    new_sim@effort[new_indices, ] <- sim@effort[2:no_t, ]
    return(list(sim=new_sim,waa=waa_save))
  }
  return(list(sim=sim,waa=waa_save,mort=mort_save))
}

###############
##
###############

project_simple_waa <- 
  function(params, 
           n = params@initial_n,
           n_pp = params@initial_n_pp,
           n_other = params@initial_n_other,
           effort = params@initial_effort,
           t = 0, dt = 0.1, steps,
           resource_dynamics_fn = get(params@resource_dynamics),
           other_dynamics_fns = lapply(params@other_dynamics, get),
           rates_fns = lapply(params@rates_funcs, get),
           waa, ...) {    
    # Handy things ----
    no_sp <- nrow(params@species_params) # number of species
    no_w <- length(params@w) # number of fish size bins
    idx <- 2:no_w
    # Hacky shortcut to access the correct element of a 2D array using 1D notation
    # This references the egg size bracket for all species, so for example
    # n[w_min_idx_array_ref] = n[,w_min_idx]
    w_min_idx_array_ref <- (params@w_min_idx - 1) * no_sp + (1:no_sp)
    # Matrices for solver
    a <- matrix(0, nrow = no_sp, ncol = no_w)
    b <- matrix(0, nrow = no_sp, ncol = no_w)
    S <- matrix(0, nrow = no_sp, ncol = no_w)
    
    # Loop over time steps ----
    for (i_time in 1:steps) {
      r <- rates_fns$Rates(
        params, n = n, n_pp = n_pp, n_other = n_other,
        t = t, effort = effort, rates_fns = rates_fns, ...)
      # * Update other components ----
      n_other_current <- n_other  # So that the resource dynamics can still 
      # use the current value
      for (component in names(params@other_dynamics)) {
        n_other[[component]] <-
          other_dynamics_fns[[component]](
            params,
            n = n,
            n_pp = n_pp,
            n_other = n_other_current,
            rates = r,
            t = t,
            dt = dt,
            component = component,
            ...
          )
      }
      
      # * Update resource ----
      n_pp <- resource_dynamics_fn(params, n = n, n_pp = n_pp,
                                   n_other = n_other_current, rates = r,
                                   t = t, dt = dt,
                                   resource_rate = params@rr_pp,
                                   resource_capacity = params@cc_pp, ...)
      
      # * Update species ----
      a[, idx] <- sweep(-r$e_growth[, idx - 1, drop = FALSE] * dt, 2,
                        params@dw[idx], "/")
      b[] <- 1 + sweep(r$e_growth * dt, 2, params@dw, "/") + r$mort * dt
      S[,idx] <- n[, idx, drop = FALSE]
      # Update first size group of n
      n[w_min_idx_array_ref] <-
        (n[w_min_idx_array_ref] + r$rdd * dt / 
           params@dw[params@w_min_idx]) /
        b[w_min_idx_array_ref]
      
      ### -- update ages
      waa <- t(sapply(1:no_sp,update_waa,grow=r$e_growth,waa=waa,dt=dt,params=params))
      
      # Update n
      # for (i in 1:no_sp){ # number of species assumed small, so no need to vectorize this loop over species
      #     for (j in (params@w_min_idx[i]+1):no_w){
      #         n[i,j] <- (S[i,j] - a[i,j]*n[i,j-1]) / b[i,j]}}
      # This is implemented via Rcpp
      n <- inner_project_loop(no_sp = no_sp, no_w = no_w, n = n,
                              A = a, B = b, S = S,
                              w_min_idx = params@w_min_idx)

      # * Update time ----
      t <- t + dt
    }
    
    return(list(n = n, n_pp = n_pp, n_other = n_other, rates = r, waa = waa))
  }

#########
##
#########

getNumbers <- function (sim){
  numbers <- sweep(sim@n, 3, sim@params@dw, 
                   "*")
  return(numbers)
}

getNumbers_mat <- function (sim){
  numbers <- sweep(sweep(sim@n, c(2, 3), sim@params@maturity, 
                               "*"), 3, sim@params@dw, "*") 
  return(numbers)
}

##########
####
##########

get_numbers_at_age <- function(input,max_ages=(dim(input[[2]])[3]-1),mature=F){
  sim <- input$sim
  waa <- input$waa
  years <- dim(waa)[1]
  nspec <- nrow(sim@params@species_params)
  ret <- array(0,dim=c(years,nspec,max_ages+2),dimnames = list(rownames(sim@n),colnames(sim@n),c(0:max_ages,"plus")))
  tmp <- rep(0,max_ages+2)
  if (mature==T){
    num_at_age <- getNumbers_mat(sim)
  }
  else{
    num_at_age <- getNumbers(sim)
  }
  for(i in 1:years){
    for(j in 1:nspec){
      for(k in 1:(max_ages+1)){
        tmp[k] <- sum(num_at_age[i,j,which(sim@params@w < waa[i,j,k])])
      }
      tmp[max_ages+2] <- sum(num_at_age[i,j,]) ## plus group
      ret[i,j,1] <- tmp[1]
      ret[i,j,-1] <- diff(tmp)
    }
  }
  return(ret)
}

###########
##
###########
getCatchNumbers <- function (sim,t_save){
  numbers <- sweep(sim@n, 3, sim@params@dw, 
                   "*")
  f <- getFMort(sim, drop = FALSE) * t_save
  yield <- f * numbers#apply(f * numbers, c(1, 2), sum)
  return(yield)
}

##########
###
#########


getCatchBiomass <- function (sim,t_save){
  biomass <- sweep(sim@n, 3,sim@params@w * sim@params@dw, 
                   "*")
  f <- getFMort(sim, drop = FALSE) * t_save
  yield <- f * biomass#apply(f * numbers, c(1, 2), sum)
  return(yield)
}

##########
##
##########

get_catch_at_age <- function(input,t_save,max_ages=(dim(input[[2]])[3]-1),biomass=F){
  sim <- input$sim
  waa <- input$waa
  years <- dim(waa)[1]
  nspec <- nrow(sim@params@species_params)
  ret <- array(0,dim=c(years,nspec,max_ages+2),dimnames = list(rownames(sim@n),colnames(sim@n),c(0:max_ages,"plus")))
  tmp <- rep(0,max_ages+2)
  if(biomass==T){
    cat_at_age <- getCatchBiomass(sim,t_save)
  }else{
    cat_at_age <- getCatchNumbers(sim,t_save)
  }
  for(i in 1:years){
    for(j in 1:nspec){
      for(k in 1:(max_ages+1)){
        tmp[k] <- sum(cat_at_age[i,j,which(sim@params@w < waa[i,j,k])])
      }
      tmp[max_ages+2] <- sum(cat_at_age[i,j,]) ## plus group
      ret[i,j,1] <- tmp[1]
      ret[i,j,-1] <- diff(tmp)
    }
  }
  return(ret)
}

##########
##
##########

getGrowthCurves1 <- 
function (object, species = NULL, max_age = 20, percentage = FALSE) 
{
  if (is(object, "MizerSim")) {
    params <- object@params
    params <- setInitialValues(params, object)
  }
  else if (is(object, "MizerParams")) {
    params <- validParams(object)
  }
  else {
    stop("The first argument to `getGrowthCurves()` must be a ", 
         "MizerParams or a MizerSim object.")
  }
  species <- valid_species_arg(params, species)
  idx <- which(params@species_params$species %in% species)
  species <- params@species_params$species[idx]
  age <- 0:max_age
  ws <- array(dim = c(length(species), length(age)), dimnames = list(Species = species, 
                                                                     Age = age))
  g <- getEGrowth(params)
  for (j in seq_along(species)) {
    i <- idx[j]
    g_fn <- stats::approxfun(c(params@w, params@species_params$w_inf[[i]]), 
                             c(g[i, ], 0))
    myodefun <- function(t, state, parameters) {
      return(list(g_fn(state)))
    }
    ws[j, ] <- deSolve::ode(y = params@w[params@w_min_idx[i]], 
                            times = age, func = myodefun)[, 2]
    if (percentage) {
      ws[j, ] <- ws[j, ]/params@species_params$w_inf[i] * 
        100
    }
  }
  return(ws)
}

########
##
########

getMF <- function(f_a_w,n_a_w,mort){
  Zs = -log(sum(n_a_w * exp(-mort))) + log(sum(n_a_w))
  ## 
  Fs <- -log(sum(n_a_w * exp(-f_a_w))) + log(sum(n_a_w))
  Ms <- Zs - Fs
  return(list(Fs=Fs,Ms=Ms))
}

########
##
########

getMortbyAge <- function(input,max_ages=(dim(input[[2]])[3]-1)){
  sim <- input$sim
  waa <- input$waa
  mort <- input$mort
  years <- dim(waa)[1]
  nspec <- nrow(sim@params@species_params)
  n_a_w <- getNumbers(sim)
  f_a_w <- getFMort(sim)
  
  ret_f <- ret_m <- array(0,dim=c(years,nspec,max_ages+2),dimnames = list(rownames(sim@n),colnames(sim@n),c(0:max_ages,"plus")))
  
  for(i in 1:years){
    for(j in 1:nspec){
      ## 1st age group
      w_tmp <- which(sim@params@w <= waa[i,j,1]) ## shouldn't need equal bit but does just incase
      tmp <- getMF(f_a_w[i,j,w_tmp],n_a_w[i,j,w_tmp],mort[i,j,w_tmp])
      ret_f[i,j,1] <- tmp$Fs ; ret_m[i,j,1] <- tmp$Ms
      ## middle groups
      for(k in 2:(max_ages+1)){
        w_tmp <- which(sim@params@w <= waa[i,j,k] & sim@params@w > waa[i,j,k-1])
        if (length(w_tmp) > 0){
          tmp <- getMF(f_a_w[i,j,w_tmp],n_a_w[i,j,w_tmp],mort[i,j,w_tmp])
          ret_f[i,j,k] <- tmp$Fs ; ret_m[i,j,k] <- tmp$Ms
        }
      }
      ## plus groups
      w_tmp <- which(sim@params@w > waa[i,j,max_ages+1] )
      if (length(w_tmp) > 0){
        tmp <- getMF(f_a_w[i,j,w_tmp],n_a_w[i,j,w_tmp],mort[i,j,w_tmp])
        ret_f[i,j,max_ages+2] <- tmp$Fs ; ret_m[i,j,max_ages+2] <- tmp$Ms
      }
    }
  }
  return(list(Fs=ret_f,Ms=ret_m))
}

#########
##
########

getAnnualRate <- function(input){
  annual_rate <- array(0,dim=c(length(unique(floor(as.numeric(dimnames(input)[[1]]) - 0.1))),dim(input)[2],dim(input)[3]),dimnames=list(unique(floor(as.numeric(dimnames(input)[[1]]) - 0.1)),dimnames(input)[[2]],dimnames(input)[[3]]))
  for (i in 1:(dim(input)[3])){
    annual_rate[,,i] <- t(simplify2array(by(input[,,i],floor(as.numeric(dimnames(input)[[1]]) - 0.1),colMeans)))
  }
  return(annual_rate)
}

#######
##
######

getAnnualsum<- function(input){
  annual_rate <- array(0,dim=c(length(unique(floor(as.numeric(dimnames(input)[[1]]) - 0.1))),dim(input)[2],dim(input)[3]),dimnames=list(unique(floor(as.numeric(dimnames(input)[[1]]) - 0.1)),dimnames(input)[[2]],dimnames(input)[[3]]))
  for (i in 1:(dim(input)[3])){
    annual_rate[,,i] <- t(simplify2array(by(input[,,i],floor(as.numeric(dimnames(input)[[1]]) - 0.1),colSums)))
  }
  return(annual_rate)
}

########
##
########

project_dyn_waa_ad_one_year <- function (object, effort, t_max = 100, dt = 0.1, t_save = 1, 
                             t_start = 0, initial_n, initial_n_pp, append = TRUE, progress_bar = TRUE, inter_mat,kappa_dyn,waa_initial,nages=20,theta,lin=F,def.sel,
                             ...) 
{
  params <- object
  kappa_dyn <- exp(theta[25]+kappa_dyn)
  initial_n_other <- params@initial_n_other
  params@initial_n[] <- initial_n
  params@initial_n_pp[] <- initial_n_pp
  params@initial_n_other <- initial_n_other
  no_sp <- length(params@w_min_idx)

  no_gears <- dim(params@catchability)[1]
  gear_names <- dimnames(params@catchability)[[1]]
  time_effort <- as.numeric(dimnames(effort)[[1]])
  
  current_effort <- rep(0,no_gears)
  names(current_effort) <- gear_names
  
  t_start <- 0 #time_effort[[1]]
  t_end <- 1 #time_effort[[length(time_effort)]] + 1
  t_max <- t_end - t_start
  if (t_max < t_save) {
    t_save <- t_max
  }
  if ((t_save < dt) || !isTRUE(all.equal((t_save - round(t_save/dt) * 
                                          dt), 0))) 
    stop("t_save must be a positive multiple of dt")
  skip <- round(t_save/dt)
  t_dimnames <- seq(t_start, t_end, by = t_save)
  # browser()
  sim <- MizerSim(params, t_dimnames = t_dimnames)
  waa_save <- array(0,dim = c(length(t_dimnames),nrow(waa_initial),ncol(waa_initial)))
  mort_save <- sim@n
  mort_save[] <- 0
  #### add names
  resource_dynamics_fn <- get(params@resource_dynamics)
  other_dynamics_fns <- lapply(params@other_dynamics, get)
  rates_fns <- lapply(params@rates_funcs, get)
  t <- t_start
  i_save_time <- 2
  i_effort <- 1
  n_list <- list(n = initial_n, n_pp = initial_n_pp, n_other = initial_n_other)
  sim@n[1, , ] <- initial_n
  sim@n_pp[1, ] <- initial_n_pp
  sim@n_other[1, ] <- initial_n_other
  
  waa <- waa_initial
  waa <- cbind(params@species_params$w_min,waa[,-(ncol(waa))]) ## it's a new year
  waa_save[1,,] <- waa
  current_effort <- get_mizer_effort(effort,waa,params@selectivity,params@w)
  #browser()
  sim@effort[1, ] <- current_effort
  mort_save[1,,] <- 0
  r_w <- (params@w_full < params@resource_params$w_pp_cutoff)
  params@interaction[,] <- inter_mat## to add the interaction matrix
  params@cc_pp[r_w] <- kappa_dyn*(params@w_full[r_w])^(-params@resource_params$lambda) ## resource dynamics
  for (i in 2:length(t_dimnames)) {
    ##
    
    ## 
    n_list <- project_simple_waa(params, n = n_list$n, n_pp = n_list$n_pp, 
                                 n_other = n_list$n_other, t = t, dt = dt, steps = skip, 
                                 effort = current_effort, resource_dynamics_fn = resource_dynamics_fn, 
                                 other_dynamics_fns = other_dynamics_fns, rates_fns = rates_fns,waa=waa)
    t <- t + t_save
    sim@n[i, , ] <- n_list$n
    sim@n_pp[i, ] <- n_list$n_pp
    sim@n_other[i, ] <- n_list$n_other
    sim@effort[i, ] <- current_effort
    waa_save[i,,] <- n_list$waa
    mort_save[i,,] <- n_list$rates$mort
    waa <- n_list$waa
    current_effort <- get_mizer_effort(effort,waa,params@selectivity,params@w,lin=lin,def.sel=def.sel)

     }
  return(list(sim=sim,waa=waa_save,mort=mort_save))
}


####
##
####
get_mizer_effort<- function(effort,waa,selectivity,w,lin=F,def.sel){
  if(lin==T){
    return(get_mizer_effort_lin(effort,waa,selectivity,w,def.sel))
  }
  return(get_mizer_effort_step(effort,waa,selectivity,w))
}

#####
##
#####

get_mizer_effort_step <- function(effort,waa,selectivity,w){
  ret <- rep(0,nrow(selectivity))
  names(ret) <- row.names(selectivity)
  age <- selectivity[1,,]
  age[,] <- ncol(waa)
  n_s <- nrow(age)
  n_w <- ncol(age)
  for(j in 1:n_s){
    for(k in (ncol(waa)-1):1){
      age[j,which(w < waa[j,k])] <- k
    }
    ret[(j-1)*n_w+(1:n_w)] <- effort[j,age[j,]]
  }
  return(ret)
}

######
##
######

get_mizer_effort_lin <- function(effort,waa,selectivity,w,def.sel){
  ret <- rep(0,nrow(selectivity))
  names(ret) <- row.names(selectivity)
  age <- selectivity[1,,]
  age[,] <- ncol(waa)
  n_s <- nrow(age)
  n_w <- ncol(age)
  for(j in 1:n_s){
    for(k in (ncol(waa)-1):1){
      age[j,which(w < waa[j,k])] <- k
    }
    
    def.sel
    ret[(j-1)*n_w+(1:n_w)] <- effort[j,age[j,]] + get_rat_sel(def.sel,age,j)
  }
  return(ifelse(ret>=0,ret,0))
}

######
##
######

get_rat_sel <- function(def.sel,age,j){
  ret <- def.sel[j,]
  unlist(tapply(ret,age[j,],function(x){x-mean(x)}) )
}

######
##
######

run_to_2019 <- function(theta,inter, dyn_kappa,NS_species,selectivity,catchability,lambda=2.05,Fs,spin_up = 100,max_age=20){
  inter <- inter[,,1:36]
  NS_species$R_max <- exp(theta[1:12])
  NS_species$erepro <- exp(theta[13:24])
  kappa <- theta[25]
  kappa_su <- kappa + theta[40]
  kappa_dyn <- exp(c(rep(kappa_su,spin_up),dyn_kappa + kappa))
  spin_upFs <- 0 * Fs[1,]
  spin_upFs[paste(NS_species$species,20,sep="")] <- theta[26:37]
  ## popes postulate
  Fs[,paste(NS_species$species[c(5,8)],20,sep="")] <- t(t(Fs[,paste(NS_species$species[c(5,8)],20,sep="")]) * theta[38:39])
  Fs <- rbind(matrix(spin_upFs,spin_up,length(spin_upFs),byrow=T),Fs)
  
  spec_overlap1 <- array(inter[,,1],dim=c(12,12,spin_up))
  spec_overlap <- abind::abind(spec_overlap1,inter)
  
  NS_species$interaction_resource <- NS_species$f0 *NS_species$h * NS_species$beta ^(2 - lambda) / ((1 - NS_species$f0)* sqrt(2 *pi) * NS_species$gamma * exp(kappa_su) * NS_species$sigma)
  if(max(NS_species$interaction_resource) > 1 |min(NS_species$interaction_resource) < 0){return(list(catch=rep(-Inf,12),survey=rep(-Inf,12)))}
  params <- newMultispeciesParams(NS_species, spec_overlap[,,1],selectivity = selectivity,catchability=catchability,kappa=exp(kappa_su),lambda=lambda)
  row.names(Fs) <- (2019 - (36 + spin_up - 1)):2019
  test <- project_dyn_waa(params,effort=Fs,kappa_dyn=kappa_dyn,inter_mat = spec_overlap,t_save=0.1,waa_initial=getGrowthCurves1(params,max_age=max_age)[,-1],dt=0.1)
  return(list(params=params,mod=test))
}

########
##
## get rates and stuff
#######

get_transfer <- function(mod){
  num_at_age <- get_numbers_at_age(mod)
  end_of_year <- num_at_age[round(as.numeric(dimnames(num_at_age)[[1]]))==(as.numeric(dimnames(num_at_age)[[1]])),,] ## numbers at age at end of year
  num_at_age_m <- get_numbers_at_age(mod,mature=T)
  end_of_year_m <- num_at_age_m[round(as.numeric(dimnames(num_at_age_m)[[1]]))==(as.numeric(dimnames(num_at_age_m)[[1]])),,]
  mort_by_age <- getMortbyAge(mod)
  break_mort <- lapply(mort_by_age,getAnnualRate)
  cat_num_aa <- get_catch_at_age(mod,t_save=0.1) ## in numbers
  cat_bio_aa <- get_catch_at_age(mod,t_save=0.1,biomass=T) ## in grams
  an_cat_num_aa <-getAnnualsum(cat_num_aa) ## in numbers
  an_cat_bio_aa <-getAnnualsum(cat_bio_aa) ## in grams
  return(list(num_at_age=end_of_year[-1,,],Ms=break_mort$Ms[-1,,],Fs=break_mort$Fs[-1,,],cat_bio_at_age=an_cat_bio_aa[-1,,],cat_num_at_age=an_cat_num_aa[-1,,],mean_waa=an_cat_bio_aa[-1,,]/an_cat_num_aa[-1,,],prop_mat=end_of_year_m[-1,,]/end_of_year[-1,,]))
}

###########
###
###########

append_mods<- function(mod_run,test,dt=0.1){
  tmp1<-max(as.numeric(rownames(mod_run$n_pp)))
  tmp <- test$sim@n[-1,,]
  tmp2 <-test$sim@n_pp[-1,]
  row.names(tmp2) <-row.names(tmp) <- paste(seq(tmp1+dt,tmp1 + (nrow(test$sim@n_pp)-1)*dt,dt))
  mod_run$n <- abind::abind(mod_run$n,tmp,along=1)
  mod_run$n_pp <- rbind(mod_run$n_pp,tmp2)
  return(mod_run)
}

########## 
##
##########

get_initial <- function(kappa_dyn,inter,num,theta,thetas,num_thetas=num,NS_species,Fs,selectivity,catchability,lambda=2.05,spin_up=50,max_age=20){
  if(missing(theta)){
    theta <-as.numeric(thetas[num_thetas,1:40])
  }
  kappa_dyn <- kappa_dyn[num,]
  testing_m <- run_to_2019(theta,inter=inter, dyn_kappa=kappa_dyn,NS_species=NS_species,selectivity=selectivity,catchability=catchability,lambda=lambda,Fs=Fs,spin_up = spin_up,max_age=max_age)
  mod_run <- list(n=testing_m$mod$sim@n[paste(seq(1984,2020,0.1)),,],n_pp=testing_m$mod$sim@n_pp[paste(seq(1984,2020,0.1)),])
  theBoy<-get_transfer(testing_m$mod)
  mod_run <- list(n=testing_m$mod$sim@n[paste(seq(1984,2020,0.1)),,],n_pp=testing_m$mod$sim@n_pp[paste(seq(1984,2020,0.1)),])
  ## edit the selection and catchability 
  params <- testing_m$params
  def.sel <- apply(params@selectivity[paste(params@species_params$species,20,sep=""),,],c(1,3),sum)
  ## change selectivity params to one per weight bin
  n_gear_names <- paste(rep(params@species_params$species,each=100),rep(round(params@w,4),12),sep="")
  params@selectivity <- array(0,dim=c(dim(selectivity)[2]*dim(selectivity)[3],dim(selectivity)[2],dim(selectivity)[3])
                              ,dimnames = list(gear=n_gear_names,sp=dimnames(selectivity)$sp,w=dimnames(selectivity)$w))
  params@catchability <- matrix(0,length(n_gear_names),12,dimnames = list(n_gear_names,dimnames(def.sel)$sp))
  #### sort catchability and selectivity now
  for (i in 1:12){
    for(j in 1:100){
      params@selectivity[paste(params@species_params$species[i],round(params@w[j],4),sep=""),i,j] <- 1
      params@catchability[paste(params@species_params$species[i],round(params@w[j],4),sep=""),i] <- 1
    }
  }
  
  return(list(mod=testing_m$mod,mod_run=mod_run,age_stuff=theBoy,params=params,kappa_dyn=kappa_dyn,inter=inter,theta=theta,def.sel=def.sel))
}

#######
##
#######

progress_one_year <- function(effort,year,mod_run,waa,kappa_dyn,inter,params,theta,prev_run,lin=F,def.sel,append=T){
  if(!missing(prev_run)){
    return(progress_one_year(effort=effort,year=year,
                      mod_run=prev_run$mod_run,waa=tail(prev_run$mod$waa,n=1)[1,,],kappa_dyn=prev_run$kappa_dyn,inter=prev_run$inter,params=prev_run$params,theta=prev_run$theta,lin=lin,def.sel=prev_run$def.sel,append=append
                        ))
  }
  #browser()
  test <- project_dyn_waa_ad_one_year(params,effort=effort,kappa_dyn=kappa_dyn[year],inter_mat = inter[,,year],t_save=0.1,waa_initial=waa,dt=0.1,initial_n=tail(mod_run$n,n=1)[1,,],initial_n_pp=tail(mod_run$n_pp,n=1)[1,],theta=theta,lin=lin,def.sel = def.sel)
  theBoy<-get_transfer(test)
  if(append==T){
    mod_run <- append_mods(mod_run,test)
  }
  return(list(mod=test,mod_run=mod_run,age_stuff=theBoy,params=params,kappa_dyn=kappa_dyn,inter=inter,theta=theta,def.sel=def.sel))
}



######################
##
## Condition FLBiol objects
## using outputs from Mizer fits
##
#####################
# name_list is a list of names from the mizer model and the ices stock names
# they encompass. e.g. "Cod" = "cod.27.47d20"
# spun_up is the burn in time for the mizer model

cond_biols_mizer <- function(name_list, mizer_fits, spin_up = 50, reduce_ages = FALSE)  {
biols <- FLBiols(lapply(seq_len(length(name_list)), function(st) {
# yr range
yrs  <- names(mizer_fits[[1]]$age_stuff$Ms[,names(name_list[st]),1] )[-c(1:spin_up)]
# age range
ages <- seq(0, length(mizer_fits[[1]]$age_stuff$num_at_age[1,names(name_list[st]),])-1,1)
## If the numbers aren't present for all years in a given age, remove these ages
if(reduce_ages) {
ages <- ages[apply(t(mizer_fits[[1]]$age_stuff$num_at_age[ac(2009:2019),names(name_list[st]),]) != 0, 1, any)]
ages <- c(min(ages)):c(max(ages[-length(ages)])+1) ## new maximum age without plus group
}
n    <- FLQuant(dimnames=list(year = yrs, age = ages, iter = seq_len(length(mizer_fits))))
m    <- FLQuant(dimnames=list(year = yrs, age = ages, iter = seq_len(length(mizer_fits))))
wt   <- FLQuant(dimnames=list(year = yrs, age = ages, iter = seq_len(length(mizer_fits))))
spwn <- FLQuant(dimnames=list(year = yrs, age = ages, iter = seq_len(length(mizer_fits))))
mat  <- FLQuant(dimnames=list(year = yrs, age = ages, iter = seq_len(length(mizer_fits))))
fs   <- FLQuant(dimnames=list(year = yrs, age = ages, iter = seq_len(length(mizer_fits))))
for(i in 1:length(mizer_fits)) {
        n[,,,,,i]   <- t(mizer_fits[[i]]$age_stuff$num_at_age[as.character(as.numeric(yrs)+1),names(name_list[st]),c(ages[-length(ages)],"plus")])/1e3
        m[,,,,,i]   <- t(mizer_fits[[i]]$age_stuff$Ms[yrs,names(name_list[st]),c(ages[-length(ages)],"plus")])
        wt[,,,,,i]  <- t(mizer_fits[[i]]$age_stuff$mean_waa[yrs,names(name_list[st]),c(ages[-length(ages)],"plus")])/1e3
        mat[,,,,,i] <- t(mizer_fits[[i]]$age_stuff$prop_mat[yrs,names(name_list[st]),c(ages[-length(ages)],"plus")])
        fs[,,,,,i]  <- t(mizer_fits[[i]]$age_stuff$Fs[yrs,names(name_list[st]),c(ages[-length(ages)],"plus")])
        spwn[]     <- 0 # ?correct
        ## Recalculate Ns so as the start of the year, not the end
          #n[,,,,,i] <- n[,,,,,i]/(exp(-(fs[,,,,,i]+ m[,,,,,i])))
          ## Can't have zeros in an age class, add a small number - but real inconsistency between n's and catches
          ## causing problems with very high q's
          n[,,,,,i][n[,,,,,i]==0] <- 10
          mat[,,,,,i][is.na(mat[,,,,,i])] <- 0
         }
      res <- FLBiol(name = names(name_list[st]),
                                  n = n,
                                  m = m,
                                  wt = wt,
                                  spwn = spwn)
      res@mat$mat[] <- mat
      res@fec$fec[] <- 1
      res@m[is.na(res@m)] <- 0.01 ## ages over those in the stock, need to check whether we should remove overhanging ages
      res@wt[is.na(res@wt)] <- 0.01 ## ages over those in the stock, need to check whether we should remove overhanging ages
      res@n[is.na(res@n)] <- 0 ## ages over those in the stock, need to check whether we should remove overhanging ages
      units(res@n) <- '1000'
      units(res@wt) <- "kg"
      range(res)[["plusgroup"]] <- max(ages)
      return(res)
      }))
  return(biols)
   }


###########################
## function to generate FLFleetsExt
## with mizer outputs and an existing FLFleetsExt
## with catch data
###########################
## OK, so method is:

# Use the fleet catches as proportional, to match the mizer fit
# Spread the plus group fleet catches across the ages above the plus group in mizer proportionally.
# Make an entirely new FLFleetsExt object to fill these with....
# See if any of that makes sense.

# name_list is a list of names from the mizer model and the ices stock names
# they encompass. e.g. "Cod" = "cod.27.47d20"

create_mizer_fleets <- function(fleets, mizer_fits, name_list, reduce_ages = FALSE) {
  
flts <- FLFleetsExt(lapply(fleets, function(fl) {
  print(fl@name)
  mets <- FLMetiersExt(lapply(fl@metiers, function(mt) {
  print(mt@name)
    catches <- FLCatchesExt(lapply(seq_len(length(name_list)), function(st) {
      stk.name <- names(name_list[st])
      stk.ices <- name_list[[st]]
      print(stk.name)
      if(stk.name=="Sandeel") {stk.ices <- stk.ices[1:4]}  ## horrible hack, but only catches in first 4 sandeel stocks
      if(!all(stk.ices %in% mt@catches@names)) {res <- FLCatchExt()} else {
        # yr range
          yrs  <- dimnames(mt@catches[[stk.ices[1]]])$year
          # age range
            ages <- seq(0, length(mizer_fits[[1]]$age_stuff$cat_bio_at_age[1,names(name_list[st]),])-1,1)
          if(reduce_ages) {
            ages <- ages[apply(t(mizer_fits[[1]]$age_stuff$num_at_age[ac(2009:2019),names(name_list[st]),]) != 0, 1, any)]
            ages <- c(min(ages)):c(max(ages[-length(ages)])+1) ## new maximum age without plus group
            }
          ## Set up empty quants
          landings.n  <- discards.n  <- landings.wt <- discards.wt <-
          landings.sel <- discards.sel <- price <- alpha <- beta <-
          FLQuant(dimnames=list(year = yrs, age = ages, iter = seq_len(length(mizer_fits))))
          ## Extract the fleet data
          if(length(stk.ices)==1) {
          price[] <-mt@catches[[stk.ices]]@price[1,]
          ## metier catch
          met_land    <- mt@catches[[stk.ices]]@landings.n ## in 1000s
          met_catch   <- mt@catches[[stk.ices]]@landings.n + mt@catches[[stk.ices]]@discards.n ## in 1000s
          tot_catch   <- catchStock(fleets, stk.ices)
          }
          ## Case of multiple stocks
          ## Need to sum the relevant variables
          if(length(stk.ices) != 1) {
          price <-lapply(stk.ices, function(x) mt@catches[[x]]@price)
          price <- Reduce("+", price)/length(stk.ices)
          price <- propagate(price, iter = length(mizer_fits), fill.iter = TRUE)
          ## metier catch
          met_land    <- Reduce("+",lapply(stk.ices, function(x) mt@catches[[x]]@landings.n))
          met_catch   <- Reduce("+",lapply(stk.ices, function(x) mt@catches[[x]]@landings.n)) +
          Reduce("+",lapply(stk.ices, function(x) mt@catches[[x]]@discards.n))  ## in 1000s
          tot_catch   <- Reduce("+", lapply(stk.ices, function(x) catchStock(fleets, x)))
          }
          ## Case where only one age
          if(length(dimnames(met_catch)$age)==1) {
          ages_fill   <- dimnames(met_catch)$age
          met_land    <- expand(met_land, age = ages)
          met_land[]  <- met_land[ages_fill, ]
          met_catch   <- expand(met_catch, age = ages)
          met_catch[] <- met_catch[ages_fill, ]
          tot_catch   <- expand(tot_catch, age = ages)
          tot_catch[] <- tot_catch[ages_fill, ]
          ls   <- met_land/met_catch
          }
          ## Case where the minimum age is larger then zero
          if(min(as.numeric(dimnames(met_catch)$age)) != 0) {
          met_land  <- expand(met_land, age = 0:max(as.numeric(dimnames(met_land)$age)))
          met_land[is.na(met_land)] <- 0
          met_catch <- expand(met_catch, age = 0:max(as.numeric(dimnames(met_catch)$age)))
          met_catch[is.na(met_catch)] <- 0
          tot_catch <- expand(tot_catch, age = 0:max(as.numeric(dimnames(tot_catch)$age)))
          tot_catch[is.na(tot_catch)] <- 0
          ls        <- met_land/met_catch
          }
          ## Catch where catches greater than our chosen ages
          surplus_ages <- as.numeric(dimnames(met_catch)$age)[!as.numeric(dimnames(met_catch)$age) %in% ages]
          if(length(surplus_ages) != 0) {
          met_land[ac(min(surplus_ages)-1)]  <- apply(met_land[ac(c(min(surplus_ages)-1,surplus_ages))],2:6, sum,na.rm=TRUE)
          met_catch[ac(min(surplus_ages)-1)] <- apply(met_catch[ac(c(min(surplus_ages)-1,surplus_ages))],2:6, sum,na.rm=TRUE)
          tot_catch[ac(min(surplus_ages)-1)] <- apply(tot_catch[ac(c(min(surplus_ages)-1,surplus_ages))],2:6, sum,na.rm=TRUE)
          met_land  <- trim(met_land, age = ages)
          met_catch <- trim(met_catch, age = ages)
          tot_catch <- trim(tot_catch, age = ages)
          }
          ## Determine the landings and discards split
          ls <- met_land / met_catch
          for(i in 1:length(mizer_fits)) {
          ## If the numbers aren't present for all years in a given age, remove these ages
          ## Final group is a plus group, which we account for elsewhere
          #ages <- ages[apply(t(mizer_fits[[1]]$age_stuff$num_at_age[-c(1:spin_up),names(name_list[st]),]) != 0, 1, any)]
          #ages <- c(min(ages)):c(max(ages[-length(ages)])) ## new maximum age without plus group
          ## Mizer catches shared by the proportion of overall catch
          mizer_catch <- t(mizer_fits[[i]]$age_stuff$cat_num_at_age[yrs,stk.name,c(ages[-length(ages)],"plus")])/1e3 ## in 1s, covert to 000s
          ## Reallocate the catch for ages outside the age range....
          mismatch_ages <- ages[!ages %in% dimnames(met_catch)$age]
          match_ages    <- ages[ages %in% dimnames(met_catch)$age]
          if(max(ages) == max(as.numeric(dimnames(met_catch)$age))) {
          match_ages    <- match_ages[-length(match_ages)]
          }
          #[matching ages (except final plus group age, and first age): split proportionally]
          match_ages <- matrix(c((met_catch/tot_catch)[ac(match_ages),]),ncol = length(dimnames(met_catch)$year)) *
          mizer_catch[rownames(mizer_catch) %in% match_ages,]
          #[older: split evenly from the plus group]
          older_ages <- c(mismatch_ages[mismatch_ages > as.numeric(max(dimnames(met_catch)$age))], "plus")
          #old_ages <- t(sapply(matrix(mizer_catch[rownames(mizer_catch) %in% older_ages,],nrow = 1),function(x) x *
          #                      matrix(met_catch[nrow(met_catch),]/tot_catch[nrow(met_catch),], ncol = length(dimnames(met_catch)$year))))
          #if(length(older_ages)==0) {
          #  old_ages <- NULL
          #} else {
          tmp <- matrix(mizer_catch[rownames(mizer_catch) %in% older_ages,],ncol = length(dimnames(met_catch)$year))
          old_ages <- t(t(tmp) * as.numeric(met_catch[nrow(met_catch),]/tot_catch[nrow(met_catch),]))
          #old_ages <- matrix(rep(old_ages, times = length(older_ages)), ncol =   length(dimnames(met_catch)$year))
          #}
          #[0 year olds: split according to the first age group]
          young_ages <- mismatch_ages[mismatch_ages < as.numeric(min(dimnames(met_catch)$age))]
          if(length(young_ages)==0) {
          age_0 <- NULL
          } else {
          tmp <- matrix(mizer_catch[rownames(mizer_catch) %in% young_ages,],ncol = length(dimnames(met_catch)$year))
          age_0 <- t(t(tmp) * as.numeric(met_catch[1,]/tot_catch[1,]))
          }
          ## Full catch-at-age
          caa <- rbind(age_0, match_ages, old_ages)
          ## Expand landings and discards split
          ## Old ages
          if(all(older_ages == "plus")) {
          ls_old <- NULL
          } else {
          ls_old <- old_ages
          if(!is.null(ls_old)) {
          ls_old[] <- matrix(rep(c(ls[nrow(ls),]),times = nrow(old_ages)), ncol = ncol(old_ages), byrow = TRUE)
          }
          }
          # first age
          ls_young   <- age_0
          if(!is.null(ls_young)) {
          ls_young[] <- matrix(rep(c(ls[nrow(ls),]),times = nrow(age_0)), ncol = ncol(age_0), byrow = TRUE)
          }
          ## Combine
          lsa <- rbind(ls_young, matrix(ls, ncol = ncol(ls)), ls_old)
          dsa <- 1-lsa
          # Fill slots
          landings.n[,,,,,i]   <- caa * lsa
          discards.n[,,,,,i]   <- caa * dsa
          landings.sel[,,,,,i] <- lsa
          discards.sel[,,,,,i] <- dsa
          landings.wt[,,,,,i]  <- t(mizer_fits[[i]]$age_stuff$mean_waa[yrs,stk.name,c(ages[-length(ages)],"plus")])/1e3 ## in kg
          discards.wt[,,,,,i] <- t(mizer_fits[[i]]$age_stuff$mean_waa[yrs,stk.name,c(ages[-length(ages)],"plus")])/1e3
          }
          alpha[] <- 1
          beta[] <- 1
          ## Expand price
          price   <- expand(price, age = ages)
          price[] <- price[1,] ## same price for each age
          res <- FLCatchExt(name = stk.name,
                                          landings.n = landings.n,
                                          landings.wt = landings.wt,
                                          landings.sel = landings.sel,
                                          discards.n = discards.n,
                                          discards.wt = discards.wt,
                                          discards.sel = discards.sel,
                                          price = price,
                                          alpha = alpha,
                                          beta = beta)
          ## calculate total landings etc..
          res@landings <- computeLandings(res)
          res@discards <- computeDiscards(res)
          ## Update units
          units(res@landings.n) <- units(res@discards.n) <-'1000'
          units(res@landings.wt) <- units(res@discards.wt) <- 'kg'
          units(res@landings) <- units(res@discards) <- 'tonnes'
          units(res@price) <- 'euros'
          ## fill in plus group details
          range(res)[["plusgroup"]] <- max(ages)
          }
      return(res)
      }))
    ## Remove any catches that are NULL
    stks <- as.character(sapply(catches, name))
    names(catches) <- stks
    stks <- stks[stks != "NA"]
    catches <- catches[stks]
    return(FLMetierExt(name = mt@name,
                              effshare = propagate(mt@effshare, iter = length(mizer_fits), fill.iter = TRUE),
                              vcost = propagate(mt@vcost, iter = length(mizer_fits), fill.iter = TRUE),
                              catches = catches))
    }))
  return(FLFleetExt(name = fl@name,
                    effort = propagate(fl@effort, iter = length(mizer_fits), fill.iter = TRUE),
                    fcost = propagate(fl@fcost, iter = length(mizer_fits), fill.iter = TRUE),
                    capacity = propagate(fl@capacity, iter = length(mizer_fits), fill.iter = TRUE),
                    crewshare = propagate(fl@crewshare, iter = length(mizer_fits), fill.iter = TRUE),
                    metiers = mets))
  }))
return(flts)
}


###################################
## Mizer age-structured population 
## growth model
##
##################################


MizerGrowth <- function(biols, SRs, fleets, year, season, covars, ...) {
  
  yr <- year
  
  ni <- dim(biols[[1]]@n)[6] ## number iterations
            
            for(i in 1:ni) {
              
              for(st in names(biols)) {
                
                ## Update numbers following Fs and Ms in mizer
                #biols[[st]]@n[,yr,,,,i] <- biols[[st]]@n[,yr-1,,,,i] * 
                #  (1-exp(-(covars$mizer_fits[[i]]$age_stuff$Fs[yr,st,] + 
                #          covars$mizer_fits[[i]]$age_stuff$Ms[yr,st,])))
                # case of first year
                if(length(dim(covars$mizer_fits[[1]]$age_stuff$num_at_age))==3) {
                biols[[st]]@n[,yr,,,,i]      <- covars$mizer_fits[[i]]$age_stuff$num_at_age[yr-1,st,]/1e3 ## in 000s
                ## weights, Ms and maturity
                biols[[st]]@wt[,yr,,,,i]      <- covars$mizer_fits[[i]]$age_stuff$mean_waa[yr,st,]/1e3 ## in kg
                biols[[st]]@m[,yr,,,,i]       <- covars$mizer_fits[[i]]$age_stuff$Ms[yr,st,]
                biols[[st]]@mat$mat[,yr,,,,i] <- covars$mizer_fits[[i]]$age_stuff$prop_mat[yr,st,] } else {
			
                # case of subsequent years``:w

                biols[[st]]@n[,yr,,,,i]      <- covars$mizer_fits[[i]]$age_stuff$num_at_age[st,]/1e3 ## in 000s
                ## weights, Ms and maturity
                biols[[st]]@wt[,yr,,,,i]      <- covars$mizer_fits[[i]]$age_stuff$mean_waa[st,]/1e3 ## in kg
                biols[[st]]@m[,yr,,,,i]       <- covars$mizer_fits[[i]]$age_stuff$Ms[st,]
                biols[[st]]@mat$mat[,yr,,,,i] <- covars$mizer_fits[[i]]$age_stuff$prop_mat[st,]   
                  
                }
                
                ## Can't have zeros in numbers so add a small value
                biols[[st]]@n[,yr,,,,i][is.na(biols[[st]]@n[,yr,,,,i])]   <- 10 ## 10,000 fish
                biols[[st]]@n[,yr,,,,i][biols[[st]]@n[,yr,,,,i]==0]       <- 10 ## 10,000 fish
                biols[[st]]@wt[,yr,,,,i][is.na(biols[[st]]@wt[,yr,,,,i])] <- 0
                biols[[st]]@m[,yr,,,,i][is.na(biols[[st]]@m[,yr,,,,i])]   <- 0
                
              }
              
            }
  
  return(list(biols = biols))
            
}


##########################
## Run Mizer
## This works in the covars.om
##
##########################

runMizer <- function(biols, fleets, covars, year) {
  
  yr <- year
  ni <- dim(biols[[1]]@n)[6]
  
 ## First get the Fs from the catches, numbers-at-age and Ms
 catchN <- lapply(biols@names, function(x) catchStock(fleets, x))
 names(catchN) <- names(biols)
 
 # Borrow the FLash function - is quick 
 F_st <- lapply(biols@names, function(st) FLash:::calcF(m(biols[[st]])[,yr],catchN[[st]][,yr],n(biols[[st]][,yr])))
 names(F_st) <- biols@names
 
 ## Now reformat for mizer  
 
 
 F_mizer <- lapply(1:ni, function(i) {
   
   F_length <- max(sapply(F_st, function(x) dim(x)[1]))
   
   F_mat <-t(sapply(names(F_st), function(st) { 
     Fs <- c(F_st[[st]][,,,,,i]) 
     Fs <- c(Fs, rep(Fs[length(Fs)],times = F_length - length(Fs)))
     }))
   
   ## reorder to the mizer inputs
   F_mat <- F_mat[mizer_fits[[1]]$params@species_params$species,]
   F_mat[is.na(F_mat)| F_mat < 0] <- 0
   
   return(F_mat)
   
 }) 
 
 covars$mizer_fits <- lapply(1:ni, function(m) {
   progress_one_year(year=yr,effort=F_mizer[[m]],prev_run=mizer_fits[[m]])
 })
 
 return(covars)
  
}

#######################
## Need perfect Obs Mizer
## function here!!
##
## This will update the biols
## with the mizer outputs
#######################


