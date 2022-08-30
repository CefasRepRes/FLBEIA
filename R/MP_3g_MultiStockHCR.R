#-------------------------------------------------------------------------------
#  MultiStock HCR based on the ICES MSY framework HCR.
#
# Biomass Based HCR.
#  Reference Points: Btrigger, Blim and Fmsy.
#   No formal porposal for any of this, usually :
#           - Btrigger = Bpa
#           - Blim  = ??? Blim?YYY?
#           - Fmsy = F0.1, Fmax....
#
#  - TAC advice depending on F in relation to BRP of ALL the stocks:
#     1. First for each stock we calculate the single stock  Ftarget depending on 
#       its status in relation to the BRPs.
#           - Ftarget = Fmsy             B >= Btrigger.
#           - Ftarget = Fmsy*B/Btrigger  B <  Btrigger.
#           - Ftarget = 0.               B <  Blim (our proposal)
#     2.Calculate the ratio between Ftarget and Fsq and calculate the maximum:
#           Fadv0[st] = lambda0*Fsq : lambda0 = max(Ftarget/Fsq)
#       (=> there is only one stock for which Fadv0 = Ftarget and for the rest Fadv0 > Ftarget)

#     3.Calculate the ratio between Fupp and Fsq and calculate the minimum:
#            x_st = Fupp[st]/Fadv0[st] for all the stocks st.
#              -  If x_st >= 1 for all the stocks: lambda1 = 1 
#              -  If exist st : x_st < 1 => lambda1 = min(Fupp[st]/Fadv0[st])
#                                           and Fadv1[st] = lambda1*Fadv0 :
#     (=> there is only one stock for which Fadv0 = Ftarget and for the rest Fadv0 > Ftarget)
#
#     4.Fadv[st] = lambda1*lambda0*Fsq => TAC[st] = C[Fadv[st]]
#
#     In relation to IcesMSY HCR this new HCR has two additional arguments:
#       * advice.ctrl[['stocksInHCR']] : A vector with the name of the stocks that are taken into account in the calculation of advice.
#       * advice.ctrl[['stknm']][['ref.pts']]: A new row in the matrix with Fupp value.
#
#
# 08/06/2016 - Created:  Dorleta Garcia.
# 18/06/2016 - Modified: Dorleta Garcia.
# 20/12/2016 - Modified: Dorleta Garcia. Incorporation of Data Limited Stocks into the HCR.
# 03/09/2018 - Modified: Dorleta Garcia. Add new settings to the HCR mean and min.
# 13/12/2018 - Modified: Dorleta Garcia. When using the min  option check that all the Fadv are above the 
#                                        lower range, if not the F-s are moved trying to put them all within
#                                        the ranges, always maintaining the proportionality
# 07/02/2022 - Modified: Matthew Pace. Implemented short term forecast to each stock 
#                                      before applying ICES rule for compatability with
#                                      single stock methods. Additional annotation of code
#
#---------------------------------------------------------------------------------------------------
#' @rdname annualTAC
#' @aliases MultiStockHCR
MultiStockHCR <- function(stocks, indices, advice, advice.ctrl, year, stknm,...){
  
  #------------------------------------------------------------------#
  # SECTION 1 - Carry out initial forecast to get SSB in TAC year    #
  #------------------------------------------------------------------#
  #
  # First, simultaneously-managed stocks are identified. A common set of 
  # intermediate year assumptions are defined. Then, for each stock (within a loop), 
  # the corresponding intermediate year F is calculated and a forward control
  # object is built.
  
  ## Extract vector of simultaneously managed stocks and management category 
  stocksInHCR    <- advice.ctrl[['stocksInHCR']]
  stocksCat      <- advice.ctrl[['stocksCategory']]
  approach       <- ifelse(is.null(advice.ctrl[['approach']]), 'max', advice.ctrl[['approach']])
  
  # Define parameters to project the stock 3 years
  # (current year, TAC year, TAC year + 1 for ssb or biomass constraints).
  # These are hard-coded because conditions must be identical for all stocks
  nyears      <- 3
  wts.nyears  <- 3
  fbar.nyears <- 3
  f.rescale   <- ifelse(is.null(advice.ctrl[[stknm]][['f.rescale']]), 
                        TRUE, 
                        advice.ctrl[[stknm]][['f.rescale']])
  
  ## Define dimension length
  iter <- dim(advice.ctrl[[stknm]]$ref.pts)[2]
  
  ## Prepare objects to store F reference points
  Ftg <- Fupp  <- Flow <- Fsq <- Cst <- matrix(NA, length(stocksInHCR), iter)
  rownames(Ftg) <- rownames(Fupp) <- rownames(Flow) <- rownames(Fsq) <- rownames(Cst) <- stocksInHCR
  
  ## loop over each simultaneously-managed stock
  for(st in stocksInHCR){
    
    #------------------------------------------#
    # SECTION 1.1 - Prepare stock for forecast #
    #------------------------------------------#
    
    # The first step is to extract data for the i'th stock and prepare the object
    # for a short-term forecast (to estimate SSB in the TAC year)
    
    ## extract the i'th stock
    stk0 <- stocks[[st]]
    
    ## Fill the 0-s and NA-s with almost 0 values to avoid problems when the 
    ## fishery is closed for example, or there is no catch...
    stk0@harvest[stk0@harvest <= 0] <- 0.00001
    
    ## infer whether age structure is present
    ageStruct <- ifelse(dim(stk0@m)[1] > 1, TRUE, FALSE)
    
    ## If current stock is category 1 (assessment includes population dynamics) then
    ## extend stock data in preparation for a three year forecast 
    ## (current year, TAC year, TAC+1 year)
    if(stocksCat[st] == 1) {
      stk0 <- stf(stk0, 
                  nyears      = nyears, 
                  wts.nyears  = wts.nyears, 
                  fbar.nyears = fbar.nyears, 
                  f.rescale   = f.rescale)
    }
    
    ## extract reference points for the i'th stock
    ref.pts_st <- advice.ctrl[[st]][['ref.pts']] # matrix[6,it]  
    
    iter     <- dim(stk0@m)[6]
    yrsnames <- dimnames(stk0@m)[[2]]
    yrsnumbs <- as.numeric(yrsnames)
    
    assyrname <- yrsnames[year]
    assyrnumb <- yrsnumbs[year]
    
    # The next step is to construct a forward control object to facilitate the
    # short-term forecast of the stock
    
    ## extract intermediate year setting - I'm  keeping this flexible for each
    ## stock rather than a single control for all stocks
    ## Use F status quo as an intermediate year condition if none is specified
    int.yr <- advice.ctrl[[st]]$intermediate.year
    int.yr <- ifelse(is.null(int.yr), 'Fsq', int.yr)
    
    ## extract Ftg (first using Fmsy)
    Ftgi <- ref.pts_st['Fmsy',]
    
    #----------------------------------------------------------#
    # SECTION 1.2 - Carry out initial forecast (if Category 1) #
    #----------------------------------------------------------#
    
    ## If stocks are Category 1
    if(stocksCat[st] == 1){
      
      ## At this point we need to loop over iterations
      for(i in 1:iter){
        
        # SECTION 1.2.1 - Build fwd.ctrl object
        #---------------------------------------#
        #
        # Having extracted stock reference points, we now generate the stock
        # forward control object. This is typically driven by intermediate year
        # assumptions of F status quo or TAC.
        
        ## Generate fwd control object
        if(int.yr == 'Fsq') {
          
          ## Extract iteration data
          stki <- iter(stk0, i)
          
          # Calculate Fsq
          if(fbar.nyears == 1 | f.rescale) {
            Fsq[st,i] <- mean(fbar(stki)[,(year-1)]) 
          } else {
            Fsq[st,i] <- mean(fbar(stki)[,(year-fbar.nyears):(year-1)])
          }
          
          fwd.ctrl <- FLash::fwdControl(data.frame(year     = c(0, 1),  
                                                   val      = c(Fsq[st,i], 
                                                                Ftgi[i]), 
                                                   quantity = c( 'f', 'f'), 
                                                   rel.year = c(NA,NA))) 
        } else {
          
          fwd.ctrl <- FLash::fwdControl(data.frame(year     = c(0, 1), 
                                                   val      = c(advice$TAC[st,
                                                                           year, 
                                                                           drop=TRUE][i], 
                                                                Ftgi[i]), 
                                                   quantity = c( 'catch', 'f')))
        }
        
        # Refresh the years in fwd!!
        fwd.ctrl@target$year     <- fwd.ctrl@target$year + assyrnumb
        fwd.ctrl@target$rel.year <- fwd.ctrl@target$rel.year + assyrnumb
        
        # if in <year 0> quantity = catch => set TAC in <year 0> in val
        if(fwd.ctrl@target[fwd.ctrl@target$year == assyrnumb,'quantity'] == 'catch'){
          k <- which(fwd.ctrl@target$year == assyrnumb)
          fwd.ctrl@target[k,'val']     <- advice$TAC[st, year,,,,i]
          fwd.ctrl@trgtArray[k, 'val',] <- advice$TAC[st ,year,,,,i]
        }
        
        # SECTION 1.2.2 - Forecast age-structured stocks
        #-----------------------------------------------#
        #
        # This next section of code is executed if stock dynamics are age-structured.
        # If so, then we must extract relevant information on stock recruitment
        # to pass to the forecast function.
        
        ## If stock is age-structured
        if(ageStruct){
          
          # First estimate/extract the SR model and params.
          sr.pars  <- advice.ctrl[[st]]$sr$params # sr parameters if specified.
          sr.model <- advice.ctrl[[st]]$sr$model  # sr model, mandatory.
          
          # if params missing => estimate the parameters using the specified years.
          if(is.null(sr.pars)){
            
            if(is.null(advice.ctrl[[st]]$sr$years)) {
              
              # yr0 missing => use all data years, except the assessment year for 
              # which rec is unknown
              sr.yrs <- which(round(quantSums(stocks[[st]]@stock.n))!=0)[1]:(year-1)
              
            } else{
              
              y.rm <- as.numeric(advice.ctrl[[st]]$sr$years['y.rm'])
              nyrs <- as.numeric(advice.ctrl[[st]]$sr$years['num.years'])
              sr.yrs <- yrsnames[(year-y.rm-nyrs + 1):(year-y.rm)]
              
            }
            
            ## Extract recruitment and SSB data to fit an SR model
            rec <- stki@stock.n[1,sr.yrs]
            ssb <- ssb(stki)[,sr.yrs]
            
            # if rec.age != 0 adjust rec and ssb.
            rec.age <- as.numeric(dimnames(rec)[[1]])
            
            if(rec.age != 0){
              rec <- rec[, -(1:rec.age),]
              ssb <- ssb[, 1:(dim(ssb)[2] - rec.age),]
            }
            
            if(sr.model != 'geomean') {
              sr.pars <- try(params(fmle(FLSR(rec   = rec, 
                                              ssb   = ssb, 
                                              model = sr.model))), 
                             silent = TRUE)
            }
            
            if(class(sr.pars) == 'try-error' | sr.model == 'geomean'){
              sr.model <- 'geomean'
              sr.pars <- 10^6*c(prod(c(rec/10^6))^(1/length(c(rec))))
              sr.pars <- FLPar(a = ifelse(is.na(sr.pars), 0, sr.pars))
            }
            
            sr1 <- sr.pars
            
          } else{ # sr.pars not null
            if(i == 1){
              sr1 <- iter(sr.pars,i)
            }
            sr1[] <-  iter(sr.pars,i)[]
            
          }
          
          ## carry out forecast
          stki <- FLash::fwd(stki, 
                             ctrl = fwd.ctrl, 
                             sr = list(model =sr.model, 
                                       params = sr1))
          
        } else {
          
          # SECTION 1.2.3 - Forecast age-aggregated stocks
          #-----------------------------------------------#
          
          # Extract the years to calculate the mean historical growth of the stock
          if(is.null(advice.ctrl[[st]]$growth.years)) {
            growth.years <- max(1,(year - 11)):(year-1)
          } else{
            y.rm <- as.numeric(advice.ctrl[[st]]$growth.years['y.rm'])
            nyrs  <- as.numeric(advice.ctrl[[st]]$growth.years['num.years'])
            growth.years <- yrsnames[(year-y.rm-nyrs + 1):(year-y.rm)]
          }
          
          stki <- fwdBD(stki, fwd.ctrl, growth.years)
          
        }
        
        ## retain prepared stock and forward control objects if this is the
        ## assessed stock
        if (st == stknm) {
          
          
          
        }
        
        ## Insert forecast outputs into stock iteration slot
        stk0[,,,,,i] <- stki
        
      }# END iteration loop
      
      # SECTION 1.2.3 - Apply ICES HCR
      #---------------------------------------#
      
      ## if the i'th stock is category 1 and age structured, then extract SSB for
      ## the final data year. If not age structured extract overall biomass for 
      ## the final data year
      
      if(ageStruct)
        b.datyr <- ssb(stk0)[,year+1,drop = TRUE] # [it]
      else
        b.datyr <- (stk0@stock.n*stk0@stock.wt)[,year+1,drop = TRUE] # [it]
      
      ## this next line of code is somewhat more complicated. Generating a 1 x iter
      ## matrix, it identifies whether the biomass is below Blim [0], between
      ## Blim and Btrigger [1] or above Btrigger [2] for each iteration
      
      b.pos <- apply(matrix(1:iter,1,iter), 2, 
                     function(i) {findInterval(b.datyr[i], 
                                               ref.pts_st[c('Blim', 'Btrigger'),i])
                     })  # [it]
      
      ## Apply HCR to find Ftarget
      Ftg[st,] <- ifelse(b.pos == 0, 
                         0, 
                         ifelse(b.pos == 1, 
                                ref.pts_st['Fmsy',]*b.datyr/ref.pts_st[ 'Btrigger',], 
                                ref.pts_st['Fmsy',]))
      
      ## Extract the age range to calculate Fbar (Why??)
      minfbar <- stocks[[st]]@range['minfbar']
      maxfbar <- stocks[[st]]@range['maxfbar']
      
      ## Calculate Fstatusquo as an average of last 3 data years
      ## if this has not been calculated in a previous step
      if(sum(is.na(Fsq[st,])) > 0) {
        Fsq[st,] <- yearMeans(fbar(stocks[[st]])[,(year-3):(year-1)])
      }

      
      ## extract upper and lower limit reference points for F
      Fupp[st,] <- ref.pts_st['Fupp',]
      Flow[st,] <- ref.pts_st['Flow',]
      
    } # END Category 1 stocks
    
    #--------------------------------------------------------------#
    # SECTION 1.3 - Calculate TAC multiplier for Category 3 stocks #
    #--------------------------------------------------------------#
    #
    # If the stock is category 3, then stock indices and reference points are
    # used to calculate a TAC multiplier that is then translated into a Catch
    # multiplier
    
    if(stocksCat[st] == 3){
      
      ## This is the mean index in the final 2 data years divided by the 
      ## mean index in the preceding 3 data years
      Brat    <- c(mean(indices[[st]][[1]]@index[,(year-2):(year-1)])/
                     mean(indices[[st]][[1]]@index[,(year-3):(year-5)])) # [it]
      
      ## The Catch status quo is the mean of last 3 data years 
      Cst[st,]    <- yearMeans(stocks[[st]]@catch[,(year-3):(year-1)])[drop=T]   # [it]
      
      ## Extract TAC, alpha, beta reference points
      tac     <- advice[['TAC']][st, year-1,drop=T]
      alpha   <- advice.ctrl[[st]][["ref.pts"]]["alpha", ]
      beta    <- advice.ctrl[[st]][["ref.pts"]]["beta", ]
      betaUp  <- advice.ctrl[[st]][["ref.pts"]]["betaUp", ]
      
      ## Apply parameters to calculate TAC multiplier
      tacUpMult <- ifelse(Brat < 1-alpha, 1 - beta*0.5,
                          ifelse(Brat < 1+alpha, 1 + beta, 1 + betaUp))
      tacMult   <- ifelse(Brat < 1-alpha, 1-beta, ifelse(Brat < 1+alpha, 1, 1 + beta))
      
      # translate the TAC multiplier to C multiplier.
      CupMult <- tacUpMult*tac/Cst[st,]
      CMult   <- tacMult*tac/Cst[st,]
      
      # For these stocks we fill the Fsq, Fupp and Ftg in terms of catch, the 
      # linearity is assumed in effort catch becasue we don't have any other 
      # information.
      Ftg[st,]  <- CMult 
      Fupp[st,] <- CupMult 
      Fsq[st,]  <- 1
      
    } # END Category 3 stocks
  } # END StocksinHCR loop
  
  #----------------------------------------------------------------------------#
  # SECTION 3 - Apply the conditions of the HCR
  #----------------------------------------------------------------------------#
  #
  # In this section, the Ftargets calculated in the previous section are first
  # scaled to eliminate overshoot of Fupper and then scaled to reduce undershoot
  # of Flower.
  # Recall that Ftg and Fsq are matrices [st, iter]. Hence lambda will be a vector
  # of length niter, and Fadv will be a matrix [st, niter].
  
  ## if Fsq is zero, impute a very low value
  Fsq[Fsq == 0] <- 1e-6
  
  ## For each column (iteration) find the ratio of Ftg/Fsq that matches the 
  ## approach (default is maximum) 
  lambda0 <- apply(Ftg/Fsq,2,approach)
  
  ## For each column (iteration) scale F status quo by lambda 
  Fadv0 <- Fadv <- sweep(Fsq,2,lambda0, "*") # [nst,it]
  
  ## The next section addresses the case where any scaled F is now higher than
  ## the upper limit reference for F (only possible if approach != min). If so,
  ## then a correction multiplier is calculated to scale down Fadv to lie within
  ## Fupp in all cases.
  
  # The checks in the second step depend on the approach
  if((approach %in% c('max', 'mean')) & any(Fadv0 > Fupp)){
    lambda1 <- apply(Fupp/Fadv0,2,min) # The F multiplier (where Fadv most exceeds Fupp)
    Fadv0 <- Fadv <- sweep(Fadv0,2,lambda1,"*")
  }
  
  ## Similarly the next chuck of code addresses the case where the scaled advice
  ## is below the lower limit reference for F
  # If for some stock, under any option, we are below Flow, we increase F as 
  # much as possible within the ranges.
  
  if(any(Fadv0 < Flow)){#if((approach == 'min') & any(Fadv0 < Flow)){
    
    ## Find the case where advice is most below Flow for each iteration - vector
    lambda2 <- apply(Flow/Fadv0,2,max) #[it]
    
    ## calculate the ratio of Fupper/Fadvice (should all be >=1) - matrix[st, iter]
    upps <- Fupp/Fadv0
    
    ## The purpose of the code below is a little complex. We divide the matrix
    ## containing the ratios of Fupper/Fadvice (all >=1) [st, iter] by the vector
    ## of max(Flower/Fadvice) [iter] - the maximum undershoot of quota.
    ##
    ## This gives a matrix [st, iter] representing the effect of increasing F
    ## to offset quota undershoot. If the 1/value > 1 then we now have overshoot
    ## for this stock. We sum across columns to check the number of stocks affected.
    ##
    ## sum across columns
    comp1 <- colSums((1/sweep(upps,2,lambda2,"/")) >1) # [it]
    
    # The correction will be always the minimum of the Fupp/Fadv0 ratios, because 
    # for NONE of the stock we can be above.
    ## We find the minimum undershoot for each iteration (most likely 1)
    pos <- apply(upps,2,which.min) # [it]
    
    # For each iteration:
    #
    # If comp1 == 0, for this iteration none of the stocks is above Fupp by 
    #                offsetting quota undershoot. Offsetting is no problem
    # If comp1 == 1, there exist one stock above Fupp
    # If comp1 > 1, there are several stocks above Fupp, so from the corresponding
    #               multipliers we need to select the lowest one.
    #
    # If comp1 >= 1 then we need to update our scaling factor using the minimum
    # ratio of Fupper/Fadvice (typically 1) per iteration
    lambda2 <- ifelse(comp1 == 0, 
                      lambda2,  
                      sapply(1:length(comp1), # apply over each iteration
                             function(x) upps[pos[x],x] ))
    
    ## Apply scaling factor to minimise quota undershoot
    Fadv <- sweep(Fadv0,2,lambda2,"*")
  }
  
  ## Extract advised F for current stock
  Fadv_st <- Fadv[stknm,]
  
  ## extract intermediate year assumptions
  int.yr <- advice.ctrl[[stknm]]$intermediate.year
  
  ## If no advice exists, impute very small F
  Fadv_st[is.na(Fadv_st)] <- 1e-6
  
  print(Fadv_st)
  
  #----------------------------------------------------------------------------#
  # SECTION 4: For 'stknm' calculate the corresponding TAC advice as done in the 
  #            rest of the HCRs 
  #----------------------------------------------------------------------------#
  #
  # In this section the F advice is used as the F target for the TAC year.
  # This is a HUUUUGELY inelegant approach, but I need to recreate the stock
  # forecast object. If I could dispense with the need to loop over iterations,
  # I could retain the necessary objects for the focal stock.
  
  ## extract data for the current stock and impute a low fishing mortality if zero
  stk <- stocks[[stknm]]
  stk@harvest[stk@harvest < 0] <- 0.00001
  
  ## infer whether age structure is present
  ageStruct <- ifelse(dim(stk@m)[1] > 1, TRUE, FALSE)
  
  ## If current stock is category 1 (assessment includes population dynamics) then
  ## extend stock data in preparation for a three year forecast (current year, TAC year, TAC+1 year)
  if(stocksCat[stknm] == 1) stk <- stf(stk, 
                                       nyears = 3, 
                                       wts.nyears = 3, 
                                       fbar.nyears = 3, 
                                       f.rescale = f.rescale)
  
  ## Extrate the basis for advice (landings or catch)
  Cadv <- ifelse(advice.ctrl[[stknm]][['AdvCatch']][year+1] == TRUE, 'catch', 'landings')
  
  iter     <- dim(stk@m)[6]
  yrsnames <- dimnames(stk@m)[[2]]
  yrsnumbs <- as.numeric(yrsnames)
  
  assyrname <- yrsnames[year]
  assyrnumb <- yrsnumbs[year]
  
  ## Loop over each iteration
  for(i in 1:iter){
    
    ## If stocks are Category 1
    if(stocksCat[stknm] == 1){
      
      ## extract i'th iteration from stock
      stki <- iter(stk, i)
      
      ## Impute a low TAC if advice is not available or is zero
      if(is.na(Fadv_st[i]) | Fadv_st[i] == 0){
        advice[['TAC']][stknm,year+1,,,,i] <- 0.1
        next
      }
      
      int.yr <- ifelse(is.null(int.yr), 'Fsq', int.yr)
      
      # Calculate Fsq
      if(fbar.nyears == 1 | f.rescale) {
        Fsq <- mean(fbar(stki)[,(year-1)]) 
      } else {
        Fsq <- mean(fbar(stki)[,(year-fbar.nyears):(year-1)])
      }
      
      ## Build forward control object
      if(int.yr == 'Fsq')
        fwd.ctrl <- FLash::fwdControl(data.frame(year = c(0, 1),  
                                                 val = c(Fsq, Fadv_st[i]), 
                                                 quantity = c( 'f', 'f'), 
                                                 rel.year = c(NA,NA)))
      else
        fwd.ctrl <- FLash::fwdControl(data.frame(year = c(0, 1),  
                                                 val = c(advice$TAC[stknm,year, 
                                                                    drop=TRUE][i], 
                                                         Fadv_st[i]), 
                                                 quantity = c( 'catch', 'f')))
      
      # Refresh the years in fwd!!
      fwd.ctrl@target$year     <- fwd.ctrl@target$year + assyrnumb
      fwd.ctrl@target$rel.year <- fwd.ctrl@target$rel.year + assyrnumb
      

      
      # if in <year 0> quantity = catch => set TAC in <year 0> in val
      if(fwd.ctrl@target[fwd.ctrl@target$year == assyrnumb,'quantity'] == 'catch'){
        k <- which(fwd.ctrl@target$year == assyrnumb)
        fwd.ctrl@target[k,'val']     <- advice$TAC[stknm,year,,,,i]
        fwd.ctrl@trgtArray[k, 'val',] <- advice$TAC[stknm,year,,,,i]
      }
      
      if(dim(stki@m)[1] > 1){
        # First estimate/extract the SR model and params.
        sr.pars  <- advice.ctrl[[stknm]]$sr$params # sr parameters if specified.
        sr.model <- advice.ctrl[[stknm]]$sr$model  # sr model, mandatory.
        if(is.null(sr.pars)){                   # if params missing => estimate the parameters using the specified years.
          if(is.null(advice.ctrl[[stknm]]$sr$years)) sr.yrs <- which(round(quantSums(stocks[[stknm]]@stock.n))!=0)[1]:(year-1)# yr0 missing => use all data years, except the assessment year for which rec is unknown
          else{
            y.rm <- as.numeric(advice.ctrl[[stknm]]$sr$years['y.rm'])
            nyrs <- as.numeric(advice.ctrl[[stknm]]$sr$years['num.years'])
            sr.yrs <- yrsnames[(year-y.rm-nyrs + 1):(year-y.rm)]
          }
          rec <- stki@stock.n[1,sr.yrs]
          ssb <- ssb(stki)[,sr.yrs]
          
          # if rec.age != 0 adjust rec and ssb.
          rec.age <- as.numeric(dimnames(rec)[[1]])
          if(rec.age != 0){
            rec <- rec[, -(1:rec.age),]
            ssb <- ssb[, 1:(dim(ssb)[2] - rec.age),]
          }
          
          if(sr.model != 'geomean') sr.pars <- try(params(fmle(FLSR(rec = rec, 
                                                                    ssb = ssb,
                                                                    model = sr.model))), 
                                                   silent = TRUE) 
          
          if(class(sr.pars) == 'try-error' | sr.model == 'geomean'){
            sr.model <- 'geomean'
            sr.pars <- c(prod(c(rec))^(1/length(c(rec))))
            sr.pars <- FLPar(a = ifelse(is.na(sr.pars), 0, sr.pars))
          }
          
          sr1 <- sr.pars
        }
        else{ # sr.pars not null
          if(i == 1){
            sr1 <- iter(sr.pars,i)
          }
          sr1[] <-  iter(sr.pars,i)[]
          
        }
        
        stki <- FLash::fwd(stki, ctrl = fwd.ctrl, sr = list(model =sr.model, params = sr1))
      }
      else{
        
        # Extract the years to calculate the mean historical growth of the stock
        if(is.null(advice.ctrl[[stknm]]$growth.years))   growth.years <- max(1,(year - 11)):(year-1)
        else{
          y.rm <- as.numeric(advice.ctrl[[stknm]]$growth.years['y.rm'])
          nyrs  <- as.numeric(advice.ctrl[[stknm]]$growth.years['num.years'])
          growth.years <- yrsnames[(year-y.rm-nyrs + 1):(year-y.rm)]
        }
        
        stki <- fwdBD(stki, fwd.ctrl, growth.years)
      }
      
      yy <- ifelse(slot(stki, Cadv)[,year+1] == 0, 1e-6, slot(stki, Cadv)[,year+1])
      yy <- ifelse(is.na(yy), 0.1, yy)
      
      advice[['TAC']][stknm,year+1,,,,i] <- yy # The TAC is given in terms of CATCH.
      
    }
    if(stocksCat[stknm] == 3){
      advice[['TAC']][stknm,year+1,,,,i] <- Cst[st]*Fadv_st
    }
    
    
  } # END iteration loop
  
  return(advice)
  
} # END function
