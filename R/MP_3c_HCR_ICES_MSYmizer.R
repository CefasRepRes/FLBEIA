#-------------------------------------------------------------------------------
# HCR proposed by ICES under MSY framework
# Biomass Based HCR.
#  Reference Points: Btrigger, Blim and Fmsy.
#   No formal porposal for any of this, usually :
#           - Btrigger = Bpa
#           - Blim  = ??? Blim?YYY?
#           - Fmsy = F0.1, Fmax....
#
#  - TAC advice depending on F in relation to BRP is:
#           - TAC[Fmsy]             B >= Btrigger.
#           - TAC[Fmsy*B/Btrigger]  B <  Btrigger.
#           - 0.                    B <  Blim (our proposal)
#
# 07/09/2011 12:20:24
# Rewritten to work with mizer
#-------------------------------------------------------------------------------
#' @rdname annualTAC
#' @aliases IcesHCR
IcesHCRmizer <- function(stocks, advice, advice.ctrl, year, stknm,...){

   
  Cadv <- ifelse(advice.ctrl[[stknm]][['AdvCatch']][year+1] == TRUE, 'catch', 'landings')
  
    stk <- stocks[[stknm]]
    iter     <- dim(stk@m)[6]
    yrsnames <- dimnames(stk@m)[[2]]
    yrsnumbs <- as.numeric(yrsnames)
    
    assyrname <- yrsnames[year]
    assyrnumb <- yrsnumbs[year]
    yr_mizer <- length(covars.ctrl$mizer$first.year:as.numeric(dimnames(biols[[1]]@n)$year[1]))-1 + year # correct ref year for mizer
    
   ## Run per iteration
    for(i in 1:iter){
    

      ref.pts <- lapply(names(biols), function(x) advice.ctrl[[x]]$ref.pts[,i]) # matrix[3,it]  rows = Blim, Btrigger, Fmsy
      names(ref.pts) <- names(biols)
        
      # For Ftg we first use Fmsy and then rerun fwd using the updated Ftg depending on the value of SSB
      Ftg   <- sapply(ref.pts, function(x) x[['Fmsy']])
      Brefs <- sapply(ref.pts, function(x) c(x[[c("Blim")]],x[[c("Btrigger")]]))
        

            # SSB in the advice year.
        
           ## run mizer forward at sq F
           mizer <- covars[["mizer"]]
           
           mizer_i <- mizer[[i]]
           
           ## Run mizer at status quo F
           F_mat <- if(length(dim(mizer_i$age_stuff$Fs))==3) {
                           mizer_i$age_stuff$Fs[nrow(mizer_i$age_stuff$Fs[,1,]),,] 
             } else {
              mizer_i$age_stuff$Fs 
               }
           
           ## Running mizer
           mizer_i <- FLBEIA:::progress_one_year(year=yr_mizer,effort=F_mat,prev_run=mizer_i) ## Need to fix year, is different reference, i.e. 37 is 2020
           
           ## Find the SSB from all the stocks as we need the F target for all the stocks
           b.datyr <- apply(mizer_i$age_stuff$cat_num_at_age * mizer_i$age_stuff$mean_waa * mizer_i$age_stuff$prop_mat,1,sum,na.rm = TRUE)/1e6 ## in tonnes
            
            # Find where the SSB (Age structured) OR Biomass (Aggregated) in relation to reference points.
            b.pos <- sapply(names(b.datyr), function(x) findInterval(b.datyr[[x]], Brefs[,x]))
            
            ## Ftarget across stocks
            Ftg <- sapply(names(Ftg), function(st) {
              ifelse(b.pos[[st]] == 0, 0, ifelse(b.pos[[st]] == 1, ref.pts[[st]][['Fmsy']]*b.datyr[[st]]/ref.pts[[st]][['Btrigger']], ref.pts[[st]][['Fmsy']]))
            })
            
            print(Ftg)
     
            ## Now fish in the TAC year at the target Fs
          
            F_mat <- mizer_i$age_stuff$Fs
            
            ## Scale F_mat to Ftg
            for(st in names(Ftg)) {
              F_mat[st,] <- (F_mat[st,]/max(F_mat[st,], na.rm = TRUE)) * Ftg[[st]]
              }
            
            mizer_ii <- FLBEIA:::progress_one_year(year=yr_mizer+1,effort=F_mat,prev_run=mizer_i)
            
            yy <- sum(mizer_ii$age_stuff$cat_bio_at_age[stknm,])/1e6 ## in tonnes
     
        advice[['TAC']][stknm,year+1,,,,i] <- yy # The TAC is given in terms of CATCH.

    }
    return(advice)
}

