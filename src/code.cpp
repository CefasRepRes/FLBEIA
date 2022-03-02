#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List condition_fleet_effort(List fl, IntegerVector dim, IntegerVector sim_yrs, IntegerVector mean_yrs)
{
  
  // Loop over fleets
  for(int f=0; f<fl.length();f++){
  
  S4 x = fl[f]; // a fleet
  
  S4 eff = x.slot("effort");              // effort slot
  NumericVector dat = eff.slot(".Data");  
  S4 cap = x.slot("capacity");              // capacity slot
  NumericVector cap_dat = cap.slot(".Data");  

// Loop through updating the variable //

  int na = dim[0];  // number ages

// loop over the years
  for(int y=sim_yrs[0]; y<sim_yrs[sim_yrs.size()-1]+1; y++) {
    
    
  for(int a=0;a<na; a++) { // loop over ages

    // Calculate the mean for the year/age combination
    NumericVector hist_yrs(mean_yrs.size());
    NumericVector cap_hist_yrs(mean_yrs.size());

    for(int i=mean_yrs[0];i<mean_yrs[mean_yrs.size()-1]+1; i++) { // loop over the years to fill the matrix
      hist_yrs[i-mean_yrs[0]] = dat[a*mean_yrs.size() + i];
      if(R_IsNA(hist_yrs[i-mean_yrs[0]])) hist_yrs[i-mean_yrs[0]] = 0;
      cap_hist_yrs[i-mean_yrs[0]] = cap_dat[a*mean_yrs.size() + i];
      if(R_IsNA(cap_hist_yrs[i-mean_yrs[0]])) cap_hist_yrs[i-mean_yrs[0]] = 0;
                  }
   
           double meanval = mean(na_omit((hist_yrs))); // calculate the mean
           double cap_meanval = mean(na_omit((cap_hist_yrs)));
           if(R_IsNA(meanval)) meanval = 0; // zeros if no value
           if(R_IsNA(cap_meanval)) cap_meanval = 0; 
           
    dat[a*sim_yrs.size() + y] = meanval;
    cap_dat[a*sim_yrs.size() + y] = cap_meanval;
                           }
                    }
   
// Return the dimensions to the object    
dat.attr("dim") = dim;
cap_dat.attr("dim") = dim;
  
// Return to FLFleet
eff.slot(".Data") = dat;
x.slot("effort") = eff;
cap.slot(".Data") = cap_dat;
x.slot("capacity")= cap;


// Metier loop

List met = x.slot("metiers");             //  metiers

for(int m=0; m<met.length();m++){
 
 S4 mt = met[m];
  
  S4 effsh = mt.slot("effshare");
  NumericVector effsh_dat = effsh.slot(".Data");  
  
  // Loop through updating the variable //

  // loop over the years
  for(int y=sim_yrs[0]; y<sim_yrs[sim_yrs.size()-1]+1; y++) {
    
    
    for(int a=0;a<na; a++) { // loop over ages
      
      // Calculate the mean for the year/age combination
      NumericVector effsh_hist_yrs(mean_yrs.size());
      
      for(int i=mean_yrs[0];i<mean_yrs[mean_yrs.size()-1]+1; i++) { // loop over the years to fill the matrix
        effsh_hist_yrs[i-mean_yrs[0]] = effsh_dat[a*mean_yrs.size() + i];
        if(R_IsNA(effsh_hist_yrs[i-mean_yrs[0]])) effsh_hist_yrs[i-mean_yrs[0]] = 0;// If it#s NA, it's really a zero - important for calculating the mean correctly
      }
      
      double effsh_meanval = mean(na_omit(effsh_hist_yrs)); // calculate the mean
      if(R_IsNA(effsh_meanval)) effsh_meanval = 0;   // zeros if no value
      
      
      effsh_dat[a*sim_yrs.size() + y] = effsh_meanval;
    }
  }
  
  // Return the dimensions to the object    
  effsh_dat.attr("dim") = dim;
  
  // Return to FLFleet
  effsh.slot(".Data") = effsh_dat;
  mt.slot("effshare") = effsh;
  
  met[m] = mt;
  
}

x.slot("metiers") = met;

fl[f] = x; // return fleet to the list

  }
  return(fl);
}

// [[Rcpp::export]]
List condition_flcatches(List fl,
                         NumericVector SLwt,
                         NumericVector SDwt,
                            NumericVector B, 
                            CharacterVector st, 
                            IntegerVector mean_yrs_q,
			    IntegerVector mean_yrs_wts,
			    IntegerVector mean_yrs_sel,
                            IntegerVector sim_yrs,
			    bool LO,
			    bool UseCWt4Discards)
{
  // Loop over fleets
  for(int f=0; f<fl.length();f++){
    
    S4 x = fl[f]; // a fleet
    
   NumericVector Ef = x.slot("effort");       // fleet effort
   
   List met = x.slot("metiers");             //  metiers
   
    // Loop over metier
   for(int m=0; m<met.length();m++){
     
     S4 mt = met[m];
     
     NumericVector mteffsh = mt.slot("effshare");
     
     S4 catches = mt.slot("catches");
     CharacterVector catchNames = catches.slot("names");
     List CL = catches.slot(".Data");   // the FLCatch Data as a list
     
     if(is_true(any(in(st, catchNames)))) {
       
      IntegerVector pos = match(st, catchNames) -1;     // this is an R function, so references from 1     
      int pos2 = pos[0];      // trick to get back to int from IntegerVector

         S4 C = CL[pos2];                               // The specific stock data
         
         S4 flq = C.slot("landings.n");                // just used to get the right dimensions
         NumericVector dim = flq.attr("dim");          // save the dimensions attribute
         
         int na = dim[0];                             // number ages for stock
         
	 NumericVector Efm = Ef * mteffsh;
         NumericVector E = rep_each(Efm, na);     // effort vector, repeated across all ages
    
         S4 L = C.slot("landings.n");     // landings.n
         NumericVector L_dat = L.slot(".Data");
         
         S4 D = C.slot("discards.n");      // discards.n
         NumericVector D_dat = D.slot(".Data");
         
         S4 alpha = C.slot("alpha");                  // alpha
         NumericVector alpha_dat = alpha.slot(".Data");
         
         S4 beta = C.slot("beta");        // beta
         NumericVector beta_dat = beta.slot(".Data");
         
         S4 Lwt = C.slot("landings.wt");  // landings.wt
         NumericVector Lwt_dat = Lwt.slot(".Data");
         
         S4 Dwt = C.slot("discards.wt");  // discards.wt
         NumericVector Dwt_dat = Dwt.slot(".Data");
         
         S4 pr = C.slot("price");         // price
         NumericVector pr_dat = pr.slot(".Data");
        
        
        //  L and D set NAs to zero
        for(int i=0; i<L_dat.size(); i++) {
          if(R_IsNA(L_dat[i])) L_dat[i] = 0;
          if(R_IsNA(D_dat[i])) D_dat[i] = 0;
        }
        
        NumericVector Cf = L_dat + D_dat;       // catch fleet
                             
        NumericVector q(Cf.size());     // a container for catchability
        
        for(int i=0; i<Cf.size(); i++) {
        q[i] = Cf[i]/(pow(E[i],alpha_dat[i])*pow(B[i],beta_dat[i]));    // catchability
                        }
        
        NumericVector Lsel = L_dat/(L_dat+D_dat);                       // landings selection in numbers

        
        // Now to condition the simulation years: 
        // catch.q, landings.sel, discards.sel, landings.wt, discards.wt
        // price, alpha, beta
      
         for(int a=0;a<na; a++) { // loop over ages
                
              // Calculate the mean for the year/age combination
              NumericVector hist_yrs_q(mean_yrs_q.size());
              NumericVector hist_yrs_s(mean_yrs_sel.size());
              NumericVector hist_yrs_lw(mean_yrs_wts.size());
              NumericVector hist_yrs_dw(mean_yrs_wts.size());
              NumericVector hist_yrs_Slw(mean_yrs_wts.size());
              NumericVector hist_yrs_Sdw(mean_yrs_wts.size()); // stock discard weights
              NumericVector hist_yrs_p(mean_yrs_wts.size());
              NumericVector hist_yrs_a(mean_yrs_q.size());
              NumericVector hist_yrs_b(mean_yrs_q.size());


	      // catchability related
              for(int i=mean_yrs_q[0];i<=mean_yrs_q[mean_yrs_q.size()-1]; i++) { // loop over the years to fill the vector
                hist_yrs_q[i-mean_yrs_q[0]] = q[na*i + a];
                if(!R_finite(hist_yrs_q[i-mean_yrs_q[0]])) hist_yrs_q[i-mean_yrs_q[0]] = 0; // remove infintes, but keep good values
                hist_yrs_a[i-mean_yrs_q[0]] = alpha_dat[na*i + a];
                hist_yrs_b[i-mean_yrs_q[0]] = beta_dat[na*i + a];
              }


	      //weight related
 for(int i=mean_yrs_wts[0];i<=mean_yrs_wts[mean_yrs_wts.size()-1]; i++) { // loop over the years to fill the vector
                hist_yrs_lw[i-mean_yrs_wts[0]] = Lwt_dat[na*i + a];
                hist_yrs_dw[i-mean_yrs_wts[0]] = Dwt_dat[na*i + a];
                hist_yrs_Slw[i-mean_yrs_wts[0]] = SLwt[na*i + a];
                hist_yrs_Sdw[i-mean_yrs_wts[0]] = SDwt[na*i + a];
                hist_yrs_p[i-mean_yrs_wts[0]] = pr_dat[na*i + a];
              }

 		// selectivity related
 for(int i=mean_yrs_sel[0];i<=mean_yrs_sel[mean_yrs_sel.size()-1]; i++) { // loop over the years to fill the vector
                hist_yrs_s[i-mean_yrs_sel[0]] = Lsel[na*i + a];
              }
             
		 // compute means
              double meanval_q = mean(na_omit((hist_yrs_q))); // calculate the mean
              double meanval_s = mean(na_omit((hist_yrs_s))); 
              double meanval_lw = mean(na_omit((hist_yrs_lw))); 
              double meanval_dw = mean(na_omit((hist_yrs_dw)));
              double meanval_Slw = mean(na_omit((hist_yrs_Slw))); 
              double meanval_Sdw = mean(na_omit((hist_yrs_Sdw)));
              double meanval_p = mean(na_omit((hist_yrs_p))); 
              double meanval_a = mean(na_omit((hist_yrs_a))); 
              double meanval_b = mean(na_omit((hist_yrs_b))); 
              
                          
              if(R_IsNA(meanval_q)) meanval_q = 0;    // zeros if no value
              if(R_IsNaN(meanval_q)) meanval_q = 0;   // zeros if NaN
              if(!R_finite(meanval_q)) meanval_q = 0; // If we have no effort but catches, get infinte q.
              if(R_IsNA(meanval_s)) meanval_s = 0; 
              if(R_IsNaN(meanval_s)) meanval_s = 0; 
              if(R_IsNA(meanval_lw)) meanval_lw = 0; 
              if(R_IsNaN(meanval_lw)) meanval_lw = 0; 
              if(meanval_lw==0) meanval_lw = meanval_Slw; // need a value where missing, use stock value
              if(R_IsNA(meanval_dw)) meanval_dw = 0; 
              if(R_IsNaN(meanval_dw)) meanval_dw = 0;     
              if(meanval_dw==0) meanval_dw = meanval_Sdw;  // need a value where missing, use stock value
	      if(meanval_dw==0) meanval_dw = meanval_Slw;  // If it's still zero, use the landings weight
	      if(R_IsNA(meanval_dw)) {meanval_dw = meanval_Slw;} // If its NA, use the landings weight
              if(R_IsNA(meanval_p)) meanval_p = 0; 
              if(R_IsNaN(meanval_p)) meanval_p = 0; 
              if(R_IsNA(meanval_a)) meanval_a = 1; 
              if(R_IsNaN(meanval_a)) meanval_a = 1; 
              if(R_IsNA(meanval_b)) meanval_b = 1; 
              if(R_IsNaN(meanval_b)) meanval_b = 1;
 
	      // the discarded weights-at-age becomes the catch weights-at-age
	      	      // In this case discards are of the full population,
	      // but landings retain the landings weight
	      if(UseCWt4Discards) {meanval_dw = (meanval_s * meanval_lw) + ((1-meanval_s) * meanval_dw);
	      }

	      // and the landings selection becomes 1
	      if(LO) {meanval_s = 1;}

              for(int y=sim_yrs[0]; y<=sim_yrs[sim_yrs.size()-1]; y++) { // loop over the years
              q[(na*y) + a] = meanval_q;    // sim year values
              Lsel[(na*y) + a] = meanval_s; 
              Lwt_dat[(na*y) + a] = meanval_lw; 
              Dwt_dat[(na*y) + a] = meanval_dw; 
              pr_dat[(na*y) + a] = meanval_p; 
              alpha_dat[(na*y) + a] = meanval_a; 
              beta_dat[(na*y) + a] = meanval_b; 
                }
          }
          
          NumericVector Dsel = 1.0 - Lsel;                    // Discards selection
          
          // Return the dimensions to the object    
          q.attr("dim") = dim;
          Lsel.attr("dim") = dim;
          Dsel.attr("dim") = dim;
          Lwt_dat.attr("dim") = dim;
          Dwt_dat.attr("dim") = dim;
          pr_dat.attr("dim") = dim;
          alpha_dat.attr("dim") = dim;
          beta_dat.attr("dim") = dim;
          
          // Return to stock //
          
          // Catchability
          S4 Q_slot = C.slot("catch.q");
          Q_slot.slot(".Data") = q;
          C.slot("catch.q") = Q_slot;
          
          S4 Lsel_slot = C.slot("landings.sel");
          Lsel_slot.slot(".Data") = Lsel;
          C.slot("landings.sel") = Lsel_slot;
          
          S4 Dsel_slot = C.slot("discards.sel");
          Dsel_slot.slot(".Data") = Dsel;
          C.slot("discards.sel") = Dsel_slot;
          
          Lwt.slot(".Data") = Lwt_dat;
          C.slot("landings.wt") = Lwt;
          
          Dwt.slot(".Data") = Dwt_dat;
          C.slot("discards.wt") = Dwt;
          
          pr.slot(".Data") = pr_dat;
          C.slot("price") = pr;
          
          alpha.slot(".Data") = alpha_dat;
          C.slot("alpha") = alpha;
          
          beta.slot(".Data") = beta_dat;
          C.slot("beta") = beta;
          
          CL[pos2] = C; // Return to CL
          catches.slot(".Data") = CL; // Return CL to catches
          
       }
        mt.slot("catches") = catches; // return FLCatches to the metier
        met[m] = mt;                  // return metier to the list of metier
   }
        x.slot("metiers") = met; // return metier to fleet
        fl[f] = x;   // return fleet to the fleet list
        
  }
  
  return(fl);
  
  
}