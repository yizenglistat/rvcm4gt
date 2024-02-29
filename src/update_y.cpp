#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]

// See McMahan's MHS_GTR_IP file for supplemental information
// N = number of individuals
// np = number of pools current individual is in
// p = current vector of probabilities, i.e. inverse link function
// Y = matrix with pool information for individuals
// G = matrix with observations and other information

// NumericVector update_y(int N,NumericVector p, NumericMatrix Y,NumericMatrix G,NumericVector W,NumericVector U,NumericVector se,NumericVector sp) {
//   int np, tmp, Gk, ck, as, ysum, id;
//   float sek, spk, pi1, pi2, pistar;

//   for(int i=0; i<N; i++){
//     pi1 = p(i);
//     pi2 = 1 - p(i);
//     np = Y(i,1);  // number of pools this individual was assigned to
//     for(int j=0; j<np; j++){
//       tmp = Y(i,(j+2)) - 1; // obtain a pool this individual was assigned to
//       Gk = G(tmp,0);    // observation for that pool
//       ck = G(tmp,1);    // number in that pool (pool size)
//       as = G(tmp,2) - 1;    // assay used
//       sek = se(as);   // sensitivity for that pool
//       spk = sp(as);   // specificity for that pool
//       ysum = 0;
//       Y(i,0) = 0;     // coding trick to 'remove' individual i, look at set Pij
//       for(int k=0; k<ck; k++){
//         id = G(tmp,(k+3))-1;    // indices of individuals in pool
//         ysum = ysum + Y(id,0);  // to decide if an individual in pool is positive
//       }
//       pi1 = pi1*(sek*Gk + (1-sek)*(1-Gk));
//       if(ysum > 0)
//         pi2 = pi2*(sek*Gk + (1-sek)*(1-Gk));
//       else
//         pi2 = pi2*((1-spk)*Gk + spk*(1-Gk));
//     }
//     pistar = (pi1 / (pi1 + pi2));
//     p(i) = pistar;
//     if(U(i) <= pistar)
//       Y(i,0) = 1;
//     else
//       Y(i,0) = 0;
//     W(i) = Y(i,0);
//   }
//   return W;
// }

NumericVector update_y(int N, NumericVector p, NumericMatrix Y, NumericMatrix G, NumericVector W, NumericVector U,NumericVector se,NumericVector sp) {
  int np, tmp, Gk, ck, as, ysum, id;
  float sek, spk, pi1, pi2, pistar;

  for(int i=0; i<N; i++){
    pi1 = p(i);
    pi2 = 1 - p(i);
    np = Y(i,1);  // number of pools this individual was assigned to
    for(int j=0; j<np; j++){
      tmp = Y(i,(j+2)) - 1; // obtain a pool this individual was assigned to
      Gk = G(tmp,0);    // observation for that pool
      ck = G(tmp,1);    // number in that pool (pool size)
      as = G(tmp,2) - 1;    // assay used
      
      // sek = se(as);   // sensitivity for that pool
      // spk = sp(as);   // specificity for that pool
      
      if(se.size()==1){
        sek = se(0);
      }else{
        sek = se(as);   // sensitivity for that pool
      }

      if(sp.size()==1){
        spk = sp(0);
      }else{
        spk = sp(as);   // specificity for that pool
      }

      ysum = 0;
      Y(i,0) = 0;     // coding trick to 'remove' individual i, look at set Pij
      for(int k=0; k<ck; k++){
        id = G(tmp,(k+3))-1;    // indices of individuals in pool
        ysum = ysum + Y(id,0);  // to decide if an individual in pool is positive
      }
      pi1 = pi1*(sek*Gk + (1-sek)*(1-Gk));
      if(ysum > 0)
        pi2 = pi2*(sek*Gk + (1-sek)*(1-Gk));
      else
        pi2 = pi2*((1-spk)*Gk + spk*(1-Gk));
    }
    pistar = (pi1 / (pi1 + pi2));
    p(i) = pistar;
    if(U(i) <= pistar)
      Y(i,0) = 1;
    else
      Y(i,0) = 0;
    W(i) = Y(i,0);
  }
  return W;
}