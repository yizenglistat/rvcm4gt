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
NumericMatrix update_sesp(int N,int npools,NumericMatrix Y,NumericMatrix Z,NumericMatrix PS,int na) {
  int Ysum, tr, cj, ida, id;
  //NumericMatrix PS(ida,4);
  for(int j=0;j<npools; j++){
    Ysum=0;
    tr=Z(j,0);
    cj=Z(j,1);
    ida=Z(j,2);
    for(int i=0;i<cj;i++){
      id=Z(j,i+3);
      Ysum=Ysum+Y(id-1,0);
    }
    if(Ysum>0){
      if(tr>0){
        PS(ida-1,0)=PS(ida-1,0)+1;
      }else{
        PS(ida-1,1)=PS(ida-1,1)+1;
      }
    }else{
      if(tr>0){
        PS(ida-1,3)=PS(ida-1,3)+1;
      }else{
        PS(ida-1,2)=PS(ida-1,2)+1;
      }
    }
  }
  return PS;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


