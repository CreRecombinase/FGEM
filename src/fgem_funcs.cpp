#include "RcppArmadillo.h"
//[[Rcpp::depends(RcppArmadillo)]]
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
NumericVector FGEM_Logit_log_lik_cpp(const arma::vec &Beta, const arma::mat &x, const arma::vec B){
  arma::vec pvec = 1/(1+arma::exp(-(x*Beta)));
  arma::vec uvec = (pvec%B)/((pvec%B)+(1-pvec));
  return(arma::sum(uvec*(arma::log(pvec)+arma::log(B))+(1-uvec)*log(1-pvec)));
}
// FGEM_Logit_fit <- function(Beta,x,B){
//   pvec  <- 1/(1+exp(-(x%*%Beta)))
//   uvec <- (pvec*B)/((pvec*B)+(1-pvec))
//   return(-sum(uvec*(log(pvec)+log(B))+(1-uvec)*log(1-pvec)))
//   return(Beta)
// }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

