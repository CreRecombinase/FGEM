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
double FGEM_Logit_log_lik_cpp(const arma::vec &Beta, const arma::mat &x, const arma::vec B){

  arma::vec pvec = 1/(1+arma::exp(-(x*Beta)));
  // Rcpp::Rcout<<"pvec:"<<pvec<<std::endl;

//  arma::vec uvec = (pvec%B)/((pvec%B)+(1-pvec));
  // Rcpp::Rcout<<"uvec:"<<uvec<<std::endl;
  return(-arma::sum(log(pvec%B+(1-pvec))));
}


