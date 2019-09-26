#include <RcppDist.h>
#include "extra_cpps.h"

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;
using namespace arma;

// 1. my_invgamma : inverse gamma generator as of conventional Wikipedia notation
// [[Rcpp::export]]
double my_invgamma(double alpha, double beta){
  return(1.0/R::rgamma(alpha,1.0/beta));
}

// 2. my_dinvgamma : inverse gamma evaluator
double my_dinvgamma(double x, double alpha, double beta){
  return(1.0/R::dgamma(x, alpha, 1.0/beta, 0));
}