#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double emds_gamma0(arma::mat dmat){
  // parameters
  int N = dmat.n_rows;
  double gamma0 = 0.0;
  double theval = 0.0;
  
  // iterate.. triplet
  for (int i=0;i<N;i++){
    for (int j=0;j<N;j++){
      for (int k=0;k<N;k++){
        theval = std::abs(dmat(i,j)+dmat(i,k)-dmat(j,k));
        if (theval > gamma0){
          gamma0 = theval;
        }
      }
    }
  }
  
  // report
  return(gamma0);
}