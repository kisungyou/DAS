#ifndef _DAS_extra_cpps_H
#define _DAS_extra_cpps_H

#define ARMA_NO_DEBUG

#include <RcppDist.h>

using namespace arma;
using namespace Rcpp;

/*
 * LIST OF FUNCTIONS
 * 
 */
double my_invgamma(double alpha, double beta);
double my_dinvgamma(double x, double alpha, double beta);

#endif