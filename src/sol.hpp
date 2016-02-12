/***********************************************************************************
 * sol.hpp
 * 
 * Interface to sol.cpp
 * 
 * 27jan2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#ifndef SOL_HPP
#define SOL_HPP

#include <RcppArmadillo.h>
#include <math.h>
#include <string.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;
//using namespace std::string ;

arma::mat euler_hat( arma::mat coeffs, arma::mat X, std::string model, 
                      int lags, List params, int n_exog, int n_endog,
                      arma::rowvec rho, arma::rowvec sig_eps, int n_integ,
                      int N, arma::rowvec upper, arma::rowvec lower, bool cheby,
                      arma::mat exog_innov_mc, bool quad, int n_nodes ) ;

#endif