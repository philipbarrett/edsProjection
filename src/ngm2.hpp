/***********************************************************************************
 * ngm2.hpp
 * 
 * Interface to ngm2.cpp
 * 
 * 10feb016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#ifndef NGM2_HPP
#define NGM2_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::rowvec integrand_ngm_2( arma::mat exog, arma::mat endog, arma::rowvec exog_lead, 
                       List params, arma::mat coeffs, int n_exog, int n_endog, 
                       int N, arma::rowvec upper, arma::rowvec lower, 
                       bool cheby ) ;
arma::rowvec euler_hat_ngm_2( arma::mat exog, arma::mat endog, arma::mat exog_innov_integ, 
                  List params, arma::mat coeffs, int n_exog, int n_endog,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights, bool print_rhs ) ;
#endif