/***********************************************************************************
 * ngm.hpp
 * 
 * Interface to ngm.cpp
 * 
 * 26jan2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#ifndef NGM_HPP
#define NGM_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

double integrand_ngm( arma::mat exog, arma::mat endog, arma::rowvec exog_lead, 
                       List params, arma::mat coeffs, int n_exog, int n_endog, 
                       int N, arma::rowvec upper, arma::rowvec lower, 
                       bool cheby ) ;

arma::rowvec euler_hat_ngm( 
                  arma::mat exog, arma::mat endog, arma::mat exog_innov_integ, 
                  List params, arma::mat coeffs, int n_exog, int n_endog,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights, 
                  bool print_rhs ) ;

#endif