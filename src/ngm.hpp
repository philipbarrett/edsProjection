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

arma::mat integrand_ngm_D( arma::mat exog, arma::mat endog, 
                arma::rowvec exog_lead, List params, arma::mat coeffs, 
                int n_exog, int n_endog, int N, arma::rowvec upper, 
                arma::rowvec lower, bool cheby ) ;

double err_ngm( arma::mat exog, arma::mat endog, arma::mat exog_innov_integ, 
                  List params, arma::mat coeffs, int n_exog, int n_endog,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::vec weights, 
                  bool print_rhs ) ;

arma::mat err_ngm_D( arma::mat exog, arma::mat endog, arma::mat exog_innov_integ, 
                     List params, arma::mat coeffs, int n_exog, int n_endog,
                     arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                     arma::rowvec lower, bool cheby, arma::vec weights ) ;

#endif