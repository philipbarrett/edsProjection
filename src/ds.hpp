/***********************************************************************************
 * ds.hpp
 * 
 * Interface to ds.cpp
 * 
 * 14aug2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#ifndef DS_HPP
#define DS_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::rowvec integrand_ds( 
                arma::rowvec endog, arma::rowvec exog_lead, 
                double gamma, arma::mat coeffs_cont, 
                int n_exog, int n_endog, int n_cont, int N, 
                arma::rowvec upper, arma::rowvec lower, bool cheby ) ;

arma::rowvec euler_hat_ds( 
                  arma::rowvec exog, arma::rowvec endog, arma::rowvec cont,
                  arma::mat exog_innov_integ, double betta,
                  double gamma, arma::mat coeffs_cont, 
                  int n_exog, int n_endog, int n_cont,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights,
                  bool print_rhs ) ;

arma::rowvec contemp_eqns_ds( 
                arma::mat exog, arma::mat endog, arma::rowvec cont, List params ) ;

// arma::rowvec irbc_reg( 
//                   arma::mat exog, arma::mat endog, arma::rowvec cont,
//                   arma::mat exog_innov_integ, 
//                   List params, arma::mat coeffs, arma::mat coeffs_cont, 
//                   int n_exog, int n_endog, int n_cont,
//                   arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
//                   arma::rowvec lower, bool cheby, arma::rowvec weights, 
//                   bool print_rhs ) ;

#endif