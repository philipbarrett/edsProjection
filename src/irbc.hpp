/***********************************************************************************
 * irbc.hpp
 * 
 * Interface to irbc.cpp
 * 
 * 14feb2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#ifndef IRBC_HPP
#define IRBC_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::rowvec integrand_irbc( 
                arma::rowvec endog, arma::rowvec exog_lead, 
                double gamma, arma::mat coeffs_cont, 
                int n_exog, int n_endog, int n_cont, int N, 
                arma::rowvec upper, arma::rowvec lower, bool cheby ) ;

arma::rowvec euler_hat_irbc( 
                  arma::rowvec exog, arma::rowvec endog, arma::rowvec cont,
                  arma::mat exog_innov_integ, double betta,
                  double gamma, arma::mat coeffs_cont, 
                  int n_exog, int n_endog, int n_cont, int n_fwd,
                  arma::mat rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights,
                  bool print_rhs ) ;

//arma::rowvec x_eqns_irbc( arma::mat exog, arma::rowvec cont, List params ) ;

arma::rowvec contemp_eqns_irbc( 
    arma::mat exog, arma::mat endog, arma::rowvec cont, List params, List extra_args ) ;

arma::rowvec irbc_reg( 
                  arma::mat exog, arma::mat endog, arma::rowvec cont,
                  arma::mat exog_innov_integ, 
                  List params, arma::mat coeffs, arma::mat coeffs_cont, 
                  int n_exog, int n_endog, int n_cont, int n_fwd,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights, 
                  List extra_args, bool print_rhs ) ;

#endif