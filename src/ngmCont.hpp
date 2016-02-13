/***********************************************************************************
 * ngmCont.hpp
 * 
 * Interface to ngm_cont.cpp
 * 
 * 12feb2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#ifndef NGMCONT_HPP
#define NGMCONT_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

double integrand_ngm_cont( 
               arma::mat exog, arma::mat endog, arma::rowvec cont, 
               arma::rowvec exog_lead, List params, arma::mat coeffs, 
               arma::mat coeffs_cont, int n_exog, int n_endog, int n_cont, 
               int N, arma::rowvec upper, arma::rowvec lower, bool cheby ) ;

arma::rowvec euler_hat_ngm_cont( 
                  arma::mat exog, arma::mat endog, arma::rowvec cont,
                  arma::mat exog_innov_integ, 
                  List params, arma::mat coeffs, arma::mat coeffs_cont, 
                  int n_exog, int n_endog, int n_cont,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights,
                  bool print_rhs ) ;

arma::rowvec con_eqns_ngm( 
  arma::mat exog, arma::mat endog, List params, arma::mat coeffs, 
  arma::mat coeffs_cont, int n_exog, int n_endog, int n_cont, int N, 
  arma::rowvec upper, arma::rowvec lower, bool cheby ) ;

arma::rowvec ngm_reg( 
                  arma::mat exog, arma::mat endog, arma::rowvec cont,
                  arma::mat exog_innov_integ, 
                  List params, arma::mat coeffs, arma::mat coeffs_cont, 
                  int n_exog, int n_endog, int n_cont,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights,
                  bool print_rhs ) ;

#endif