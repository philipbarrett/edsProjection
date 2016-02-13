/***********************************************************************************
 * ngmCont.hpp
 * 
 * Interface to ngm_cont.cpp
 * 
 * 13feb2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#ifndef NGMCONT2_HPP
#define NGMCONT2_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::rovec integrand_ngm_cont_2( 
               arma::mat exog, arma::mat endog, arma::rowvec cont, 
               arma::rowvec exog_lead, List params, arma::mat coeffs, 
               arma::mat coeffs_cont, int n_exog, int n_endog, int n_cont, 
               int N, arma::rowvec upper, arma::rowvec lower, bool cheby ) ;

arma::rowvec euler_hat_ngm_cont_2( 
                  arma::mat exog, arma::mat endog, arma::rowvec cont,
                  arma::mat exog_innov_integ, 
                  List params, arma::mat coeffs, arma::mat coeffs_cont, 
                  int n_exog, int n_endog, int n_cont,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights,
                  bool print_rhs ) ;

arma::rowvec con_eqns_ngm_2( 
  arma::mat exog, arma::mat endog, List params, arma::mat coeffs, 
  arma::mat coeffs_cont, int n_exog, int n_endog, int n_cont, int N, 
  arma::rowvec upper, arma::rowvec lower, bool cheby ) ;

arma::rowvec ngm_reg_2( 
                  arma::mat exog, arma::mat endog, arma::rowvec cont,
                  arma::mat exog_innov_integ, 
                  List params, arma::mat coeffs, arma::mat coeffs_cont, 
                  int n_exog, int n_endog, int n_cont,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights,
                  bool print_rhs ) ;

#endif