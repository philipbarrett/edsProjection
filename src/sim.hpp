/***********************************************************************************
 * sim.hpp
 * 
 * Interface to sim.cpp
 * 
 * 25jan2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#ifndef SIM_HPP
#define SIM_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::vec ar1_sim( int n_pds, double rho, double sig_eps, 
            bool init_flag, double init ) ;
    // Creates an AR(1) simulation
arma::rowvec endog_update( arma::rowvec exog, arma::rowvec endog_old, arma::mat coeffs, 
                            int n_exog, int n_endog, int N,
                            arma::rowvec upper, arma::rowvec lower, bool cheby ) ;
    // Uses the updating rule defined by coeffs, upper, lower, and cheby to
    // update the endogenous states
arma::mat endog_sim( int n_out, arma::mat exog_sim, arma::mat coeffs, int N,
                      arma::rowvec upper, arma::rowvec lower, arma::rowvec endog_init, 
                      bool cheby, int kappa, int burn, bool lag ) ;
    // Returns a matrix of simulated exogenous and endogenous variables
arma::mat cont_sim( arma::mat xn_sim, arma::mat coeffs_cont, int N,
                      int n_endog, int n_exog, int n_cont, 
                      arma::rowvec upper, arma::rowvec lower, bool cheby ) ;
    // Creates a simulated path of the control variables from a simulation of 
arma::mat irf_create( int pds, int n_sim, int N, int shk_idx,
                      arma::rowvec rho, arma::rowvec sig_eps, 
                      arma::mat coeffs, arma::rowvec upper, arma::rowvec lower, 
                      arma::rowvec init, int n_endog, int n_exog, 
                      double shk, bool cheby ) ;
    // Creates the impulse responses

#endif