/***********************************************************************************
 * irbcSol.hpp
 * 
 * Interface to irbcSol.cpp
 * 
 * 22feb2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#ifndef IRBCSOL_HPP
#define IRBCSOL_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::mat euler_hat_grid( 
            arma::mat coeffs, arma::mat coeffs_cont, 
            arma::mat X, int lags, List params, 
            int n_exog, int n_endog, int n_cont, int n_fwd,
            arma::rowvec rho, arma::rowvec sig_eps, int n_integ,
            int N, arma::rowvec upper, arma::rowvec lower, bool cheby,
            arma::mat exog_innov_mc, bool quad, int n_nodes ) ;

//arma::mat x_eqns_irbc_grid( arma::mat X, int lags, List params,
//                            int n_exog, int n_endog, int n_cont ) ;

arma::mat contemp_eqns_irbc_grid( arma::mat X, int lags, List params,
          int n_exog, int n_endog, int n_cont, List extra_args, std::string model ) ;

#endif