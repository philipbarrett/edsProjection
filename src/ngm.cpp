/***********************************************************************************
 * ngm.cpp
 * 
 * Code to compute the errors on the neoclassical growth model
 * 
 * 26jan2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#include "ngm.hpp"
#include "sim.hpp"

// [[Rcpp::export]]
double integrand_ngm( arma::mat exog, arma::mat endog, arma::rowvec exog_lead, 
                       List params, arma::mat coeffs, int n_exog, int n_endog, 
                       int N, arma::rowvec upper, arma::rowvec lower, 
                       bool cheby=false ){
// Computes the integrad for a particular realization of the exogenous states. 
// Inputs are:
//   - exog   A matrix of lags of the exogenous states.  Eg. exog(1,1) is the 
//            first lag of the second exogenou state, and exog(0,1) is the 
//            current value of the second exogenous state.
//   - endog  A matrix of lags of the endogenous states.
//   - exog_lead  A vector of realizations of the exogenous variable in the 
//            next period
//   - params A list of parameters
//   - coeffs A matrix of parameters of the approximation
  
  double A = params["A"] ;
  double alpha = params["alpha"] ;
  double delta = params["delta"] ;
  double gamma = params["gamma"] ;
  
  rowvec endog_lead = endog_update( exog_lead, endog.row(0), coeffs, n_exog, 
                                    n_endog, N, upper, lower, cheby ) ;
      // Create the leads of the endogenous variables
  double c_t = ( 1 - delta ) * endog( 1, 0 ) - endog( 0, 0 ) + 
                  exp( exog( 0, 0 ) ) * alpha * A * pow( endog( 1, 0 ), alpha ) ;
      // Current period consumption
  double c_t1 = ( 1 - delta ) * endog( 0, 0 ) - endog_lead( 0 ) + 
                  exp( exog_lead( 0 ) ) * alpha * A * pow( endog( 0, 0 ), alpha ) ;
      // Next period consumption
  double integrand = pow( c_t / c_t1, - gamma ) * 
                ( 1 - delta * exp( exog(0) ) * alpha * A * ( endog( 0, 0 ) ) ) ;
      // Calculate the integrand
  return integrand ;
}

// [[Rcpp::export]]
double err_ngm_mc( arma::mat exog, arma::mat endog, arma::mat exog_innov, 
                  List params, arma::mat coeffs, int n_exog, int n_endog,
                  arma::rowvec rho, arma::rowvec sig_eps, int n_mc,
                  int N, arma::rowvec upper, arma::rowvec lower, bool cheby=false ){
// Computes the single-period error on the neocassical growth model equilibrium 
// condition using a Monte Carlo approach.  NB: THE INNOVATIONS EXOG_INNOV MUST
// BE STANDARD, IE. MEAN 0, VARIANCE 1
  
  double betta = params["betta"] ;
      // Extract beta
  mat exog_lead = zeros( n_mc, n_exog ) ;
      // Initalize the draws of the exogenous variables in the next period
  for( int i = 0 ; i < n_mc ; i++ ){
    exog_lead.row(i) = rho % exog.row(0) + sig_eps % exog_innov.row(i) ;
        // Multiply the most recent exogenous draw by the appropriate rho and
        // add the innovation
  }
  
  double rhs = 0 ;
      // Initialize the righ hand side
  for( int i = 0 ; i < n_mc ; i++ ){
    rhs = rhs + betta / n_mc * 
                integrand_ngm( exog, endog, exog_lead.row(i), params, coeffs, 
                                n_exog, n_endog, N, upper, lower, cheby ) ;
  }   // Compute the integral via monte carlo
  return ( 1 - rhs ) ;
      // Because the target equation is: 1 = beta * E( integrand )
}




