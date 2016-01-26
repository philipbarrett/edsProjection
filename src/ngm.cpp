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
  double betta = params["betta"] ;
  
      Rcout << "Bing" << std::endl ;
  
  rowvec endog_lead = endog_update( exog_lead, endog.row(0), coeffs, n_exog, 
                                    n_endog, N, upper, lower, cheby ) ;
      // Create the leads of the endogenous variables
  
      Rcout << "Bong" << std::endl ;
  
  double c_t = ( 1 - delta ) * endog( 1, 0 ) - endog( 0, 0 ) + 
                  exp( exog( 0, 0 ) ) * alpha * A * pow( endog( 1, 0 ), alpha ) ;
      // Current period consumption
  double c_t1 = ( 1 - delta ) * endog( 0, 0 ) - endog_lead( 0 ) + 
                  exp( exog_lead( 0 ) ) * alpha * A * pow( endog( 0, 0 ), alpha ) ;
      // Next period consumption
      
      Rcout << "Bang" << std::endl ;
      
  double integrand = pow( c_t / c_t1, - gamma ) * 
                ( 1 - delta * exp( exog(0) ) * alpha * A * ( endog( 0, 0 ) ) ) ;
      // Calculate the integrand
  return integrand ;
}






