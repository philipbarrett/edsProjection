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
double err_ngm( arma::mat exog, arma::mat endog, arma::mat exog_innov_integ, 
                  List params, arma::mat coeffs, int n_exog, int n_endog,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::vec weights ){
// Computes the single-period error on the neocassical growth model equilibrium 
// condition using a Monte Carlo approach.  NB: THE INNOVATIONS exog_innov_integ
// MUST ALREADY BE SCALED (IE. HAVE THE APPROPRIATE VARIANCE)
  
  double betta = params["betta"] ;
      // Extract beta
  mat exog_lead = zeros( n_integ, n_exog ) ;
      // Initalize the draws of the exogenous variables in the next period
  for( int i = 0 ; i < n_integ ; i++ ){
    exog_lead.row(i) = rho % exog.row(0) + exog_innov_integ.row(i) ;
        // Multiply the most recent exogenous draw by the appropriate rho and
        // add the innovation
  }
  
  double rhs = 0 ;
  double err = 0 ;
      // Initialize the right hand side
  for( int i = 0 ; i < n_integ ; i++ ){
    err = betta * integrand_ngm( exog, endog, exog_lead.row(i), params, coeffs, 
                                n_exog, n_endog, N, upper, lower, cheby ) ;
    rhs = rhs + weights(i) * err ;
                                
//      Rcout << "i = " << i << std::endl ;
//      Rcout << "weights(i) = " << weights(i) << std::endl ;
//      Rcout << "err = " << err << std::endl ;
//      Rcout << "rhs = " << rhs << std::endl ;
//      Rcout << "exog_lead.row(i):\n" << exog_lead.row(i) << std::endl ;
      
  }   // Compute the integral
  return fabs( ( 1.0 - rhs ) / 1.0 ) ;
      // Because the target equation is: 1 = beta * E( integrand )
      // Division just to make clear that this is a relative error
}




