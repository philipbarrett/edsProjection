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
                  exp( exog( 0, 0 ) ) * A * pow( endog( 1, 0 ), alpha ) ;
      // Current period consumption
  double c_t1 = ( 1 - delta ) * endog( 0, 0 ) - endog_lead( 0 ) + 
                  exp( exog_lead( 0 ) ) * A * pow( endog( 0, 0 ), alpha ) ;
      // Next period consumption
  double integrand = pow( c_t1 / c_t, - gamma ) * 
          ( 1 - delta + 
            exp( exog(0, 0) ) * alpha * A * pow( endog( 0, 0 ), alpha - 1 ) ) ;
      // Calculate the integrand
      
//      Rcout << "endog:\n" << endog <<std::endl ;
//      Rcout << "exog:\n" << exog <<std::endl ;
//      Rcout << "exog_lead:" << exog_lead <<std::endl ;
//      Rcout << "c_t = " << c_t <<std::endl ;
//      Rcout << "c_t1 = " << c_t1 <<std::endl ;
//      Rcout << "integrand = " << integrand <<std::endl ;
      
  return integrand ;
}

// [[Rcpp::export]]
arma::mat integrand_ngm_D( arma::mat exog, arma::mat endog, 
                arma::rowvec exog_lead, List params, arma::mat coeffs, 
                int n_exog, int n_endog, int N, arma::rowvec upper, 
                arma::rowvec lower, bool cheby=false ){
// Computes the derivative of the integrand_ngm w.r.t the polynomial
// coefficients, computed as:
//    d/d(coeffs)[integrand] = gamma / c_t1 * integrand * basis_vector

  double A = params["A"] ;
  double alpha = params["alpha"] ;
  double delta = params["delta"] ;
  double gamma = params["gamma"] ;
  
  rowvec endog_lead = endog_update( exog_lead, endog.row(0), coeffs, n_exog, 
                                    n_endog, N, upper, lower, cheby ) ;
      // Create the leads of the endogenous variables
  double c_t = ( 1 - delta ) * endog( 1, 0 ) - endog( 0, 0 ) + 
                  exp( exog( 0, 0 ) ) * A * pow( endog( 1, 0 ), alpha ) ;
      // Current period consumption
  double c_t1 = ( 1 - delta ) * endog( 0, 0 ) - endog_lead( 0 ) + 
                  exp( exog_lead( 0 ) ) * A * pow( endog( 0, 0 ), alpha ) ;
      // Next period consumption
  double integrand = pow( c_t1 / c_t, - gamma ) * 
          ( 1 - delta + 
            exp( exog(0, 0) ) * alpha * A * pow( endog( 0, 0 ), alpha - 1 ) ) ;
      // Calculate the integrand
  double common = gamma / c_t1 * integrand ;
      // The common part of the integral
  
  mat basis = zeros<mat>( coeffs.n_rows, coeffs.n_cols ) ;
      // Initialize the matrix of basis polynomials
  mat basis_coeffs = basis ;
      // The matrix of unit coefficients (to be filled in as we go)
  
  basis_coeffs(0,0) = 1 ;
  vec temp = endog_update( exog_lead, endog.row(0), basis_coeffs, n_exog, 
                                    n_endog, N, upper, lower, cheby ) ;
  basis(0) = temp(0) ;
      // Extract the contribution from the first basis term
  for( int i = 1 ; i < coeffs.n_elem ; i++ ){
      basis_coeffs(i-1) = 0 ;
      basis_coeffs(i) = 1 ;
      temp = endog_update( exog_lead, endog.row(0), basis_coeffs, n_exog, 
                                    n_endog, N, upper, lower, cheby ) ;
      basis(i) = temp(0) ;
      // Iterate over the other basis terms
  }
  
  mat out = common * basis ;
  return out ;
}

// [[Rcpp::export]]
double err_ngm( arma::mat exog, arma::mat endog, arma::mat exog_innov_integ, 
                  List params, arma::mat coeffs, int n_exog, int n_endog,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::vec weights,
                  bool print_rhs=false ){
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
  
  vec rhs = zeros(1) ;
  rowvec err = zeros<rowvec>(n_integ) ;
      // Initialize the right hand side
      
  for( int i = 0 ; i < n_integ ; i++ ){
    err(i) = betta * integrand_ngm( exog, endog, exog_lead.row(i), params, coeffs, 
                                n_exog, n_endog, N, upper, lower, cheby ) ;
                                
  }   // Compute the integral

  rhs = err * weights ;
  
    if( print_rhs ){
      Rcout << "err: " << err << std::endl ;
      Rcout << "weights:\n" << weights << std::endl ;
      Rcout << "rhs: " << rhs << std::endl ;
      Rcout << "rhs(0) - 1 = " << rhs(0) - 1 << std::endl ;
    }
  
  double rel_err = ( 1.0 - rhs(0) )  / 1.0 ;
  return pow( rel_err, 2.0 ) ;
//  return fabs( ( 1.0 - rhs ) / 1.0 ) ;
      // Because the target equation is: 1 = beta * E( integrand )
      // Division just to make clear that this is a relative error
}

// [[Rcpp::export]]
arma::mat err_ngm_D( arma::mat exog, arma::mat endog, arma::mat exog_innov_integ, 
                     List params, arma::mat coeffs, int n_exog, int n_endog,
                     arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                     arma::rowvec lower, bool cheby, arma::vec weights ){
// Computes the derivative of the intergated error
  
  double betta = params["betta"] ;
      // Extract beta
  mat exog_lead = zeros( n_integ, n_exog ) ;
      // Initalize the draws of the exogenous variables in the next period
  for( int i = 0 ; i < n_integ ; i++ ){
    exog_lead.row(i) = rho % exog.row(0) + exog_innov_integ.row(i) ;
        // Multiply the most recent exogenous draw by the appropriate rho and
        // add the innovation
  }
  
  mat rhs_D = zeros(coeffs.n_rows, coeffs.n_cols) ;
  mat err_D = rhs_D ;
      // Initialize the right hand side
  double rhs = 0 ;
  double err = 0 ;
      // Initialize the right hand side
  for( int i = 0 ; i < n_integ ; i++ ){
    err_D = betta * integrand_ngm_D( exog, endog, exog_lead.row(i), params, coeffs, 
                                    n_exog, n_endog, N, upper, lower, cheby ) ;
    rhs_D = rhs_D + weights(i) * err_D ;
        // The derivative
    err = betta * integrand_ngm( exog, endog, exog_lead.row(i), params, coeffs, 
                                n_exog, n_endog, N, upper, lower, cheby ) ;
    rhs = rhs + weights(i) * err ;
        // The level.  Needed to make sure that the sign is ok.  Otherwise the 
        // absolute part causes innaccuracies. There is surely a better way to
        // do this.
  }   // Compute the integral
  
  double rel_err = ( 1.0 - rhs )  / 1.0 ;
  return - 2 * rhs_D * rel_err ;
      // Because d/d(coeff)[ .5 * rel_err ^ 2 ] = - d/d(coeff)[rhs] * rel_err
  
//  double sign = ( 1.0 - rhs >= 0 ) ? 1 : -1 ;
//      // The sign of the non-absolute error
//  return - rhs_D * sign ;
}


