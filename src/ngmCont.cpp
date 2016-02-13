/***********************************************************************************
 * ngmCont.cpp
 * 
 * Code to compute the errors on the neoclassical growth model with controls
 * 
 * 12feb2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#include "ngmCont.hpp"
#include "sim.hpp"

// [[Rcpp::export]]
double integrand_ngm_cont( 
               arma::mat exog, arma::mat endog, arma::rowvec cont, 
               arma::rowvec exog_lead, List params, arma::mat coeffs, 
               arma::mat coeffs_cont, int n_exog, int n_endog, int n_cont, 
               int N, arma::rowvec upper, arma::rowvec lower, bool cheby=false ){
// Computes the integrad for a particular realization of the exogenous states. 
// Inputs are:
//   - exog   A matrix of lags of the exogenous states.  Eg. exog(1,1) is the 
//            first lag of the second exogenou state, and exog(0,1) is the 
//            current value of the second exogenous state.
//   - endog  A matrix of lags of the endogenous states.
//   - exog_lead  A vector of realizations of the exogenous variable in the 
//            next period
//   - params A list of parameters
//   - coeffs A matrix of parameters of the approximation of the endogenous 
//            states
//   - coeffs_cont The equivalent matrix of coefficients for the controls
  
  double A = params["A"] ;
  double alpha = params["alpha"] ;
  double delta = params["delta"] ;
  double gamma = params["gamma"] ;
  
  rowvec cont_lead = endog_update( exog_lead, endog.row(0), coeffs_cont, n_exog, 
                                    n_endog, N, upper, lower, cheby ) ;
      // Create the next-period control.  Remember, the current-period state is
      // an end-of-period variable, so the control depends directly on its lag.
      
//    Rcout << "cont: " << cont << std::endl ;
//    Rcout << "endog_lead: " << endog_lead << std::endl ;
//    Rcout << "cont_lead: " << cont_lead << std::endl ;
      
  double integrand = pow( cont_lead(0) / cont(0), - gamma ) * 
          ( 1 - delta + 
            exp( exog_lead( 0 ) ) * alpha * A * pow( endog( 0, 0 ), alpha - 1 ) ) ;
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
arma::rowvec euler_hat_ngm_cont( 
                  arma::mat exog, arma::mat endog, arma::rowvec cont,
                  arma::mat exog_innov_integ, 
                  List params, arma::mat coeffs, arma::mat coeffs_cont, 
                  int n_exog, int n_endog, int n_cont,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights,
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
  
  vec rhs = zeros(n_endog) ;
  vec err = zeros(n_integ) ;
      // Initialize the right hand side
      
  for( int i = 0 ; i < n_integ ; i++ ){
    err(i) = betta * 
              integrand_ngm_cont( exog, endog, cont, exog_lead.row(i), params, 
                                  coeffs, coeffs_cont, n_exog, n_endog, 
                                  n_cont, N, upper, lower, cheby ) ;
  }   // Compute the integral

  rhs = weights * err ;
  
    if( print_rhs ){
      Rcout << "err: \n" << err << std::endl ;
      Rcout << "weights:" << weights << std::endl ;
      Rcout << "rhs: " << rhs << std::endl ;
      Rcout << "rhs(0) - 1 = " << rhs(0) - 1 << std::endl ;
    }
  
  rowvec endog_hat = rhs % endog.row(0) ;
  return endog_hat ;
      // Return k * the Euler Error
  
}

// [[Rcpp::export]]
arma::rowvec con_eqns_ngm( 
  arma::mat exog, arma::mat endog, List params, arma::mat coeffs, 
  arma::mat coeffs_cont, int n_exog, int n_endog, int n_cont, int N, 
  arma::rowvec upper, arma::rowvec lower, bool cheby=false ){
// Computes the predicted controls.  In general can depend on the controls. 
// Happens not to in this case.

  // Extract parameters
  double A = params["A"] ;
  double alpha = params["alpha"] ;
  double delta = params["delta"] ;

  rowvec out = zeros<rowvec>(n_cont) ;
      // Initialize the output vector
  out(0) = ( 1 - delta ) * endog( 1, 0 ) - endog( 0, 0 ) + 
                  exp( exog( 0, 0 ) ) * A * pow( endog( 1, 0 ), alpha ) ;
  return( out ) ;
}

// [[Rcpp::export]]
arma::rowvec ngm_reg( 
                  arma::mat exog, arma::mat endog, arma::rowvec cont,
                  arma::mat exog_innov_integ, 
                  List params, arma::mat coeffs, arma::mat coeffs_cont, 
                  int n_exog, int n_endog, int n_cont,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights, 
                  bool print_rhs=false ){
// Computes the dependent variables for the regression problem.  Aggregates both
// the states and controls.
  rowvec out = zeros<rowvec>( n_endog + n_cont ) ;
      // Initialize the output
  out.head(n_endog) = 
        euler_hat_ngm_cont( exog, endog, cont, exog_innov_integ, params, 
                            coeffs, coeffs_cont, n_exog, n_endog, n_cont,
                            rho, n_integ, N, upper, lower, cheby, weights,
                            print_rhs ) ;
      // The endogenous states
  out.tail(n_cont) = 
        con_eqns_ngm( exog, endog, params, coeffs, coeffs_cont, n_exog, 
                      n_endog, n_cont, N, upper, lower, cheby ) ;
      // The controls
  return out ;
}