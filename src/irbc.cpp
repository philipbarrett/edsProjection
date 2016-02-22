/***********************************************************************************
 * irbc.cpp
 * 
 * Code to compute the errors on the international RBC model with asset market 
 * restrictions
 * 
 * 14feb2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#include "irbc.hpp"
#include "sim.hpp"

// [[Rcpp::export]]
arma::rowvec integrand_irbc( 
                arma::mat exog, arma::mat endog, arma::rowvec cont, 
                arma::rowvec exog_lead, List params, arma::mat coeffs, 
                arma::mat coeffs_cont, int n_exog, int n_endog, int n_cont, 
                int N, arma::rowvec upper, arma::rowvec lower, bool cheby=false ){

  double gamma = params["gamma"] ;
  
  rowvec cont_lead = endog_update( exog_lead, endog.row(0), coeffs_cont, n_exog, 
                                    n_endog, N, upper, lower, cheby ) ;
      // Create the next-period control.  Remember, the current-period state is
      // an end-of-period variable, so the control depends directly on its lag.

  rowvec integrand(4) ;
  
  double R1 = cont(0) ;
  double R2 = cont(1) ;
      // The nominal interest rates
  double cons_r1 = cont_lead(2) / cont(2) ;
  double cons_r2 = cont_lead(3) / cont(3) ;
      // The consumption ratios
  double p_r1 = cont_lead(4) / cont(4) ;
  double p_r2 = cont_lead(5) / cont(5) ;
      // The price ratios (inflation rates)
  double e_r12 = cont_lead(6) / cont(6) ;
      // The ratio of exchange rates (appreciation/depreciation)
  
  integrand(0) = R1 * pow( cons_r1, - gamma ) / p_r1 ;
  integrand(1) = R2 * pow( cons_r2, - gamma ) / p_r2 ;
  integrand(2) = pow( cons_r2, - gamma ) / ( p_r2 * e_r12 ) ;
  integrand(3) = pow( cons_r1, - gamma ) * e_r12 / p_r1 ;
      // Two countries have same consumption but diefferent capital stocks and
      // technologies
      
  return integrand ;
}

// [[Rcpp::export]]
arma::rowvec euler_hat_irbc( 
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
  
  rowvec integral = zeros<rowvec>( n_endog + 2 ) ;
  mat integrand = zeros( n_integ, n_endog + 2 ) ;
      // Initialize the right hand side.
      // Add the extra two columns for the extra equations for the prices R1 and R2
      
  for( int i = 0 ; i < n_integ ; i++ ){
    integrand.row(i) = integrand_irbc( exog, endog, cont, exog_lead.row(i), params, 
                                      coeffs, coeffs_cont, n_exog, n_endog, 
                                      n_cont, N, upper, lower, cheby ) ;
  }   // Compute the integral

  integral = weights * integrand ;
  
    if( print_rhs ){
      Rcout << "err: \n" << integrand << std::endl ;
      Rcout << "weights:" << weights << std::endl ;
      Rcout << "integral: " << integral << std::endl ;
      Rcout << "integral(0) - 1 = " << integral(0) - 1 << std::endl ;
    }
  
  rowvec endog_hat(4) ;
  endog_hat(0) = endog(0,0) * betta * integral(0) ;
  endog_hat(1) = endog(0,1) * betta * integral(1) ;
  endog_hat(2) = 1 / ( betta * integral(2) ) ;
  endog_hat(3) = 1 / ( betta * integral(3) ) ;
      // The predictors
  return endog_hat ;
      // Return predictors for (B11, B22, R1, R2)
  
}


// [[Rcpp::export]]
arma::rowvec cont_eqns_irbc( 
  arma::mat exog, arma::mat endog, arma::rowvec cont, List params, arma::mat coeffs, 
  arma::mat coeffs_cont, int n_exog, int n_endog, int n_cont, int N, 
  arma::rowvec upper, arma::rowvec lower, bool cheby=false ){
// Computes the predicted controls.

  // Extract parameters
  double alpha = params["alpha"] ;
  double P1_bar = params["P1.bar"] ;
  double P2_bar = params["P2.bar"] ;

  rowvec out = zeros<rowvec>(n_cont - 2) ;
      // Initialize the output vector.  Lose two equations to the Euler part
  
  double R_1 = cont(0) ;
  double R_2 = cont(1) ;
  double C_1 = cont(2) ;
  double C_2 = cont(3) ;
  double P_1 = cont(4) ;
  double P_2 = cont(5) ;
  double e_12 = cont(6) ;
  double X_11 = cont(7) ;
  double X_22 = cont(8) ;
  double X_12 = cont(9) ;
      // Extract the controls
  rowvec A = exp( exog.row(0) ) ;
  double A_1 = A(0) ;
  double A_2 = A(1) ;
  double B_11 = endog(0,0) ;
  double B_22 = endog(0,1) ;
  double B_11_lag = endog(1,0) ;
  double B_22_lag = endog(1,1) ;
      // Extract the states
  
  out(0) = pow( X_11, alpha ) * pow( X_12, 1 - alpha ) ;
  out(1) = pow( X_22, alpha ) * pow( A_1 - X_11, 1 - alpha ) ;
      // Prodcution technologies
  out(2) = ( A_1 + B_11_lag - e_12 * B_22_lag - B_11 / R_1 + 
                                      e_12 * B_22 / R_2 ) / C_1 ;
  out(3) = ( A_2 + B_22_lag - B_11_lag / e_12  - B_22 / R_2 + 
                                      B_11 / R_1 / e_12 ) / C_2 ;
      // Budget constraints
//  out(4) = ( 1 - alpha ) * C_1 / ( A_2 - X_22 ) * P_1 / P2_bar ;
  out(4) = pow( C_1 / C_2 * P_1 / P_2 * (A_1-X_11) / X_12, .5 ) ;
      // Definitiion of nominal exchange rate
  out(5) = alpha * C_1 * P_1 / P1_bar ;
  out(6) = alpha * C_2 * P_2 / P2_bar ;
      // Optimal intermediate inputs
  out(7) = pow( 1 / alpha - 1, 2 ) * X_11 * X_22 / ( A_1 - X_11 ) ;
      // The law of one price
      
  // NB: Walras Law => the GMC cond for country 2's production is redundant 
  
  return( out ) ;
}

// [[Rcpp::export]]
arma::rowvec irbc_reg( 
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
  out.head(n_endog+2) = 
        euler_hat_irbc( exog, endog, cont, exog_innov_integ, params, 
                            coeffs, coeffs_cont, n_exog, n_endog, n_cont,
                            rho, n_integ, N, upper, lower, cheby, weights,
                            print_rhs ) ;
      // The endogenous states
  out.tail(n_cont-2) = 
        cont_eqns_irbc( exog, endog, cont, params, coeffs, coeffs_cont, n_exog, 
                      n_endog, n_cont, N, upper, lower, cheby ) ;
      // The controls
  return out ;
}