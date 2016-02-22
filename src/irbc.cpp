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
                arma::mat exog, arma::mat endog, arma::rowvec exog_lead, 
                List params, arma::mat coeffs_cont, 
                int n_exog, int n_endog, int n_cont, int N, 
                arma::rowvec upper, arma::rowvec lower, bool cheby=false ){

  double gamma = params["gamma"] ;
  
  rowvec cont_lead = endog_update( exog_lead, endog.row(0), coeffs_cont, n_exog, 
                                    n_endog, N, upper, lower, cheby ) ;
      // Create the next-period control.  Remember, the current-period state is
      // an end-of-period variable, so the control depends directly on its lag.

  rowvec log_integrand(4) ;
  
  double c1_lead = cont_lead(0) ;
  double c2_lead = cont_lead(1) ;
      // The log consumption leads
      // The consumption ratios
  double p1_lead = cont_lead(9) ;
  double p2_lead = cont_lead(10) ;
      // The price ratios (inflation rates)
  double e12_lead = cont_lead(13) ;
      // The log exchange rates lead
  
  log_integrand(0) = - ( gamma * c1_lead + p1_lead ) ;
  log_integrand(1) = - ( gamma * c2_lead + p2_lead ) ;
  log_integrand(2) = - ( e12_lead + gamma * c2_lead + p2_lead ) ;
  log_integrand(3) = e12_lead - gamma * c1_lead - p1_lead ;
      // The integrands in the log Euler equations
      
  return exp( log_integrand ) ;
}

// [[Rcpp::export]]
arma::rowvec euler_hat_irbc( 
                  arma::mat exog, arma::mat endog, arma::rowvec cont,
                  arma::mat exog_innov_integ, 
                  List params, arma::mat coeffs_cont, 
                  int n_exog, int n_endog, int n_cont,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights,
                  bool print_rhs=false ){
// Computes the single-period error on the neocassical growth model equilibrium 
// condition using a Monte Carlo approach.  NB: THE INNOVATIONS exog_innov_integ
// MUST ALREADY BE SCALED (IE. HAVE THE APPROPRIATE VARIANCE)
  
  double betta = params["betta"] ;
  double gamma = params["gamma"] ;
  double rho_pref = - std::log( betta ) ;
      // Extract coefficients
  double B_11 = endog( 0, 0 ) ;
  double B_22 = endog( 0, 1 ) ;
      // Extract the endogenous states
  double c1 = cont(0) ;
  double c2 = cont(1) ;
  double r1 = cont(2) ;
  double r2 = cont(3) ;
  double p1 = cont(9) ;
  double p2 = cont(10) ;
  double e12 = cont(13) ;
      // Extract controls
  
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
    integrand.row(i) = integrand_irbc( exog, endog, exog_lead.row(i), params, 
                                      coeffs_cont, n_exog, n_endog, 
                                      n_cont, N, upper, lower, cheby ) ;
  }   // Compute the integral

  integral = weights * integrand ;
  
    if( print_rhs ){
      Rcout << "err: \n" << integrand << std::endl ;
      Rcout << "weights:" << weights << std::endl ;
      Rcout << "integral: " << integral << std::endl ;
      Rcout << "integral(0) - 1 = " << integral(0) - 1 << std::endl ;
    }
  
  rowvec out(4) ;
  out(0) = B_11 * std::exp( ( c1 - 1 / gamma * ( rho_pref - r1 - p1 
                                              - std::log( integral(0) ) ) ) ) ;
  out(1) = B_22 * std::exp( ( c2 - 1 / gamma * ( rho_pref - r2 - p2 
                                              - std::log( integral(1) ) ) ) ) ;
  out(2) = rho_pref - e12 - gamma * c2 - p2 - std::log( integral(2) ) ;
  out(3) = rho_pref + e12 - gamma * c1 - p1 - std::log( integral(3) ) ;
      // The predictors.  Set up B_11, B_22 s.t. if current consumption is too 
      // high for the Euler equation to hold then B_{t+1} increases.  This will
      // pull down on consumption in the contemporaneous block later.
  return out ;
      // Return predictors for (B11, B22, r1, r2)
  
}


// [[Rcpp::export]]
arma::rowvec x_eqns_irbc( arma::mat exog, arma::rowvec cont, List params ){
// Computes the error on the production function market clearing conditions
// given consumption.  Predicts x_i^i(t) from:
//    x_i^i(t) = log( e^a_i(t) - e^x_j^i(t) )
//    x_j^i(t) = ( alpha * x_j^j(t) - c_j(t) ) / ( 1 - alpha )
  
  double alpha = params["alpha"] ;
      // Extract alpha from params
  double a_1 = exog(0,0) ;
  double a_2 = exog(0,1) ;
      // Extract the exogenous technology processes
  double c_1 = cont(0) ;
  double c_2 = cont(1) ;
      // Consumption
  double x_11_in = cont(5) ;
  double x_22_in = cont(6) ;
      // Own-country intermediate shares
  double x_12 = ( alpha * x_11_in - c_1 ) / ( 1 - alpha ) ;
  double x_21 = ( alpha * x_22_in - c_2 ) / ( 1 - alpha ) ;
      // Cross-country intermediate shares
  double x_11_out = std::log( std::exp( a_1 ) - std::exp( x_21 ) ) ;
  double x_22_out = std::log( std::exp( a_2 ) - std::exp( x_12 ) ) ;
      // The implied own-country shares
      // NB: Do I need a >0 guard here in the log???
  rowvec out(2) ;
  out << x_11_out << x_22_out << endr ;
      // Theb output vector
  return out ;
}


// [[Rcpp::export]]
arma::rowvec contemp_eqns_irbc( 
  arma::mat exog, arma::mat endog, arma::rowvec cont, List params ){
// Computes the block of contemporaneous equations.

  // Extract parameters
  double alpha = params["alpha"] ;
  double eta = params["eta"] ;
  double log_alpha = std::log( alpha ) ;
  double log_1_alpha = std::log( 1 - alpha ) ;
  double P1_bar = params["P1.bar"] ;
  double P2_bar = params["P2.bar"] ;
  double p1_bar = std::log( P1_bar ) ;
  double p2_bar = std::log( P2_bar ) ;

  rowvec out( 9 ) ;
      // Initialize the output vector.  Defines the equations for:
      //   c_1, c_2, x_12, x_21, p_1, p_2, p_12, p_21, e_12

  rowvec A = exp( exog.row(0) ) ;
  double A_1 = A(0) ;
  double A_2 = A(1) ;
  double B_11 = endog(0,0) ;
  double B_22 = endog(0,1) ;
  double B_11_lag = endog(1,0) ;
  double B_22_lag = endog(1,1) ;
      // Extract the states
  double c_1 = cont(0) ;
  double c_2 = cont(1) ;
  double r_1 = cont(2) ;
  double r_2 = cont(3) ;
  double x_11 = cont(4) ;
  double x_22 = cont(5) ;
      // Extract the controls
  
  double x_12 = ( alpha * x_11 - c_1 ) / ( 1 - alpha ) ;
  double x_21 = ( alpha * x_22 - c_2 ) / ( 1 - alpha ) ;
      // Cross-country intermediates
  double p_1 = p1_bar - log_alpha + ( x_11 - c_1 ) / eta ;
  double p_2 = p2_bar - log_alpha + ( x_22 - c_2 ) / eta ;
      // The aggregate price levels
  double p_12 = log_1_alpha - log_alpha + p1_bar + 
                      ( x_11 - c_1 ) / ( eta * ( 1 - alpha ) ) ;
  double p_21 = log_1_alpha - log_alpha + p2_bar + 
                      ( x_22 - c_2 ) / ( eta * ( 1 - alpha ) ) ;
      // The cross-country price levels
  double e_12 = .5 * ( p1_bar - p2_bar + p_1 - p_2 ) ;
      // the nominal exchange rate
  double c_1_new = - p_1 + 
        std::log( A_1 * std::exp( p_1 ) + B_11_lag - B_22_lag * std::exp( e_12 ) 
                  - B_11 * std::exp( - r_1 ) - B_22 * std::exp( e_12 - r_2 ) ) ;
  double c_2_new = - p_2 + 
        std::log( A_2 * std::exp( p_2 ) + B_22_lag - B_11_lag * std::exp( - e_12 ) 
                  - B_22 * std::exp( - r_2 ) - B_11 * std::exp( - e_12 - r_1 ) ) ;
      // The consumption levels impled by the budget constraints
  out << c_1_new << c_2_new << x_12 << x_21 << p_1 << p_2 
                 << p_12 << p_21 << e_12 << endr ;
      // The output vector    
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
                            coeffs_cont, n_exog, n_endog, n_cont,
                            rho, n_integ, N, upper, lower, cheby, weights,
                            print_rhs ) ;
      // The endogenous states
  out.tail(n_cont-2) = 
        contemp_eqns_irbc( exog, endog, cont, params ) ;
      // The controls
  return out ;
}