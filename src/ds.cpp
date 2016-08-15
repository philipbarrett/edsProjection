/***********************************************************************************
 * ds.cpp
 * 
 * Code to compute the errors on the Devreux-Sutherland/dynare version of the 
 * IRBC model
 * 
 * 14aug2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#include "ds.hpp"
#include "sim.hpp"

// [[Rcpp::export]]
arma::rowvec integrand_ds( 
                arma::rowvec endog, arma::rowvec exog_lead, 
                double gamma, arma::mat coeffs_cont, 
                int n_exog, int n_endog, int n_cont, int N, 
                arma::rowvec upper, arma::rowvec lower, bool cheby=false ){

  rowvec cont_lead = endog_update( exog_lead, endog, coeffs_cont, n_exog, 
                                    n_endog, N, upper, lower, cheby ) ;
      // Create the next-period control.  Remember, the current-period state is
      // an end-of-period variable, so the control depends directly on its lag.

  rowvec log_integrand(4) ;
  
  double c1_lead = cont_lead(0) ;
  double c2_lead = cont_lead(1) ;
      // The log consumption leads
      // The consumption ratios
  double rb1_lead = cont_lead(2) ;
  double rb2_lead = cont_lead(3) ;
      // The price ratios (inflation rates)
  double q_lead = cont_lead(15) ;
      // The log real exchange rate lead (#16 in R-style counting)
  
  log_integrand(0) = - gamma * c1_lead + rb1_lead ;
  log_integrand(1) = - gamma * c1_lead + rb2_lead ;
  log_integrand(2) = - gamma * c2_lead + rb1_lead - q_lead ;
  log_integrand(3) = - gamma * c2_lead + rb2_lead - q_lead ;
      // The integrands in the log Euler equations
      
  return exp( log_integrand ) ;
}

// [[Rcpp::export]]
arma::rowvec euler_hat_ds( 
                  arma::rowvec exog, arma::rowvec endog, arma::rowvec cont,
                  arma::mat exog_innov_integ, double betta, 
                  double gamma, arma::mat coeffs_cont, 
                  int n_exog, int n_endog, int n_cont,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights,
                  bool print_rhs=false ){
// Computes the single-period error on the Euler equations
  
//  double betta = params["betta"] ;
//  double gamma = params["gamma"] ;
  double rho_pref = - std::log( betta ) ;
      // Extract coefficients
  double NFA = endog( 0 ) ;
  double z1 = endog( 1 ) ;
  double z2 = endog( 2 ) ;
      // Extract the endogenous states
  double c1 = cont(0) ;
  double c2 = cont(1) ;
  double q = cont(15) ;
  double af1 = cont(16) ;
      // Extract controls
  
  mat exog_lead = zeros( n_integ, n_exog ) ;
      // Initalize the draws of the exogenous variables in the next period
  exog_lead = ones(n_integ) * ( rho % exog ) + exog_innov_integ ;
        // Multiply the most recent exogenous draw by the appropriate rho and
        // add the innovation
  
  rowvec integral = zeros<rowvec>( n_endog + 2 ) ;
  mat integrand = zeros( n_integ, n_endog + 2 ) ;
      // Initialize the right hand side.
      // Add the extra two columns for the extra equations for the prices 
      // r1 and r2
      
  for( int i = 0 ; i < n_integ ; i++ ){
    integrand.row(i) = integrand_irbc( endog, exog_lead.row(i), gamma, 
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
  out(0) = z1  - ( gamma * c1 - rho_pref + std::log( integral(0) ) ) ;
  out(1) = af1 - ( gamma * c1 - rho_pref + std::log( integral(1) ) ) ;
  out(2) = z2  - ( gamma * c2 - rho_pref + std::log( integral(2) ) ) ;
  out(2) = q   - ( gamma * c2 - rho_pref + std::log( integral(3) ) ) ;
      // The predictors.  Set up z1, Z_2 s.t. if current consumption is too high
      // for the Euler equation to hold then z1 decreases.  This will pull down 
      // on consumption in the contemporaneous block later.  Similarly, for af1:
      // if the expected return on asset 2 is too high then the investment share
      // in asset 1 declines.  Finally, if q(+1) is too high
  return out ;
      // Return predictors for (B11, B22, r1, r2)
  
}


// [[Rcpp::export]]
arma::rowvec contemp_eqns_ds( 
  arma::mat exog, arma::mat endog, arma::rowvec cont, List params, 
  double y_1_ss=0 ){
// Updates the controls in the contemporaneous block of the model

// Takes as given: z1, z2, af1, q, and guesses for x12, x21
// Outputs values of the contemporaneous controls consistent with these, as well
// as predictors (i.e. new guesses) for x12 and x21


  // Extract parameters
  double alpha = params["alpha"] ;
  double log_alpha = std::log(alpha) ;
  double log_1_alpha = std::log(1-alpha) ;
  double alphahat = params["alphahat"] ;
  double betta = params["betta"] ;
  double eta = params["eta"] ;
  double log_alphahat = std::log( alphahat ) ;
  double log_1_alphahat = std::log( 1 - alphahat ) ;
  double P1_bar = params["P1.bar"] ;
  double P2_bar = params["P2.bar"] ;
  double p1_bar = std::log( P1_bar ) ;
  double p2_bar = std::log( P2_bar ) ;

  rowvec out( cont.n_elem ) ;
      // Initialize the output vector.  Defines the equations for:
      //   c_1, c_2, rb1, rb2, x_11, x_22, x_12, x_21, p_1, p_2, 
      //   p_11, p_22, p_12, p_21, e_12, q, y_1, y_2, cd, cg

  rowvec A = exp( exog.row(0) ) ;
  double A_1 = A(0) ;
  double A_2 = A(1) ;
  double NFA_lag = endog(1,0) ;
  double z1_lag = endog(1,1) ;
  double z2_lag = endog(1,2) ;
      // Extract the states
      
  double q = cont(15) ;
  double af1 = cont(16) ;
  double x_12 = cont(6) ;
  double x_21 = cont(7) ;
      // Extract the controls
  
  double x_11= std::log( std::max( A_1 - std::exp( x_21 ), 1e-08 ) ) ;
  double x_22 = std::log( std::max( A_2 - std::exp( x_12 ), 1e-08 ) ) ;
      // Goods market clearing.  Guards to make sure that logs don't fail here.
  double c_1, c_2 ;
  if( eta == 1.0 ){
    c_1 = alpha * x_11 + ( 1 - alphahat ) * x_12 ;
    c_2 = alpha * x_22 + ( 1 - alphahat ) * x_21 ;
  }else{
    c_1 = eta / ( eta - 1 ) * std::log( alphahat * std::exp( ( 1 - 1 / eta ) * x_11 ) +
              ( 1 - alphahat ) * std::exp( ( 1 - 1 / eta ) * x_12 ) ) ;
    c_2 = eta / ( eta - 1 ) * std::log( alphahat * std::exp( ( 1 - 1 / eta ) * x_22 ) +
              ( 1 - alphahat ) * std::exp( ( 1 - 1 / eta ) * x_21 ) ) ;
  }
  
      // Consumption aggregators

  double p_11 = p1_bar ;
  double p_22 = p2_bar ;
      // Producer prices
  double p_1 = p_11 - 1 / eta * ( log_alpha + c_1 - x_11 ) ;
  double p_2 = p_22 - 1 / eta * ( log_alpha + c_2 - x_22 ) ;
      // Aggregate price levels
  double e = q + p_1 - p_2 ;
      // The nominal exchange rate
  double p_12 = e + p_22 ;
  double p_21 = e + p_11 ;
      // Imported goods prices
  double rb1 = - p_1 - z1_lag ;
  double rb2 = e - p_1 - z2_lag ;
      // Realized returns
  double y_1 = A_1 - p_1 ;
  double y_2 = A_2 - p_2 ;
      // Real incomes
  double NFA = std::exp( rb2 ) * NFA_lag + std::exp( y_1 ) - std::exp( c_1 ) +
                  std::exp( y1_ss ) * betta * af1 * std::exp( rb1 - rb2 ) ;
      // Net foreign assets
  
  /** NEED TO REWRITE FROM HERE **/
  
  // NFA = exp(rb2)*NFA(-1) + exp(Y1) - exp(C1) + exp(STEADY_STATE(Y1))*( BT*af1*(exp(rb1) - exp(rb2)) + zeta );
      
  double x_11_new = c_1 + eta * ( p_1 - p1_bar + log_alphahat ) ;
  double x_22_new = c_2 + eta * ( p_2 - p2_bar + log_alphahat ) ;
      // The resulting factor demands from the remaining optimality condition
  out << c_1 << c_2 << r_1 << r_2 << x_11_new << x_22_new << x_12 << x_21 
             << p_1 << p_2 << p_12 << p_21 << e_12 << endr ;
      // The output vector    
  return( out ) ;
}

// // [[Rcpp::export]]
// arma::rowvec irbc_reg( 
//                   arma::mat exog, arma::mat endog, arma::rowvec cont,
//                   arma::mat exog_innov_integ, 
//                   List params, arma::mat coeffs, arma::mat coeffs_cont, 
//                   int n_exog, int n_endog, int n_cont,
//                   arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
//                   arma::rowvec lower, bool cheby, arma::rowvec weights, 
//                   bool print_rhs=false ){
// // Computes the dependent variables for the regression problem.  Aggregates both
// // the states and controls.
//   double betta = params["betta"] ;
//   double gamma = params["gamma"] ;
// 
//   rowvec out = zeros<rowvec>( n_endog + n_cont ) ;
//       // Initialize the output
//   out.head(n_endog+2) = 
//         euler_hat_irbc( exog, endog, cont, exog_innov_integ, betta, gamma, 
//                             coeffs_cont, n_exog, n_endog, n_cont,
//                             rho, n_integ, N, upper, lower, cheby, weights,
//                             print_rhs ) ;
//       // The endogenous states
//   out.tail(n_cont-2) = 
//         contemp_eqns_irbc( exog, endog, cont, params ) ;
//       // The controls
//   return out ;
// }