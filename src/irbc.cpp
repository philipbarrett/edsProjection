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
  double p1_lead = cont_lead(8) ;
  double p2_lead = cont_lead(9) ;
      // The price ratios (inflation rates)
  double e12_lead = cont_lead(12) ;
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
                  arma::rowvec exog, arma::rowvec endog, arma::rowvec cont,
                  arma::mat exog_innov_integ, double betta, double theta,
                  double gamma, arma::mat coeffs_cont, 
                  int n_exog, int n_endog, int n_cont, int n_fwd,
                  arma::mat rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights,
                  bool print_rhs=false ){
// Computes the single-period error on the neocassical growth model equilibrium 
// condition using a Monte Carlo approach.  NB: THE INNOVATIONS exog_innov_integ
// MUST ALREADY BE SCALED (IE. HAVE THE APPROPRIATE VARIANCE)
  
//  double betta = params["betta"] ;
//  double gamma = params["gamma"] ;
  double rho_pref = - std::log( betta ) ;
      // Extract coefficients
  double B_11 = endog( 0 ) ;
  double B_22 = endog( 1 ) ;
      // Extract the endogenous states
  double c1 = cont(0) ;
  double c2 = cont(1) ;
  double r1 = cont(2) ;
  double r2 = cont(3) ;
  double p1 = cont(8) ;
  double p2 = cont(9) ;
  double e12 = cont(12) ;
      // Extract controls

  mat exog_lead = zeros( n_integ, n_exog ) ;
      // Initalize the draws of the exogenous variables in the next period
  exog_lead = ones(n_integ) * ( exog * rho ) + exog_innov_integ ;
        // Multiply the most recent exogenous draw by the appropriate rho and
        // add the innovation
        // Would be easy to convert to VAR(1) form here. DONE!
  rowvec integral = zeros<rowvec>( n_fwd ) ;
  mat integrand = zeros( n_integ, n_fwd ) ;
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
//  out(0) = B_11 * std::exp( ( c1 - 1 / gamma * ( rho_pref - r1 - p1 
//                                              - std::log( integral(0) ) ) ) ) ;
//  out(1) = B_22 * std::exp( ( c2 - 1 / gamma * ( rho_pref - r2 - p2 
//                                              - std::log( integral(1) ) ) ) ) ;
//  out(2) = rho_pref - e12 - gamma * c2 - p2 - std::log( integral(2) ) ;
//  out(3) = rho_pref + e12 - gamma * c1 - p1 - std::log( integral(3) ) ;
  out(0) = B_11 -  
              ( ( gamma - theta ) * c2 - rho_pref + r1 + e12 + p2 + std::log( integral(2) ) ) ;
  out(1) = ( gamma - theta ) * c1 - rho_pref + r2 + p1 + std::log( integral(3) ) ;
  out(2) = rho_pref - ( gamma - theta ) * c1 - p1 - std::log( integral(0) ) ;
  out(3) = rho_pref - ( gamma - theta ) * c2 - p2 - std::log( integral(1) ) ;
      // The predictors, for B11, e_12, R1 and R2.  Set up B_11 s.t. if current consumption is too 
      // high for the Euler equation to hold then B_{t+1} increases.  This will
      // pull down on consumption in the contemporaneous block later.
  return out ;
      // Return predictors for (B11, e_12, r1, r2)
  
}


// [[Rcpp::export]]
arma::rowvec contemp_eqns_irbc( 
  arma::mat exog, arma::mat endog, arma::rowvec cont, List params, List extra_args ){
// Updates the controls in the contemporaneous block of the model
// Takes as given: B11, B22, r1, r2, and guesses for x11, x22
// Outputs values of the contemporaneous controls consistent with these, as well
// as predictors (i.e. new guesses) for x11 and x22


  // Extract parameters
  double alpha = params["share"] ;
  double log_alpha = std::log(alpha) ;
  double log_1_alpha = std::log(1-alpha) ;
  double eta = params["eta"] ;
  double mu = params["mu"] ;
  double log_mu = std::log(mu) ;
  double log_1_mu = std::log(1-mu) ;
  double xi = params["xi"] ;
  double alphahat = pow( alpha, 1 / eta )  ;
  double alphahat_1 = pow( 1 - alpha, 1 / eta )  ;
  double mu_hat = pow( mu, 1 / xi )  ;
  double mu_1_hat = pow( 1 - mu, 1 / xi )  ;

  rowvec out( cont.n_elem + endog.n_cols ) ;
      // Initialize the output vector.  Defines the equations for:
      //   c_1, c_2, r_1, r_2, x_11, x_22, x_12, x_21, p_1, p_2, p_12, p_21, e_12

  double A_1 = mu * exp( exog(0,0) ) ;
  double A_2 = mu * exp( exog(0,1) ) ;
      // Because the shock is mean-zero
  double log_A_1 = std::log( A_1 ) ;
  double log_A_2 = std::log( A_2 ) ; 
      // The log of the level of the shock
  double p_11 = exog(0,2) ;
  double p_22 = exog(0,3) ;
      
  double B_11 = endog(0,0) ;
      // The edogenous state solved from the forward-looking equations
  double B_11_lag = endog(1,0) ;
  double B_22_lag = endog(1,1) ;
      // The lagged states
  double r_1 = cont(2) ;
  double r_2 = cont(3) ;
  double e_12 = cont(12) ;
      // The controls in the foward-looking part of the model
  double x_12 = cont(6) ;
  double x_21 = cont(7) ;
      // The hyopthesized controls (drives the conetmporaneous block)
  
  double cn_1 = std::log( 1 - mu ) ;
  double cn_2 = std::log( 1 - mu ) ;
  double x_11= std::log( std::max( A_1 - std::exp( x_21 ), 1e-08 ) ) ;
  double x_22 = std::log( std::max( A_2 - std::exp( x_12 ), 1e-08 ) ) ;
      // Goods market clearing.  Guards to make sure that logs don't fail here.
  double ct_1, ct_2 ;
  if( eta == 1.0 ){
    ct_1 = alphahat * x_11 + alphahat_1 * x_12 ;
    ct_2 = alphahat * x_22 + alphahat_1 * x_21 ;
  }else{
    ct_1 = eta / ( eta - 1 ) * std::log( alphahat * std::exp( ( 1 - 1 / eta ) * x_11 ) +
      alphahat_1 * std::exp( ( 1 - 1 / eta ) * x_12 ) ) ;
    ct_2 = eta / ( eta - 1 ) * std::log( alphahat * std::exp( ( 1 - 1 / eta ) * x_22 ) +
      alphahat_1 * std::exp( ( 1 - 1 / eta ) * x_21 ) ) ;
  }
      // Tradeable goods consumption aggregators
  double c_1, c_2 ;
  if( xi == 1.0 ){
    c_1 = mu_hat * x_11 + mu_1_hat * x_12 ;
    c_2 = mu_hat * x_22 + mu_1_hat * x_21 ;
  }else{
    c_1 = xi / ( xi - 1 ) * std::log( mu_hat * std::exp( ( 1 - 1 / xi ) * ct_1 ) +
      mu_1_hat * std::exp( ( 1 - 1 / xi ) * cn_1 ) ) ;
    c_2 = xi / ( xi - 1 ) * std::log( mu_hat * std::exp( ( 1 - 1 / xi ) * ct_2 ) +
      mu_1_hat * std::exp( ( 1 - 1 / xi ) * cn_2 ) ) ;
  }
      // Final goods consumption aggregators

  double pt_1 = p_11 - 1 / eta * ( log_alpha + ct_1 - x_11 ) ;
  double pt_2 = p_22 - 1 / eta * ( log_alpha + ct_2 - x_22 ) ;
  double p_1 = pt_1 - 1 / xi * ( log_mu + c_1 - ct_1 ) ;
  double p_2 = pt_2 - 1 / xi * ( log_mu + c_2 - ct_2 ) ;
  double pn_1 = p_1 + 1 / xi * ( log_1_mu + c_1 - cn_1 ) ;
  double pn_2 = p_2 + 1 / xi * ( log_1_mu + c_2 - cn_2 ) ;
      // Tradeables price levels, from factor demands
  double q_12 = e_12 - p_1 + p_2 ;
      // The real exchange rate
  double p_12 = p_22 + e_12 ;
  double p_21 = p_11 - e_12 ;
      // Imported goods prices
  double B_22 = std::exp( r_2 ) * ( 
            + B_22_lag - std::exp( - e_12 ) * ( A_1 * exp( p_11 ) + exp( cn_1 + pn_1 ) - 
              B_11 * std::exp( - r_1 ) + B_11_lag - std::exp( c_1 + p_1 ) ) ) ;
      // Debt in country 2
  double x_12_new = ct_1 + log_1_alpha - eta * ( p_12 - pt_1 ) ;
  double x_21_new = ct_2 + log_1_alpha - eta * ( p_21 - pt_2 ) ;
      // The resulting factor demands from the remaining optimality condition
  out << B_11 << B_22 <<
          c_1 << c_2 << r_1 << r_2 << x_11 << x_22 << x_12_new << x_21_new <<
          p_1 << p_2 << p_12 << p_21 << e_12 << q_12 << log_A_1 << log_A_2 <<
          pn_1 << pn_2 << pt_1 << pt_2 << cn_1 << cn_2 << ct_1 << ct_2 << endr ;
      // The output vector    
      
  return( out ) ;
}

// [[Rcpp::export]]
arma::rowvec irbc_reg( 
                  arma::mat exog, arma::mat endog, arma::rowvec cont,
                  arma::mat exog_innov_integ, 
                  List params, arma::mat coeffs, arma::mat coeffs_cont, 
                  int n_exog, int n_endog, int n_cont, int n_fwd,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights, 
                  List extra_args, bool print_rhs=false ){
// Computes the dependent variables for the regression problem.  Aggregates both
// the states and controls.
  double betta = params["betta"] ;
  double gamma = params["gamma"] ;
  double theta = params["theta"] ;

  rowvec out = zeros<rowvec>( n_endog + n_cont ) ;
      // Initialize the output
  out.head(n_endog+2) = 
        euler_hat_irbc( exog, endog, cont, exog_innov_integ, betta, theta, gamma, 
                            coeffs_cont, n_exog, n_endog, n_cont, n_fwd,
                            rho, n_integ, N, upper, lower, cheby, weights,
                            print_rhs ) ;
      // The endogenous states
  out.tail(n_cont-2) = 
        contemp_eqns_irbc( exog, endog, cont, params, extra_args ) ;
      // The controls
  return out ;
}