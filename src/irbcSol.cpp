/***********************************************************************************
 * irbcSol.cpp
 * 
 * Code to compute the errors on the international RBC model with asset market 
 * restrictions on a matrix of state grid points
 * 
 * 22feb2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#include "irbcSol.hpp"
#include "irbc.hpp"
#include "ds.hpp"
#include "sim.hpp"
#include "quad.hpp"
#include "mono.hpp"

// [[Rcpp::export]]
arma::mat euler_hat_grid( 
            arma::mat coeffs, arma::mat coeffs_cont, 
            arma::mat X, int lags, List params, 
            int n_exog, int n_endog, int n_cont, int n_fwd,
            arma::rowvec rho, arma::rowvec sig_eps, int n_integ,
            int N, arma::rowvec upper, arma::rowvec lower, bool cheby,
            arma::mat exog_innov_mc, bool quad=true, int n_nodes=0,
            std::string model="irbc", std::string mono="none" ){
// Creates the vector of errors on the Euler equations

  int n_pts = X.n_rows ;
      // The number of points at which the error is assessed
  rowvec exog = zeros<rowvec>( n_exog ) ;
  rowvec endog = zeros<rowvec>( n_endog ) ;
  rowvec cont = zeros<rowvec>( std::max( n_cont, 1 ) ) ;
      // Temporary containers used in the loop.  Make cont bigger than size 0
      // here - just passing a useless empty container
  mat err = zeros(n_pts, n_fwd ) ;
      // Becuase n_fwd is the number of forward-looking expectational equations
  double betta = params["betta"] ;
  double gamma = params["gamma"] ;

  /** Create the integration nodes and weights **/
  n_integ = quad ? pow( n_nodes, n_exog ) : n_integ ;
      // Update the number of points if using quadrature
  if( mono == "m1" ) n_integ = 2 * n_exog ;
      // Monomial rule
  rowvec weights( n_integ ) ;
  mat nodes( n_integ, n_exog ) ;
  vec v_sig_eps = conv_to<vec>::from( sig_eps ) ;
      // The weights and integration nodes
      
  if( quad ){
    if(mono == "m1"){
      mat m_sig2 = eye( n_exog, n_exog ) ;
          // Initiate covariance matrix
      for ( int i = 0 ; i < n_exog ; i++ ){
        m_sig2(i,i) = pow( v_sig_eps(i), 2.0 ) ;
      }
          // Create the variance-covariance matrix in the uncorrelated case
      mat m_quad = M1_nodes_weights_mat( zeros(n_exog), m_sig2 ) ;
      vec temp = m_quad.col( n_exog ) ;
      weights = conv_to<rowvec>::from( temp ) ;
      nodes = m_quad.head_cols(n_exog) ;
          // Monomial rule
    }
    else
    {
      mat m_quad = quad_nodes_weights_mat( n_nodes, n_exog, 
                          v_sig_eps, zeros(n_exog) ) ;
      vec temp = m_quad.col( n_exog ) ;
      weights = conv_to<rowvec>::from( temp ) ;
      nodes = m_quad.head_cols(n_exog) ;
          // Quadrature
    }
  }
  else
  {
    weights = ones<rowvec>( n_pts ) / n_integ ;
    nodes = exog_innov_mc ;
        // Monte Carlo integration
  }
  
  /** Define the model function **/
  rowvec (*euler_hat_fn)( 
             arma::rowvec exog, arma::rowvec endog, arma::rowvec cont,
                  arma::mat exog_innov_integ, double betta, 
                  double gamma, arma::mat coeffs_cont, 
                  int n_exog, int n_endog, int n_cont, int n_fwd,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights,
                  bool print_rhs ) ;
      // Pointer to model evaluation function
  
  if( model == "irbc" )
    euler_hat_fn = euler_hat_irbc ;
        // The first attempt at the Adams-Barrett model
  if( model == "ds" )
    euler_hat_fn = euler_hat_ds ;
        // The Devreux-Sutherland version
  
  
  /** Now compute the model errors **/
  for( int i = 0 ; i < n_pts ; i++ ){
  // Loop over the evaluation points
  
    exog = X.row(i).head(n_exog) ;
    endog = X.row(i).subvec( n_exog, n_exog + n_endog - 1 ) ;
        // Fill in the endogenous and exogenous matrices
    if( n_cont > 0 )
      cont = X.row(i).tail( n_cont ) ;
        // The controls
    err.row(i) = euler_hat_fn(
                    exog, endog, cont, nodes, betta, gamma,
                    coeffs_cont, n_exog, n_endog, n_cont, n_fwd, rho, 
                    n_integ, N, upper, lower, cheby, weights, false ) ;
  }   // The error on the states according to the Euler equations
  
  return err ;
}

// [[Rcpp::export]]
arma::mat contemp_eqns_irbc_grid( arma::mat X, int lags, List params,
                          int n_exog, int n_endog, int n_cont, List extra_args,
                          std::string model="irbc" ){
// Computes contemp_eqns_irbc on the state/control grid
 
  int n_fwd = extra_args["n.fwd"] ;
      // Number of forward-looking equations
  int n_pts = X.n_rows ;
      // The number of points at which the error is assessed
  mat exog = zeros( 1 + lags, n_exog ) ;
  mat endog = zeros( 1 + lags, n_endog ) ;
  rowvec cont = zeros<rowvec>( std::max( n_cont, 1 ) ) ;
      // Temporary containers used in the loop.  Make cont bigger than size 0
      // here - just passing a useless empty container
  mat err = zeros(n_pts, n_cont + n_endog ) ;
      // The number of static equations to be solved is the number of
      // controls plus one for each state that cannot be paired with a control.
  
    /** Define the model function **/
  rowvec (*contemp_eqns_fn)( 
             arma::mat exog, arma::mat endog, arma::rowvec cont, List params, 
             List extra_args ) ;
      // Pointer to model evaluation function
  
  if( model == "irbc" )
    contemp_eqns_fn = contemp_eqns_irbc ;
        // The first attempt at the Adams-Barrett model
  if( model == "ds" )
    contemp_eqns_fn = contemp_eqns_ds ;
        // The Devreux-Sutherland version
  
  
  /** Now compute the model errors **/
  for( int i = 0 ; i < n_pts ; i++ ){
  // Loop over the evaluation points
    for( int j = 0 ; j < 1 + lags ; j++ ){
    // Loop over the lags
      exog.row(j) = X.row(i).subvec( j*(n_exog+n_endog), 
                                        j*(n_exog+n_endog) + n_exog - 1 ) ;
      endog.row(j) = X.row(i).subvec( j*(n_exog+n_endog) + n_exog, 
                                        (j+1)*(n_exog+n_endog) - 1 ) ;
    }   // Fill in the endogenous and exogenous matrices
    if( n_cont > 0 )
      cont = X.row(i).tail( n_cont ) ;
        // The controls
    err.row(i) = contemp_eqns_fn( exog, endog, cont, params, extra_args ) ;
  }   // The error on the contemporaneous block
  return err ;
}

