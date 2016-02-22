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
#include "sim.hpp"
#include "quad.hpp"

// [[Rcpp::export]]
arma::mat euler_hat_grid( 
            arma::mat coeffs, arma::mat coeffs_cont, 
            arma::mat X, int lags, List params, 
            int n_exog, int n_endog, int n_cont,
            arma::rowvec rho, arma::rowvec sig_eps, int n_integ,
            int N, arma::rowvec upper, arma::rowvec lower, bool cheby,
            arma::mat exog_innov_mc, bool quad=true, int n_nodes=0 ){
// Creates the vector of errors on the Euler equations

  int n_pts = X.n_rows ;
      // The number of points at which the error is assessed
  mat exog = zeros( 1 + lags, n_exog ) ;
  mat endog = zeros( 1 + lags, n_endog ) ;
  rowvec cont = zeros<rowvec>( std::max( n_cont, 1 ) ) ;
      // Temporary containers used in the loop.  Make cont bigger than size 0
      // here - just passing a useless empty container
  mat err = zeros(n_pts, 4 ) ;
      // Becuase there are four Euler equations

  /** Create the integration nodes and weights **/
  n_integ = quad ? pow( n_nodes, n_exog ) : n_integ ;
      // Update the number of points if using quadrature
  rowvec weights( n_integ ) ;
  mat nodes( n_integ, n_exog ) ;
  vec v_sig_eps = conv_to<vec>::from( sig_eps ) ;
      // The weights and integration nodes
      
  if( quad ){
    mat m_quad = quad_nodes_weights_mat( n_nodes, n_exog, 
                        v_sig_eps, zeros(n_exog) ) ;
    vec temp = m_quad.col( n_exog ) ;
    weights = conv_to<rowvec>::from( temp ) ;
    nodes = m_quad.head_cols(n_exog) ;
        // Quadrature
  }
  else
  {
    weights = ones<rowvec>( n_pts ) / n_integ ;
    nodes = exog_innov_mc ;
        // Monte Carlo integration
  }
  
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
    err.row(i) = euler_hat_irbc(
                    exog, endog, cont, nodes, params,
                    coeffs_cont, n_exog, n_endog, n_cont, rho, 
                    n_integ, N, upper, lower, cheby, weights, false ) ;
  }   // The error on the states according to the Euler equations
  
  return err ;
}

arma::mat x_eqns_irbc_grid( arma::mat X, int lags, List params,
                              int n_exog, int n_endog, int n_cont ){
// Compute the errors on the intermediates
 
   int n_pts = X.n_rows ;
      // The number of points at which the error is assessed
  mat exog = zeros( 1 + lags, n_exog ) ;
  rowvec cont = zeros<rowvec>( std::max( n_cont, 1 ) ) ;
      // Temporary containers used in the loop.  Make cont bigger than size 0
      // here - just passing a useless empty container
  mat err = zeros(n_pts, 2 ) ;
      // Becuase there are as many equations as endogenous variablesfour Euler
      // equations
  
  /** Now compute the model errors **/
  for( int i = 0 ; i < n_pts ; i++ ){
  // Loop over the evaluation points
    for( int j = 0 ; j < 1 + lags ; j++ ){
    // Loop over the lags
      exog.row(j) = X.row(i).subvec( j*(n_exog+n_endog), 
                                        j*(n_exog+n_endog) + n_exog - 1 ) ;
    }   // Fill in the exogenous matrices
    
    if( n_cont > 0 )
      cont = X.row(i).tail( n_cont ) ;
        // The controls
    err.row(i) = x_eqns_irbc( exog, cont, params ) ;
  }   // The error on the intermediate goods equations
  return err ;
}

arma::mat contemp_eqns_irbc_grid( arma::mat X, int lags, List params,
                                  int n_exog, int n_endog, int n_cont ){
// Compute the errors on the intermediates
 
   int n_pts = X.n_rows ;
      // The number of points at which the error is assessed
  mat exog = zeros( 1 + lags, n_exog ) ;
  mat endog = zeros( 1 + lags, n_endog ) ;
  rowvec cont = zeros<rowvec>( std::max( n_cont, 1 ) ) ;
      // Temporary containers used in the loop.  Make cont bigger than size 0
      // here - just passing a useless empty container
  mat err = zeros(n_pts, 2 ) ;
      // Becuase there are as many equations as endogenous variablesfour Euler
      // equations
  
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
    err.row(i) = contemp_eqns_irbc( exog, endog, cont, params ) ;
  }   // The error on the contemporaneous block
  return err ;
}