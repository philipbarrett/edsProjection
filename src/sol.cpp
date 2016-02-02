/***********************************************************************************
 * sol.cpp
 * 
 * Generic solution algorithm
 * 
 * 27jan2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#include "sol.hpp"
#include "ngm.hpp"
#include "quad.hpp"

// [[Rcpp::export]]
double eval_err( arma::mat coeffs, arma::mat X, std::string model, 
                  int lags, List params, int n_exog, int n_endog,
                  arma::rowvec rho, arma::rowvec sig_eps, int n_integ,
                  int N, arma::rowvec upper, arma::rowvec lower, bool cheby,
                  arma::mat exog_innov_mc, bool quad=true, int n_nodes=0 ){
// Evaluates the error on the model equation at the points in X

  int n_pts = X.n_rows ;
      // The number of points at which the error is assessed
  double err = 0 ;
      // The error
  mat exog = zeros( 1 + lags, n_exog ) ;
  mat endog = zeros( 1 + lags, n_endog ) ;
      // Temporary containers used in the loop

  /** Create the integration nodes and weights **/
  n_integ = quad ? pow( n_nodes, n_exog ) : n_integ ;
      // Update the number of points if using quadrature
  vec weights( n_integ ) ;
  mat nodes( n_integ, n_exog ) ;
      // The weights and integration nodes
  if( quad ){
    mat quad = quad_nodes_weights_mat( pow( n_nodes, n_exog ), n_exog, 
                        sig_eps, zeros(n_exog) ) ;
    weights = quad.col(n_exog) ;
    nodes = quad.head_cols(n_exog) ;
        // Quadrature
  }
  else
  {
    vec weights = ones( n_pts ) / n_integ ;
    mat nodes = exog_innov_mc ;
        // Monte Carlo integration
  }
  
  /** Now compute the model errors **/
  if( model == "ngm" ){
  // The neoclassical growth model
    for( int i = 0 ; i < n_pts ; i++ ){
    // Loop over the evaluation points
      for( int j = 0 ; j < 1 + lags ; j++ ){
      // Loop over the lags
        exog.row(j) = X.row(i).subvec( j*(n_exog+n_endog), 
                                          j*(n_exog+n_endog) + n_exog - 1 ) ;
        endog.row(j) = X.row(i).subvec( j*(n_exog+n_endog) + n_exog, 
                                          (j+1)*(n_exog+n_endog) - 1 ) ;
      }   // Fill in the endogenous and exogenous matrices
      err += err_ngm( exog, endog, nodes, params, coeffs, 
                          n_exog, n_endog, rho, n_integ, N, 
                          upper, lower, cheby, weights, false ) / n_pts ;
    }   // The average absolute relative error
  }
      
  return err ;
}

// [[Rcpp::export]]
arma::vec eval_err_D( arma::mat coeffs, arma::mat X, std::string model, 
                      int lags, List params, int n_exog, int n_endog,
                      arma::rowvec rho, arma::rowvec sig_eps, int n_integ,
                      int N, arma::rowvec upper, arma::rowvec lower, bool cheby,
                      arma::mat exog_innov_mc, bool quad=true, int n_nodes=0 ){
// The error on the model equation at the points in X

  int n_pts = X.n_rows ;
      // The number of points at which the error is assessed
  vec err_D = zeros( coeffs.n_elem ) ;
      // The error
  mat exog = zeros( 1 + lags, n_exog ) ;
  mat endog = zeros( 1 + lags, n_endog ) ;
      // Temporary containers used in the loop

  /** Create the integration nodes and weights **/
  n_integ = quad ? pow( n_nodes, n_exog ) : n_integ ;
      // Update the number of points if using quadrature
  vec weights( n_integ ) ;
  mat nodes( n_integ, n_exog ) ;
      // The weights and integration nodes
  if( quad ){
    mat quad = quad_nodes_weights_mat( pow( n_nodes, n_exog ), n_exog, 
                        sig_eps, zeros(n_exog) ) ;
    weights = quad.col(n_exog) ;
    nodes = quad.head_cols(n_exog) ;
        // Quadrature
  }
  else
  {
    vec weights = ones( n_pts ) / n_integ ;
    mat nodes = exog_innov_mc ;
        // Monte Carlo integration
  }
  
  /** Now compute the model errors **/
  if( model == "ngm" ){
  // The neoclassical growth model
    for( int i = 0 ; i < n_pts ; i++ ){
    // Loop over the evaluation points
      for( int j = 0 ; j < 1 + lags ; j++ ){
      // Loop over the lags
        exog.row(j) = X.row(i).subvec( j*(n_exog+n_endog), 
                                          j*(n_exog+n_endog) + n_exog - 1 ) ;
        endog.row(j) = X.row(i).subvec( j*(n_exog+n_endog) + n_exog, 
                                          (j+1)*(n_exog+n_endog) - 1 ) ;
      }   // Fill in the endogenous and exogenous matrices
      err_D += err_ngm_D( exog, endog, nodes, params, coeffs, 
                          n_exog, n_endog, rho, n_integ, N, 
                          upper, lower, cheby, weights ) / n_pts ;
    }   // The average absolute relative error
  }
      
  return err_D ;
}