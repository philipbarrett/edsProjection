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
arma::vec euler_hat( arma::mat coeffs, arma::mat X, std::string model, 
                      int lags, List params, int n_exog, int n_endog,
                      arma::rowvec rho, arma::rowvec sig_eps, int n_integ,
                      int N, arma::rowvec upper, arma::rowvec lower, bool cheby,
                      arma::mat exog_innov_mc, bool quad=true, int n_nodes=0 ){
// Creates the vactor of "hatted" endogenous states by multiplying through by the Euler scale error

  int n_pts = X.n_rows ;
      // The number of points at which the error is assessed
  mat exog = zeros( 1 + lags, n_exog ) ;
  mat endog = zeros( 1 + lags, n_endog ) ;
      // Temporary containers used in the loop
  mat err = zeros(n_pts, n_endog) ;
      // Becuase there are as many equations as endogenous variables

  /** Create the integration nodes and weights **/
  n_integ = quad ? pow( n_nodes, n_exog ) : n_integ ;
      // Update the number of points if using quadrature
  vec weights( n_integ ) ;
  mat nodes( n_integ, n_exog ) ;
      // The weights and integration nodes
  if( quad ){
    mat m_quad = quad_nodes_weights_mat( pow( n_nodes, n_exog ), n_exog, 
                        sig_eps, zeros(n_exog) ) ;
    weights = m_quad.col(n_exog) ;
    nodes = m_quad.head_cols(n_exog) ;
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
      
      err.row(i) = euler_hat_ngm( exog, endog, nodes, params, coeffs, 
                              n_exog, n_endog, rho, n_integ, N, 
                              upper, lower, cheby, weights, false ) ; // / n_pts ;
    }   // The average absolute relative error
  }
  
//    Rcout << "n_integ = " << n_integ << std::endl ;
//    Rcout << "weights:\n" << weights << std::endl ;
//    Rcout << "nodes:\n" << nodes << std::endl ;
//    Rcout << "err:\n" << err << std::endl ;
  
  return err ;
}
