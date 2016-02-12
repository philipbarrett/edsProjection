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
#include "ngm2.hpp"
#include "ngmCont.hpp"
#include "quad.hpp"



// [[Rcpp::export]]
arma::mat euler_hat( arma::mat coeffs, arma::mat coeffs_cont, 
                      arma::mat X, std::string model, int lags, List params, 
                      int n_exog, int n_endog, int n_cont,
                      arma::rowvec rho, arma::rowvec sig_eps, int n_integ,
                      int N, arma::rowvec upper, arma::rowvec lower, bool cheby,
                      arma::mat exog_innov_mc, bool quad=true, int n_nodes=0 ){
// Creates the vactor of "hatted" endogenous states by multiplying through by the Euler scale error

  int n_pts = X.n_rows ;
      // The number of points at which the error is assessed
  mat exog = zeros( 1 + lags, n_exog ) ;
  mat endog = zeros( 1 + lags, n_endog ) ;
  rowvec cont = zeros( n_cont ) ;
      // Temporary containers used in the loop
  mat err = zeros(n_pts, n_endog + n_cont ) ;
      // Becuase there are as many equations as endogenous variables

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

  /** Define the model function **/
  rowvec (*euler_hat_fn) ( 
                  arma::mat exog, arma::mat endog, arma::rowvec cont,
                  arma::mat exog_innov_integ, 
                  List params, arma::mat coeffs, arma::mat coeffs_cont, 
                  int n_exog, int n_endog, int n_cont,
                  arma::rowvec rho, int n_integ, int N, arma::rowvec upper, 
                  arma::rowvec lower, bool cheby, arma::rowvec weights,
                  bool print_rhs ) ;
      // The pointer to model evaluation function
  if( model == "ngm" )
    euler_hat_fn = euler_hat_ngm ;
        // The one-country neoclassical growth model
  if( model == "ngm.cont" )
    euler_hat_fn = ngm_reg ;
        // The one-country neoclassical growth model
  if( model == "ngm2" )
    euler_hat_fn = euler_hat_ngm_2 ;
        // The two-country neoclassical growth model
  
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
    err.row(i) = euler_hat_fn( 
                    exog, endog, cont, nodes, params, coeffs, 
                    coeffs_cont, n_exog, n_endog, n_cont, rho, 
                    n_integ, N, upper, lower, cheby, weights, false ) ;
  }   // The error on the states according to the Euler equations
  
//    Rcout << "n_integ = " << n_integ << std::endl ;
//    Rcout << "weights:\n" << weights << std::endl ;
//    Rcout << "nodes:\n" << nodes << std::endl ;
//    Rcout << "err:\n" << err << std::endl ;
  
  return err ;
}
