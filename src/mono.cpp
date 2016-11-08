/***********************************************************************************
 * mono.cpp
 * 
 * Code to produce monomial integration nodes and weights
 * 
 * 08nov2016
 * Philip Barrett << Chicago
 * 
 ***********************************************************************************/

#include "mono.hpp"

// [[Rcpp::export]]
arma::mat M1_nodes_weights_mat( arma::vec mu, arma::mat sigma ){
// Computes the M1 rule for weights and nodes
  int n_dim = mu.n_elem ;
      // Number of dimensions
  mat out( 2 * n_dim, n_dim + 1 ) ;
      // Initialize the output matrix
  mat R = sqrt(n_dim) * chol(sigma).t() ;
      // Because JMM diefinition of Cholesky decomp of A is LTri B s.t. A=BB',
      // whereas Armadillo returns UTri C s.t. A=C'C.  B=C' as A=A'.
  mat M = ones( n_dim ) * mu.t() ;
      // The matrix of means
  out( span( 0, n_dim - 1 ), span( 0, n_dim - 1 ) ) = M + R ;
  out( span( n_dim, 2 * n_dim - 1), span( 0, n_dim - 1 ) ) = M - R ;
  out.col(n_dim) = ones(2 * n_dim) / ( 2 * n_dim ) ;
      // Fill in the matrix of nodes and weights
  return out ;
}