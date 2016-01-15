/***********************************************************************************
 * eds.cpp
 * 
 * Code to compute the epsilon-distinguishable sets as described in 
 * Maliars & Maliars (2015), Quantitative Economics
 * 
 * 14jan2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

// [[Rcpp::export]]
arma::mat pcd_scaling( arma::mat & X ) {
// Reutrns the mean-zero, spherical unit-variance equivalent of a cloud of
// points via PC decomposition

  rowvec X_mean = mean( X, 0 ) ;
      // The vector of column means
  mat X_center = ( X - ones( X.n_rows ) * X_mean ) ;
      // Centering the X matrix on its mean
  mat X_pc = princomp( X_center ) ;
      // The principal components of X
  mat X_pcd = X_center * X_pc ;
      // The principal component projection of X
  rowvec X_sd = stddev( X_pcd, 0, 0 ) ;
      // The vector of standard deviations
  mat X_norm = X_pcd / ( ones( X.n_rows ) * X_sd ) ;
      // Rescale to make spherical unit variance
      
  return X_norm ;

}

// [[Rcpp::export]]
arma::mat p_eps( arma::mat & X, double eps ) {
// The Basic EDS algorithm for fixed distance epsilon
  
  mat X_norm = pcd_scaling( X ) ;
      // The PCD rescaling of X (removes units)
  
  /** 1. Calculate the rows to retain **/
  int n_rows = X.n_rows ;
      // The number of rows of X
  vec keep = ones( n_rows ) ;
      // Initiate the loop indicator
  double dist = 0 ;
      // The Euclidean distance
  for ( int i = 0 ; i < n_rows - 1 ; i++ ){
    // Loop over rows of X_norm
    if( keep(i) == 1 ){
      // If the point is not a keeper, no need to looki for nearest neighbours
      for( int j = i + 1 ; j < n_rows ; j++ ){
        // Loop over other rows
        dist = norm( X_norm.row(i) - X_norm.row(j), 2 ) ;
            // Compute the Euclidean distance
        if( dist < eps ) keep(j) = 0 ;
            // If the distance is small set the keeping indicator to zero
      }
    }
  }
  
  /** 2. Construct a return matrix with those rows **/
  int n_z = sum( keep ) ;
      // The number of points to keep
  mat Z = zeros( n_z, X.n_cols ) ;
      // Initialize the return matrix
  int i_z = 0 ;
      // Counter for rows of z
  for( int i= 0 ; i < n_rows ; i++ ){
    // Loops over rows of X
    if( keep(i) == 1 ){
      // If i is a keeper, then keep it!
      Z.row( i_z ) = X.row( i ) ;
          // Add to Z
      i_z++ ;
          // Increment Z
    }
  }
  
  return Z ;
}
