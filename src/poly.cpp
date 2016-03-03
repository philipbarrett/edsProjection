/***********************************************************************************
 * poly.cpp
 * 
 * Code to calculate fast polynomial approximations
 * 
 * 23jan2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/


#include "poly.hpp"
#include "eds.hpp"

// [[Rcpp::export]]
arma::urowvec idx_increment( arma::urowvec idx, int N, int K ){
// Increments a row vector subject to the constraint that the sum of the
// coefficients is less than n.
//   idx: Starting vector
//   N:   Maximum sum of entries
//   K:   Vector length
  
  if( idx(0) == N ){
  // If we already are at the last case, return to the start
    return zeros<urowvec>( K ) ;
  }
  
  if( sum( idx ) < N ){
      // The "usual" case, where only the last element needs incrementing
    urowvec out = idx ;
    out( K - 1 ) ++ ;
    return out ;
        // Increment the final index
  }
  
  bool complete = false ;
      // I
  int k = K - 1 ;
      // The counter for each vector element
  urowvec idx_new = idx ;
      // Initialize the updated vector
  while( !complete & k > 0 ){
    idx_new( k ) = 0 ;
    idx_new( k - 1 ) ++ ;
        // Set the current elemt to zero and increment the preceding one
    complete = ( sum( idx_new ) <= N ) ;
        // Stop if the sum of terms is less than the order
    k -- ;
        // Decrement counter
  }
      
  return idx_new ;
}

// [[Rcpp::export]]
int idx_count( int N, int K ){
// Counts the number of distinct indices 
  urowvec init = zeros<urowvec>( K ) ;
      // The initial indexation vector
  int num_idx = 1 ;
  urowvec inc = init ;
      // Initialize loop variables
  while( inc(0) < N ){
    inc = idx_increment( inc, N, K ) ;
    num_idx ++ ;
        // Increment the index and the counter
  }
  return num_idx ;
}

// [[Rcpp::export]]
arma::umat idx_create( int N, int K ){
// Creates the matrix of all polynomial powers for an order N approximation of a
// K-dimensional space
  int n_idx = idx_count( N, K ) ;
      // Numer of rows
  umat out = zeros<umat>( n_idx, K ) ;
      // Initialize the output matrix
  for( int i = 1 ; i < n_idx ; i++ ){
    out.row(i) = idx_increment( out.row(i-1), N, K ) ;
        // Add the new row
  }
  return out ;
}

// [[Rcpp::export]]
arma::mat ordinary_create( arma::rowvec x, int N, int K ){
// Creates the order N ordinary polynomials for each of the K dimensions of a
// single vector x (NB: x should be of length K)
  mat out = ones( N + 1, K ) ;
      // Initialize output matrix
  for( int n = 1 ; n < N + 1 ; n++ ){
    out.row(n) = out.row(n-1) % x ;
  }
  return out ;
}

// [[Rcpp::export]]
arma::mat cheby_create( arma::rowvec x, int N, int K ){
// Creates the N chebychev polynomials for the elemnts of a K-vector
  mat out = ones( N + 1, K ) ;
      // Initialize output matrix
  out.row(1) = x ;
      // The linear term
  if( N == 1 ) return out ;
      // If N=2 we are done
  for( int n = 2 ; n < N + 1 ; n++ ){
    out.row(n) = 2 * out.row(n-1) % x - out.row( n - 2) ;
  }
  return out ;
}

// [[Rcpp::export]]
arma::cube basis_cube( arma::mat X, int N, int K, bool cheby=false ){
// Computes the cube of polynomial bases for a matrix of multi-dimensional
// points X
  int M = X.n_rows ;
  cube out = ones( N + 1, K, M ) ;
      // Initialize the output
  for( int m = 0 ; m < M ; m ++ ){
    if( cheby ){
      out.slice(m) = cheby_create( X.row(m), N, K ) ; 
    }else{
      out.slice(m) = ordinary_create( X.row(m), N, K ) ; 
    }
  }
  return out ;
}

// [[Rcpp::export]]
arma::mat X_rescale( arma::mat X_in, int K, bool pc_rescale=false ){
// Rescales the evaluation points into a sphere in [-1,1]^K

  /** NEED TO RESCALE ON A FIXED INTERVAL **/

  mat X_sphere = X_in ;
  
  /** PC-scaling: Do not use yet **/
  pc_rescale = false ;
  if( pc_rescale ) X_sphere = pcd_scaling( X_in ) ;
      // Create unit-variance sphere of points
  rowvec X_max = max( X_sphere ) ;
  rowvec X_min = min( X_sphere ) ;
      // The max and min of X_sphere
  rowvec scale = ones<rowvec>( K ) ;
      // Initialize the vector of scale factors
  for( int k = 0 ; k < K ; k++ ){
    scale(k) = std::max( X_max(k), - X_min(k) ) ;
        // The scaling factor is the maximum element-wise distance from the origin
  }
  mat X_out = X_sphere / ( ones( X_sphere.n_rows ) * scale ) ;
      // The output matrix.
  return X_out ;
}

// [[Rcpp::export]]
arma::mat X_limits( arma::mat X_in, arma::rowvec lower, arma::rowvec upper, 
                      int M, int K ){
// Rescales the evaluation points into [-1,1]^K based on predefined end-points

  mat X = 2 * ( X_in - ones(M) * lower ) / ( ones( M ) * ( upper - lower ) ) - 
              ones( M, K ) ;
      // The output matrix.
  return X ;
}

// [[Rcpp::export]]
arma::vec poly_eval( arma::vec a, arma::mat X_in, int N, 
                      arma::rowvec lower, arma::rowvec upper,
                      bool cheby=false ){
// Computes the order-N polynomial approximation defined by the vector of
// coefficients a evaluated at the cloud of points X_in
  
  int K = X_in.n_cols ;
  int M = X_in.n_rows ;
      // The number of dimensions and points of X
//  mat X = X_rescale( X_in, K, pc_rescale ) ;
  mat X = X_limits( X_in, lower, upper, M, K ) ; 
      // Rescale to the unit sphere if required
  cube basis = basis_cube( X, N, K, cheby ) ;
      // The cube of basis polynomials
  umat indices = idx_create( N, K ) ;
      // The matrix of orders
  int n_terms = indices.n_rows ;
      // The number of terms in the expansion
  vec out = zeros(M) ;
      // Initialize the output matrix
  vec v_basis = zeros(K) ;
      // Holding vector for basis terms in each dimension
  for( int m = 0 ; m < M ; m++ ){
      // Loop over points
//    for( int i = 0 ; i < n_terms ; i++){
//        // Loop over terms in the expansion
//      for( int k = 0 ; k < K ; k++ ){
//          // Loop over dimensions
//        v_basis(k) = basis( indices( i, k ), k, m ) ;
//            // Extract the appropriate combination  of basis elements
//      }
//      out(m) = out(m) + a(i) * prod( v_basis ) ;
//    }
    out(m) = poly_eval_core( a, basis.slice(m), indices, n_terms, K, M ) ;
  }
  
  return out ;
}

// [[Rcpp::export]]
double poly_eval_core( arma::vec a, arma::mat m_basis, 
                          arma::umat indices, 
                          int n_terms, int K, int M){
// The central part of polynomial evaluation at just one point.  
// Needs as inputs:
//   - a: coefficients
//   - basis: rescaled basis values
//   - indices: the coefficient indices
//   - n_terms: the number of terms in the approximation
//   - K: the number of dimensions
//   - M: the number of points at which to evaluate
  double out = 0 ;
  
  vec v_basis = zeros(K) ;
    // Holding vector for basis terms in each dimension
  for( int i = 0 ; i < n_terms ; i++){
      // Loop over terms in the expansion
    for( int k = 0 ; k < K ; k++ ){
        // Loop over dimensions
      v_basis(k) = m_basis( indices( i, k ), k ) ;
          // Extract the appropriate combination  of basis elements
    }
    out = out + a(i) * prod( v_basis ) ;
  }
  return out ;
}

// [[Rcpp::export]]
arma::vec coeff_reg( arma::vec y, arma::mat X_in, int N, 
                      arma::rowvec lower, arma::rowvec upper,
                      bool cheby=false ){
// Computes the coefficients of the polynomial approximating y on the data X

  int K = X_in.n_cols ;
  int M = X_in.n_rows ;
      // The number of dimensions and points of X
//  mat X = X_rescale( X_in, K, pc_rescale ) ;
  mat X = X_limits( X_in, lower, upper, M, K ) ; 
      // Rescale to the unit sphere if required
      
//      Rcout << "X.row(0) = " << X.row(0) << std::endl ;
//      Rcout << "N = " << N << std::endl ;
//      Rcout << "K = " << K << std::endl ;
//      Rcout << "cheby = " << cheby << std::endl ;
      
  cube basis = basis_cube( X, N, K, cheby ) ;
      // The cube of basis polynomials
  umat indices = idx_create( N, K ) ;
      // The matrix of orders
  int n_terms = indices.n_rows ;
      // The number of terms in the expansion
  vec v_basis = zeros(K) ;
      // Holding vector for basis terms in each dimension
  mat X_reg = zeros( M, n_terms ) ;
      // The dependent variables for the regression
  for( int m = 0 ; m < M ; m++ ){
      // Loop over points
    for( int i = 0 ; i < n_terms ; i++){
        // Loop over terms in the expansion
      for( int k = 0 ; k < K ; k++ ){
          // Loop over dimensions
        v_basis(k) = basis( indices( i, k ), k, m ) ;
            // Extract the appropriate combination  of basis elements
      }
      X_reg(m,i) = prod( v_basis ) ;
    }
  }
  
//    Rcout << "X_in:\n" << X_in << std::endl ;
//    Rcout << "X:\n" << X << std::endl ;
//    Rcout << "basis:\n" << basis << std::endl ;
//    Rcout << "indices:\n" << indices << std::endl ;
//    Rcout << "X_reg:\n" << X_reg << std::endl ;
  
  vec coeff = solve( X_reg, y ) ;
  return coeff ;
      // The regression coefficinets
}


