/***********************************************************************************
 * poly.hpp
 * 
 * Interface to poly.cpp
 * 
 * 23jan2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/


#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::urowvec idx_increment( arma::urowvec idx, int N, int K ) ;
  // Increments a row vector subject to the constraint that the sum of the
  // coefficients is less than n.
  //   idx: Starting vector
  //   N:   Maximum sum of entries
  //   K:   Vector length
int idx_count( int N, int K ) ;
  // Counts the number of distinct indices
arma::umat idx_create( int N, int K ) ;
  // Creates the matrix of all polynomial powers for an order N approximation of a
  // K-dimensional space
arma::mat ordinary_create( arma::rowvec x, int N, int K ) ;
  // Creates the order N ordinary polynomials for each of the K dimensions of a
  // single vector x (NB: x should be of length K)
arma::mat cheby_create( arma::rowvec x, int N, int K ) ;
  // Creates the N chebychev polynomials for the elemnts of a K-vector