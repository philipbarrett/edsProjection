/***********************************************************************************
 * poly.hpp
 * 
 * Interface to poly.cpp
 * 
 * 23jan2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#ifndef POLY_HPP
#define POLY_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

/** TO DO: Convert all the passing by value to by reference **/

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
arma::cube basis_cube( arma::mat X, int N, int K, bool cheby ) ;
  // Computes the cube of polynomial bases for a matrix of multi-dimensional
  // points X
arma::mat X_rescale( arma::mat X_in, int K, bool pc_rescale ) ;
  // Rescales the evaluation points into a sphere in [-1,1]^K
arma::mat X_limits( arma::mat X_in, arma::rowvec lower, 
                    arma::rowvec upper, int M, int K ) ;
  // Rescales the evaluation points into [-1,1]^K based on predefined end-points
arma::vec poly_eval( arma::vec a, arma::mat X_in, int N, 
                      arma::rowvec lower, arma::rowvec upper, bool cheby ) ;
  // Computes the order-N polynomial approximation defined by the vector of
  // coefficients a evaluated at the cloud of points X_in
double poly_eval_core( arma::vec a, arma::mat m_basis, arma::umat indices, 
                          int n_terms, int K, int M) ;
  // The core of the polynomial evaluation.  Assumes scaling and basis creation
  // is already complete.
arma::vec coeff_reg( arma::vec y, arma::mat X_in, int N, 
                       arma::rowvec lower, arma::rowvec upper,
                       bool cheby ) ;
#endif