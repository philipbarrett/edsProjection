/***********************************************************************************
 * quad.hpp
 * 
 * Interface to quad.cpp
 * 
 * 28jan2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#ifndef QUAD_HPP
#define QUAD_HPP

#include <RcppArmadillo.h>
#include <math.h>
#include <string.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;
//using namespace std::string ;

arma::vec quad_nodes_1d( int nodes ) ;
arma::vec quad_weights_1d( int nodes ) ;
List quad_nodes_weights( int n_nodes, int n_dim, 
                          arma::vec sig, arma::vec mu ) ;
arma::mat quad_nodes_weights_mat( int n_nodes, int n_dim, 
                          arma::vec sig, arma::vec mu ) ;

#endif