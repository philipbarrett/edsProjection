/***********************************************************************************
 * mono.hpp
 * 
 * Interface to mono.cpp
 * 
 * 08nov2016
 * Philip Barrett, Washington DC
 * 
 ***********************************************************************************/

#ifndef MONO_HPP
#define MONO_HPP

#include <RcppArmadillo.h>
#include <math.h>
#include <string.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;
//using namespace std::string ;

arma::mat M1_nodes_weights_mat( arma::vec mu, arma::mat sigma ) ;

#endif