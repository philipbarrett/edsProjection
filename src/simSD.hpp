/***********************************************************************************
 * simSD.hpp
 * 
 * Interface to simSD.cpp
 * 
 * 12aug2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#ifndef SIMSD_HPP
#define SIMSD_HPP

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;
using namespace arma ;

arma::mat simDS( arma::mat u, List params ) ;
    // Returns the simulated DS model given the matrix of shocks
    
// IRF
// Long simulation

#endif