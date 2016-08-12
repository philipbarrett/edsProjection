/***********************************************************************************
 * simSD.cpp
 * 
 * Code simulate the Devreux-Sutherland style model
 * 
 * 12aug2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#include "simSD.hpp"

// [[Rcpp::export]]
arma::mat simDS( arma::mat u, List mod ){
// Returns the simulated DS model given the matrix of shocks
  
  /* 0. unpack parameters */
  int pds = u.n_rows ;
      // Number of simulation periods
  int nvar = mod["nvar"] ;
  int nsk = mod["nsk"] ;
  int nx = mod["nx"] ;
      // Number of total, static, and predetermined variables
  rowvec gamma = mod["gamma"] ;
  rowvec delta = mod["delta"] ;
      // The responses for asset holdings and excess returns
  mat ghu = mod["ghu"] ;
      // The shock impact
  mat ghx = mod["ghx"] ;
      // The shock propagation
    
  /* 1. Create the innovations */
  mat BB = trans(ghu.cols(0,1)) ;
      // Drop the effect of zeta
  mat innovs = u * BB ;
      // The variable innovations
  
  /* 2. Compute the error propagation */
  mat out = zeros(pds, nvar) ;
      // Initialize the output matrix
  mat AA = trans(ghx) ;
      // The propagation matrix
  out.row(0) = innovs.row(0) ;
      // The impact of the shock
  for ( int i = 1 ; i < pds ; i++ ){
    out.row(i) = out.row(i-1).subvec(nsk,nsk+nx-1) * AA + innovs.row(i) ;
        // The propagated shock
  }
  
  // TO ADD: ASSET SHARES AND EXCESS RETURNS //
  //         AND COLUMN NAMES //
  
  return out ;
}