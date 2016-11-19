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
arma::mat simDS( arma::mat u, List mod, double betta, int y1_idx=0 ){
// Returns the simulated DS model given the matrix of shocks.
// Default: Y1 is the first variable.  **Must change if varo[1] != 'Y1 '**
  
  /* 0. unpack parameters */
  int pds = u.n_rows ;
      // Number of simulation periods
  int nshk = u.n_cols ;
      // Number of shocks
  int na = mod["na"] ;
  int nvar = mod["nvar"] ;
  int nsk = mod["nsk"] ;
  int nx = mod["nx"] ;
      // Number of assets, and total, static, and predetermined variables
  vec gamma = mod["gamma"] ;
  vec delta = mod["delta"] ;
      // The responses for asset holdings and excess returns
  rowvec ys = mod["ys"] ;
  rowvec alpha_tilde = mod["alpha.tilde"] ;
      // Steady states
  mat ghu = mod["ghu"] ;
      // The shock impact
  mat ghx = mod["ghx"] ;
      // The shock propagation
    
  /* 1. Create the innovations */
  
  mat BB = trans(ghu.cols(0,nshk-1)) ;
      // Drop the effect of zeta
  mat innovs = u * BB ;
      // The variable innovations
  
  /* 2. Compute the error propagation */
  mat sim = zeros(pds, nvar) ;
      // Initialize the output matrix
  mat AA = trans(ghx) ;
      // The propagation matrix
  sim.row(0) = innovs.row(0) ;
      // The impact of the shock
  for ( int i = 1 ; i < pds ; i++ ){
    sim.row(i) = sim.row(i-1).subvec(nsk,nsk+nx-1) * AA + innovs.row(i) ;
        // The propagated shock
  }

  
  /* 2. Asset shares and excess returns */
  vec shares = sim.cols(nsk, nsk+nx-1) * gamma ;
//  mat rx = sim.cols(nsk, nsk+nx-1) * delta ;
      // Asset shares and excess returns
  double y1ss = ys[y1_idx] ;
      // Steady state Y1
      
  /* 3. Scaling */      
  sim = sim + ones(pds) * ys ;
      // Add on the steady states
  shares = ( shares + ones(pds) * alpha_tilde ) * betta * exp( y1ss );
      // Because the shares need scaling too
  mat out = join_horiz( sim, shares ) ;
      // The output matrix
  
  return out ;
}

// [[Rcpp::export]]
arma::mat stoch_simDS( List mod, arma::mat sigma,
                      double betta, int pds=1e6, int burn=1e4, int y1_idx=0 ){
// Computes a long stochastic simulation for the DS-style model
  
  int nshk = sigma.n_rows ;
  mat u = randn(pds+burn,nshk) * chol(sigma) ;
      // The matrix of shocks
  mat sim = simDS( u, mod, betta, y1_idx) ;
      // Create the simulation
  return sim.rows( burn, pds+burn-1 ) ;
      // Drop the burn periods
}
