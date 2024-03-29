/***********************************************************************************
 * sim.cpp
 * 
 * Code to run fast stochastic simulations
 * 
 * 25jan2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#include "sim.hpp"
#include "poly.hpp"

// [[Rcpp::export]]
arma::vec ar1_sim( int n_pds, double rho, double sig_eps, 
            bool init_flag=false, double init=0 ){
// Creates an AR(1) simulation
  vec innov = sig_eps * randn( n_pds ) ;
      //  The innovation vector
  vec out = zeros(n_pds) ;
      //  Initialize the output
  if( init_flag ){
    out(0) = init ;
        // Initialize the output
  }
  else
  {
    double sig_uncond = sig_eps / std::sqrt( 1 - std::pow( rho, 2.0 ) ) ;
        // The unconditional variance
    vec v_init = sig_uncond * randn( 1 ) ;
    out(0) = v_init(0) ;
        // Fill the first element of the output with a random draw from the
        // unconditional distribution
  }
  for( int i = 1 ; i < n_pds ; i++ ){
    out(i) = rho * out(i-1) + innov(i) ;
        // Fill the output with the AR(1) process
  }
  return out ;
}
/** TO DO: Add correlated shocks in a VAR **/

// [[Rcpp::export]]
arma::mat var1_sim( int n_pds, arma::mat A, arma::mat sig_eps, 
                   bool init_flag=false, double init=0 ){
// Creates a VAR(1) simulation X(+1) = X * A + e, where e ~ N( 0, sigma )
// where X is a row vector.
  int n_dim = A.n_cols ;
      // Number of dimensions
  mat innov_1 = randn( n_pds, n_dim ) ;
      //  The uncorrelated innovation matrix
  mat omega = chol( sig_eps ) ;
      // The cholesky decomposition of the innovation matrix
  mat innov = innov_1 * omega ;
      // The correlated innovation matrix
  mat out = zeros(n_pds, n_dim) ;
      //  Initialize the output
  if( init_flag ){
    out(0) = init ;
        // Initialize the output
  }
  for( int i = 1 ; i < n_pds ; i++ ){
    out.row(i) = out.row(i-1) * A + innov.row(i) ;
        // Fill the output with the VAR(1) process
  }
  return out ;
}


// [[Rcpp::export]]
arma::rowvec endog_update_slow( arma::rowvec exog, arma::rowvec endog_old, arma::mat coeffs, 
                            int n_exog, int n_endog, int N,
                            arma::rowvec upper, arma::rowvec lower, bool cheby=false ){
// Uses the updating rule defined by coeffs, upper, lower, and cheby to update
// the endogenous states
  mat X = zeros( 1, n_exog + n_endog ) ;
      // The point to evaluate at
  X.row(0).head(n_exog) = exog ;
  X.row(0).tail(n_endog) = endog_old ;
      // Fill in the evaluation point
  vec temp(1) ;
      // Temporary container
  int n_out = coeffs.n_cols ;
      // The number of points to output
  rowvec out = zeros<rowvec>( n_out ) ;
      // Intialize output vector
  for( int i=0; i < n_out ; i++ ){
    temp = poly_eval( coeffs.col(i), X, N, lower, upper, cheby ) ;
    out(i) = temp(0) ;
        // Update each of the endogenous states
  }
  return out ;
}


// [[Rcpp::export]]
arma::rowvec endog_update( arma::rowvec exog, arma::rowvec endog_old, 
                            arma::mat coeffs, int n_exog, int n_endog, int N,
                            arma::rowvec upper, arma::rowvec lower, bool cheby=false ){
// Uses the updating rule defined by coeffs, upper, lower, and cheby to update
// the endogenous states
  
  int K = n_exog + n_endog ;
  int M = 1 ;
      // The number of dimensions and points of X
  mat X_in = zeros( M, K ) ;
      // The point to evaluate at
  X_in.row(0).head(n_exog) = exog ;
  X_in.row(0).tail(n_endog) = endog_old ;
      // Fill in the evaluation point
  int n_out = coeffs.n_cols ;
      // The number of points to output
  mat X = X_limits( X_in, lower, upper, M, K ) ; 
      // Rescale to the unit sphere if required
  umat indices = idx_create( N, K ) ;
      // The matrix of orders
  mat basis = ones( N + 1, K ) ;
  basis = cheby ? cheby_create( X.row(0), N, K ) : ordinary_create( X.row(0), N, K ) ;
      // The matrix of basis polynomials
  int n_terms = indices.n_rows ;
      // The number of terms in the expansion
  
  rowvec out = zeros<rowvec>( n_out ) ;
      // Intialize output vector
      
  for( int i=0; i < n_out ; i++ ){
    out(i) =  poly_eval_core( coeffs.col(i), basis, indices, n_terms, K, M ) ;
        // Update each of the endogenous states
  }
  return out ;
}


// [[Rcpp::export]]
arma::mat endog_sim( int n_out, arma::mat exog_sim, arma::mat coeffs, int N,
                      arma::rowvec upper, arma::rowvec lower, arma::rowvec endog_init, 
                      bool cheby=false, int kappa=1, int burn=0, bool lag=false ){
// Returns a matrix of simulated exogenous and endogenous variables

  int n_exog = exog_sim.n_cols ;
  int n_endog = coeffs.n_cols ;
      // Dimensions of the problem
  rowvec endog_old = endog_init ;
  rowvec endog_new = endog_init ;
      // The endogenous and exogenous states
  int n_pds = kappa * n_out ;
      // Skipping every kappa periods in the record
//  mat out = zeros<mat>( n_out, ( 1 + leads ) * ( n_endog + n_exog ) ) ;
  int n_col = !lag ? ( n_endog + n_exog ) : ( 2 * ( n_endog + n_exog ) ) ;
      // Number of columns
  mat out = zeros<mat>( n_out, n_col ) ;
      // The vector of outputs

  if( burn > 0 ){
    for( int i = 0 ; i < burn ; i ++ ){
      endog_old = endog_new ;
      endog_new = endog_update( exog_sim.row(i), endog_old, coeffs, n_exog, n_endog, 
                                N, upper, lower, cheby ) ;
    }
  }     // Create the inital point via simulation
  
  
  int out_row = 0 ;
      // The counter for the row of the output matrix
  for( int i = 0 ; i < n_pds ; i ++ ){
    endog_old = endog_new ;
    endog_new = endog_update( exog_sim.row(i+burn), endog_old, coeffs, n_exog, n_endog, 
                              N, upper, lower, cheby ) ;
                              
    if( i == out_row * kappa ){
      out.row(out_row).subvec(0,n_exog-1) = exog_sim.row(i+burn) ;
          // Record the exogenous variable
      out.row(out_row).subvec(n_exog,n_endog+n_exog-1) = endog_new ;
          // Record the endogenous variable
      if( lag ){
        if( i == 0 & burn == 0 ){
          out.row(out_row).subvec(n_endog+n_exog,n_endog+2*n_exog-1) = zeros<rowvec>( n_exog ) ;
        }else{
          out.row(out_row).subvec(n_endog+n_exog,n_endog+2*n_exog-1) = exog_sim.row(i+burn-1) ;          
        }
        out.row(out_row).subvec(2*n_exog+n_endog,2*(n_endog+n_exog)-1) = endog_old ;
      }   // Save the lagged variables where required
      out_row ++ ;
          // Increment the output row
    }
  }     // Create and store the rest of the simulation
  return out ;
}

// [[Rcpp::export]]
arma::mat cont_sim( arma::mat xn_sim, arma::mat coeffs_cont, int N,
                      int n_endog, int n_exog,
                      arma::rowvec upper, arma::rowvec lower, bool cheby=false ){
// Creates a simulated path of the control variables from a simulation of states.
// NB: Can actually be used to generate endogenous states as well on a fixed grid.

  int n_out = xn_sim.n_rows ;
  mat out = zeros<mat>( n_out, coeffs_cont.n_cols ) ;
      // The vector of outputs

  for( int i = 0 ; i < n_out ; i ++ ){
    out.row(i) = endog_update( 
          xn_sim.row(i).head(n_exog), 
          xn_sim.row(i).tail(n_endog ),
                // Because the control depends on the **lagged** endogenous
                // state (at least as written) because the endogenous state is
                // an end-of-period variable
          coeffs_cont, n_exog, n_endog, N, upper, lower, cheby ) ;
  }     // Create and store the simulated controls
  return out ;
}

// [[Rcpp::export]]
arma::mat irf_create( int pds, int n_sim, int N, int shk_idx,
                      arma::rowvec rho, arma::rowvec sig_eps, 
                      arma::mat coeffs, arma::rowvec upper, arma::rowvec lower, 
                      arma::rowvec init, int n_endog, int n_exog, 
                      double shk = 0, bool cheby=false ){
// Computes an IRF via simulation

  cube exog_base = zeros( pds, n_exog, n_sim ) ;
  cube exog_alt = zeros( pds, n_exog, n_sim ) ;
  cube endog_base = zeros( pds, n_endog, n_sim ) ;
  cube endog_alt = zeros( pds, n_endog, n_sim ) ;
      // Initialize the different simulation matrices
  rowvec impulse_0 = zeros<rowvec>( n_exog ) ;
  impulse_0( shk_idx ) = shk ;
      // The initial impulse
  if( shk == 0 ){
    impulse_0( shk_idx ) = sig_eps( shk_idx ) / sqrt( 1 - pow( rho( shk_idx ), 2 ) ) ;
      // The impulse is a one standard-deviation shock if not specified
  }
  
  mat rho_decline = ones( pds, n_exog ) ;
  rho_decline.row(0) = zeros<rowvec>( n_exog ) ;
      // The initial period has no shock
  for( int i = 2 ; i < pds ; i++ ){
    rho_decline.row(i) = rho % rho_decline.row(i-1) ;
  }   // The unscaled vector of decreasing shocks for the exogenous variables
  mat impulse = ( ones(pds) * impulse_0 ) % rho_decline ;
      // The matrix of impulses
  
  for( int j = 0 ; j < n_sim ; j++ ){
    for( int i = 0 ; i < n_exog ; i++ ){
      exog_base.slice(j).col(i) = ar1_sim( pds, rho(i), sig_eps(i), true, init(i) ) ;
    }
    exog_alt.slice(j) = impulse + exog_base.slice(j) ;
  } // Create the exogenous simulations
  
  mat temp = zeros( n_exog + n_endog, pds ) ;
      // Temporary storage of the output simulation (endog+exog)
  for( int j = 0 ; j < n_sim ; j++ ){
    for( int i = 0 ; i < n_exog ; i++ ){
      temp = endog_sim( pds, exog_base.slice(j), coeffs, N, upper, lower, 
                        init.tail(n_endog), cheby, 1, 0 ) ;
      endog_base.slice(j) = temp.cols( n_exog, n_exog + n_endog - 1 ) ;
          // The base endogenous simulation
      temp = endog_sim( pds, exog_alt.slice(j), coeffs, N, upper, lower, 
                        init.tail(n_endog), cheby, 1, 0 ) ;
      endog_alt.slice(j) = temp.cols( n_exog, n_exog + n_endog - 1 ) ;
          // The alternative endogenous simulation
    }
  }
  mat out = zeros( pds, n_exog + n_endog ) ;
      // Initialize output matrix
  for( int i = 0 ; i < n_sim ; i++ ){
    out.cols( 0, n_exog - 1 ) = out.cols( 0, n_exog - 1 ) +
          ( exog_alt.slice(i) - exog_base.slice(i) ) / n_sim ;
    out.cols( n_exog, n_exog + n_endog - 1 ) = 
                out.cols( n_exog, n_exog + n_endog - 1 ) +
                  ( endog_alt.slice(i) / endog_base.slice(i) - 1 ) / n_sim ;
        // Create the IRFs
  }
  return out ;
}







