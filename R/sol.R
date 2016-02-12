#########################################################################
# sol.R
#
# Generic solution algorithm for the EDS-projection algorithm
# Philip Barrett, Chicago
# 02feb2016
#########################################################################

sol.iterate <- function( coeff.init, coeff.cont.init, opt, params, debug.flag=FALSE ){
# The main solution iteration loop
  
  # Extract settings
  model <- opt$model
  lags <- opt$lags
  n.exog <- opt$n.exog
  n.endog <- opt$n.endog
  n.cont <- opt$n.cont
  N <- opt$N
  cheby <- opt$cheby
  upper <- opt$upper
  lower <- opt$lower
  quad <- opt$quad
  n.quad <- opt$n.quad
  diff.tol <-opt$diff.tol 
  n.iter <-opt$n.iter
  burn <- opt$burn
  kappa <- opt$kappa
  n.sim <- opt$n.sim
  eps <- opt$eps
  delta <- opt$delta
  h <- opt$h
  endog.init <- opt$endog.init
  gain <- opt$gain
  sr <- opt$sr
  adapt.gain <- opt$adapt.gain
  adapt.exp <- opt$adapt.exp
  
  # Extract parameters
  rho <- params$rho
  sig.eps <- params$sig.eps
  
  n.terms <- idx_count( N, n.exog + n.endog )
  
  set.seed(1234)
  exog.sim <- sapply( 1:n.exog, function(i) ar1_sim( n.sim * kappa + burn, 
                                                     rho[i], sig.eps[i] ) )
      # Create the exogenous simulation
  n.state <- n.endog + n.exog
  coeff.new <- coeff <- coeff.init
  coeff.cont.new <- coeff.cont <- coeff.cont.init
  i.iter <- 0
  diff <- 2 * diff.tol
      # Initiate loop variables
  
  state.select <- 1:(n.exog)
  if( lags > 0 ) for( j in 1:lags ) state.select <- c( state.select, j*n.state + n.exog + 1:n.endog )
      # The states selected for the EDS algorithm
  
  while( i.iter < n.iter & diff > diff.tol ){
    
    if(debug.flag) browser()
    
    coeff.old <- coeff
        # Store the old coefficient
    
    message('\n\n\n********************************************************')
    message('Iteration ', i.iter + 1)
    message('  Simulating...')
    
    endog.sim <- endog_sim( n.sim, exog.sim, coeff, N, upper, lower, 
                            endog.init, cheby, kappa, burn, (lags>0) )
        # The simulation
    
    if( max(abs(endog.sim)) > 1e05 )
      stop('Candidate solution explodes.  Try reducing the gain.')
    
    message('  ...complete\n  State reduction...')
    
    if( sr ){
      idces <- p_eps_cheap_const_idx( endog.sim[,state.select], eps, delta )
          # Need to include the lagged state in the evaluation set for the
          # equilibrium condition
      X <- endog.sim[ idces == 1, ]
          # The restricted simulation
    }else{
      X <- endog.sim
    }
    
    message('  ...complete\n  State reduced to ', nrow(X), ' points.')
    
    if( n.cont > 0 ){
      message('Computing controls...' )
      cont.sim <- cont_sim( X, coeff.cont, N, n.endog, n.exog, n.cont, upper, lower, cheby )
      X.cont <- cbind( X, cont.sim )
      message('...complete' )
    }else{
      X.cont <- X
    }
    
    message('Computing new solution...' )
    
    k.hat <- euler_hat( coeff, coeff.cont, X.cont, model, lags, params, n.exog, 
                        n.endog, n.cont, rho, sig.eps, 0, N, upper, lower, cheby, 
                        matrix(0,1,1,), TRUE, n.quad )
    for( i in 1:n.endog )
      coeff.new[,i] <- coeff_reg( k.hat[,i], X[,state.select], N, lower, upper, cheby )
    if( n.cont > 0 ){
      for( i in 1:n.cont )
        coeff.cont.new[,i] <- coeff_reg( k.hat[,n.endog+i], X[,state.select], N, lower, upper, cheby )      
    }
    
    message('...complete' )
    
#     arma::mat coeffs, arma::mat coeffs_cont, 
#     arma::mat X, std::string model, int lags, List params, 
#     int n_exog, int n_endog, int n_cont,
#     arma::rowvec rho, arma::rowvec sig_eps, int n_integ,
#     int N, arma::rowvec upper, arma::rowvec lower, bool cheby,
#     arma::mat exog_innov_mc, bool quad, int n_nodes

    diff <- max( abs( coeff.old - coeff.new ) / abs( coeff.old ) )
        ## The difference to the new estimate
    if( adapt.gain ){
      gain <- max( exp( - adapt.exp * diff ), gain )
      message('  Adaptive gain = ', round( gain, 5 ) )
    } # Adaptive gain
    coeff <- ( 1 - gain ) * coeff + gain * 
      matrix( coeff.new, nrow=n.terms, ncol=n.endog )
    coeff.cont <- ( 1 - gain ) * coeff.cont + gain * 
      matrix( coeff.cont.new, nrow=n.terms, ncol=n.endog )
    
    endog.init <- apply( matrix( X[, n.exog + 1:n.endog ], ncol=n.endog), 2, mean )
        # Housekeeping
    i.iter <- i.iter + 1
    
    par(mfrow=c(2,1))
    x.vals <- ( 1:nrow(coeff) - 1 ) * ( 1 + n.endog ) + .5
    plot.coeffs( coeff.new, main='New' )
    for( i in 1: n.endog ) points( x.vals + i, coeff[,i], col='blue', pch=16 )
    plot.coeffs( coeff.old, main='Old' )
    for( i in 1: n.endog ) points( x.vals + i, coeff[,i], col='blue', pch=16 )
    par(mfrow=c(1,1))
    
    message('  Maximum normalized difference = ', round( diff, 5 ) , "\n" )
  }
  
  
  # Checks
  #  1. Error on a new set of shocks
  #  2. Bounds
  
  return( coeff )
}

sol.check <- function( sol, opt, params ){
# Computes a high-precision error on the Euler equation and checks that the bounds are observed
  
  # Extract settings
  model <- opt$model
  lags <- opt$lags
  n.exog <- opt$n.exog
  n.endog <- opt$n.endog
  N <- opt$N
  cheby <- opt$cheby
  upper <- opt$upper
  lower <- opt$lower
  quad <- opt$quad
  n.quad <- opt$n.quad
  diff.tol <-opt$diff.tol 
  n.iter <-opt$n.iter
  burn <- opt$burn
  kappa <- opt$kappa
  n.sim <- opt$n.sim
  eps <- opt$eps
  delta <- opt$delta
  h <- opt$h
  endog.init <- opt$endog.init
  gain <- opt$gain
  reg.method <- opt$reg.method
  sr <- opt$sr
  adapt.gain <- opt$adapt.gain
  adapt.exp <- opt$adapt.exp
  
  rho <- params$rho
  sig.eps <- params$sig.eps
  
  # Simulation
  set.seed(4321)
  exog.sim <- sapply( 1:n.exog, function(i) ar1_sim( n.sim * kappa + burn, 
                                                     rho[i], sig.eps[i] ) )
  endog.sim <- endog_sim( n.sim, exog.sim, sol, N, upper, lower, 
                          endog.init, cheby, kappa, burn, (lags>0) )
  k.hat <- euler_hat( sol, endog.sim, model, lags, params, n.exog, n.endog, 
                      rho, sig.eps, 0, N, upper, lower, cheby, 
                      matrix(0,1,1,), TRUE, n_nodes=10 )
      # Always increase the number of quadrature nodes to the maximum
  state.select <- n.exog + 1:(n.endog)
  rel.err <- abs( k.hat / endog.sim[, state.select] - 1 ) * 100
  upper.check <- sapply( 1:(n.exog+n.endog), function( i ) all( endog.sim[,i] < upper[i] ) )
  lower.check <- sapply( 1:(n.exog+n.endog), function( i ) all( endog.sim[,i] < upper[i] ) )
  
  return( list( max=apply(rel.err, 2, max), ave=apply(rel.err, 2, mean),
                upper.check=upper.check, lower.check=lower.check ) )
}