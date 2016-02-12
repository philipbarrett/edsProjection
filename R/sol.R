#########################################################################
# sol.R
#
# Generic solution algorithm for the EDS-projection algorithm
# Philip Barrett, Chicago
# 02feb2016
#########################################################################

reg.update <- function(){
## Computes the updated guess of the endogeunous states using the regression
## approach
  
}

err.min <- function( coeffs.init, X, model, lags, params, n.exog, 
                     n.endog, rho, sig.eps, n.integ, N, upper, lower, cheby,
                     exog.innov.mc, quad, n.nodes, use.case=FALSE ){
## Minimizes the error on the model equations given a grid X
  
  n.terms <- idx_count( N, n.exog + n.endog )
      # Number of terms in the polynomial approximation
  f <- function(x) eval_err( matrix( x, nrow=n.terms, ncol=n.endog ), 
                             X, model, lags, params, n.exog, n.endog,
                             rho, sig.eps, n.integ, N, upper, lower, cheby,
                             exog.innov.mc, quad, n.nodes )
      # The error evaluated on the restricted grid X
  f.grad <- function(x) eval_err_D( matrix( x, nrow=n.terms, ncol=n.endog ), 
                                    X, model, lags, params, n.exog, n.endog,
                                    rho, sig.eps, n.integ, N, upper, lower, 
                                    cheby,exog.innov.mc, quad, n.nodes )
      # The gradient of the error
  
  if( use.case ){
    return( list( f=f(coeffs.init), f.grad=f.grad(coeffs.init) ) )
  }
  
    ### TEMPORARY HACK
#     g <- function(x) c( x[2] - 1, - ( x[2] - 1 ) )
#     g.grad <- function(x) rbind( c(0,1,0), c(0,-1,0) )
        # Make the coefficinet on k less than one
  
  opts <- list("algorithm"="NLOPT_LD_LBFGS", maxeval=400,
               "xtol_rel"=1.0e-8, print_level=1,
               "check_derivatives" = TRUE,
               "check_derivatives_print" = "all")
      # Set some options
  sol <- nloptr::nloptr( x0=coeffs.init, eval_f = f, eval_grad_f = f.grad, opts=opts )
#   sol <- nloptr::nloptr( x0=coeffs.init, eval_f = f, eval_grad_f = f.grad, 
#                          eval_g_ineq = g, eval_jac_g_ineq = g.grad,
#                          opts=opts )
      # Find the minimum
  
  if( max( abs( f.grad( sol$solution ) ) < 1e-07 ) ){
    message("  Minimum error found, avg relative error = ", 
            round( sqrt( sol$objective ), 4 ) )
    message("                             max gradient = ", 
            max( abs( f.grad( sol$solution ) ) ) )
  }
      # Reporting
  return( list( err=sqrt( sol$objective ), coeff=sol$solution ) )
}


sol.iterate <- function( coeff.init, opt, params, debug.flag=FALSE ){
# The main solution iteration loop
  
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
    message('Computing new solution...' )
    
    if( reg.method ){
      k.hat <- euler_hat( coeff, X, model, lags, params, n.exog, n.endog, 
                          rho, sig.eps, 0, N, upper, lower, cheby, 
                          matrix(0,1,1,), TRUE, n.quad )
      for( i in 1:n.endog )
        coeff.new[,i] <- coeff_reg( k.hat[,i], X[,state.select], N, lower, upper, cheby )
          ## NEED TO CHANGE THE REGRESSION STAGE HERE
    }else{
      update <- err.min(  coeff, X, model, lags, params, n.exog, n.endog, rho, sig.eps, 
                          0, N, upper, lower, cheby, matrix(0,1,1), TRUE, n.quad )
      coeff.new <- update$coeff
    }
    
    message('...complete' )
    

    diff <- max( abs( coeff.old - coeff.new ) / abs( coeff.old ) )
        ## The difference to the new estimate
    if( adapt.gain ){
      gain <- max( exp( - adapt.exp * diff ), gain )
      coeff <- ( 1 - gain ) * coeff + gain * 
        matrix( coeff.new, nrow=n.terms, ncol=n.endog )
      message('  Adaptive gain = ', round( gain, 5 ) )
    }else{
      coeff <- ( 1 - gain ) * coeff + gain * 
        matrix( coeff.new, nrow=n.terms, ncol=n.endog )
    } # Adaptive gain
    
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