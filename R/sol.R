#########################################################################
# sol.R
#
# Generic solution algorithm for the EDS-projection algorithm
# Philip Barrett, Chicago
# 02feb2016
#########################################################################

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
  err.tol <-opt$err.tol 
  n.iter <-opt$n.iter
  burn <- opt$burn
  kappa <- opt$kappa
  n.sim <- opt$n.sim
  eps <- opt$eps
  delta <- opt$delta
  h <- opt$h
  endog.init <- opt$endog.init
  gain <- opt$gain
  
  # Extract parameters
  rho <- params$rho
  sig.eps <- params$sig.eps
  
  n.terms <- idx_count( N, n.exog + n.endog )
  
  set.seed(1234)
  exog.sim <- ar1_sim( n.sim * kappa + burn, rho, sig.eps )
      # Create the exogenous simulation
  n.state <- n.endog + n.exog
  coeff <- coeff.init
  i.iter <- 0
  err <- 2 * err.tol
      # Initiate loop variables
  
#     browser()
  
  while( i.iter < n.iter & err > err.tol ){
    
    if(debug.flag) browser()
    
    message('\n\n\n********************************************************')
    message('Iteration ', i.iter)
    message('  Simulating...')
    
    endog.sim <- endog_sim( n.sim, exog.sim, coeff, N, upper, lower, 
                            endog.init, cheby, kappa, burn, (lags>0) )
        # The simulation
    
    message('  ...complete\n  State reduction...')
    
    state.select <- 1:(n.state)
    if( lags > 0 ) for( j in 1:lags ) state.select <- c( state.select, j*n.state + n.exog + 1:n.endog )
            # The states selected for the EDS algorithm
    idces <- p_eps_cheap_const_idx( endog.sim[,state.select], eps, delta )
        # Need to include the lagged state in the evaluation set for the
        # equilibrium condition
    X <- endog.sim[ idces == 1, ]
        # The restricted simulation

    message('  ...complete\n  State reduced to ', nrow(X), ' points.')
    message('Minimizing error...' )
    
    update <- err.min(  coeff, X, 'ngm', lags, params, n.exog, n.endog, rho, sig.eps, 
                        0, N, upper, lower, cheby, matrix(0,1,1), TRUE, n.quad )
    
    message('...complete' )
    
    err <- update$err
    coeff <- ( 1 - gain ) * coeff + gain * matrix( update$coeff, nrow=n.terms, ncol=n.endog )
    endog.init <- apply( matrix( X[, n.exog + 1:n.endog ], ncol=n.endog), 2, mean )
        # Housekeeping
    i.iter <- i.iter + 1
  }
  
  
  # Checks
  #  1. Error on a new set of shocks
  #  2. Bounds
  
  return( coeff )
}
