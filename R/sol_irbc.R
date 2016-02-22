#########################################################################
# sol_irbc.R
#
# Solution algorithm for the EDS-projection algorithm for the IRBC model
# Philip Barrett, Chicago
# 22feb2016
#########################################################################

sol.irbc.iterate <- function( coeff.init, opt, params, coeff.cont.init, debug.flag=FALSE ){
  # The main solution iteration loop
  
  # Extract settings
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
  x.gain <- opt$x.gain
  c.gain <- opt$c.gain
  sr <- opt$sr
  adapt.gain <- opt$adapt.gain
  adapt.exp <- opt$adapt.exp
  image <- opt$image
  x.iter <- opt$x.iter
  x.tol <- opt$x.tol
  c.iter <- opt$cons.iter
  c.tol <- opt$cons.tol
  sym.reg <- opt$sym.reg
  l.sym.ave <- opt$l.sym.ave
  l.pairs <- opt$l.pairs
  
  # Extract parameters
  rho <- params$rho
  sig.eps <- params$sig.eps
  n.terms <- idx_count( N, n.exog + n.endog )
  
  # Create the exogenous simulation
  set.seed(1234)
  exog.sim <- sapply( 1:n.exog, function(i) ar1_sim( n.sim * kappa + burn, 
                                                     rho[i], sig.eps[i] ) )

  # Initiate loop variables
  n.state <- n.endog + n.exog
  if( is.null( coeff.cont.init ) ) coeff.cont.init <- 0 * coeff.init
  coeff <- coeff.init
  coeff.cont.new <- coeff.cont <- coeff.cont.init
  i.iter <- 0
  diff <- 2 * diff.tol
  state.select <- 1:(n.exog)
  if( lags > 0 ) for( j in 1:lags ) state.select <- c( state.select, j*n.state + n.exog + 1:n.endog )
      # The states selected for the EDS algorithm
  
  while( i.iter < n.iter & diff > diff.tol ){
    
    if(debug.flag) browser()
    
    coeff.old <- coeff
    # Store the old coefficient
    
    #### 0. APPROXIMATE THE STATE SPACE BY SIMULATION AND REDUCTION ####
    
    message('\n\n\n********************************************************')
    message('Iteration ', i.iter + 1)
    message('  Simulating...')
    
    endog.sim <- endog_sim( n.sim, exog.sim, coeff, N, upper, lower, 
                            endog.init, cheby, kappa, burn, (lags>0) )
        # Simulate the endogenous states
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
    
    #### 1. ITERATE OVER RULES FOR X AND C TO SOLVE CONTEMPORANEOUS EQUATIONS ####
    
    message('  Solving for the consumption rule:' )
    
    # Within iteration: Solve for intermediates, then consumption, updating the
    # grid of these controls at each step
    
    # After iteration: Update the rules for the other controls, update the grid
    # of all controls
    
    #### 2. ITERATE OVER RULES FOR B AND R TO SOLVE EULER EQUATIONS ####
    
    message('  Solving for the state rules:' )
    
    
    
    
    
    message('  Inner loop:' )
    j <- 1
    inner.diff <- 2* inner.tol
    coeff.inner.new <- coeff.inner <- coeff
    # Initialize loop variables
    while( j < inner.iter + 1 & inner.diff > inner.tol ){
      # The inner loop on a fixed set
      message( '    Iteration ', j )
      
      if( j > 1 ){
        endog.new <- cont_sim( X, coeff.inner, N, n.endog, n.exog, upper, lower, cheby )
        X[, n.exog + 1:n.endog] <- endog.new
        # Update the matrix of the model variables on the grid
      }
      
      if( n.cont > 0 ){
        message('       Calculating controls on set...' )
        cont.sim <- cont_sim( X, coeff.cont, N, n.endog, n.exog, upper, lower, cheby )
        X.cont <- cbind( X, cont.sim )
        message('       ...complete' )
      }else{
        X.cont <- X
      }
      
      message('       Computing new inner solution...' )
      k.hat <- euler_hat( coeff.inner, coeff.cont, X.cont, model, lags, params, n.exog, 
                          n.endog, n.cont, rho, sig.eps, 0, N, upper, lower, cheby, 
                          matrix(0,1,1,), TRUE, n.quad )
      
      for( i in 1:n.endog )
        coeff.inner.new[,i] <- coeff_reg( k.hat[,i], X[,state.select], N, lower, upper, cheby )
      if( sym.reg ) 
        coeff.inner.new <- m.sym.ave.pair( coeff.inner.new, l.sym.ave, l.pairs )
      # Symmetry regularization  
      
      
      if( n.cont > 0 ){
        for( i in 1:n.cont )
          coeff.cont.new[,i] <- coeff_reg( k.hat[,n.endog+i], X[, state.select], 
                                           N, lower, upper, cheby )
        if( sym.reg ) 
          coeff.cont.new <- m.sym.ave.pair( coeff.cont.new, l.sym.ave, l.pairs.cont )
        # Symmetry regularization  
      }
      message('       ...complete' )
      
      inner.diff.old <- inner.diff
      inner.diff <- max( abs( coeff.inner.new / coeff.inner - 1 ) )
      # The inner difference
      
      message('       Inner normalized difference = ', round( inner.diff, 5 ) )
      # Screen updating
      coeff.inner <- gain * coeff.inner.new + ( 1 - gain ) * coeff.inner
      coeff.cont <- gain * coeff.cont.new + ( 1 - gain ) * coeff.cont
      j <- j + 1
      # Updating 
      #       if( j > 2 & inner.diff > inner.diff.old ){
      #         j <- inner.iter * 2
      #         message( '     Inner loop divergence detected.  Aborting.' )
      #       }
      
    }
    
    diff <- max( abs( coeff.old - coeff.inner.new ) / abs( coeff.old ) )
    ## The difference to the new estimate
    if( adapt.gain ){
      gain <- max( exp( - adapt.exp * diff ), gain )
      message('  Adaptive gain = ', round( gain, 5 ) )
    } # Adaptive gain
    message('  Maximum normalized difference = ', round( diff, 5 ) , "\n" )
    
    #     coeff <- ( 1 - gain ) * coeff + gain * 
    #       matrix( coeff.inner.new, nrow=n.terms, ncol=n.endog )
    #     coeff.cont <- coeff.cont.new
    #     coeff.cont <- ( 1 - gain ) * coeff.cont + gain * 
    #       matrix( coeff.cont.new, nrow=n.terms, ncol=n.endog )
    
    endog.init <- apply( matrix( X[, n.exog + 1:n.endog ], ncol=n.endog), 2, mean )
    i.iter <- i.iter + 1
    # Housekeeping
    
    if(image){
      par(mfrow=c(2,1))
      x.vals <- ( 1:nrow(coeff) - 1 ) * ( 1 + n.endog ) + .5
      plot.coeffs( coeff.inner.new, main='New' )
      for( i in 1:n.endog ) points( x.vals + i, coeff.inner[,i], col='blue', pch=16 )
      plot.coeffs( coeff.old, main='Old' )
      for( i in 1: n.endog ) points( x.vals + i, coeff.inner[,i], col='blue', pch=16 )
      par(mfrow=c(1,1))
    }
    
    coeff <- coeff.inner
  }
  
  
  # Checks
  #  1. Error on a new set of shocks
  #  2. Bounds
  
  out <- list( coeff=coeff )
  if( n.cont > 0 ) out$coeff.cont <- coeff.cont
  # Set up the output
  
  return( out )
}
