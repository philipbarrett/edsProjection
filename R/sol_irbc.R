#########################################################################
# sol_irbc.R
#
# Solution algorithm for the EDS-projection algorithm for the IRBC model
# Philip Barrett, Chicago
# 22feb2016
#########################################################################

sol.irbc.iterate <- function( coeff.init, opt, params, coeff.cont.init, debug.flag=FALSE ){
  # The main solution iteration loop
  
  #### Extract settings ####
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
  burn <- opt$burn
  kappa <- opt$kappa
  n.sim <- opt$n.sim
  eps <- opt$eps
  delta <- opt$delta
  h <- opt$h
  endog.init <- opt$endog.init
  sr <- opt$sr
  adapt.gain <- opt$adapt.gain
  adapt.exp <- opt$adapt.exp
  image <- opt$image
  iter <- opt$iter
  tol <- opt$tol
  n.iter <- opt$n.iter
  n.tol <- opt$n.tol
  n.gain <- opt$n.gain
  c.iter <- opt$c.iter
  c.tol <- opt$c.tol
  c.gain <- opt$c.gain
  k.iter <- opt$k.iter
  k.tol <- opt$k.tol
  k.gain <- opt$k.gain
  sym.reg <- opt$sym.reg
  l.sym.ave <- opt$l.sym.ave
  l.pairs <- opt$l.pairs
  l.pairs.cont <- opt$l.pairs.cont
  n.space <- opt$n.space
  
  #### Extract parameters ####
  rho <- params$rho
  sig.eps <- params$sig.eps
  n.terms <- idx_count( N, n.exog + n.endog )
  
  ##### Create the exogenous simulation ####
  set.seed(1234)
  exog.sim <- sapply( 1:n.exog, function(i) ar1_sim( n.sim * kappa + burn, 
                                                     rho[i], sig.eps[i] ) )

  #####  Initiate loop variables   ##### 
  n.state <- n.endog + n.exog
  coeff.old <- coeff <- coeff.init
  coeff.cont.new <- coeff.cont <- coeff.cont.init
  i.iter <- 0
  diff <- 2 * tol
  state.select <- 1:(n.exog)
  if( lags > 0 ) for( j in 1:lags ) state.select <- c( state.select, j*n.state + n.exog + 1:n.endog )
      # The states selected for the EDS algorithm
  
  while( i.iter < iter & diff > tol ){
  # Main outer loop
    
    ######## 1. APPROXIMATE THE STATE SPACE BY SIMULATION AND REDUCTION ########
    
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
    
    
    ######## 2. ITERATE OVER RULES FOR X TO SOLVE CONTEMPORANEOUS EQUATIONS ########
    
    if( debug.flag ) browser()
    
    message('  Solving for the controls:' )
    i.c <- 1
    c.diff <- 2 * c.tol
    c.gain <- opt$c.gain
        # The loop variables for c
    
    ### 2.1 Add the controls to the grid ###
    cont.sim <- cont_sim( X, coeff.cont, N, n.endog, n.exog, upper, lower, cheby )
    X.cont <- cbind( X, cont.sim )
        # Compute the new simulation for the controls
    
        #### ADD TOTAL EQUATION ERRORS HERE ####
    
    while( i.c <= c.iter & c.diff > c.tol ){
      
      ### 2.3 Compute the **values** of the controls satisfying the
      ### contemporaneous block
      cont.hat <- contemp_eqns_irbc_grid( X.cont, lags, params, n.exog, n.endog, n.cont )
          # The predictors in the contemporaneous block
      c.diff <- max( abs( cont.hat - cont.sim ) )
      i.c <- i.c + 1
          # Update the loop controls. Level diff because in logs
      cont.sim <- c.gain * cont.hat + ( 1 - c.gain ) * cont.sim
          # Compute the new coefficients on the consumption rule from the consumption predictors
      X.cont <- cbind( X, cont.sim )
          # Update the coefficients and the grid controls
      if( adapt.gain ) c.gain <- max( exp( - adapt.exp * c.diff ), c.gain )
          # Update the gain where required
    }
    
    message('  Control rules complete.\n    Iterations = ', i.c,  
            '\n    Difference = ', round( c.diff, 5 ),
            '\n    Adaptive gain = ', round( c.gain, 5 ) )
    
    ### 2.4 Update the control rules and grid values ###
    for(i in 1:n.cont) coeff.cont[,i] <- coeff_reg( cont.sim[,i], X[, state.select], 
                                                  N, lower, upper, cheby )
        # Update the coefficients on the other controls in accordance with the
        # solutions for consumption and intermediates.
        # NB: Re-computing the r rules here, but doesn't matter
    if( sym.reg ) coeff.cont <- m.sym.ave.pair( coeff.cont, l.sym.ave, l.pairs.cont )
        # Symmetry regularization
    cont.sim <- cont_sim( X, coeff.cont, N, n.endog, n.exog, upper, lower, cheby )
    X.cont <- cbind( X, cont.sim )
        # Compute the new simulation for the controls
    
    
    
    ######## 3. ITERATE OVER RULES FOR B AND R TO SOLVE EULER EQUATIONS ########
    if( debug.flag ) browser()
    
    message('  Solving for the forward-looking rules:' )
    
    i.n <- 1
    n.diff <- 2 * n.tol
    n.gain <- opt$n.gain
    coeff.n <- coeff.n.new <- coeff
    coeff.r <- coeff.r.new <- coeff.cont[ , 3:4 ]
        # The loop variables for the rule for B and r

    while( i.n <= n.iter & n.diff > n.tol ){

      message('    Iteration ', i.n )
      message('      Solving the Euler equations...' )
      
      k.diff <- 2 * opt$k.tol
      k.gain <- opt$k.gain
      i.k <- 0
          # Intialize the loop variables
      b.r.idces <- c(n.exog+1:2,(1+lags)*(n.endog+n.exog)+3:4)
          # Indices in X.cont
      while( k.diff > k.tol & i.k < k.iter ){
        k.hat <- euler_hat_grid( coeff.n, coeff.cont, X.cont, lags, params, n.exog, 
                                 n.endog, n.cont, rho, sig.eps, 0, N, upper, lower, cheby, 
                                 matrix(0,1,1,), TRUE, n.quad )
            # The forward-looking variable predictors
        v.k.diff <- apply( abs( k.hat - X.cont[, b.r.idces] ), 2, mean )
        k.diff <- max( v.k.diff )
        i.k <- i.k + 1
            # Update loop variables
        X.cont[,b.r.idces] <- k.gain * k.hat + ( 1 - k.gain ) * X.cont[,b.r.idces]
            # Update the simulations for B and r
#         if( adapt.gain ) k.gain <- max( exp( - adapt.exp * k.diff ), k.gain )
            # Update the gain where required
        if( i.k %% 5 == 0 ) 
          message('        i.k = ', i.k, ',  v.k.diff = (', 
                  round(v.k.diff[1],4), ', ', round(v.k.diff[2],4),  ', ',
                  round(v.k.diff[3],4), ', ', round(v.k.diff[4],4), ' )' )
      }
      
      message('      ...complete.\n        Iterations = ', i.k,  
              '\n        Difference = ', round( k.diff, 5 ),
              '\n        Adaptive gain = ', round( k.gain, 5 ) )
      message('      Updating rules...')
      
      
      ### 3.3 Compute the error-minimizing coefficients for B and r ###
      for( i in 1:n.endog )
        coeff.n.new[,i] <- coeff_reg( k.hat[,i], X[,state.select], N, lower, upper, cheby )
      for( i in 1:2 )
        coeff.r.new[,i] <- coeff_reg( k.hat[,2+i], X[,state.select], N, lower, upper, cheby )
          # The updated coefficients

      ### 3.4 Damp this part of the loop ###
      n.diff <- max( abs( cbind( coeff.n, coeff.r ) - cbind( coeff.n.new, coeff.r.new ) ) )
      i.n <- i.n + 1
          # Update the loop controls
      coeff.n <- n.gain * coeff.n.new + ( 1 - n.gain ) * coeff.n
      coeff.r <- n.gain * coeff.r.new + ( 1 - n.gain ) * coeff.r
      coeff.cont[, 3:4 ] <- coeff.r
          # Update the rules
      if( adapt.gain ) n.gain <- max( exp( - adapt.exp * n.diff ), n.gain )
          # Update the adaptive gain
      if( sym.reg ){
        coeff.n <- m.sym.ave.pair( coeff.n, l.sym.ave, l.pairs )
        coeff.cont <- m.sym.ave.pair( coeff.cont, l.sym.ave, l.pairs.cont )
      }   # Symmetry regularization
      
      ### 3.5 Evaluate the rules on the state grid ###
      X.cont[, n.exog + 1:n.endog] <- cont_sim( X, coeff.n, N, n.endog, n.exog, upper, lower, cheby )
      X.cont[, ( 1 + lags ) * ( n.exog + n.endog ) + 3:4 ] <- 
            cont_sim( X, coeff.cont[,3:4], N, n.endog, n.exog, upper, lower, cheby )
      
      message('      ...complete.\n        Difference    = ', round( n.diff, 5 ),
              '\n        Adaptive gain = ', round( n.gain, 5 ) )
    }
    

    ###### 4. MEASURE THE CHANGES IN THE STATE VARIABLE COEFFICIENTS ######
    diff <- max( abs( coeff.old - coeff.n.new ) )
        ## The difference to the new estimate
    message('  Outer maximum normalized difference = ', round( diff, 5 ) , "\n" )
        # The overall change
    eq.err.n <- euler_hat_grid( coeff.n, coeff.cont, X.cont, lags, params, n.exog, 
                                 n.endog, n.cont, rho, sig.eps, 0, N, upper, lower, cheby, 
                                 matrix(0,1,1,), TRUE, n.quad ) - 
                    X.cont[, c( n.exog + 1:n.endog, 
                                ( 1 + lags ) * ( n.exog + n.endog ) + 3:4 ) ]
    eq.err.cont <- contemp_eqns_irbc_grid( X.cont, lags, params, n.exog, n.endog, n.cont ) - 
          X.cont[, ( 1 + lags ) * ( n.exog + n.endog ) + 1:n.cont ]
    max.eq.err <- max( apply( abs( eq.err.n ), 2, mean ), 
                       apply( abs( eq.err.cont ), 2, mean ) )
    message('  Maximum average absolute equation error = ', round(max.eq.err, 4) )
  
    endog.init <- apply( matrix( X[, n.exog + 1:n.endog ], ncol=n.endog), 2, mean )
    i.iter <- i.iter + 1
        # Housekeeping
    
    if(image){
      par(mfrow=c(2,1))
      x.vals <- ( 1:nrow(coeff) - 1 ) * ( 1 + n.endog ) + .5
      plot.coeffs( coeff.n.new, main='New' )
      for( i in 1:n.endog ) points( x.vals + i, coeff.n[,i], col='blue', pch=16 )
      plot.coeffs( coeff.old, main='Old' )
      for( i in 1: n.endog ) points( x.vals + i, coeff.n[,i], col='blue', pch=16 )
      par(mfrow=c(1,1))
    }
        # Charting
    coeff <- coeff.n
    coeff.old <- coeff
        # Update coefficients
  }
  
  out <- list( coeff=coeff )
  if( n.cont > 0 ) out$coeff.cont <- coeff.cont
  out$X.cont <- X.cont
  out$opt <- opt
  out$params <- params
      # Set up the output
  return( out )
}

sol.irbc.check <- function( sol, params=NULL, opt=NULL ){
# Checks the model equation errors on the reduced state space X.cont
  
  if( is.null(params) ) params <- sol$params
  if( is.null(opt) ) opt <- sol$opt
      # Assign the options and params when required
  
  n.check <- 10000
  n.burn <- 1000
  n.quad <- 8
  kappa <- 101
      # The check and burn numbers.  Also do high-precision integration.
  
  n.exog <- opt$n.exog
  n.endog <- opt$n.endog
  n.cont <- opt$n.cont
  upper <- opt$upper
  lower <- opt$lower
  N <- opt$N
  cheby <- opt$cheby
  lags <- opt$lags
      # Copy from options
  
  rho <- params$rho
  sig.eps <- params$sig.eps
      # Copy from parameters
  
  exog.sim <- sapply( 1:n.exog, function(i) ar1_sim( kappa * n.check + n.burn, 
                                                     rho[i], sig.eps[i] ) )
      # The exogenous simulation
  endog.sim <- endog_sim( n.check, exog.sim, sol$coeff, N, upper, lower, 
                          c(0,0), cheby, kappa, n.burn, (lags>0) )
      # The endogenous simulation (Here set kappa=1)
  cont.sim <- cont_sim( endog.sim, sol$coeff.cont, N, n.endog, n.exog, upper, lower, cheby )
      # The controls
  all.sim <- cbind( endog.sim, cont.sim )
      # The combined simulation
  
  k.hat <- euler_hat_grid( sol$coeff, sol$coeff.cont, all.sim, lags, params, n.exog, 
                  n.endog, n.cont, rho, sig.eps, 0, N, upper, lower, cheby, 
                  matrix(0,1,1,), TRUE, n.quad )
      # The predictors for B & r
  cont.hat <- contemp_eqns_irbc_grid( all.sim, lags, params, n.exog, n.endog, n.cont )
      # The predictors for the controls
  cont.hat[, 3:4] <- k.hat[, 3:4]
      # Replace the predictors for r
  
  err <- cbind( k.hat[,1:2] - endog.sim[, 3:4], cont.hat - cont.sim )
      # The errors
  out <- list( endog.sim=endog.sim, cont.sim=cont.sim, err=err )
  
  return( out )
}



