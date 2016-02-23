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
  x.gain <- opt$x.gain
  c.gain <- opt$c.gain
  n.gain <- opt$n.gain
  sr <- opt$sr
  adapt.gain <- opt$adapt.gain
  adapt.exp <- opt$adapt.exp
  image <- opt$image
  iter <- opt$iter
  tol <- opt$tol
  n.iter <- opt$n.iter
  n.tol <- opt$n.tol
  c.iter <- opt$c.iter
  c.tol <- opt$c.tol
  x.iter <- opt$x.iter
  x.tol <- opt$x.tol
  sym.reg <- opt$sym.reg
  l.sym.ave <- opt$l.sym.ave
  l.pairs <- opt$l.pairs
  l.pairs.cont <- opt$l.pairs.cont
  
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
  
  while( i.iter < n.iter & diff > tol ){
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
    
    
    
    ########2. ITERATE OVER RULES FOR X AND C TO SOLVE CONTEMPORANEOUS EQUATIONS ########
    
    message('  Solving for the consumption rule:' )
    i.c <- 1
    c.diff <- 2 * c.tol
    c.gain <- opt$c.gain
    coeff.c.new <- coeff.cont[ , 1:2 ]
        # The loop variables for c
    
    ### 2.1 Add the controls to the grid ###
    cont.sim <- cont_sim( X, coeff.cont, N, n.endog, n.exog, upper, lower, cheby )
    X.cont <- cbind( X, cont.sim )
        # Compute the new simulation for the controls
    
    while( i.c <= c.iter & c.diff > c.tol ){
      
      ### 2.2 Iterate to solve for intermediates given consumption ###
      i.x <- 1
      x.diff <- 2 * x.tol
      x.gain <- opt$x.gain
      coeff.x.new <- coeff.cont[ , 5:6 ]
          # The loop variables for x
      
      message('    Intermediates ... ')
      
      if(debug.flag) browser()
      
      ## Solve for the **values** of x_11 and x_22 given everything else in the simulation ##
      while( i.x <= x.iter & x.diff > x.tol ){
        x.hat <- x_eqns_irbc_grid( X.cont, lags, params, n.exog, n.endog, n.cont )
            # The predictors for the x values
        cont.sim[ , 5:6 ] <- x.gain * x.hat + ( 1 - x.gain ) * cont.sim[, 5:6 ]
            # Update the simulation on the grid
        x.diff <- max( abs( x.hat - cont.sim[, 5:6] ) )
        i.x <- i.x + 1
            # Update the loop controls
        X.cont <- cbind( X, cont.sim )
            # Update the coeffients and the simulation on the grid
        if( adapt.gain ) x.gain <- max( exp( - adapt.exp * x.diff ), x.gain )
            # Update the gain where required
      }
      
      message('    ... complete.\n      Iterations = ', i.x,  
              '\n      Difference = ', round( x.diff, 5 ),
              '\n      Adaptive gain = ', round( x.gain, 5 ) )
      
      ### 2.3 Compute the consumption errors and update the rule
      cont.hat <- contemp_eqns_irbc_grid( X.cont, lags, params, n.exog, n.endog, n.cont )
          # The predictors in the contemporaneous block
      for(i in 1:2) coeff.c.new[,i] <- coeff_reg( cont.hat[,i], X[, state.select], 
                                                  N, lower, upper, cheby )
          # Compute the new coefficients on the consumption rule from the consumption predictors
      c.diff <- max( abs( coeff.c.new / coeff.cont[, 1:2] - 1 ) )
      i.c <- i.c + 1
          # Update the loop controls
      coeff.cont[ 1:2, ] <- x.gain * coeff.c.new + ( 1 - x.gain ) * coeff.cont[ , 1:2 ]
      cont.sim[ 1:2, ] <- cont_sim( X, coeff.cont[ , 1:2 ], N, n.endog, n.exog, upper, lower, cheby )
      X.cont <- cbind( X, cont.sim )
          # Update the coefficients and the grid controls
      if( adapt.gain ) c.gain <- max( exp( - adapt.exp * c.diff ), c.gain )
          # Update the gain where required
    }
    
    message('  Consumption rule complete.\n    Iterations = ', i.c,  
            '\n    Difference = ', round( c.diff, 5 ),
            '\n    Adaptive gain = ', round( c.gain, 5 ) )
    
    ### 2.4 Update the control rules and grid values ###
    for(i in 1:7) coeff.cont[,6+i] <- coeff_reg( cont.hat[,2+i], X[, state.select], 
                                                  N, lower, upper, cheby )
        # Update the coefficients on the other controls in accordance with the
        # solutions for consumption and intermediates
    if( sym.reg ) coeff.cont <- m.sym.ave.pair( coeff.cont, l.sym.ave, l.pairs.cont )
        # Symmetry regularization
    cont.sim <- cont_sim( X, coeff.cont, N, n.endog, n.exog, upper, lower, cheby )
    X.cont <- cbind( X, cont.sim )
        # Compute the new simulation for the controls
    
    
    
    
    ######## 3. ITERATE OVER RULES FOR B AND R TO SOLVE EULER EQUATIONS ########

    message('  Solving for the forward-looking rules:' )
    
    
    ### 3.1 Compute the endogenous variables on X ###
    X[, n.exog + 1:n.endog] <- cont_sim( X, coeff, N, n.endog, n.exog, upper, lower, cheby )
        # Update the endogenous states on the grid points (the controls will not change)
    
    i.n <- 1
    n.diff <- 2 * n.tol
    n.gain <- opt$n.gain
    coeff.n <- coeff.n.new <- coeff
    coeff.r <- coeff.r.new <- coeff.cont[ , 3:4 ]
        # The loop variables for the rule for B and r
    
    
    while( i.n <= n.iter & n.diff > n.tol ){

      ### 3.2 Create the new predictors ###
      k.hat <- euler_hat_grid( coeff.n, coeff.cont, X.cont, lags, params, n.exog, 
                                n.endog, n.cont, rho, sig.eps, 0, N, upper, lower, cheby, 
                                matrix(0,1,1,), TRUE, n.quad )
          # The forward-looking variable predictors
      
      ### 3.3 Compute the error-minimizing coefficients for B and r ###
      for( i in 1:n.endog )
        coeff.n.new[,i] <- coeff_reg( k.hat[,i], X[,state.select], N, lower, upper, cheby )
      for( i in 1:2 )
        coeff.r.new[,2+i] <- coeff_reg( k.hat[,2+i], X[,state.select], N, lower, upper, cheby )
          # The updated coefficients

      ### 3.4 Damp this part of the loop ###
      n.diff <- max( abs( cbind( coeff.n, coeff.r ) - cbind( coeff.n.new, coeff.r.new ) / 
                            cbind( coeff.n, coeff.r ) ) )
      i.n <- i.n + 1
          # Update the loop controls
      coeff.n <- n.gain * coeff.n.new + ( 1 - gain ) * coeff.n
      coeff.r <- n.gain * coeff.r.new + ( 1 - gain ) * coeff.r
      coeff.cont[, 3:4 ] <- coeff.r
          # Update the rules
      if( adapt.gain ) n.gain <- max( exp( - adapt.exp * n.diff ), n.gain )
          # Update the adaptive gain
      if( sym.reg ){
        coeff.n <- m.sym.ave.pair( coeff.n, l.sym.ave, l.pairs )
        coeff.cont <- m.sym.ave.pair( coeff.cont, l.sym.ave, l.pairs.cont )
      }   # Symmetry regularization
      
      ### 3.5 Evaluate the rules on the state grid ###
      X[, n.exog + 1:n.endog] <- cont_sim( X, coeff.n, N, n.endog, n.exog, upper, lower, cheby )
      X[, ( 1 + lags ) * ( n.exog + n.endog ) ] <- 
            cont_sim( X, coeff.cont, N, n.endog, n.exog, upper, lower, cheby )
      
    }
    
    message('  B and r rules complete.\n    Iterations = ', i.n,  
            '\n    Difference    = ', round( n.diff, 5 ),
            '\n    Adaptive gain = ', round( n.gain, 5 ) )

    
    
    ###### 4. MEASURE THE CHANGES IN THE STATE VARIABLE COEFFICIENTS ######
    diff <- max( abs( coeff.old - coeff.n.new ) / abs( coeff.old ) )
        ## The difference to the new estimate
    message('  Outer maximum normalized difference = ', round( diff, 5 ) , "\n" )
        # The overall change
    
    endog.init <- apply( matrix( X[, n.exog + 1:n.endog ], ncol=n.endog), 2, mean )
    i.iter <- i.iter + 1
        # Housekeeping
    
    if(image){
      par(mfrow=c(2,1))
      x.vals <- ( 1:nrow(coeff) - 1 ) * ( 1 + n.endog ) + .5
      plot.coeffs( coeff.new, main='New' )
      for( i in 1:n.endog ) points( x.vals + i, coeff[,i], col='blue', pch=16 )
      plot.coeffs( coeff.old, main='Old' )
      for( i in 1: n.endog ) points( x.vals + i, coeff[,i], col='blue', pch=16 )
      par(mfrow=c(1,1))
    }
        # Charting
    coeff <- coeff.n
    coeff.old <- coeff
        # Update coefficients
  }
  
  out <- list( coeff=coeff )
  if( n.cont > 0 ) out$coeff.cont <- coeff.cont
      # Set up the output
  return( out )
}
