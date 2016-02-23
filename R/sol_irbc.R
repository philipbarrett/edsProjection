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
        # The loop variables for c
    
    ### 2.1 Add the controls to the grid ###
    cont.sim <- cont_sim( X, coeff.cont, N, n.endog, n.exog, upper, lower, cheby )
    X.cont <- cbind( X, cont.sim )
        # Compute the new simulation for the controls
    
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
    
    message('  Consumption rule complete.\n    Iterations = ', i.c,  
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

    message('  Solving for the forward-looking rules:' )
    
    i.n <- 1
    n.diff <- 2 * n.tol
    n.gain <- opt$n.gain
    coeff.n <- coeff.n.new <- coeff
    coeff.r <- coeff.r.new <- coeff.cont[ , 3:4 ]
        # The loop variables for the rule for B and r

    while( i.n <= n.iter & n.diff > n.tol ){

      
      ###################################################################
      ## TO DO: Make this part of the solution compute the B and r     ##
      ##        values that solve these equations exactly, and *then*  ##
      ##        approximate.                                           ##
      ## Idea:  Make B act on the cross-country asset holdings?        ##
      ###################################################################
      
      
      ### 3.2 Create the new predictors ###
      k.hat <- euler_hat_grid( coeff.n, coeff.cont, X.cont, lags, params, n.exog, 
                                n.endog, n.cont, rho, sig.eps, 0, N, upper, lower, cheby, 
                                matrix(0,1,1,), TRUE, n.quad )
          # The forward-looking variable predictors
      
      ### 3.3 Compute the error-minimizing coefficients for B and r ###
      for( i in 1:n.endog )
        coeff.n.new[,i] <- coeff_reg( k.hat[,i], X[,state.select], N, lower, upper, cheby )
      for( i in 1:2 )
        coeff.r.new[,i] <- coeff_reg( k.hat[,2+i], X[,state.select], N, lower, upper, cheby )
          # The updated coefficients

      ### 3.4 Damp this part of the loop ###
      n.diff <- max( abs( cbind( coeff.n, coeff.r ) - cbind( coeff.n.new, coeff.r.new ) / 
                            cbind( coeff.n, coeff.r ) ) )
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
      X.cont[, ( 1 + lags ) * ( n.exog + n.endog ) + 1:n.cont ] <- 
            cont_sim( X, coeff.cont, N, n.endog, n.exog, upper, lower, cheby )
          # Again, could make just r for extra speed
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
      plot.coeffs( coeff.n.new, main='New' )
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
