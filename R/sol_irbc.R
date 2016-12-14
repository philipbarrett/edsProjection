#########################################################################
# sol_irbc.R
#
# Solution algorithm for the EDS-projection algorithm for the IRBC model
# Philip Barrett, Chicago
# 22feb2016
#########################################################################

sol.irbc.iterate <- function( coeff.init, opt, params, coeff.cont.init, 
                              debug.flag=FALSE, c.first=TRUE ){
  # The main solution iteration loop
  
  #### Extract settings ####
  lags <- opt$lags
  n.exog <- opt$n.exog
  n.endog <- opt$n.endog
  n.cont <- opt$n.cont
  n.fwd <- opt$n.fwd
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
  exog.names <- opt$exog.names
  endog.names <- opt$endog.names
  cont.names <- opt$cont.names
  fwd.vars <- opt$fwd.vars
  model <- opt$model
  mono <- if( is.null(opt$mono) ) "none" else opt$mono
  
  #### Extract parameters ####
  rho <- params$rho
  sig.eps <- params$sig.eps
  params$alphahat <- params$share ^ params$eta
  n.terms <- idx_count( N, n.exog + n.endog )
  
  if( debug.flag ) browser()
  
  ##### Create the exogenous simulation ####
  set.seed(1234)
  # exog.sim <- sapply( 1:n.exog, function(i) ar1_sim( n.sim * kappa + burn, 
  #                                                    rho[i], sig.eps[i] ) )
  exog.sim <- var1_sim( kappa * n.sim + burn, rho, sig.eps )

  #####  Initiate loop variables   ##### 
  n.state <- n.endog + n.exog
  coeff.old <- coeff <- coeff.init
  coeff.cont.new <- coeff.cont <- coeff.cont.init
  i.iter <- 0
  diff <- 2 * tol
  state.select <- 1:(n.exog)
  if( lags > 0 ) for( j in 1:lags ) state.select <- c( state.select, j*n.state + n.exog + 1:n.endog )
      # The states selected for the EDS algorithm
  contemp.select <- c( n.exog+1:n.endog, (1+lags)*(n.endog+n.exog)+1:n.cont )
      # The indices of the contemporaenous variables
  extra.args <- list( n.fwd=opt$n.fwd, y1.ss=opt$ys['Y1'] )
  
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
    
    if( c.first || i.iter > 0 ){
      while( i.c <= c.iter & c.diff > c.tol ){
        
        ### 2.3 Compute the **values** of the controls satisfying the
        ### contemporaneous block
        contemp <- contemp_eqns_irbc_grid( X.cont, lags, params, n.exog, n.endog, 
                                            n.cont, extra.args, opt$model )
            # The predictors in the contemporaneous block
        c.diff <- max( abs( contemp - X.cont[,contemp.select] ) )
        i.c <- i.c + 1
            # Update the loop controls. Level diff because in logs
        X.cont[,contemp.select] <- 
                        c.gain * contemp + ( 1 - c.gain ) * X.cont[,contemp.select]
            # Compute the new coefficients on the consumption rule from the consumption predictors
        if( adapt.gain ) c.gain <- max( exp( - adapt.exp * c.diff ), c.gain )
            # Update the gain where required
      }
      
      message('  Contemporaneous block complete.\n    Iterations = ', i.c,  
              '\n    Difference = ', round( c.diff, 5 ),
              '\n    c.gain = ', round( c.gain, 5 ) )
      
      
      ### 2.4 Update the comtemporaneous rules and grid values ###
      for(i in 1:n.endog){
        if( !( endog.names[i] %in% fwd.vars ) ){
          # coeff[,i] <- (1-n.gain) * coeff[,i] + n.gain *
          #                   coeff_reg( contemp[,i], X[, state.select],
          #                                       N, lower, upper, cheby )
          coeff[,i] <- coeff_reg( contemp[,i], X[, state.select],
                                  N, lower, upper, cheby )
        }
      }   # Update the endogenous state coefficients which are not forward-looking
      for(i in 1:n.cont){
        if( !( cont.names[i] %in% fwd.vars ) ){
          # coeff.cont[,i] <- (1-n.gain) * coeff.cont[,i] + n.gain *
          #                     coeff_reg( contemp[,n.endog+i], X[, state.select],
          #                                           N, lower, upper, cheby )
          coeff.cont[,i] <- coeff_reg( contemp[,n.endog+i], X[, state.select],
                                          N, lower, upper, cheby )
        }
      }
          # Update the coefficients on the non-forward-looking variables in line
          # with the predictors from the conemporanneous block
      # TEST: Enforcing 100% update here
      
      if( sym.reg ) coeff.cont <- m.sym.ave.pair( coeff.cont, l.sym.ave, l.pairs.cont )
          # Symmetry regularization
    }
    sim <- cont_sim( X, cbind( coeff, coeff.cont ), N, n.endog, n.exog, upper, lower, cheby )
        # The new simulated values
    X.cont <- cbind( X[,1:n.exog], sim[,1:n.endog], 
                     X[,n.exog+1:(n.exog+n.endog)], sim[,n.endog+1:n.cont] )
        # Save the new simulation for the current-period values
    
    
    ######## 3. ITERATE OVER RULES FOR B AND R TO SOLVE EULER EQUATIONS ########
    if( debug.flag ) browser()
    
    message('  Solving for the forward-looking rules:' )
    
    i.n <- 1
    n.diff <- 2 * n.tol
    n.gain <- opt$n.gain
    coeff.n <- coeff.n.new <- coeff
    coeff.c <- coeff.c.new <- coeff.cont
        # The loop variables for the endogenous and control variable rules

    while( i.n <= n.iter & n.diff > n.tol ){

      message('    Iteration ', i.n )
      message('      Solving the Euler equations...' )
      
      k.diff <- 2 * opt$k.tol
      k.diff.min <- k.diff
      k.gain <- opt$k.gain
      i.k <- 0
          # Intialize the loop variables
      fwd.idces <- sapply( fwd.vars, 
               function(nn) 
                 if(nn %in% endog.names){
                   n.exog + which(endog.names==nn)
                 }else{
                   (1+lags)*(n.exog+n.endog)+which(cont.names==nn)
                 } )
      # b.r.idces <- c(n.exog+1:2,(1+lags)*(n.endog+n.exog)+3:4)
          # Indices in X.cont
      
      while( k.diff > k.tol & i.k < k.iter & k.diff <= k.diff.min ){
        k.hat <- euler_hat_grid( coeff.n, coeff.c, X.cont, lags, params, 
                                 n.exog, n.endog, n.cont, n.fwd, rho, sig.eps, 
                                 0, N, upper, lower, cheby, matrix(0,1,1,), 
                                 TRUE, n.quad, model, mono )
            # The forward-looking variable predictors
#         v.k.diff <- apply( abs( k.hat - X.cont[, fwd.idces] ), 2, max )
        v.k.diff <- apply( abs( k.hat - X.cont[, fwd.idces] ), 2, mean )
        k.diff <- max( v.k.diff )
        if(i.k==0) k.diff.min <- k.diff
        if( k.diff > k.diff.min ){
          break
        }else{
          k.diff.min <- min( k.diff, k.diff.min )
        }
            # Check if solution is converging and exit the iteration if not
        i.k <- i.k + 1
            # Update loop variables
        X.cont[,fwd.idces] <- k.gain * k.hat + ( 1 - k.gain ) * X.cont[,fwd.idces]
            # Update the simulations for forward looking variables
        if( adapt.gain ) k.gain <- max( exp( - adapt.exp * k.diff ), k.gain )
            # Update the gain where required
        if( i.k %% 1 == 0 ) 
          message('        i.k = ', i.k, ',  v.k.diff = (', 
                  round(v.k.diff[1],5), ', ', round(v.k.diff[2],5),  ', ',
                  round(v.k.diff[3],5), ', ', round(v.k.diff[4],5), ' )' )
      }
      
      message('      ...complete.\n        Iterations = ', i.k,  
              '\n        Difference = ', round( k.diff, 5 ),
              '\n        k.gain = ', round( k.gain, 5 ) )
      message('      Updating rules...')
      
      ### 3.3 Compute the error-minimizing coefficients for the forward-looking variables ###
      for(i in 1:n.endog){
        if( endog.names[i] %in% fwd.vars ){
          coeff.n.new[,i] <- coeff_reg( k.hat[,which(fwd.vars==endog.names[i])], 
                                              X[, state.select], N,
                                              lower, upper, cheby )
        }
      } 
      for(i in 1:n.cont){
        if( cont.names[i] %in% fwd.vars ){
          coeff.c.new[,i] <- coeff_reg( k.hat[,which(fwd.vars==cont.names[i])], 
                                        X[, state.select], N,
                                        lower, upper, cheby )
        }
      }
          # Update the coefficients on the forward-looking variables in line
          # with the predictors from the conemporanneous block

      ### 3.4 Damp this part of the loop ###
      n.diff <- max( abs( cbind( coeff.n, coeff.c ) - cbind( coeff.n.new, coeff.c.new ) ) )
      i.n <- i.n + 1
          # Update the loop controls
      coeff.n <- n.gain * coeff.n.new + ( 1 - n.gain ) * coeff.n
      coeff.c <- n.gain * coeff.c.new + ( 1 - n.gain ) * coeff.c
          # Update the rules
      if( adapt.gain ) n.gain <- max( exp( - 100 * adapt.exp * n.diff ), n.gain )
          # Update the adaptive gain
#          # Turn off for this.
      if( sym.reg ){
        coeff.n <- m.sym.ave.pair( coeff.n, l.sym.ave, l.pairs )
        coeff.c <- m.sym.ave.pair( coeff.c, l.sym.ave, l.pairs.cont )
      }   # Symmetry regularization
      
      ### 3.5 Evaluate the rules on the state grid ###
      sim <- cont_sim( X, cbind( coeff.n, coeff.c ), N, n.endog, n.exog, upper, lower, cheby )
          # The new simulated values
      X.cont <- cbind( X[,1:n.exog], sim[,1:n.endog], 
                       X[,n.exog+n.endog+1:(n.exog+n.endog)], sim[,n.endog+1:n.cont] )
          # Re-create X      
      
      message('      ...complete.\n        Difference    = ', round( n.diff, 5 ),
              '\n        n.gain = ', round( n.gain, 5 ) )
    }

    ###### 4. MEASURE THE CHANGES IN THE STATE VARIABLE COEFFICIENTS ######
    message('  Outer maximum normalized difference = ', round( diff, 5 ) , "\n" )
        # The overall change
    pred <- contemp_eqns_irbc_grid( X.cont, lags, params, n.exog, n.endog, 
                                    n.cont, extra.args, opt$model )
        # The contemporaneous predictors
    colnames(pred) <- c( endog.names, cont.names )
        # Rename the columns
    pred[,fwd.vars] <- euler_hat_grid( coeff.n, coeff.c, X.cont, lags, params,
                                        n.exog, n.endog, n.cont, n.fwd, rho, 
                                        sig.eps, 0, N, upper, lower, cheby,
                                        matrix(0,1,1), TRUE, n.quad, model, mono )
        # The forward-looking variables
    err <- pred - X.cont[,contemp.select]
    bias.old <- if(exists('bias')) bias else Inf
    aad.old <- if(exists('aad')) aad else Inf
    bias <- mean(err)
    aad <- mean(abs(err))
    max.abs.var.err <- apply( abs( err ), 2, mean )
        # Error measures
    message('  Bias = ', round(max(bias), 5) )
    message('  Ave abs deviation = ', round(aad, 5) )
    message('  Max variable-average abs error = ', round(max(max.abs.var.err), 5) )
    message('  Variable with maximum mean absolute equation error is ', 
            c( endog.names, cont.names )[which.max(max.abs.var.err)] )  
#     message('  Maximum absolute mean equation error = ', round(max(abs(max.eq.err)), 5) )
        # Maximum equation error
    diff <- max(max.abs.var.err)
  
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
    
    if( abs(bias) > abs(bias.old) & aad > aad.old ){
      message('Increasing absolute bias and aad.  Aborting iteration (try reducing n.gain?).')
      break
    }else{
      coeff.cont <- coeff.c
      coeff <- coeff.n
      coeff.old <- coeff
          # Update coefficients
    }
  }
  
  out <- list( coeff=coeff )
  if( n.cont > 0 ) out$coeff.cont <- coeff.cont
  out$X.cont <- X.cont
  out$opt <- opt
  out$params <- params
  out$err <- err
      # Set up the output
  return( out )
}

sol.irbc.check <- function( sol, params=NULL, opt=NULL ){
# Checks the model equation errors on the reduced state space X.cont
  
  if( is.null(params) ) params <- sol$params
  if( is.null(opt) ) opt <- sol$opt
      # Assign the options and params when required
  
  n.check <- 20000
  n.burn <- 1000
  n.quad <- 4
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
  endog.names <- opt$endog.names
  cont.names <- opt$cont.names
  fwd.vars <- opt$fwd.vars
  n.fwd <- opt$n.fwd
  model <- opt$model
      # Copy from options
  extra.args <- list( n.fwd=opt$n.fwd, y1.ss=opt$ys['Y1'] )
  
  rho <- params$rho
  sig.eps <- params$sig.eps
      # Copy from parameters
  
  # exog.sim <- sapply( 1:n.exog, function(i) ar1_sim( kappa * n.check + n.burn, 
  #                                                    rho[i], sig.eps[i] ) )
  exog.sim <- var1_sim( kappa * n.check + n.burn, rho, sig.eps )
      # The exogenous simulation
  

  
  endog.sim <- endog_sim( n.check, exog.sim, sol$coeff, N, upper, lower, 
                          rep(0,n.endog), cheby, kappa, n.burn, (lags>0) )
      # The endogenous simulation (Here set kappa=1)
  cont.sim <- cont_sim( endog.sim, sol$coeff.cont, N, n.endog, n.exog, upper, lower, cheby )
      # The controls
  all.sim <- cbind( endog.sim, cont.sim )
      # The combined simulation
  
  pred <- contemp_eqns_irbc_grid( all.sim, lags, params, n.exog, n.endog, 
                                  n.cont, extra.args, opt$model )
      # The contemporaneous predictors
  colnames(pred) <- c( endog.names, cont.names )
      # Rename the columns
  pred[,fwd.vars] <- euler_hat_grid( sol$coeff, sol$coeff.cont, all.sim, lags, params,
                                     n.exog, n.endog, n.cont, n.fwd, rho, 
                                     sig.eps, 0, N, upper, lower, cheby,
                                     matrix(0,1,1), TRUE, n.quad, model )
      # The forward-looking predictors  
  contemp.select <- c( n.exog+1:n.endog, (1+lags)*(n.endog+n.exog)+1:n.cont )
  err <- cbind( pred - all.sim[,contemp.select])
                               # The indices of the contemporaenous variables] )
      # The errors
  out <- list( endog.sim=endog.sim, cont.sim=cont.sim, err=err )
  
  return( out )
}



