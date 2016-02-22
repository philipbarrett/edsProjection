#########################################################################
# sol.R
#
# Generic solution algorithm for the EDS-projection algorithm
# Philip Barrett, Chicago
# 02feb2016
#########################################################################

sol.iterate <- function( coeff.init, opt, params, coeff.cont.init = NULL, debug.flag=FALSE ){
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
  image <- opt$image
  inner.iter <- opt$inner.iter
  inner.tol <- opt$inner.tol
  sym.reg <- opt$sym.reg
  l.sym.ave <- opt$l.sym.ave
  l.pairs <- opt$l.pairs
  
  # Extract parameters
  rho <- params$rho
  sig.eps <- params$sig.eps
  
  n.terms <- idx_count( N, n.exog + n.endog )
  
  set.seed(1234)
  exog.sim <- sapply( 1:n.exog, function(i) ar1_sim( n.sim * kappa + burn, 
                                                     rho[i], sig.eps[i] ) )
      # Create the exogenous simulation
  n.state <- n.endog + n.exog
  if( is.null( coeff.cont.init ) ) coeff.cont.init <- 0 * coeff.init
  coeff <- coeff.init
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
      

      
      ########################################################################
      # HACK ALERT!! HACK ALERT!! HACK ALERT!!
      #  >n<  >n<  >n<  >n<  >n<  >n<  >n<  >n<
      ########################################################################
      
#       coeff.inner.new[c(3,5,6),] <- 1e-09
      
      ########################################################################
      # HACK PASSED!! HACK PASSED!! HACK PASSED!! 
      #  >n<  >n<  >n<  >n<  >n<  >n<  >n<  >n<
      ########################################################################
      
      
      
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

sol.check <- function( sol, opt, params ){
# Computes a high-precision error on the Euler equation and checks that the bounds are observed
  
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
  image <- opt$image
  
  # Extract parameters
  rho <- params$rho
  sig.eps <- params$sig.eps
  
  coeff <- sol$coeff
  coeff.cont <- if( n.cont > 0 ) sol$coeff.cont else 0 * coeff
  
  # Simulation
  set.seed(4321)
  exog.sim <- sapply( 1:n.exog, function(i) ar1_sim( n.sim * kappa + burn, 
                                                     rho[i], sig.eps[i] ) )
  endog.sim <- endog_sim( n.sim, exog.sim, coeff, N, upper, lower, 
                          endog.init, cheby, kappa, burn, (lags>0) )
  if( n.cont > 0 ){
    cont.sim <- cont_sim( endog.sim, coeff.cont, N, n.endog, n.exog, upper, lower, cheby )
    X.cont <- cbind( endog.sim, cont.sim )
  }else{
    X.cont <- endog.sim
  }
  
  k.hat <- euler_hat( sol$coeff, coeff.cont, X.cont, model, lags, params, n.exog, 
                      n.endog, n.cont, rho, sig.eps, 0, N, upper, lower, cheby, 
                      matrix(0,1,1), TRUE, 10 )
      # Always increase the number of quadrature nodes to the maximum
  state.select <- n.exog + 1:(n.endog)
  if( n.cont > 0 ) state.select <- 
                c( state.select,
                    ( 1 + lags ) * ( n.exog + n.endog ) + 1:n.cont )
      # Also select the controls
  rel.err <- matrix( abs( k.hat / X.cont[, state.select] - 1 ) * 100, ncol=n.endog+n.cont)
  upper.check <- sapply( 1:(n.exog+n.endog), function( i ) all( endog.sim[,i] < upper[i] ) )
  lower.check <- sapply( 1:(n.exog+n.endog), function( i ) all( endog.sim[,i] < upper[i] ) )
  
  return( list( max=apply(rel.err, 2, max), ave=apply(rel.err, 2, mean),
                upper.check=upper.check, lower.check=lower.check ) )
}

sym.pair <- function( x ){
# Regularizes two pairs of coeffenients by imposing symmetry on them.
  t.l <- mean( c( x[1,1], x[2,2] ) )
  b.l <- mean( c( x[2,1], x[1,2] ) )
      # The top left and bottom left means
  return( matrix( c( t.l, b.l, b.l, t.l ), 2, 2 ) )
}

ave.pair <- function( x ){
# Regularizes a pair of coefficients by asserting that they are the same
  return( rep( mean(x), 2 ) )
}

sym.ave.pair <- function( X, l.sym.ave ){
# Regularizes a 2-column matrix of coefficients X by imposing symmetry and 
# equality on the row dictated in the list l.sym.ave with members $sym and $ave
  out <- X
  for( idx in l.sym.ave$sym )
    out[ idx, ] <- sym.pair( X[ idx, ] )
  for( idx in l.sym.ave$ave )
    out[ idx, ] <- ave.pair( X[ idx, ] )
  return( out )
}

m.sym.ave.pair <- function( X, l.sym.ave, l.pairs ){
# Applies sym.ave.pair to the list of pairs in l.pairs
  out <- X
  for( idx in l.pairs )
    out[, idx] <- sym.ave.pair( X[, idx], l.sym.ave )
  return( out )
}