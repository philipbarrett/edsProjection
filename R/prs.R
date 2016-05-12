###############################################################################
# Solving the model with perfect risk-sharing
###############################################################################

## To do: Turn this into a function

prs.sol.sim <- function( n.pds, n.sample, n.burn, params ){
# Produce a simulation from the perfect risk-sharing problem
  
  alpha <- params$share
  eta <- params$eta
  betta <- params$betta
  p1.bar <- log( params$P1.bar )
  p2.bar <- log( params$P2.bar )
      # Extract parameters
  n.sim <- n.pds * n.sample + n.burn
  
  a.1 <- ar1_sim( n.sim, params$rho[1], params$sig.eps[1] )
  a.2 <- ar1_sim( n.sim, params$rho[2], params$sig.eps[2] )
      # Full technology simulation
  a <- matrix( NA, n.pds, 2 )
  idx <- n.burn + 1:n.pds * n.sample
      # Indices of the subsample
  a[,1] <- a.1[idx]
  a[,2] <- a.2[idx]
#   a <- cbind( ar1_sim( n.pds + n.burn, params$rho[1], params$sig.eps[1] ),
#               ar1_sim( n.pds + n.burn, params$rho[2], params$sig.eps[2] ) )
      # The vector of (log) shocks
  x.11 <- log( alpha ) + a[,1]
  x.12 <- log( 1 - alpha ) + a[,2]
  x.21 <- log( 1 - alpha ) + a[,1]
  x.22 <- log( alpha ) + a[,2]
      # The log intermediates
  c.1 <- alpha * x.11 + ( 1 - alpha ) * x.12
  c.2 <- alpha * x.22 + ( 1 - alpha ) * x.21
      # The log consumptions
  p.12 <- log( ( 1 - alpha ) * params$P1.bar / alpha ) + ( x.11 - x.12 ) / eta
  p.21 <- log( ( 1 - alpha ) * params$P2.bar / alpha ) + ( x.22 - x.21 ) / eta
      # Trade prices
  p.1 <- - log( 1 - alpha ) + p.12 - c.1 + x.12
  p.2 <- - log( 1 - alpha ) + p.21 - c.2 + x.21
      # Aggregate prices
  e.12 <- .5 * ( p1.bar - p2.bar + p.12 - p.21 )
      # The nominal exchange rate
  r.1 <- c( -log(betta) + params$gamma * ( c.1[-1] - c.1[-length(c.1)] ) + 
              p.1[-1] - p.1[-length(p.1)], NA )
  r.2 <- c( -log(betta) + params$gamma * ( c.2[-1] - c.2[-length(c.2)] ) + 
              p.2[-1] - p.2[-length(p.1)], NA )
      # Nominal interest rates
  states <- cbind( a, a.1[idx-1], a.2[idx-1] )
      # Forming the states with lags
  cont <- cbind( c.1, c.2, r.1, r.2, x.11, x.22, x.12, x.21, 
                 p.1, p.2, p.12, p.21, e.12 )
      # The controls
  return( list( states=states, cont=cont ) )
}

prs.sol.cont.a <- function( sol, opt ){
# Computes the coefficients on a for a solution for the rs model in terms of the
# technology shocks a
  n.terms <- idx_count( opt$N, opt$n.exog )
  coeff.cont <- matrix( NA, n.terms, opt$n.cont )
      # Initialize the matrix of control coefficients
  for( i in 1:opt$n.cont )
    coeff.cont[,i] <- coeff_reg( sol$cont[,i], sol$states[,1:opt$n.exog], opt$N, 
                                 opt$lower[1:opt$n.exog], opt$upper[1:opt$n.exog], 
                                 opt$cheby )
  for( i in 3:4 )
    coeff.cont[,i] <- coeff_reg( sol$cont[-nrow(sol$cont),i], 
                                 sol$states[-nrow(sol$cont),1:opt$n.exog], opt$N, 
                                 opt$lower[1:opt$n.exog], opt$upper[1:opt$n.exog], 
                                 opt$cheby )  
  return( coeff.cont )
}

prs.sol.sim.coeff <- function( states, coeff.cont, opt ){
# Simulates the controls for a given set of states using the coefficients
# (only really useful for r, but a good check anyway)
  out <- cont_sim( states, coeff.cont, opt$N, opt$n.exog, 0, 
                   opt$upper[1:opt$n.exog], opt$lower[1:opt$n.exog], opt$cheby )
  return(out)
}

prs.sol.sim.B <- function( states, cont ){
# Simulates the sequence of debt levels B consistent with the states A and the
# controls
  n.sim <- nrow(cont)
  B.11 <- rep( 0, n.sim )
  B.22 <- rep( 0, n.sim )
  B.11[1] <- B.22[1] <- .5
      # Initialize the simulation for B
  a <- states[,1:2]
  c.1 <- cont[,1]
  c.2 <- cont[,2]
  r.1 <- cont[,3]
  r.2 <- cont[,4]
  p.1 <- cont[,9]
  p.2 <- cont[,10]
  e.12 <- cont[,13]
      # Extract the controls
  for( i in 2:(n.sim)){
    B.11[i] <- ( - exp(a[i,1]) + ( 1 + exp( e.12[i] ) ) * B.11[i-1] + 
                   exp( p.1[i] + c.1[i] -  exp( e.12[i] ) * ( 1 - exp( - r.2[i] ) ) ) ) / 
                      ( exp( - r.1[i] ) - exp( e.12[i] - r.2[i] ) )
    B.22[i] <- ( - exp(a[i,2]) + ( 1 + exp( - e.12[i] ) ) * B.22[i-1] + 
                   exp( p.2[i] + c.2[i] ) -  exp( - e.12[i] ) * ( 1 - exp( - r.1[i] ) ) ) /
                  ( exp( - r.2[i] ) - exp( - e.12[i] - r.1[i] ) )
  }
      # Compute debt holdings from the budget constraint (and B.11=B.22)
  return(cbind(B.11, B.22))
}

