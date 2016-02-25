##########################################################################################
# Trying to figure out why the equations for r are so sensitive
#
# 24feb2016
# Philip Barrett, Chicago
##########################################################################################

params <- list( alpha = .85, gamma = 5, P1.bar=1, P2.bar=1, betta=.99,
                rho=c(.95,.95), sig.eps=c(.01,.01), eta=1 )
# Parameters
lower <- sd.x<- sqrt( params$sig.eps / ( 1 - params$rho ^ 2 ) )
upper <- c(  3 * sd.x, rep( 1, 2 ) )
lower <- -upper
# Bounds

l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                   ave=c(1,5,14) )
l.pairs <- list( c(1,2) )
l.pairs.cont <- list( c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12) )
# Symmetry constraints

opt <- list( lags=1, n.exog=2, n.endog=2, n.cont=13, N=2, cheby=FALSE,
             upper = upper, lower=lower, quad=TRUE, n.quad=3,  burn=1000,
             kappa=25, n.sim=10000, eps = .7, delta=.05, endog.init=c(0, 0), 
             c.iter=400, c.tol=1e-07, c.gain=.25,
             k.iter=0, k.tol=1e-05, k.gain=.1,
             n.iter=0, n.tol=1e-05, n.gain=.05, tol=1e-05, iter=1,
             sr=TRUE, adapt.gain=TRUE, adapt.exp=20, image=FALSE,
             sym.reg=FALSE, l.sym.ave=l.sym.ave, l.pairs=l.pairs, 
             l.pairs.cont=l.pairs.cont )
    ## NB: No iterations for the state rules.  That's what I want to look at
    ## here.  Just solve for consistent control rules.

coeff.init <- matrix( 0, 15, 2 )
coeff.init[2,] <- c(-.1,.4)
coeff.init[4,] <- c(.4,-.1)
coeff.init[7,] <- c(.2,.5)
coeff.init[11,] <- c(.5,.2)

coeff.cont.init <- matrix( 0, 15, opt$n.cont )
c.ss <- params$alpha ^ params$alpha * ( 1 - params$alpha ) ^ ( 1 - params$alpha )
p.ss <- 1 / c.ss
coeff.cont.init[1, ] <- log( c( c.ss, c.ss, 1 / params$betta, 1 / params$betta, 
                                params$alpha, params$alpha, 1 - params$alpha, 1 - params$alpha,
                                p.ss, p.ss, 1, 1, 1 ) )
coeff.cont.init[2, ] <- c( -.1, .5, -.2, .5, rep( .1, 9 ) )
coeff.cont.init[4, ] <- c(  .5,-.1,  .5,-.2, rep( .1, 9 ) )
coeff.cont.init[7, ] <- coeff.cont.init[2, ]
coeff.cont.init[11, ] <- coeff.cont.init[4, ]

sol <- sol.irbc.iterate( coeff.init, opt, params, coeff.cont.init )
    ## Create a solution for just the controls

test <- sol$X.cont[,-(1:8)] - contemp_eqns_irbc_grid( sol$X.cont, opt$lags, params, opt$n.exog, 
                                                      opt$n.endog, opt$n.cont )
max(abs(test))


### SOLVE THE EULER EQUATIONS BY SIMULTANEOUS ITERATION ###
gain <- .5
n.loop <- 60
    # Number of iterations of the informal loop to compute B and r
X.cont <- sol$X.cont
    # Initiate the grid
B.11 <- X.cont[,3]
B.22 <- X.cont[,4]
r.1 <- X.cont[,8+3]
r.2 <- X.cont[,8+4]

for( i in 1:n.loop){

  k.hat <- euler_hat_grid( sol$coeff, sol$coeff.cont, X.cont, opt$lags, params, opt$n.exog, 
                           opt$n.endog, opt$n.cont, params$rho, params$sig.eps, 0, opt$N, 
                           opt$upper, opt$lower, opt$cheby, matrix(0,1,1,), TRUE, opt$n.quad )
      # Create new predictors for B and r
  B.11 <- cbind( B.11, k.hat[,1] )
  B.22 <- cbind( B.22, k.hat[,2] )
  r.1 <- cbind( r.1, k.hat[,3] )
  r.2 <- cbind( r.2, k.hat[,4] )
      # Update the guesses
  diff <- apply( abs( cbind( k.hat[,1:2] - X.cont[ , 3:4],
                             k.hat[,3:4] - X.cont[ , 11:12] ) ), 2, max )
  if( i %% 5 == 0 ){
    message( "i = ", i)
    message(  ' diff = ( ', round(diff[1],4), ', ', round(diff[2],4),  ', ',
                            round(diff[3],4), ', ', round(diff[4],4), ' )')
  }
  X.cont[ , c(3,4,11,12)] <- gain * k.hat + ( 1 - gain ) * X.cont[ , c(3,4,11,12)]
}
    # Min max error almost zero


### SOLVE FOR THE INTEREST RATE ONLY ###
n.loop <- 5
gain <- 1
X.cont <- sol$X.cont
    # Initiate the grid
r.1 <- X.cont[,8+3]
r.2 <- X.cont[,8+4]

for( i in 1:n.loop){
  message( "i = ", i)
  k.hat <- euler_hat_grid( sol$coeff, sol$coeff.cont, X.cont, opt$lags, params, opt$n.exog, 
                           opt$n.endog, opt$n.cont, params$rho, params$sig.eps, 0, opt$N, 
                           opt$upper, opt$lower, opt$cheby, matrix(0,1,1,), TRUE, opt$n.quad )
  # Create new predictors for B and r
  r.1 <- cbind( r.1, k.hat[,3] )
  r.2 <- cbind( r.2, k.hat[,4] )
  # Update the guesses
  diff <- apply( abs( k.hat[,3:4] - X.cont[ , c(11,12)] ), 2, mean )
  message(  ' diff = ( ', round(diff[1],4), ', ', round(diff[2],4), ' )')
  X.cont[ , c(11,12)] <- gain * k.hat[,3:4] + ( 1 - gain ) * X.cont[ , c(11,12)]
}
    # This is just one step (as it should be)


### Now solve B, then r, then B, then r etc... ###
n.inner <- 10
n.outer <- 40
gain.inner <- .05
gain.outer <- .5
X.cont <- sol$X.cont
    # Initiate the grid
r.1 <- X.cont[,8+3]
r.2 <- X.cont[,8+4]

for( i in 1:n.outer){
  message( "i = ", i)
  
  for( j in 1:n.inner){
    k.hat <- euler_hat_grid( sol$coeff, sol$coeff.cont, X.cont, opt$lags, params, opt$n.exog, 
                             opt$n.endog, opt$n.cont, params$rho, params$sig.eps, 0, opt$N, 
                             opt$upper, opt$lower, opt$cheby, matrix(0,1,1,), TRUE, opt$n.quad )
    inner.diff <- apply( abs( k.hat[,1:2] - X.cont[ , 3:4] ), 2, mean )
    message(  ' inner.diff = ( ', round(inner.diff[1],4), ', ', round(inner.diff[2],4), ' )')
    X.cont[ , c(3,4)] <- gain.inner * k.hat[,1:2] + ( 1 - gain.inner ) * X.cont[ , c(3,4)]
  }
  
  k.hat <- euler_hat_grid( sol$coeff, sol$coeff.cont, X.cont, opt$lags, params, opt$n.exog, 
                           opt$n.endog, opt$n.cont, params$rho, params$sig.eps, 0, opt$N, 
                           opt$upper, opt$lower, opt$cheby, matrix(0,1,1,), TRUE, opt$n.quad )
      # Create new predictors for B and r
  r.1 <- cbind( r.1, k.hat[,3] )
  r.2 <- cbind( r.2, k.hat[,4] )
      # Update the guesses
  outer.diff <- apply( abs( k.hat[,3:4] - X.cont[ , c(11,12)] ), 2, mean )
  message(  'outer.diff = ( ', round(outer.diff[1],4), ', ', round(outer.diff[2],4), ' )\n')
  X.cont[ , c(11,12)] <- gain.outer * k.hat[,3:4] + ( 1 - gain.outer ) * X.cont[ , c(11,12)]
}



### TRY PLOTTING THE SOLUTION IN B ALONE ###

exog <- matrix( 0, 2, 2 )
endog <- exog
cont <- log( c( c.ss, c.ss, 1 / params$betta, 1 / params$betta, 
                params$alpha, params$alpha, 1 - params$alpha, 1 - params$alpha, 
                p.ss, p.ss, 1, 1, 1 ) )

node.wts <- quad_nodes_weights( opt$n.quad, opt$n.exog, params$sig.eps, c(0,0) )

test <- euler_hat_irbc( exog, endog, cont, node.wts$nodes, params, sol$coeff.cont, opt$n.exog, 
                        opt$n.endog, opt$n.cont, params$rho, opt$n.quad ^ opt$n.exog, opt$N, 
                        opt$upper, opt$lower, opt$cheby, node.wts$weights, TRUE )
test - c( 0, 0, cont[3:4] )

k.hat <- euler_hat_grid( sol$coeff, sol$coeff.cont, sol$X.cont, opt$lags, params, opt$n.exog, 
                         opt$n.endog, opt$n.cont, params$rho, params$sig.eps, 0, opt$N, 
                         opt$upper, opt$lower, opt$cheby, matrix(0,1,1,), TRUE, opt$n.quad )

i.row <- 50
just.B.11 <- function(B.11){
  this.endog <- matrix( c( B.11, sol$X.cont[i.row, c(4,7:8)] ), 2, 2, byrow=T )
  this.exog <- matrix( sol$X.cont[i.row, c(1:2,5:6)], 2, 2, byrow=T )
  this.cont <- sol$X.cont[i.row, -(1:8)]
  return( euler_hat_irbc( this.exog, this.endog, this.cont, node.wts$nodes, params, sol$coeff.cont, opt$n.exog, 
                          opt$n.endog, opt$n.cont, params$rho, opt$n.quad ^ opt$n.exog, opt$N, 
                          opt$upper, opt$lower, opt$cheby, node.wts$weights, FALSE )[1] )
}

B.11.seq <- seq( sol$X.cont[i.row,3] - .01, sol$X.cont[i.row,3] + .01, length.out=101)
plot( B.11.seq, sapply( B.11.seq, just.B.11), type='l', lwd=2 )
abline(0,1)
abline( v=sol$X.cont[i.row,3], lty=2 )
abline( h=k.hat[i.row,1], lty=2, col=2 )

######  AHAHAHAHAHAAAAA!!!!  The smoking gun.  The equations for B are not 
######  stable under repeated substitution (presumably this is the reason we 
######  have uniqueness, right? It's due to divergence of other asset 
######  holdings.). So I'm going to need to do repeated bisection (or some such)
######  for the B solutions.
######  FIXED!!! Just subtract the Euler equation instead of adding it!

just.B.22 <- function(B.22){
  this.endog <- matrix( c( sol$X.cont[i.row, 3], B.22, sol$X.cont[i.row, 7:8] ), 2, 2, byrow=T )
  this.exog <- matrix( sol$X.cont[i.row, c(1:2,5:6)], 2, 2, byrow=T )
  this.cont <- sol$X.cont[i.row, -(1:8)]
  return( euler_hat_irbc( this.exog, this.endog, this.cont, node.wts$nodes, params, sol$coeff.cont, opt$n.exog, 
                          opt$n.endog, opt$n.cont, params$rho, opt$n.quad ^ opt$n.exog, opt$N, 
                          opt$upper, opt$lower, opt$cheby, node.wts$weights, FALSE )[2] )
}

B.22.seq <- seq( sol$X.cont[i.row,4] - .01, sol$X.cont[i.row,4] + .01, length.out=101)
plot( B.22.seq, sapply( B.22.seq, just.B.22), type='l', lwd=2 )
abline(0,1)
abline( v=sol$X.cont[i.row,4], lty=2 )
abline( h=k.hat[i.row,2], lty=2, col=2 )

