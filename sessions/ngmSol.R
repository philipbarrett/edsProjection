params.full.dep <- list( A=1, alpha=.3, delta=1, gamma=1, betta=.99, 
                         rho=0, sig.eps=.01 )
k.ss <- ( ( 1 / params.full.dep$betta - 1 + params.full.dep$delta ) / 
            ( params.full.dep$A * params.full.dep$alpha )  ) ^ ( 1 / ( params.full.dep$alpha - 1 ) )
# The steady state
sd.x<- sqrt( params.full.dep$sig.eps / ( 1 - params.full.dep$rho ^ 2 ) )
upper <- c(  3 * sd.x, k.ss * 1.4 )
lower <- c( -3 * sd.x, k.ss * 0.6 )
# Create the bounds
opt <- list( model='ngm', lags=1, n.exog=1, n.endog=1, N=1, cheby=FALSE, 
             upper = upper, lower=lower, quad=TRUE, n.quad=7,
             err.tol=1e-05, n.iter=10, burn=1000, kappa=25, n.sim=10000,
             eps = .3, delta=.02, endog.init=k.ss, gain=.1 )
# Solution options
a.b.range <- upper - lower
coeff.init <- matrix( c( k.ss, 
                         params.full.dep$alpha * a.b.range[2] / 2, 
                         k.ss * a.b.range[1] / 2 ), 3, 1)
# The linear-approx solution evaluated @ the steady state:
#   k' = alpha * betta * exp(x) * A * k ^ alpha
#     ~= k.ss + k.ss * x + alpha * ( k - k.ss )
#   NB: The range is already centered on k.ss, so can just use 
#       these coeffs

k.ss.prime <- endog_update( exog = c(0), k.ss, coeff.init, opt$n.exog,
                            opt$n.endog, opt$N, upper, lower, opt$cheby )
expect_true( k.ss == k.ss.prime )

k.ss.euler <- 
  integrand_ngm( matrix( 0, 1, 1), matrix( k.ss, nrow=2, ncol=1), 0, 
                 params.full.dep, coeff.init, 1, 1, 1, upper, lower, FALSE )
expect_equal( k.ss.euler * params.full.dep$betta, 1 )

nodes <- quad_nodes_1d( opt$n.quad, params.full.dep$sig.eps, 0 )
wts <- quad_weights_1d( opt$n.quad )

single.integral <- 
  err_ngm( matrix( 0, 1, 1), matrix( k.ss, nrow=2, ncol=1), nodes, 
                  params.full.dep, coeff.init, opt$n.exog, opt$n.endog, 
                  params.full.dep$rho, opt$n.quad, opt$N, upper, lower, FALSE, wts, TRUE )

many.pts <- 
  eval_err( coeff.init, matrix( c( 0, k.ss, 0, k.ss ), nrow=1 ), 'ngm', opt$lags, 
            params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
            params.full.dep$sig.eps,  0, opt$N, upper, lower, FALSE, 
            matrix(0,1,1,), TRUE, opt$n.quad )


# sol <- sol.iterate( coeff.init, opt, params.full.dep, FALSE )

exact.sol <- function( k, x, params ){
# The exact solution
  return( params$alpha * params$betta * exp(x) * params$A * k ^ params$alpha )
}

KK <- seq( lower[2], upper[2], length.out=201 )
XX <- seq( lower[1], upper[1], length.out=5 )
pol.sol.exact <- sapply( XX, function(x) 
  sapply(KK, function(k) exact.sol(k,x,params.full.dep) ) )
plot( range(KK), range(pol.sol.exact), typ='n' )
for( i in 1:5 ) lines( KK, pol.sol.exact[, i], col=i, lwd=2 )
abline( v=k.ss )
abline( a=0, b=1 )
for( i in 1:5 ) abline( a=(1-params.full.dep$alpha+XX[i])*k.ss, 
                        b=params.full.dep$alpha, col=i, lwd=2, lty=2 )

## Evaluate the method on a small grid near k.ss
exog.0 <- matrix( rep( 0, 20 ), ncol=1 )
endog.sim.l <- endog_sim( 20, exog.0, coeff.init, opt$N, opt$upper,
                          opt$lower, .8*k.ss, FALSE, 1, 0, opt$lags )
endog.sim.u <- endog_sim( 20, exog.0, coeff.init, opt$N, opt$upper,
                          opt$lower, 1.2*k.ss, FALSE, 1, 0, opt$lags )
plot( c(1,20), range( endog.sim.l[,4], endog.sim.u[,4] ), type='n' )
lines( 1:20, endog.sim.l[,4], lwd=2 )
lines( 1:20, endog.sim.u[,4], lwd=2 )
endog.sim <- rbind( endog.sim.u, endog.sim.l )

eval_err( coeff.init, endog.sim[15:20,], 'ngm', opt$lags, 
          params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
          params.full.dep$sig.eps,  0, opt$N, opt$upper, opt$lower, FALSE, 
          matrix(0,1,1,), TRUE, opt$n.quad )
    # Error for very near ss
eval_err( coeff.init, endog.sim, 'ngm', opt$lags, 
          params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
          params.full.dep$sig.eps,  0, opt$N, opt$upper, opt$lower, FALSE, 
          matrix(0,1,1,), TRUE, opt$n.quad )
    # Error further from the ss

eval_err_D( coeff.init, endog.sim[15:20,], 'ngm', opt$lags, 
          params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
          params.full.dep$sig.eps,  0, opt$N, opt$upper, opt$lower, FALSE, 
          matrix(0,1,1,), TRUE, opt$n.quad )
    # Error for very near ss
eval_err_D( coeff.init, endog.sim, 'ngm', opt$lags, 
          params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
          params.full.dep$sig.eps,  0, opt$N, opt$upper, opt$lower, FALSE, 
          matrix(0,1,1,), TRUE, opt$n.quad )
    # Error further from the ss

coeff.new <- matrix(
  err.min( coeff.init, endog.sim, 'ngm', opt$lags, 
           params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
           params.full.dep$sig.eps,  0, opt$N, opt$upper, opt$lower, FALSE, 
           matrix(0,1,1,), TRUE, opt$n.quad )$coeff, ncol=1 )

k.prime.apx <- sapply( KK, function(k) 
                  endog_update( 0, k, coeff.init, opt$n.exog, 
                          opt$n.endog, opt$N, opt$upper, opt$lower, FALSE ) )
k.prime.apx.new <- sapply( KK, function(k) 
  endog_update( 0, k, coeff.new, opt$n.exog, 
                opt$n.endog, opt$N, opt$upper, opt$lower, FALSE ) )

plot( range(KK), range(pol.sol.exact), typ='n' )
lines( KK, pol.sol.exact[, 3], col=3, lwd=2 )
abline( v=k.ss )
abline( a=0, b=1 )
lines( KK, k.prime.apx, col=3, lwd=2, lty=2 )
lines( KK, k.prime.apx.new, col=3, lwd=2, lty=3 )
abline(v=endog.sim[,4], lty=3)

## CHECKING THE ERROR EVALUATION
err_ngm( matrix(endog.sim[1,1],1,1), matrix(endog.sim[1,c(2,4)],2,1), nodes, 
         params.full.dep, coeff.init, opt$n.exog, opt$n.endog, params.full.dep$rho, opt$n.quad, 
         opt$N, opt$upper, opt$lower, FALSE, wts, TRUE )

err_ngm( matrix(endog.sim[1,1],1,1), matrix(endog.sim[1,c(4,2)],2,1), nodes, 
         params.full.dep, coeff.init, opt$n.exog, opt$n.endog, params.full.dep$rho, opt$n.quad, 
         opt$N, opt$upper, opt$lower, FALSE, wts, TRUE )

eval_err( coeff.init, matrix( endog.sim[1,], nrow=1), 'ngm', opt$lags, 
          params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
          params.full.dep$sig.eps,  0, opt$N, upper, lower, FALSE, 
          matrix(0,1,1), TRUE, opt$n.quad )

eval_err_D( coeff.init, matrix( endog.sim[1,], nrow=1), 'ngm', opt$lags, 
            params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
            params.full.dep$sig.eps,  0, opt$N, upper, lower, FALSE, 
            matrix(0,1,1), TRUE, opt$n.quad )

eval_err_D( coeff.init, matrix( endog.sim[19:20,], nrow=2), 'ngm', opt$lags, 
            params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
            params.full.dep$sig.eps,  0, opt$N, upper, lower, FALSE, 
            matrix(0,1,1), TRUE, opt$n.quad )
    ## WHY IS THIS NOT ZERO!?!?!?!?

err.min( coeff.init, matrix( endog.sim[1,], nrow=1), 'ngm', opt$lags, params.full.dep,
         opt$n.exog, opt$n.endog,
         params.full.dep$rho, params.full.dep$sig.eps, 0, opt$N, upper, lower, FALSE, 
         matrix(0,1,1), TRUE, opt$n.quad, TRUE )

err.min( coeff.init, endog.sim, 'ngm', opt$lags, params.full.dep,
         opt$n.exog, opt$n.endog,
         params.full.dep$rho, params.full.dep$sig.eps, 0, opt$N, upper, lower, FALSE, 
         matrix(0,1,1), TRUE, opt$n.quad, TRUE )

# sol <- sol.iterate( coeff.init, opt, params.full.dep, TRUE )

coeff.k.inc <- seq( -coeff.init[2], coeff.init[2], length.out=201 )
err.seq <- sapply( coeff.k.inc, function(coeff.k) 
  eval_err( coeff.init + c(0,coeff.k,0), endog.sim[19:20,], 'ngm', opt$lags, 
#   eval_err( coeff.init + c(0,coeff.k,0), endog.sim, 'ngm', opt$lags, 
            params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
            params.full.dep$sig.eps,  0, opt$N, upper, lower, FALSE, 
            matrix(0,1,1), TRUE, opt$n.quad ) )
plot( coeff.init[2] + coeff.k.inc, err.seq, type='l' )
abline( v=coeff.init[2], lty=2 )
abline( v=coeff.new[2], lty=2, col=2 )




## Iterating over successivve solutions on a grid with x=0
coeff.new <- coeff.init
par( mfrow=c(3,1) )
for( i in 1:2 ){
  coeff.old <- coeff.new
  message( "i = ", i )
  message( "coeff.old = (", coeff.old[1], " ", 
           coeff.old[2], " ", coeff.old[3], " )" )
  endog.sim.l <- endog_sim( 10, exog.0, coeff.new, opt$N, opt$upper,
                            opt$lower, .8*k.ss, FALSE, 1, 0, opt$lags )
  endog.sim.u <- endog_sim( 10, exog.0, coeff.new, opt$N, opt$upper,
                            opt$lower, 1.2*k.ss, FALSE, 1, 0, opt$lags )
  endog.sim <- rbind( endog.sim.u, endog.sim.l )
  plot( c(1,10), range( endog.sim.l[,4], endog.sim.u[,4] ), type='n', 
        main=paste( "i =", i ) )
  lines( 1:10, endog.sim.l[,4], lwd=2 )
  lines( 1:10, endog.sim.u[,4], lwd=2 )
  plot(endog.sim[,c(4,2)])
  coeff.new <- matrix(
    err.min( coeff.new, endog.sim, 'ngm', opt$lags, 
             params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
             params.full.dep$sig.eps,  0, opt$N, opt$upper, opt$lower, FALSE, 
             matrix(0,1,1,), TRUE, opt$n.quad )$coeff, ncol=1 )
  k.prime.apx.new <- sapply( KK, function(k) 
    endog_update( 0, k, coeff.new, opt$n.exog, 
                  opt$n.endog, opt$N, opt$upper, opt$lower, FALSE ) )
  
  plot( range(KK), range(pol.sol.exact), typ='n' )
  lines( KK, pol.sol.exact[, 3], col=3, lwd=2 )
  abline( v=k.ss )
  abline( a=0, b=1 )
  lines( KK, k.prime.apx, col=3, lwd=2, lty=2 )
  lines( KK, k.prime.apx.new, col=3, lwd=2, lty=3 )
  abline(v=endog.sim[,4], lty=3)
}
par( mfrow=c(1,1) )

# Use the second-iteration coefficients to examine the Euler equation errors
k.ss.prime.2 <- endog_update( exog = c(0), k.ss, coeff.old, opt$n.exog,
                            opt$n.endog, opt$N, upper, lower, opt$cheby )
expect_true( k.ss == k.ss.prime )

k.ss.euler.2 <- 
  integrand_ngm( matrix( 0, 1, 1), matrix( k.ss, nrow=2, ncol=1), 0, 
                 params.full.dep, coeff.old, 1, 1, 1, upper, lower, FALSE )
# expect_equal( k.ss.euler.2 * params.full.dep$betta, 1 )
    # THIS FAILS.  OF COURSE IT SHOULD

# single.integral.2 <- 
#   err_ngm( matrix( 0, 1, 1), matrix( k.ss, nrow=2, ncol=1), nodes, 
#            params.full.dep, coeff.old, opt$n.exog, opt$n.endog, 
#            params.full.dep$rho, opt$n.quad, opt$N, upper, lower, FALSE, wts, TRUE )
# 
# many.pts <- 
#   eval_err( coeff.old, matrix( c( 0, k.ss, 0, k.ss ), nrow=1 ), 'ngm', opt$lags, 
#             params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
#             params.full.dep$sig.eps,  0, opt$N, upper, lower, FALSE, 
#             matrix(0,1,1,), TRUE, opt$n.quad )
# 
# err.on.set.old <- 
#   eval_err( coeff.old, endog.sim, 'ngm', opt$lags, 
#             params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
#             params.full.dep$sig.eps,  0, opt$N, upper, lower, FALSE, 
#             matrix(0,1,1,), TRUE, opt$n.quad )
# 
# err.on.set.new <- 
#   eval_err( coeff.new, endog.sim, 'ngm', opt$lags, 
#             params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
#             params.full.dep$sig.eps,  0, opt$N, upper, lower, FALSE, 
#             matrix(0,1,1,), TRUE, opt$n.quad )
#     # **THIS** is the crazy thing.  That the error on the alternate coeffieicnt 
#     # is smaller.  Is this a legit (just explosive?) solution.  Can I constrain 
#     # the lag on k to be positive (or better, just to be non-explosive - this is
#     # more general)

## SOMEWHERE I AM SETTING N_INTEG = 0


## Where is the problem? It must be somewhere in the minimization routine.  It 
## simply isn't giving the same coefficients as are being input.  So it has to 
## be this.  Well, why is that?  It isn't the integrand - that appears to be 
## working (see the check above).  Then what's left? 1) The integration (can I 
## try MC instead - or is this just too slow?) 2) The gradient (that looks ok) 
## 3) The set of X is too small (doesn't seem right - it's pretty big - or could
## it be to do with needless weighting of the outliers?  Can I fix that 
## somehow?) 4) [MOST LIKELY??] It is the scaling.  Need to check that these 
## aren't getting switched somewhere, because the simulations look good and the 
## regression coefficients appear fine.  This might entail a full check through 
## the polyomial code.  But that's probably a good idea anyway.  Check that the
## upper & lower arguments are always the same way round too.


## TO DO:
## Add the burn
## Extend to higher-order solutions
## Check how delta is working in the namings
## Extend to the non-full-dep case