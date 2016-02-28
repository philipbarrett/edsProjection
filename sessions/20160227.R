rm(list=ls())
load( '24feb2016.RData' )
    # Load the initial guess
params <- list( alpha = .85, gamma = 5, P1.bar=1, P2.bar=1, betta=.99,
                rho=c(.95,.95), sig.eps=c(.01,.01), eta=1 )
    # Define parameters
opt$iter <- 125
opt$k.iter <- 15
opt$n.gain <- .02

## 1. BASELINE SOLUTION ##
sol.irbc <- sol.irbc.iterate( sol.3$coeff, opt, params, sol.3$coeff.cont )
    # The solution
check <- sol.irbc.check( sol.irbc )

print( summary( abs( check$err ) ) )
plot( 1:(opt$n.cont+opt$n.endog), log( apply( abs( check$err ), 2, mean ), 10 ) )
plot( 1:(opt$n.cont+opt$n.endog), log( apply( abs( check$err ), 2, max ), 10 ) )

expect <- e_cont( sol.irbc$coeff.cont, check$endog.sim, opt$n.exog, opt$n.endog, 
                  opt$n.cont, sol.irbc$params$rho, sol.irbc$params$sig.eps, 0, opt$N, 
                  opt$upper, opt$lower, opt$cheby, matrix(0,1,1), TRUE, 5 )

n.irf <- 10
sim.irf <- 1000
irf <- irf_create( n.irf+1, sim.irf, opt$N, 1, sol.irbc$params$rho, 
                        sol.irbc$params$sig.eps, sol.irbc$coeff, opt$upper, opt$lower, 
                        c(0,0,0,0), opt$n.endog, opt$n.exog, .01, FALSE )
plot( 0:n.irf, irf[,3], type='l', lwd=2, col='blue' )


## 2. DECREASE SHOCK PERSISTENCE ##

params$rho <- c(.5,.5)
opt$iter <- 75
opt$n.gain <- .02
opt$k.gain <- .25
sol.irbc.2 <- sol.irbc.iterate( sol.irbc$coeff, opt, params, sol.irbc$coeff.cont )

check.2 <- sol.irbc.check( sol.irbc.2 )
print( summary( abs( check.2$err ) ) )
plot( 1:(opt$n.cont+opt$n.endog), log( apply( abs( check.2$err ), 2, mean ), 10 ) )
plot( 1:(opt$n.cont+opt$n.endog), log( apply( abs( check.2$err ), 2, max ), 10 ) )

expect.2 <- e_cont( sol.irbc.2$coeff.cont, check.2$endog.sim, opt$n.exog, opt$n.endog, 
                    opt$n.cont, sol.irbc.2$params$rho, sol.irbc.2$params$sig.eps, 0, opt$N, 
                    opt$upper, opt$lower, opt$cheby, matrix(0,1,1), TRUE, 5 )

plot( sol.irbc$X.cont[,1], sol.irbc$X.cont[,9], col=2, xlab='a1', ylab='c1' )
points( sol.irbc.2$X.cont[,1], sol.irbc.2$X.cont[,9] )
plot( sol.irbc$X.cont[,1], sol.irbc$X.cont[,21], col=2, xlab='a1', ylab='e12' )
points( sol.irbc.2$X.cont[,1], sol.irbc.2$X.cont[,21] )
plot( sol.irbc$X.cont[,9], sol.irbc$X.cont[,21], col=2, xlab='c1', ylab='e12' )
points( sol.irbc.2$X.cont[,9], sol.irbc.2$X.cont[,21] )


# ## 3. INCREASE RISK AVERSION ##
# 
# params$gamma <- 10
# opt$iter <- 90
# opt$n.gain <- .05
# sol.irbc.3 <- sol.irbc.iterate( sol.irbc.2$coeff, opt, params, sol.irbc.2$coeff.cont )
# 
# check.3 <- sol.irbc.check( sol.irbc.3 )
# print( summary( abs( check.3$err ) ) )
# plot( 1:(opt$n.cont+opt$n.endog), log( apply( abs( check.3$err ), 2, mean ), 10 ) )
# plot( 1:(opt$n.cont+opt$n.endog), log( apply( abs( check.3$err ), 2, max ), 10 ) )
# 
# plot( sol.irbc.2$X.cont[,1], sol.irbc.2$X.cont[,9], col=2, xlab='a1', ylab='c1' )
# points( sol.irbc.3$X.cont[,1], sol.irbc.3$X.cont[,9] )
# plot( sol.irbc.2$X.cont[,1], sol.irbc.2$X.cont[,21], col=2, xlab='a1', ylab='e12' )
# points( sol.irbc.3$X.cont[,1], sol.irbc.3$X.cont[,21] )
# plot( sol.irbc.2$X.cont[,9], sol.irbc.2$X.cont[,21], col=2, xlab='c1', ylab='e12' )
# points( sol.irbc.3$X.cont[,9], sol.irbc.3$X.cont[,21] )
# 
# 
# ## 3. DIAL DOWN RHO SOME MORE ##
# 
# params$gamma <- 5
# params$rho <- c( .1, .1 )
# opt$iter <- 90
# opt$n.gain <- .05
# sol.irbc.4 <- sol.irbc.iterate( sol.irbc.2$coeff, opt, params, sol.irbc.2$coeff.cont )
# 
# check.4 <- sol.irbc.check( sol.irbc.4 )
# print( summary( abs( check.4$err ) ) )
# plot( 1:(opt$n.cont+opt$n.endog), log( apply( abs( check.4$err ), 2, mean ), 10 ) )
# plot( 1:(opt$n.cont+opt$n.endog), log( apply( abs( check.4$err ), 2, max ), 10 ) )
# 
# plot( sol.irbc.2$X.cont[,1], sol.irbc.2$X.cont[,9], col=2, xlab='a1', ylab='c1' )
# points( sol.irbc.4$X.cont[,1], sol.irbc.4$X.cont[,9] )
# plot( sol.irbc.2$X.cont[,1], sol.irbc.2$X.cont[,21], col=2, xlab='a1', ylab='e12' )
# points( sol.irbc.4$X.cont[,1], sol.irbc.4$X.cont[,21] )
# plot( sol.irbc.2$X.cont[,9], sol.irbc.2$X.cont[,21], col=2, xlab='c1', ylab='e12' )
# points( sol.irbc.4$X.cont[,9], sol.irbc.4$X.cont[,21] )
