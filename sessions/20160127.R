rm(list=ls())
load( '24feb2016.RData' )
    # Load the initial guess
params <- list( alpha = .85, gamma = 5, P1.bar=1, P2.bar=1, betta=.99,
                rho=c(.95,.95), sig.eps=c(.01,.01), eta=1 )
    # Define parameters
opt$iter <- 2

## 1. BASELINE SOLUTION ##
sol.irbc <- sol.irbc.iterate( sol.3$coeff, opt, params, sol.3$coeff.cont )
    # The solution

check <- sol.irbc.check( sol.irbc )

print( summary( abs( check$err ) ) )
plot( 1:(opt$n.cont+opt$n.endog), log( apply( abs( check$err ), 2, mean ), 10 ) )
plot( 1:(opt$n.cont+opt$n.endog), log( apply( abs( check$err ), 2, max ), 10 ) )


## 2. DECREASE SHOCK PERSISTENCE ##

params$rho <- c(.5,.5)
opt$iter <- 50
opt$n.gain <- .1
sol.irbc.2 <- sol.irbc.iterate( sol.irbc$coeff, opt, params, sol.irbc$coeff.cont )

check.2 <- sol.irbc.check( sol.irbc.2 )
print( summary( abs( check.2$err ) ) )
plot( 1:(opt$n.cont+opt$n.endog), log( apply( abs( check.2$err ), 2, mean ), 10 ) )
plot( 1:(opt$n.cont+opt$n.endog), log( apply( abs( check.2$err ), 2, max ), 10 ) )

plot( sol.irbc$X.cont[,1], sol.irbc$X.cont[,9], col=2, xlab='a1', ylab='c1' )
points( sol.irbc.2$X.cont[,1], sol.irbc.2$X.cont[,9] )
plot( sol.irbc$X.cont[,1], sol.irbc$X.cont[,21], col=2, xlab='a1', ylab='e12' )
points( sol.irbc.2$X.cont[,1], sol.irbc.2$X.cont[,21] )
plot( sol.irbc$X.cont[,9], sol.irbc$X.cont[,21], col=2, xlab='c1', ylab='e12' )
points( sol.irbc.2$X.cont[,9], sol.irbc.2$X.cont[,21] )


## 3. INCREASE RISK AVERSION ##

params$gamma <- 10
opt$iter <- 90
opt$n.gain <- .05
sol.irbc.3 <- sol.irbc.iterate( sol.irbc.2$coeff, opt, params, sol.irbc.2$coeff.cont )

check.3 <- sol.irbc.check( sol.irbc.3 )
print( summary( abs( check.3$err ) ) )
plot( 1:(opt$n.cont+opt$n.endog), log( apply( abs( check.3$err ), 2, mean ), 10 ) )
plot( 1:(opt$n.cont+opt$n.endog), log( apply( abs( check.3$err ), 2, max ), 10 ) )

plot( sol.irbc.2$X.cont[,1], sol.irbc.2$X.cont[,9], col=2, xlab='a1', ylab='c1' )
points( sol.irbc.3$X.cont[,1], sol.irbc.3$X.cont[,9] )
plot( sol.irbc.2$X.cont[,1], sol.irbc.2$X.cont[,21], col=2, xlab='a1', ylab='e12' )
points( sol.irbc.3$X.cont[,1], sol.irbc.3$X.cont[,21] )
plot( sol.irbc.2$X.cont[,9], sol.irbc.2$X.cont[,21], col=2, xlab='c1', ylab='e12' )
points( sol.irbc.3$X.cont[,9], sol.irbc.3$X.cont[,21] )
