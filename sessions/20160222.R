##########################################################################################
# Code to solve an international RBC with restricted nominal asset markets
#
# 22feb2016
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

opt <- list( lags=1, n.exog=2, n.endog=2, n.cont=13, N=1, cheby=FALSE,
             upper = upper, lower=lower, quad=TRUE, n.quad=5,  burn=1000,
             kappa=25, n.sim=10000, eps = .7, delta=.05, endog.init=c(0, 0), 
             x.iter=20, x.tol=1e-05, x.gain=.25, c.iter=20, c.tol=1e-05, c.gain=.05,
             n.iter=15, n.tol=1e-05, n.gain=.05, tol=1e-05, iter=4,
             sr=TRUE, adapt.gain=TRUE, adapt.exp=20,image=TRUE,
             sym.reg=TRUE, l.sym.ave=l.sym.ave, l.pairs=l.pairs, 
             l.pairs.cont=l.pairs.cont )


coeff.init <- matrix( 0, 5, 2 )
coeff.init[2,] <- c(-.1,.4)
coeff.init[3,] <- c(.4,-.1)
coeff.init[4,] <- c(.2,.5)
coeff.init[5,] <- c(.5,.2)

coeff.cont.init <- matrix( 0, 5, opt$n.cont )
c.ss <- params$alpha ^ params$alpha * ( 1 - params$alpha ) ^ ( 1 - params$alpha )
p.ss <- 1 / c.ss
coeff.cont.init[1, ] <- log( c( c.ss, c.ss, 1 / params$betta, 1 / params$betta, 
                                params$alpha, params$alpha, 1 - params$alpha, 1 - params$alpha,
                                p.ss, p.ss, 1, 1, 1 ) )
coeff.cont.init[2, ] <- c( -.1, .5, -.2, .5, rep( .1, 9 ) )
coeff.cont.init[3, ] <- c(  .5,-.1,  .5,-.2, rep( .1, 9 ) )
coeff.cont.init[4, ] <- coeff.cont.init[2, ]
coeff.cont.init[5, ] <- coeff.cont.init[3, ]

sol <- sol.irbc.iterate( coeff.init, opt, params, coeff.cont.init, debug.flag = TRUE )