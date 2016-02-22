##########################################################################################
# Code to solve an international RBC with restricted nominal asset markets
#
# 14feb2016
# Philip Barrett, Chicago
##########################################################################################

params <- list( alpha = .85, gamma = 5, P1.bar=1, P2.bar=1, betta=.99,
                rho=c(.95,.95), sig.eps=c(.01,.01) )
    # Parameters
lower <- sd.x<- sqrt( params$sig.eps / ( 1 - params$rho ^ 2 ) )
upper <- c(  3 * sd.x, rep( 1, 2 ) )
lower <- -upper
    # Bounds

l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                   ave=c(1,5,14) )
l.pairs <- list( c(1,2) )
l.pairs.cont <- list( c(1,2), c(3,4), c(5,6), c(8,9) )
    # Symmetry constraints

opt <- list( model='irbc', lags=1, n.exog=2, n.endog=2, N=2, cheby=FALSE,
             upper = upper, lower=lower, quad=TRUE, n.quad=5, diff.tol=1e-05,
             n.iter=4, burn=1000, kappa=25, n.sim=10000, eps = .7, delta=.05,
             endog.init=c(0, 0), gain=.02, reg.method=TRUE, sr=TRUE, adapt.gain=TRUE,
             adapt.exp=20, n.cont=10, inner.iter=20, inner.tol=1e-03, image=TRUE,
             sym.reg=TRUE, l.sym.ave=l.sym.ave, l.pairs=l.pairs, 
             l.pairs.cont=l.pairs.cont )


coeff.init <- matrix( 0, 15, 2 )
# coeff.init[2,] <- c(-.1,.4)
# coeff.init[4,] <- c(.4,-.1)
coeff.init[7,] <- c(.2,.5)
coeff.init[11,] <- c(.5,.2)

coeff.cont.init <- matrix( 0, 15, opt$n.cont )
c.ss <- params$alpha ^ params$alpha * ( 1 - params$alpha ) ^ ( 1 - params$alpha )
p.ss <- 1 / c.ss
coeff.cont.init[1, ] <- c( 1 / params$betta, 1 / params$betta, 
                           c.ss, c.ss, p.ss, p.ss, 1, 
                           params$alpha, params$alpha, 1 - params$alpha )
coeff.cont.init[2, ] <- c( -.1, .5, -.2, .5, rep( .1, 6 ) )
coeff.cont.init[4, ] <- c(  .5,-.1,  .5,-.2, rep( .1, 6 ) )
coeff.cont.init[7, ] <- coeff.cont.init[2, ]
coeff.cont.init[11, ] <- coeff.cont.init[4, ]
  
sol <- sol.iterate( coeff.init, opt, params, coeff.cont.init ) #, debug.flag = T )
# sol.2 <- sol.iterate( sol$coeff, opt, params, sol$coeff.cont ) #, debug.flag = T )

# opt.2 <- opt
# opt.2$inner.iter <- 50
# opt.2$gain <- .05
# sol.2 <- sol.iterate( sol$coeff, opt.2, params, sol$coeff.cont )


sol.check(sol, opt, params)

exog.sim <- cbind( ar1_sim( 11000, params$rho[1], params$sig.eps[2] ),
                   ar1_sim( 11000, params$rho[2], params$sig.eps[2] ) )
endog.sim <- endog_sim( 11000, exog.sim, sol$coeff, opt$N, opt$upper, opt$lower, c(0,0) ) 
cont.sim <- cont_sim( endog.sim[-(1:1000), ], sol$coeff.cont, opt$N, opt$n.endog, opt$n.exog, opt$upper, opt$lower )

app <- log( cont.sim[-1,7] ) - log( cont.sim[-nrow(cont.sim),7] )
lm( app ~ log( cont.sim[-10000,1] / cont.sim[-10000,2] ) )

Q.12 <- cont.sim[,7] * cont.sim[,6] / cont.sim[,5]

cor( diff( log( cont.sim[,4]) ), diff( - log( Q.12 ) ) )

mean( endog.sim[-(1:1000),4] / ( endog.sim[-(1:1000),4] - endog.sim[-(1:1000),3] / cont.sim[,7] ) )
