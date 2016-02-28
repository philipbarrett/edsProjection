
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
             k.iter=50, k.tol=2e-04, k.gain=.5,
             n.iter=1, n.tol=1e-05, n.gain=.02, 
             tol=1e-05, iter=40,
             sr=TRUE, adapt.gain=TRUE, adapt.exp=15, image=TRUE,
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
coeff.cont.init[11, ] <- coeff.cont.init[3, ]

sol <- sol.irbc.iterate( coeff.init, opt, params, coeff.cont.init ) #, debug.flag = TRUE )
    # Solve the model
err <- euler_hat_grid( sol$coeff, sol$coeff.cont, sol$X.cont, opt$lags, params, opt$n.exog, 
                       opt$n.endog, opt$n.cont, params$rho, params$sig.eps, 0, opt$N, 
                       opt$upper, opt$lower, opt$cheby, matrix(0,1,1,), TRUE, opt$n.quad ) - 
  sol$X.cont[, c(2:3,11:12)]

n.irf <- 10
sim.irf <- 100000
irf.init <- irf_create( n.irf+1, sim.irf, opt$N, 1, params$rho, 
                        params$sig.eps, sol.3$coeff, opt$upper, opt$lower, 
                        c(0,0,0,0), opt$n.endog, opt$n.exog, .01, FALSE )
plot( 0:n.irf, irf.init[,3], type='l', lwd=2, col='blue' )
