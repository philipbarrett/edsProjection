## ITT: More pushing of the calibration with high gamma, eta

library(xtable)
library(scales)

### A Very simple parameter setting
params <- list( share = .75, gamma = 2.5, P1.bar=1, P2.bar=1, betta=.95,
                rho=c(.75,.75), sig.eps=c(.01,.01), eta=1.5 )

sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
upper <- c(  3 * sd.x, rep( .5, 2 ) )
lower <- -upper

opt <- list( lags=1, n.exog=2, n.endog=2, n.cont=13, N=1, cheby=FALSE,
             upper = upper, lower=lower, quad=TRUE, n.quad=3,  burn=1000,
             kappa=25, n.sim=10000, eps = .7, delta=.05, endog.init=c(0, 0), 
             c.iter=400, c.tol=1e-07, c.gain=.25,
             k.iter=50, k.tol=2e-04, k.gain=.2,
             n.iter=1, n.tol=1e-05, n.gain=.25, 
             tol=1e-05, iter=40, model='irbc',
             sr=TRUE, adapt.gain=TRUE, adapt.exp=15, image=TRUE,
             sym.reg=TRUE )

coeff.init <- matrix( .02, 5, 2 )
coeff.init[2,] <- c( -.01,.02)
coeff.init[3,] <- c( .02,-.01)
coeff.init[4,] <- c( -.02,.01)
coeff.init[5,] <- c( .01,-.02)

coeff.cont.init <- matrix( 0, 5, opt$n.cont )
c.ss <- params$share ^ params$share * ( 1 - params$share ) ^ ( 1 - params$share )
p.ss <- 1 / c.ss
coeff.cont.init[1, ] <- log( c( c.ss, c.ss, 1 / params$betta, 1 / params$betta, 
                                params$share, params$share, 1 - params$share, 1 - params$share,
                                p.ss, p.ss, 1, 1, 1 ) )

coeff.cont.init[2, ] <- c( -.1, .5, -.2, .5, rep( .1, 9 ) ) / 10
coeff.cont.init[3, ] <- c(  .5,-.1,  .5,-.2, rep( .1, 9 ) ) / 10
coeff.cont.init[4, ] <- coeff.cont.init[2, ]
coeff.cont.init[5, ] <- coeff.cont.init[3, ]


opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$l.pairs <- list( c(1,2) )
opt$l.pairs.cont <- list( c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12) )
opt$sym.reg <- TRUE


#### 1. BASELINE ####

# Rough solution
opt$iter <- 6
sol.irbc.1.1 <- sol.irbc.iterate( coeff.init, opt, params, coeff.cont.init, debug.flag = F )

# Finer solution
opt$n.gain <- .05
opt$iter <- 12
sol.irbc.1.1 <- sol.irbc.iterate( sol.irbc.1.1$coeff, opt, params, sol.irbc.1.1$coeff.cont )
rep.irbc.1.1 <- report.data( sol.irbc.1.1 )
report.create( sol.irbc.1.1, rep.irbc.1.1 )

## Now: Find a nonlinear solution ##
opt$N <- 2
opt$n.gain <- .1
coeff.init <- matrix( 0, 15, 2 )
coeff.init[ c(1,2,4,7,11), ] <- sol.irbc.1.1$coeff
coeff.cont.init <- matrix( 0, 15, 13 )
coeff.cont.init[ c(1,2,4,7,11), ] <- sol.irbc.1.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 50
    # Try a nonlinear solution.
    # Need many iterations for convergence
sol.irbc.1.2 <- sol.irbc.iterate( coeff.init, opt, params, coeff.cont.init )
rep.irbc.1.2 <- report.data( sol.irbc.1.2 )
report.create( sol.irbc.1.2, rep.irbc.1.2 )
print( paste0( "err = ", round( max(apply( abs( rep.irbc.1.2$err ), 2, mean )) * 100, 4), "pp" ) )
