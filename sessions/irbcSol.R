##########################################################################################
# Code to solve an international RBC with restricted nominal asset markets
#
# 14feb2016
# Philip Barrett, Chicago
##########################################################################################

params <- list( alpha = .75, gamma = 1, P1.bar=1, P2.bar=2, betta=.99, 
                rho=c(.5,.5), sig.eps=c(.01,.01) )
    # Parameters
lower <- sd.x<- sqrt( params$sig.eps / ( 1 - params$rho ^ 2 ) )
upper <- c(  3 * sd.x, rep( 1, 2 ) )
lower <- -upper
    # Bounds

opt <- list( model='irbc', lags=1, n.exog=2, n.endog=2, N=1, cheby=FALSE,
             upper = upper, lower=lower, quad=TRUE, n.quad=5, diff.tol=1e-05,
             n.iter=3, burn=1000, kappa=25, n.sim=10000, eps = .8, delta=.05,
             endog.init=c(0, 0), gain=.1, reg.method=TRUE, sr=TRUE, adapt.gain=TRUE,
             adapt.exp=40, n.cont=11, inner.iter=10, inner.tol=1e-04, image=TRUE )

coeff.init <- matrix( c(  0,  0, 
                          0, .5,
                         .5,  0,
                          0, .5,
                         .5,  0 ), 5, 2, byrow=TRUE ) +  1e-05 * rnorm( 10 )

coeff.cont.init <- matrix( 0, 5, 11 )
coeff.cont.init[1, ] <- c( 1 / params$betta, 1 / params$betta, 1, 1, 1, 1, 1, 
                           params$alpha, params$alpha, 1, 1 )
coeff.cont.init <- coeff.cont.init + 1e-06 * rnorm( 55 )

sol <- sol.iterate( coeff.init, opt, params, coeff.cont.init, debug.flag = F )



