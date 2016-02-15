###################################################################
# Code to run the two-country NGM model with controls
###################################################################


#### 0. PREPARATORY WORK ####
params.1 <- list( A=1, alpha=.3, delta=1, gamma=1, betta=.99, 
                  rho=0.5, sig.eps=.01 )
# NB: With i.i.d. shocks there is not a unique solution for the coefficients.
params.1$A <- ( 1 / params.1$betta - 1 + params.1$delta ) / params.1$alpha
k.ss <- ( ( 1 / params.1$betta - 1 + params.1$delta ) / 
            ( params.1$A * params.1$alpha )  ) ^ ( 1 / ( params.1$alpha - 1 ) )
# The steady state
sd.x<- sqrt( params.1$sig.eps / ( 1 - params.1$rho ^ 2 ) )
upper <- c(  3 * sd.x, k.ss * 1.4 )
lower <- c( -3 * sd.x, k.ss * 0.6 )
# Create the bounds
opt <- list( model='ngm', lags=1, n.exog=1, n.endog=1, N=1, cheby=FALSE, 
             upper = upper, lower=lower, quad=TRUE, n.quad=5, diff.tol=1e-05, 
             n.iter=20, burn=1000, kappa=25, n.sim=10000, eps = .3, delta=.02, 
             endog.init=k.ss, gain=1, reg.method=TRUE, sr=FALSE, adapt.gain=TRUE,
             n.cont=0, adapt.exp=10, inner.iter=10, inner.tol=1e-04, image=FALSE )
# Solution options
a.b.range <- upper - lower
coeff.init <- matrix( c( k.ss, params.1$alpha * a.b.range[2] / 2, 
                         k.ss * a.b.range[1] / 2 ), 3, 1)
# The linear-approx solution evaluated @ the steady state:
#   k' = alpha * betta * exp(x) * A * k ^ alpha
#     ~= k.ss + k.ss * x + alpha * ( k - k.ss )
#   NB: The range is already centered on k.ss, so can just use 
#       these coeffs

sol.1 <- sol.iterate( coeff.init, opt, params.1 )

#### 1. TWO-COUNTRY VERSION: PERFECT DEPRECIATION, NO CONTROLS ####
params.2 <- list( A=params.1$A, alpha=.3, delta=1, gamma=1, betta=.99, 
                  rho=c(.5,.5), sig.eps=c(.01,.01) )
upper.2 <- c(  3 * sd.x,  3 * sd.x, k.ss * 1.4, k.ss * 1.4 )
lower.2 <- c( -3 * sd.x, -3 * sd.x, k.ss * 0.6, k.ss * 0.6 )
a.b.range.2 <- upper.2 - lower.2
opt.2 <- list( model='ngm2', lags=1, n.exog=2, n.endog=2, N=1, cheby=FALSE, 
               upper = upper.2, lower=lower.2, quad=TRUE, n.quad=5, diff.tol=1e-05, 
               n.iter=30, burn=1000, kappa=25, n.sim=10000, eps = .5, delta=.05, 
               endog.init=c(k.ss,k.ss), gain=1, reg.method=TRUE, sr=TRUE, adapt.gain=TRUE,
               n.cont=0, adapt.exp=10, inner.iter=10, inner.tol=1e-04, image=TRUE )

# > idx_create( 1, 4 )
# [,1] [,2] [,3] [,4]
# [1,]    0    0    0    0
# [2,]    0    0    0    1
# [3,]    0    0    1    0
# [4,]    0    1    0    0
# [5,]    1    0    0    0

coeff.init.2 <- matrix( c( sol.1$coeff[1], sol.1$coeff[1],
                           0, sol.1$coeff[2],
                           sol.1$coeff[2], 0,
                           0, sol.1$coeff[3],
                           sol.1$coeff[3], 0 ), 5, 2, byrow=TRUE )

sol.2 <- sol.iterate( coeff.init.2, opt.2, params.2 )
    # The solution
sol.check( sol.2, opt.2, params.2 )

n.irf <- 10
sim.irf <- 1000
irf.init <- irf_create( n.irf+1, sim.irf, opt$N, 0, params.2$rho, 
                        params.2$sig.eps, coeff.init.2, opt.2$upper, opt.2$lower, 
                        c(0,0,k.ss,k.ss), opt.2$n.endog, opt.2$n.exog, .01, FALSE )
irf <- irf_create( n.irf+1, sim.irf, opt$N, 0, params.2$rho, 
                   params.2$sig.eps, sol.2$coeff, opt.2$upper, opt.2$lower, 
                   c(0,0,k.ss,k.ss), opt.2$n.endog, opt.2$n.exog, .01, FALSE )
plot( 0:n.irf, irf.init[,3], type='l', lwd=2, col='blue' )
lines( 0:n.irf, irf.init[,4], type='l', lwd=2, col='blue', lty=2 )
lines( 0:n.irf, irf[,3], type='l', lwd=2, col='red' )
lines( 0:n.irf, irf[,4], type='l', lwd=2, col='red', lty=2 )
lines( 0:n.irf, irf[,1], lwd=2 )

#### 1. TWO-COUNTRY VERSION: PERFECT DEPRECIATION, WITH CONTROLS ####
c.ss <- params.2$A * k.ss ^ params.2$alpha - params.2$delta * k.ss
    # Steady state consumption
coeff.cont.init.2 <- matrix(
  c( c.ss, 
     .5 * ( 1 / params.2$betta * a.b.range.2[4] / 2 - sol.2$coeff[2] ),
     .5 * ( 1 / params.2$betta * a.b.range.2[3] / 2 - sol.2$coeff[2] ),
     .5 * ( params.1$A * k.ss ^ params.1$alpha * a.b.range[2] / 2 - sol.2$coeff[3] ),
     .5 * ( params.1$A * k.ss ^ params.1$alpha * a.b.range[1] / 2 - sol.2$coeff[3] ) ), 5, 1)
opt.2.cont <- list( model='ngm2.cont', lags=1, n.exog=2, n.endog=2, N=1, cheby=FALSE, 
               upper = upper.2, lower=lower.2, quad=TRUE, n.quad=5, diff.tol=1e-05, 
               n.iter=100, burn=1000, kappa=25, n.sim=10000, eps = .8, delta=.05, 
               endog.init=c(k.ss,k.ss), gain=.1, reg.method=TRUE, sr=TRUE, adapt.gain=TRUE,
               adapt.exp=15, n.cont=1, inner.iter=10, inner.tol=1e-04, image=TRUE )
sol.2.cont <- sol.iterate( sol.2$coeff, opt.2.cont, params.2, coeff.cont.init.2 )

#### 2. TWO-COUNTRY VERSION: PERFECT DEPRECIATION, WITH CONTROLS, NONLINEAR VERSION ####
init.3 <- matrix( 0, 15, 2 )
init.3[1, ] <- sol.2.cont$coeff[1, ]
init.3[2, ] <- sol.2.cont$coeff[2, ]
init.3[4, ] <- sol.2.cont$coeff[3, ]
init.3[7, ] <- sol.2.cont$coeff[4, ]
init.3[11, ] <- sol.2.cont$coeff[5, ]
init.cont.3 <- matrix( 0, 15, 1 )
init.cont.3[1,] <- sol.2.cont$coeff.cont[1,]
init.cont.3[2,] <- sol.2.cont$coeff.cont[2,]
init.cont.3[4,] <- sol.2.cont$coeff.cont[3,]
init.cont.3[7,] <- sol.2.cont$coeff.cont[4,]
init.cont.3[11,] <- sol.2.cont$coeff.cont[5,]
opt.3.cont <- opt.2.cont
opt.3.cont$N <- 2
opt.3.cont$inner.iter <- 30
sol.3.cont <- sol.iterate( init.3, opt.3.cont, params.2, init.cont.3 )
sol.check( sol.3.cont, opt.3.cont, params.2 )
