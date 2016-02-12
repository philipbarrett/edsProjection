

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
             endog.init=k.ss, gain=1, reg.method=TRUE, sr=FALSE, adapt.gain=TRUE )
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

#### 1. TWO-COUNTRY VERSION: PERFECT DEPRECIATION ####
params.2 <- list( A=params.1$A, alpha=.3, delta=1, gamma=1, betta=.99, rho=c(.5,.5), sig.eps=c(.01,.01) )
upper.2 <- c(  3 * sd.x,  3 * sd.x, k.ss * 1.4, k.ss * 1.4 )
lower.2 <- c( -3 * sd.x, -3 * sd.x, k.ss * 0.6, k.ss * 0.6 )
opt.2 <- list( model='ngm2', lags=1, n.exog=2, n.endog=2, N=1, cheby=FALSE, 
                upper = upper.2, lower=lower.2, quad=TRUE, n.quad=5, diff.tol=1e-05, 
                n.iter=30, burn=1000, kappa=25, n.sim=10000, eps = .3, delta=.02, 
                endog.init=c(k.ss,k.ss), gain=1, reg.method=TRUE, sr=FALSE, adapt.gain=TRUE )

    # > idx_create( 1, 4 )
    # [,1] [,2] [,3] [,4]
    # [1,]    0    0    0    0
    # [2,]    0    0    0    1
    # [3,]    0    0    1    0
    # [4,]    0    1    0    0
    # [5,]    1    0    0    0

coeff.init.2 <- matrix( c( sol.1[1], sol.1[1],
                           0, sol.1[2],
                           sol.1[2], 0,
                           0, sol.1[3],
                           sol.1[3], 0 ), 5, 2, byrow=TRUE )

exog <- matrix( 0, 2, 2 )
endog <- matrix( k.ss, 2, 2 )
exog.lead <- rep( 0, 2 )

endog_update( exog.lead, endog[1,], coeff.init.2, opt.2$n.exog, opt.2$n.endog, 
              opt.2$N, opt.2$upper, opt.2$lower, opt.2$cheby )
    # Just checking that  this works
integrand_ngm_2( exog, endog, exog.lead, params.2, coeff.init.2, opt.2$n.exog, opt.2$n.endog, 
                opt.2$N, opt.2$upper, opt.2$lower )
nodes.wts <- quad_nodes_weights( opt.2$n.quad, opt.2$n.exog, params.2$sig.eps, c( 0, 0 ) )

euler_hat_ngm_2( exog, endog, nodes.wts$nodes, params.2, coeff.init.2, opt.2$n.exog, 
                 opt.2$n.endog, params.2$rho, opt.2$n.quad ^ opt.2$n.exog, opt.2$N, opt.2$upper, opt.2$lower, 
                 opt.2$cheby, nodes.wts$weights, TRUE )

microbenchmark( euler_hat_ngm_2( exog, endog, nodes.wts$nodes, params.2, coeff.init.2, opt.2$n.exog, 
                                 opt.2$n.endog, params.2$rho, opt.2$n.quad ^ opt.2$n.exog, opt.2$N, opt.2$upper, opt.2$lower, 
                                 opt.2$cheby, nodes.wts$weights ) )

sol.2 <- sol.iterate( coeff.init.2, opt.2, params.2 )
    # The solution

n.irf <- 10
sim.irf <- 1000
irf.init <- irf_create( n.irf+1, sim.irf, opt$N, 0, params.2$rho, 
                        params.2$sig.eps, coeff.init.2, opt.2$upper, opt.2$lower, 
                        c(0,0,k.ss,k.ss), opt.2$n.endog, opt.2$n.exog, .01, FALSE )
irf <- irf_create( n.irf+1, sim.irf, opt$N, 0, params.2$rho, 
                   params.2$sig.eps, sol.2, opt.2$upper, opt.2$lower, 
                   c(0,0,k.ss,k.ss), opt.2$n.endog, opt.2$n.exog, .01, FALSE )
plot( 0:n.irf, irf.init[,3], type='l', lwd=2, col='blue' )
lines( 0:n.irf, irf.init[,4], type='l', lwd=2, col='blue', lty=2 )
lines( 0:n.irf, irf[,3], type='l', lwd=2, col='red' )
lines( 0:n.irf, irf[,4], type='l', lwd=2, col='red', lty=2 )
lines( 0:n.irf, irf[,1], lwd=2 )

#### 2. TWO-COUNTRY VERSION: PARTIAL DEPRECIATION ####
params.2 <- list( A=1, alpha=.3, delta=.6, gamma=1, betta=.99, rho=c(.95,.95), 
                  sig.eps=c(.01,.01) )
params.2$A <- ( 1 / params.2$betta - 1 + params.2$delta ) / params.2$alpha
sd.x <- sqrt( params.2$sig.eps / ( 1 - params.2$rho ^ 2 ) )
upper.2 <- c(  3 * sd.x, k.ss * 1.4, k.ss * 1.4 )
lower.2 <- c( -3 * sd.x, k.ss * 0.6, k.ss * 0.6 )
opt.2 <- list( model='ngm2', lags=1, n.exog=2, n.endog=2, N=1, cheby=FALSE, 
               upper = upper.2, lower=lower.2, quad=TRUE, n.quad=5, diff.tol=1e-05, 
               n.iter=100, burn=1000, kappa=25, n.sim=10000, eps = .3, delta=.02, 
               endog.init=c(k.ss,k.ss), gain=.1, reg.method=TRUE, sr=FALSE, adapt.gain=TRUE,
               adapt.exp=10 )

sol.del.1 <- sol.iterate( sol.2, opt.2, params.2 )
# sol.del.1 <- sol.iterate( sol.del.1, opt.2, params.2 )
err.del.1 <- sol.check( sol.del.1, opt.2, params.2 )

#### 3. TWO-COUNTRY VERSION: VERY PARTIAL DEPRECIATION ####
params.3 <- list( A=1, alpha=.3, delta=.025, gamma=1, betta=.99, rho=c(.95,.95), 
                  sig.eps=c(.01,.01) )
params.3$A <- ( 1 / params.3$betta - 1 + params.3$delta ) / params.3$alpha
opt.3 <- opt.2
opt.3$adapt.exp <- 50
opt.3$diff.tol <- 2e-04
sol.del.2 <- sol.iterate( sol.del.1, opt.3, params.3 )
err.del.2 <- sol.check( sol.del.2, opt.3, params.3 )


#### 4. TWO-COUNTRY VERSION: NONLINEAR SOLUTION ####
opt.2.nl <- opt.2
opt.2.nl$N <- 2
opt.2.nl$n.iter <- 10
opt.2.nl$sr <- TRUE
opt.2.nl$eps <- .5
opt.2.nl$delta <- .05
opt.2.nl$gain <- .75
opt.2.nl$n.iter <- 40
init.del.1.nl <- matrix( 0, idx_count( opt.2.nl$N, opt.2.nl$n.endog + opt.2.nl$n.exog ), opt.2.nl$n.endog )
init.del.1.nl[1,] <- sol.del.1[1,]
init.del.1.nl[2,] <- sol.del.1[2,]
init.del.1.nl[4,] <- sol.del.1[3,]
init.del.1.nl[7,] <- sol.del.1[4,]
init.del.1.nl[11,] <- sol.del.1[5,]
sol.del.1.nl <- sol.iterate( init.del.1.nl, opt.2.nl, params.2 )
err.del.1.nl <- sol.check( sol.del.1.nl, opt.2.nl, params.2 )

n.irf <- 40
irf.del.1 <- irf_create( n.irf+1, sim.irf, opt.2$N, 0, params.2$rho, 
                        params.2$sig.eps, sol.del.1, opt.2$upper, opt.2$lower, 
                        c(0,0,k.ss,k.ss), opt.2$n.endog, opt.2$n.exog, .05, FALSE )
irf.del.1.nl <- irf_create( n.irf+1, sim.irf, opt.2.nl$N, 0, params.2$rho, 
                         params.2$sig.eps, sol.del.1.nl, opt.2$upper, opt.2$lower, 
                         c(0,0,k.ss,k.ss), opt.2$n.endog, opt.2$n.exog, .05, FALSE )
plot( c(0, n.irf), range( irf.del.1, irf.del.1.nl ), type='n' )
lines( 0:n.irf, irf.del.1[,3], lwd=2, col='blue' )
lines( 0:n.irf, irf.del.1[,4], type='l', lwd=2, col='blue', lty=2 )
lines( 0:n.irf, irf.del.1.nl[,3], type='l', lwd=2, col='red' )
lines( 0:n.irf, irf.del.1.nl[,4], type='l', lwd=2, col='red', lty=2 )
lines( 0:n.irf, irf.del.1[,1], lwd=2 )
abline( h=0 )
