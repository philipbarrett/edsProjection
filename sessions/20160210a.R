

#### 0. PREPARATORY WORK ####
params.1 <- list( A=1, alpha=.3, delta=1, gamma=1, betta=.99, 
                         rho=0, sig.eps=.01 )
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
             endog.init=k.ss, gain=1, reg.method=TRUE, sr=FALSE )
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

#### 1. PREPARATORY WORK ####
params.2 <- list( A=params.1$A, alpha=.3, delta=1, gamma=1, betta=.99, rho=c(0,0), sig.eps=c(.01,.01) )
upper.2 <- c(  3 * sd.x,  3 * sd.x, k.ss * 1.4, k.ss * 1.4 )
lower.2 <- c( -3 * sd.x, -3 * sd.x, k.ss * 0.6, k.ss * 0.6 )
opt.2 <- list( model='ngm2', lags=1, n.exog=2, n.endog=2, N=1, cheby=FALSE, 
                upper = upper.2, lower=lower.2, quad=TRUE, n.quad=5, diff.tol=1e-05, 
                n.iter=20, burn=1000, kappa=25, n.sim=10000, eps = .3, delta=.02, 
                endog.init=c(k.ss,k.ss), gain=1, reg.method=TRUE, sr=FALSE )

    # > idx_create( 1, 4 )
    # [,1] [,2] [,3] [,4]
    # [1,]    0    0    0    0
    # [2,]    0    0    0    1
    # [3,]    0    0    1    0
    # [4,]    0    1    0    0
    # [5,]    1    0    0    0

coeff.init.2 <- cbind( c( sol.1[1], 0, sol.1[2], 0, sol.1[3] ), 
                       c( sol.1[1], sol.1[2], 0, sol.1[3], 0 ) )
    # Set up the linear approximation to the two-dimensional problem.  See the 
    # output of the idx_create command above for an explanation of the weird
    # ordering.

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
