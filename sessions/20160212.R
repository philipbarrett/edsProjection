###################################################################
# Code to run the NGM model with controls
###################################################################


### FULL DEPRECIATION CASE
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
opt <- list( model='ngm.cont', lags=1, n.exog=1, n.endog=1, N=1, cheby=FALSE, 
             upper = upper, lower=lower, quad=TRUE, n.quad=5, diff.tol=1e-05, 
             n.iter=20, burn=1000, kappa=25, n.sim=10000, eps = .3, delta=.02, 
             n.cont=1, endog.init=k.ss, gain=1, reg.method=TRUE, sr=FALSE, 
             adapt.gain=TRUE, adapt.exp = 10 )
# Solution options
a.b.range <- upper - lower
coeff.init <- matrix( c( k.ss, params.1$alpha * a.b.range[2] / 2, 
                         k.ss * a.b.range[1] / 2 ), 3, 1)

c.ss <- params.1$A * k.ss ^ params.1$alpha - params.1$delta * k.ss
    # Steady state consumption
coeff.cont.init <- matrix( 
  c( c.ss, 
     1 / params.1$betta * a.b.range[2] / 2 - coeff.init[2], 
     params.1$A * k.ss ^ params.1$alpha * a.b.range[1] / 2 - coeff.init[3] ), 3, 1)
# So c ~= c.ss + dc/dk * ( k - k.ss ) + dc/dx * x
# And:
#     dc/dk|{ss} = 1 / betta - dk'/dk
#     dc/dx|{ss} = A * k.ss ^ alpha - dk'/dx

exog <- 0
endog <- 1
endog.lag <- matrix( 1, 2, 1 )
exog.lead <- 0
cont <- c.ss

integrand_ngm_cont( exog, endog, cont, exog.lead, params.1, coeff.init, coeff.cont.init, 
                    opt$n.exog, opt$n.endog, opt$n.cont, opt$N, opt$upper, opt$lower )
nodes.wts <- quad_nodes_weights( opt$n.quad, opt$n.exog, params.1$sig.eps, 0 )
euler_hat_ngm_cont( exog, endog, cont, nodes.wts$nodes, params.1, coeff.init, 
                    coeff.cont.init, opt$n.exog, opt$n.endog, opt$n.cont, 
                    params.1$rho, opt$n.quad ^ opt$n.exog, opt$N, opt$upper, 
                    opt$lower, opt$cheby, nodes.wts$weights, TRUE )
con_eqns_ngm( exog, endog.lag, params.1, coeff.init, coeff.cont.init, 
              opt$n.exog, opt$n.endog, opt$n.cont, opt$N, opt$upper, opt$lower )
    # Should be c.ss.  It is.

X <- matrix( c( 0, k.ss, 0, k.ss, c.ss ), nrow=1 )
euler_hat( coeff.init, coeff.cont.init, X, 'ngm.cont', 1, params.1, opt$n.exog, 
           opt$n.endog, opt$n.cont, params.1$rho, params.1$sig.eps, 0, opt$N, 
           opt$upper, opt$lower, opt$cheby, matrix(0,1,1), opt$quad, opt$n.quad )
    # The Euler equation

sol <- sol.iterate( coeff.init, coeff.cont.init, opt, params.1, F )

opt.k <- list( model='ngm', lags=1, n.exog=1, n.endog=1, N=1, cheby=FALSE, 
             upper = upper, lower=lower, quad=TRUE, n.quad=5, diff.tol=1e-05, 
             n.iter=20, burn=1000, kappa=25, n.sim=10000, eps = .3, delta=.02, 
             n.cont=0, endog.init=k.ss, gain=1, reg.method=TRUE, sr=FALSE, 
             adapt.gain=TRUE, adapt.exp = 10 )
sol.k <- sol.iterate( coeff.init, 0 * coeff.init, opt.k, params.1, F )
