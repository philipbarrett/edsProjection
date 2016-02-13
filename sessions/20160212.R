###################################################################
# Code to run the NGM model with controls
###################################################################


### CHECK FIRST THAT IT STILL WORKS WITHOUT CONTROLS
params.full <- list( A=1, alpha=.3, delta=1, gamma=1, betta=.99, 
                  rho=0.5, sig.eps=.01 )
    # NB: With i.i.d. shocks there is not a unique solution for the coefficients.
params.full$A <- ( 1 / params.full$betta - 1 + params.full$delta ) / params.full$alpha
k.ss <- ( ( 1 / params.full$betta - 1 + params.full$delta ) / 
            ( params.full$A * params.full$alpha )  ) ^ ( 1 / ( params.full$alpha - 1 ) )
# The steady state
sd.x<- sqrt( params.full$sig.eps / ( 1 - params.full$rho ^ 2 ) )
upper <- c(  3 * sd.x, k.ss * 1.4 )
lower <- c( -3 * sd.x, k.ss * 0.6 )
# Create the bounds
opt <- list( model='ngm', lags=1, n.exog=1, n.endog=1, N=1, cheby=FALSE, 
             upper = upper, lower=lower, quad=TRUE, n.quad=5, diff.tol=1e-05, 
             n.iter=20, burn=1000, kappa=25, n.sim=10000, eps = .3, delta=.02, 
             endog.init=k.ss, gain=1, reg.method=TRUE, sr=FALSE, adapt.gain=TRUE,
             n.cont=0, adapt.exp=10, image=FALSE )
# Solution options
a.b.range <- upper - lower
coeff.init <- matrix( c( k.ss, params.full$alpha * a.b.range[2] / 2, 
                         k.ss * a.b.range[1] / 2 ), 3, 1)
sol <- sol.iterate( coeff.init, opt, params.full )
sol.check( sol, opt, params.full )


### FULL DEPRECIATION CASE
params.1 <- list( A=1, alpha=.3, delta=1, gamma=1, betta=.99, 
                  rho=0.5, sig.eps=.01 )
    # NB: With i.i.d. shocks there is not a unique solution for the coefficients.
params.1$A <- ( 1 / params.1$betta - 1 + params.1$delta ) / params.1$alpha
k.ss <- ( ( 1 / params.1$betta - 1 + params.1$delta ) / 
            ( params.1$A * params.1$alpha )  ) ^ ( 1 / ( params.1$alpha - 1 ) )
    # The steady state
sd.x <- sqrt( params.1$sig.eps / ( 1 - params.1$rho ^ 2 ) )
upper <- c(  3 * sd.x, k.ss * 1.4 )
lower <- c( -3 * sd.x, k.ss * 0.6 )
    # Create the bounds
opt.1 <- list( model='ngm.cont', lags=1, n.exog=1, n.endog=1, N=1, cheby=FALSE, 
               upper = upper, lower=lower, quad=TRUE, n.quad=5, diff.tol=1e-05, 
               n.iter=10, burn=1000, kappa=25, n.sim=10000, eps = .3, delta=.02, 
               n.cont=1, endog.init=k.ss, gain=1, reg.method=TRUE, sr=TRUE, 
               adapt.gain=TRUE, adapt.exp = 10, image=FALSE, 
               cont.tol=1e-03, cont.iter=10 )
# Solution options
a.b.range <- upper - lower
coeff.init.1 <- sol$coeff

c.ss <- params.1$A * k.ss ^ params.1$alpha - params.1$delta * k.ss
    # Steady state consumption
coeff.cont.init. <- matrix( 
  c( c.ss, 
     1 / params.1$betta * a.b.range[2] / 2 - coeff.init[2], 
     params.1$A * k.ss ^ params.1$alpha * a.b.range[1] / 2 - coeff.init[3] ), 3, 1)
# So c ~= c.ss + dc/dk * ( k - k.ss ) + dc/dx * x
# And:
#     dc/dk|{ss} = 1 / betta - dk'/dk
#     dc/dx|{ss} = A * k.ss ^ alpha - dk'/dx

exog <- matrix( 0, 1, 1 )
endog <-  matrix( k.ss, 1, 1 )
endog.lag <- matrix( k.ss, 2, 1 )
exog.lead <- 0
cont <- c.ss

integrand_ngm_cont( exog, endog, cont, exog.lead, params.1, coeff.init.1, 
                    coeff.cont.init, opt.1$n.exog, opt.1$n.endog, 
                    opt.1$n.cont, opt.1$N, opt.1$upper, opt.1$lower )
nodes.wts <- quad_nodes_weights( opt.1$n.quad, opt.1$n.exog, params.1$sig.eps, 0 )
euler_hat_ngm_cont( exog, endog.lag, cont, nodes.wts$nodes, params.1, coeff.init.1, 
                    coeff.cont.init, opt.1$n.exog, opt.1$n.endog, opt.1$n.cont, 
                    params.1$rho, opt.1$n.quad ^ opt.1$n.exog, opt.1$N, opt.1$upper, 
                    opt.1$lower, opt.1$cheby, nodes.wts$weights, TRUE )
con_eqns_ngm( exog, endog.lag, params.1, coeff.init, coeff.init.1, 
              opt.1$n.exog, opt.1$n.endog, opt.1$n.cont, opt.1$N, opt.1$upper, opt.1$lower )
    # Should be c.ss.  It is.

X <- matrix( c( 0, k.ss, 0, k.ss, c.ss ), nrow=1 )
euler_hat( coeff.init.1, coeff.cont.init, X, 'ngm.cont', 1, params.1, opt.1$n.exog, 
           opt.1$n.endog, opt.1$n.cont, params.1$rho, params.1$sig.eps, 0, opt.1$N, 
           opt.1$upper, opt.1$lower, opt.1$cheby, matrix(0,1,1), opt.1$quad, opt.1$n.quad )
    # The Euler equation

X <- matrix( c( .01, k.ss, 0, k.ss, c.ss ), nrow=1 )
euler_hat( coeff.init.1, coeff.cont.init, X, 'ngm.cont', 1, params.1, opt.1$n.exog, 
           opt.1$n.endog, opt.1$n.cont, params.1$rho, params.1$sig.eps, 0, opt.1$N, 
           opt.1$upper, opt.1$lower, opt.1$cheby, matrix(0,1,1), opt.1$quad, opt.1$n.quad )

euler_hat( coeff.init.1, coeff.cont.init, X, 'ngm', 1, params.1, opt.1$n.exog, 
           opt.1$n.endog, 0, params.1$rho, params.1$sig.eps, 0, opt.1$N, 
           opt.1$upper, opt.1$lower, opt.1$cheby, matrix(0,1,1), opt.1$quad, opt.1$n.quad )

c.t <- con_eqns_ngm( exog + .01, endog.lag, params.1, coeff.init, coeff.init.1, 
              opt.1$n.exog, opt.1$n.endog, opt.1$n.cont, opt.1$N, opt.1$upper, opt.1$lower )
integrand_ngm_cont( exog + .01, endog, c.t, exog.lead, params.1, coeff.init.1, 
                    coeff.cont.init, opt.1$n.exog, opt.1$n.endog, 
                    opt.1$n.cont, opt.1$N, opt.1$upper, opt.1$lower )
integrand_ngm( exog + .01, endog.lag, exog.lead, params.1, coeff.init.1, 
                    opt.1$n.exog, opt.1$n.endog, 
                    opt.1$N, opt.1$upper, opt.1$lower )
    ## So this looks good for 

sol.1 <- sol.iterate( coeff.init.1, opt.1, params.1, coeff.cont.init )
sol.check( sol.1, opt.1, params.1 )
sol.check( sol, opt, params.full )

# ### TO DO ### 
# ### DONE: TICK! ###
# # Create sims for the two versions and then calculate the two
# # versions of k.hat that occur as a result.
# n.sim <- 400
# set.seed(333)
# exog.sim <- ar1_sim( n.sim, params.1$rho, params.1$sig.eps )
# endog.sim <- endog_sim( n.sim, exog.sim, sol, 1, upper, lower, k.ss, lag=1 )
# endog.sim.cont <- endog_sim( n.sim, exog.sim, coeff.init.1, 1, upper, lower, k.ss, lag=1 )
#     # Create the two endogenous simulations.  
# max( abs( endog.sim - endog.sim.cont ) )
#     # These are the same
# cont.sim <- cont_sim( endog.sim, coeff.cont.init, 1, 1, 1, 1, upper, lower, FALSE )
#     # The continuation simulation
# cont.sim.exact <- (1-params.1$delta) * endog.sim[,4] - endog.sim[,2] + 
#         exp( endog.sim[,1] ) * params.1$A * endog.sim[,4] ^ params.1$alpha
#     # The exact 
# coeff.cont.exact <- coeff_reg( cont.sim.exact, endog.sim[,c(1,4)], 1, lower, upper )


#### NONLINEAR SOLUTION ###
opt.2 <- opt.1
opt.2$N <- 2
opt.2$n.iter <- 30
coeff.init.2 <- matrix( c( sol.1$coeff[1:2], 0, sol.1$coeff[3], 0, 0 ), 6,1 )
coeff.cont.init.2 <- matrix( c( sol.1$coeff.cont[1:2], 0, sol.1$coeff.cont[3], 0, 0 ), 6,1 )
sol.2 <- sol.iterate( coeff.init.2, opt.2, params.1, coeff.cont.init.2 )
sol.check( sol.2, opt.2, params.1 )


#### NONLINEAR SOLUTION: PARTIAL DEPRECIATION ###
params.full <- list( A=1, alpha=.3, delta=.025, gamma=1, betta=.99, 
                rho=0.95, sig.eps=.01 )
params.full$A <- ( 1 / params.full$betta - 1 + params.full$delta ) / params.full$alpha
k.ss <- ( ( 1 / params.full$betta - 1 + params.full$delta ) / 
            ( params.full$A * params.full$alpha )  ) ^ ( 1 / ( params.full$alpha - 1 ) )
sd.x<- sqrt( params.full$sig.eps / ( 1 - params.full$rho ^ 2 ) )
upper <- c(  3 * sd.x, k.ss * 1.4 )
lower <- c( -3 * sd.x, k.ss * 0.6 )
opt.full <- list( model='ngm', lags=1, n.exog=1, n.endog=1, N=2, cheby=FALSE, 
                   upper = upper, lower=lower, quad=TRUE, n.quad=5, diff.tol=1e-04, 
                   n.iter=30, burn=1000, kappa=25, n.sim=10000, eps = .3, delta=.02, 
                   n.cont=0, endog.init=k.ss, gain=.1, reg.method=TRUE, sr=FALSE, 
                   adapt.gain=TRUE, adapt.exp = .5, image=FALSE, cont.tol=1e-03, cont.iter=10 )
a.b.range <- upper - lower
coeff.init.full <- matrix( c( k.ss, 
                              params.full$alpha * a.b.range[2] / 2, 0,
                              k.ss * a.b.range[1] / 2, 0, 0 ), 6, 1)
sol.full <- sol.iterate( coeff.init.full, opt.full, params.full )
    # Use the non-cont model to get an initial guess
sol.check( sol.full, opt.full, params.full )
    # Check the solution

opt.full$model <- 'ngm.cont'
opt.full$n.cont <- 1
opt.full$n.iter <- 50
opt.full$image <- TRUE
    # Use the model with controls
c.ss <- params.full$A * k.ss ^ params.full$alpha - params.full$delta * k.ss
    # Steady state consumption
coeff.cont.init.full <-   matrix( c( c.ss, 
                                     1 / params.full$betta * a.b.range[2] / 2 - sol.full$coeff[2], 0,
                                     params.full$A * k.ss ^ params.full$alpha * a.b.range[1] / 2 -
                                     sol.full$coeff[4], 0, 0 ), 6, 1)
    # Initialize with the full-depreciation solution
opt.full$adapt.exp <- 5
    # Damp overshooting
sol.full.cont <- sol.iterate( sol.full$coeff, opt.full, params.full, coeff.cont.init.full )
    # Use the full rule for the states and the linearized one for the controls
    # from the no-control solution
sol.full.cont <- sol.iterate( coeff.init.full, opt.full, params.full, coeff.cont.init.full )
    # Now use the bad starting point for the states and the no-control solution
    # guess for the control rule.
sol.full.hard <- sol.iterate( coeff.init.2, opt.full, params.full, coeff.cont.init.2 )
    # Last, use both rules from just the full-depreciation solution
sol.check( sol.full.hard, opt.full, params.full )

