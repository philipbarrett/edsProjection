


params <- list( share = .75, gamma = 2, P1.bar=1, P2.bar=1, betta=.95,
                rho=c(.9,.9), sig.eps=c(.01,.01), eta=6 )
    # Model parameters
l.coeffs <- mod.gen(params)
    # The dynare solution
opt <- list( lags=1, n.exog=2, n.endog=2, n.fwd=4, n.cont=21, N=1, cheby=FALSE,
             upper = l.coeffs$alt.sol$upper, lower=l.coeffs$alt.sol$lower, quad=TRUE, 
             n.quad=3,  burn=1000, kappa=25, n.sim=10000, eps = .7, delta=.05, 
             endog.init=l.coeffs$mod$ys[c('NFA','af1')], 
             c.iter=400, c.tol=1e-07, c.gain=.25,
             k.iter=50, k.tol=2e-04, k.gain=.2,
             n.iter=1, n.tol=1e-05, n.gain=.25, 
             tol=1e-05, iter=40, model='irbc',
             sr=TRUE, adapt.gain=TRUE, adapt.exp=15, image=TRUE,
             sym.reg=TRUE )
    # Solution options

