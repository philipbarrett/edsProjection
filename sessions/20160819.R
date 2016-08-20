params <- list( share = .75, gamma = 2, P1.bar=1, P2.bar=1, betta=.95,
                rho=c(.9,.9), sig.eps=c(.01,.01), eta=5 )
# Model parameters
l.coeffs <- mod.gen(params)
# The dynare solution
exog.names <- c('A1','A2')
endog.names <- c( 'NFA', 'Z1', 'Z2' )
cont.names <- c( 'C1', 'C2', 'rb1', 'rb2', 'X11', 'X22', 'X12', 'X21', 
                 'P1', 'P2', 'P11', 'P22', 'P12', 'P21', 'E', 'Q', 'af1',
                 'Y1', 'Y2', 'cd', 'cg' )
# Variable names
opt <- list( lags=1, n.exog=2, n.endog=3, n.fwd=4, n.cont=21, N=1, cheby=FALSE,
             upper = l.coeffs$ds.sol$upper, lower=l.coeffs$ds.sol$lower, quad=TRUE, 
             n.quad=3,  burn=1000, kappa=25, n.sim=10000, eps = 1, delta=.05, 
             endog.init=l.coeffs$mod$ys[c('NFA','Z1', 'Z2')], 
             fwd.vars=c('Z1','Z2','af1','Q'),
             exog.names=exog.names, endog.names=endog.names, cont.names=cont.names,
             c.iter=400, c.tol=1e-07, c.gain=.8,
             k.iter=50, k.tol=2e-04, k.gain=.2,
             n.iter=1, n.tol=1e-05, n.gain=.2, 
             tol=1e-04, iter=1, model='ds',
             sr=TRUE, adapt.gain=TRUE, adapt.exp=15, image=TRUE,
             sym.reg=FALSE, ys=l.coeffs$mod$ys )
# Solution options
opt$N <- 1
dyn.check <- sol.irbc.check(l.coeffs$ds.sol, params, opt)

opt$N <- 2
idx.2 <- idx_create( opt$N, opt$n.endog+opt$n.exog )
coeff.init <- matrix( 0, nrow(idx.2), opt$n.endog )
coeff.init[ which( apply( idx.2, 1, sum ) < 2 ), ] <- l.coeffs$ds.sol$coeff
coeff.init.cont <- matrix( 0,  nrow(idx.2), opt$n.cont )
coeff.init.cont[ which( apply( idx.2, 1, sum ) < 2 ), ] <- l.coeffs$ds.sol$coeff.cont
opt$tol <- 1e-05
opt$iter <- 4
opt$n.gain <- .005

sol.2 <- sol.irbc.iterate( coeff.init, opt, params,  coeff.init.cont, FALSE )

nl.check <- sol.irbc.check(sol.2)


### SIMULATE ###
exog.sim <- sapply( 1:opt$n.exog, function(i) ar1_sim( opt$n.sim, params$rho[i], 
                                                       params$sig.eps[i] ) )
nl.sim <- endog_sim( opt$n.sim, exog.sim, sol.2$coeff, 2, opt$upper,
                     opt$lower, opt$endog.init, FALSE, 1, 0, TRUE )
nl.cont.sim <- cont_sim( nl.sim, sol.2$coeff.cont, 2, opt$n.endog, opt$n.exog, 
                         sol.2$opt$upper, sol.2$opt$lower, FALSE )
colnames(nl.cont.sim) <- cont.names
    # The nonlinear simulation
dyn.sim <- endog_sim( opt$n.sim, exog.sim, l.coeffs$ds.sol$coeff, 1, opt$upper,
                     opt$lower, opt$endog.init, FALSE, 1, 0, TRUE )
dyn.cont.sim <- cont_sim( dyn.sim, l.coeffs$ds.sol$coeff.cont, 1, opt$n.endog, opt$n.exog, 
                         sol.2$opt$upper, sol.2$opt$lower, FALSE )
colnames(dyn.cont.sim) <- cont.names


### IMPLUSE RESPONSES ###
U <- matrix(0,100,2)
U[,1] <- .01 * (params$rho[1] ^ (0:99))
    # The shock
irf.nl <- endog_sim( 100, U, sol.2$coeff, 2, l.coeffs$ds.sol$upper, 
                     l.coeffs$ds.sol$lower, l.coeffs$mod$ys[c('NFA','Z1','Z2')], 
                      FALSE, 1, 0 )
irf.nl.cont <- cont_sim( irf.nl, sol.2$coeff.cont, 2, opt$n.endog, opt$n.exog, 
                         sol.2$opt$upper, sol.2$opt$lower, FALSE )
colnames(irf.nl.cont) <- cont.names
    # Nonlinear IRFs
irf.dyn <- endog_sim( 100, U, l.coeffs$ds.sol$coeff, 1, l.coeffs$ds.sol$upper, 
                      l.coeffs$ds.sol$lower, l.coeffs$mod$ys[c('NFA','Z1','Z2')], 
                      FALSE, 1, 0 )
irf.dyn.cont <- cont_sim( irf.dyn, l.coeffs$ds.sol$coeff.cont, 1, opt$n.endog, opt$n.exog, 
                         sol.2$opt$upper, sol.2$opt$lower, FALSE )
colnames(irf.dyn.cont) <- cont.names


### Create chart with asset distributions