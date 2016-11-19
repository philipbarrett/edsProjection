## Dynare solution
params <- list( share = .86, gamma = 8, P1.bar=1, P2.bar=1, betta=.99,
                rho=diag(rep(.9,4)), sig.eps=diag(c(.01,.01, .005, .005)^2), 
                eta=2, theta=.05 )
params$rho[2,1] <- params$rho[1,2] <- .025
params$sig.eps[2,1] <- params$sig.eps[1,2] <- .005 ^ 2
    # Add non-sphericality
baseline <- mod.gen(params, check=FALSE)

# Extract things from baseline
upper <- baseline$ds.sol$upper
lower <- baseline$ds.sol$lower
n.endog <- ncol(baseline$ds.sol$coeff)
n.exog <- nrow(baseline$ds.sol$coeff) - n.endog - 1
n.cont <- ncol(baseline$ds.sol$coeff.cont)
endog.init <- tail(baseline$ds.sol$ys, n.endog)

exog.names <- c('A1','A2', 'P11', 'P22')
endog.names <- c( 'B11', 'B22' )
cont.names <- c( 'C1', 'C2', 'R_1', 'R_2', 'X11', 'X22', 'X12', 'X21', 
                 'P1', 'P2', 'P12', 'P21', 'E', 'Q' )
fwd.vars <- c('B11', 'E', 'R_1', 'R_2')

## Set global solution options
opt <- list( lags=1, n.exog=n.exog, n.endog=n.endog, n.cont=n.cont, N=1, cheby=FALSE,
             upper = upper, lower=lower, quad=TRUE, n.quad=3,  burn=1000,
             kappa=25, n.sim=10000, eps = 1.0, delta=.02, endog.init=endog.init, 
             c.iter=100, c.tol=1e-07, c.gain=.8,
             k.iter=20, k.tol=1e-07, k.gain=.4,
             n.iter=2, n.tol=1e-05, n.gain=.3, 
             tol=1e-05, iter=5, model='irbc',
             sr=TRUE, adapt.gain=TRUE, adapt.exp=15, image=TRUE,
             exog.names=exog.names, endog.names=endog.names, 
             cont.names=cont.names, fwd.vars=fwd.vars,
             sym.reg=FALSE, n.fwd=4, ys=baseline$ds.sol$ys, mono="m1" )
baseline.sol <- list( params=params, opt=opt, 
                      coeff=baseline$ds.sol$coeff, 
                      coeff.cont=baseline$ds.sol$coeff.cont )

## Linear solution
l.sym.ave <- list( ave=c(1), sym=list(c(2,3),c(4,5),c(6,7) ) )
l.pairs <- list( c(1,2) )
l.pairs.cont <- list( c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12) )
sol.1 <- sol.irbc.iterate( baseline$ds.sol$coeff, opt, params, 
                           baseline$ds.sol$coeff.cont )
sol.1.sym <- list(
  coeff=m.sym.ave.pair( sol.1$coeff, l.sym.ave, l.pairs ),
  coeff.cont=m.sym.ave.pair( sol.1$coeff.cont, l.sym.ave, l.pairs.cont ),
  opt=opt, params=params )
    # Ex post regularization

## Simulations
n.sim <- 100000
baseline.sol <- list( params=params, opt=sol.1$opt, 
                      coeff=baseline$ds.sol$coeff, 
                      coeff.cont=baseline$ds.sol$coeff.cont )
sim.baseline <- sim.sol( baseline.sol, n.sim )
sim.exog <- sim.baseline[,1:n.exog]
sim.sol.1 <- sim.sol( sol.1, n.sim, sim.exog )
sim.sol.1.sym <- sim.sol( sol.1.sym, n.sim, sim.exog )

## Errors
extra.args <- list( n.fwd=opt$n.fwd, y1.ss=opt$ys['Y1'] )
baseline.err <- sim.err( sim.baseline, baseline.sol, extra.args )
err.sol.1 <- sim.err( sim.sol.1, sol.1, extra.args )
err.sol.1.sym <- sim.err( sim.sol.1.sym, sol.1.sym, extra.args )

## Diagnsotics
par.init <- par()
par(mfrow=c(2,1), mai = c(.5, 0.5, 0.3, 0.1))
    # Settings
l.err <- list( Local=baseline.err, Global1=err.sol.1 )
barplot(t(sapply( l.err, function(x) apply(abs(x), 2, mean) )), 
        beside=TRUE, main="Ave abs err")
barplot(t(sapply( l.err, function(x) apply(x, 2, mean) )), 
        beside=TRUE, main="Ave err")
    # Baseline + global linear
# barplot(apply(abs(err.sol.1.sym), 2, mean), main="Ave abs err")
# barplot(apply((err.sol.1.sym), 2, mean), main="Ave err")
#     # Symmetric gloabl linear - it pretty bad
par(par.init)
