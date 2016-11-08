params <- list( share = .86, gamma = 2, P1.bar=1, P2.bar=1, betta=.95,
                rho=c(.9,.9, .9, .9), sig.eps=c(.01,.01, .0025, .0025), eta=2, 
                theta=.025 )

baseline <- mod.gen(params, check=FALSE)

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

opt <- list( lags=1, n.exog=n.exog, n.endog=n.endog, n.cont=n.cont, N=1, cheby=FALSE,
             upper = upper, lower=lower, quad=TRUE, n.quad=2,  burn=1000,
             kappa=25, n.sim=10000, eps = 1.0, delta=.05, endog.init=endog.init, 
             c.iter=100, c.tol=1e-07, c.gain=.6,
             k.iter=20, k.tol=1e-07, k.gain=.6,
             n.iter=1, n.tol=1e-05, n.gain=.75, 
             tol=1e-05, iter=6, model='irbc',
             sr=TRUE, adapt.gain=TRUE, adapt.exp=15, image=TRUE,
             exog.names=exog.names, endog.names=endog.names, 
             cont.names=cont.names, fwd.vars=fwd.vars,
             sym.reg=FALSE, n.fwd=4, ys=baseline$ds.sol$ys, mono="m1" )

sol.1 <- sol.irbc.iterate( baseline$ds.sol$coeff, opt, params, 
                                baseline$ds.sol$coeff.cont ) #, debug.flag = T )

euler.bias <- rbind( global=apply(sol.1$err,2,mean)[fwd.vars],
                     local=baseline$l.err$bias[fwd.vars] )
euler.aad <- log( rbind( global=apply(abs(sol.1$err),2,mean)[fwd.vars],
                         local=baseline$l.err$aad[fwd.vars] ), 10 )

## NOW: GENERATE THE SECOND ORDER SOLUTION





