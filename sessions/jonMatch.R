## Parameters
rho <- diag(c(.94,.94,.9,.9))
sig <- diag(c(.000016, .000016, .000031, .000031))
params <- list( share = .86, gamma = 2, P1.bar=1, P2.bar=1, betta=.99,
                rho=rho, sig.eps=sig, eta=2.5, theta=0.05, mu=.55, xi=.44 )
baseline <- mod.gen(params, check=FALSE ) #, err.deets = TRUE )

## Convert to R-style solution object
# Extract things from baseline
upper <- baseline$ds.sol$upper
lower <- baseline$ds.sol$lower
n.endog <- ncol(baseline$ds.sol$coeff)
n.exog <- nrow(baseline$ds.sol$coeff) - n.endog - 1
n.cont <- ncol(baseline$ds.sol$coeff.cont)
endog.init <- tail(baseline$ds.sol$ys, n.endog)

exog.names <- c('shk1','shk2', 'P11', 'P22')
endog.names <- c( 'B11', 'B22' )
cont.names <- c( 'C1', 'C2', 'R_1', 'R_2', 'X11', 'X22', 'X12', 'X21', 
                 'P1', 'P2', 'P12', 'P21', 'E', 'Q', 'A1', 'A2', 'PN1', 'PN2', 
                 'PT1', 'PT2', 'CN1', 'CN2', 'CT1', 'CT2' )
fwd.vars <- c('B11', 'E', 'R_1', 'R_2')

# Set global solution options
opt <- list( lags=1, n.exog=n.exog, n.endog=n.endog, n.cont=n.cont, N=1, cheby=FALSE,
             upper = upper, lower=lower, quad=TRUE, n.quad=3,  burn=1000,
             kappa=25, n.sim=40000, eps = 1.0, delta=.02, endog.init=endog.init, 
             c.iter=100, c.tol=1e-07, c.gain=.8,
             k.iter=20, k.tol=1e-07, k.gain=.4,
             n.iter=1, n.tol=1e-05, n.gain=.4, 
             tol=1e-05, iter=10, model='irbc',
             sr=TRUE, adapt.gain=TRUE, adapt.exp=15, image=FALSE,
             exog.names=exog.names, endog.names=endog.names, 
             cont.names=cont.names, fwd.vars=fwd.vars,
             sym.reg=FALSE, n.fwd=4, ys=baseline$ds.sol$ys, mono="m1" )
baseline.sol <- list( params=params, opt=opt, 
                      coeff=baseline$ds.sol$coeff, 
                      coeff.cont=baseline$ds.sol$coeff.cont )
sol.1 <- sol.irbc.iterate( baseline$ds.sol$coeff, opt, params, 
                           baseline$ds.sol$coeff.cont ) #, debug.flag = TRUE )

## Create simulation
n.sim <- 100000
sim.baseline <- sim.sol( baseline.sol, n.sim )
sim.exog <- sim.baseline[,1:n.exog]
sim.sol.1 <- sim.sol( sol.1, n.sim, sim.exog )
extra.args <- list( n.fwd=opt$n.fwd, y1.ss=opt$ys['Y1'] )
baseline.err <- sim.err( sim.baseline, baseline.sol, extra.args )
sol.1.err <- sim.err( sim.sol.1, sol.1, extra.args )

l.sim <- list( sim.baseline, sim.sol.1)
bs <- sapply( l.sim, function(sim) cor( diff(sim[,'C1']-sim[,'C2']), diff(sim[,'Q']) ) )
bs.level <- sapply( l.sim, function(sim) cor( exp(sim[,'C1'])-exp(sim[,'C2']), exp(sim[,'Q']) ) )
uip.coeff <- sapply( l.sim, function(sim) lm( diff(sim[,'E']) ~ (sim[,'R_1']-sim[,'R_2'])[-n.sim] )$coeff[2] )

l.err <- list( baseline.err, sol.1.err)
err.sumy <- lapply( l.err, function(err) summary( err ) )
err.abs.sumy <- lapply( l.err, function(err) summary( abs(err) ) )
err.mean <- sapply( l.err, function(err) apply( err,2, mean ) )
err.abs <- sapply( l.err, function(err) apply( abs(err),2, mean ) )

