## Parameters
rho <- diag(rep(.975,4))
sig <- diag(c(.01,.01,.005,.005)^2)
params <- list( share = .9, gamma = 8, P1.bar=1, P2.bar=1, betta=.99,
                rho=rho, sig.eps=sig, eta=3, theta=0.025 )
params.0 <- params
params.0$theta <- 0
baseline <- mod.gen(params, check=FALSE ) #, err.deets = TRUE )
baseline.0 <- mod.gen(params.0, check=FALSE ) #, err.deets = TRUE )

## Convert to R-style solution object
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

# Set global solution options
opt <- list( lags=1, n.exog=n.exog, n.endog=n.endog, n.cont=n.cont, N=1, cheby=FALSE,
             upper = upper, lower=lower, quad=TRUE, n.quad=3,  burn=1000,
             kappa=25, n.sim=10000, eps = 1.0, delta=.02, endog.init=endog.init, 
             c.iter=100, c.tol=1e-07, c.gain=.8,
             k.iter=20, k.tol=1e-07, k.gain=.4,
             n.iter=2, n.tol=1e-05, n.gain=.1, 
             tol=1e-05, iter=4, model='irbc',
             sr=TRUE, adapt.gain=TRUE, adapt.exp=15, image=TRUE,
             exog.names=exog.names, endog.names=endog.names, 
             cont.names=cont.names, fwd.vars=fwd.vars,
             sym.reg=FALSE, n.fwd=4, ys=baseline$ds.sol$ys, mono="m1" )
baseline.sol <- list( params=params, opt=opt, 
                      coeff=baseline$ds.sol$coeff, 
                      coeff.cont=baseline$ds.sol$coeff.cont )

## Linear solution
sol.1 <- sol.irbc.iterate( baseline$ds.sol$coeff, opt, params, 
                           baseline$ds.sol$coeff.cont )
opt$iter <- 2
sol.1.0 <- sol.irbc.iterate( sol.1$coeff, opt, params.0, sol.1$coeff.cont )

## Quadratic solution
opt$N <- 2
opt$iter <- 2
n.coeff <- idx_count( opt$N, n.exog + n.endog )
idx.coeff <- apply( idx_create(opt$N, n.exog + n.endog), 
                    1, function(x) sum(x) <=1 )
coeff.init <- matrix(0,n.coeff,n.endog)
coeff.init.cont <- matrix(0,n.coeff,n.cont)
coeff.init[ idx.coeff, ] <- sol.1$coeff
coeff.init.cont[ idx.coeff, ] <- sol.1$coeff.cont
    # Set up the new initial guess
opt$iter <- 2
sol.2 <- sol.irbc.iterate( coeff.init, opt, params, coeff.init.cont )
coeff.init[ idx.coeff, ] <- sol.1.0$coeff
coeff.init.cont[ idx.coeff, ] <- sol.1.0$coeff.cont
sol.2.0 <- sol.irbc.iterate( coeff.init, opt, params.0, coeff.init.cont )


## Simulations
n.sim <- 100000
sim.baseline <- sim.sol( baseline.sol, n.sim )
sim.exog <- sim.baseline[,1:n.exog]
sim.sol.1 <- sim.sol( sol.1, n.sim, sim.exog )
sim.sol.1.0 <- sim.sol( sol.1.0, n.sim, sim.exog )
sim.sol.2 <- sim.sol( sol.2, n.sim, sim.exog )
sim.sol.2.0 <- sim.sol( sol.2.0, n.sim, sim.exog )

## Errors
extra.args <- list( n.fwd=opt$n.fwd, y1.ss=opt$ys['Y1'] )
baseline.err <- sim.err( sim.baseline, baseline.sol, extra.args )
err.sol.1 <- sim.err( sim.sol.1, sol.1, extra.args )
err.sol.1.0 <- sim.err( sim.sol.1.0, sol.1.0, extra.args )
err.sol.2 <- sim.err( sim.sol.2, sol.2, extra.args )
err.sol.2.0 <- sim.err( sim.sol.2.0, sol.2.0, extra.args )

# Outputs
l.sim <- list( local=sim.baseline, global.1=sim.sol.1 , global.2=sim.sol.2,
               global.1.0=sim.sol.1.0, global.2.0=sim.sol.2.0 )
l.err <- list( local=baseline.err, global.1=err.sol.1 , global.2=err.sol.2,
               global.1.0=err.sol.1.0, global.2.0=err.sol.2.0 )
bs.log <- sapply( l.sim, 
                  function(sim) cor( diff(sim[,'C1']-sim[,'C2']), diff(sim[,'Q']) ) )
uip.coeff <- sapply( l.sim, 
                     function(sim) lm( diff(sim[,'E']) ~ (sim[,'R_1']-sim[,'R_2'])[-n.sim] )$coeff[2] )
names(uip.coeff) <- names(bs.log)
err.ave.sumy <- sapply( l.err, function(x) apply(x, 2, mean) )
err.abs.sumy <- sapply( l.err, function(x) apply(abs(x), 2, mean) )



### NOW
# 1. PACKAGE AS FUNCTION
# 2. SET TO LOOP  (By Sunday)
# 3. Start to write up something qulaitative

