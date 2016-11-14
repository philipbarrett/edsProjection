## Dynare solution
params <- list( share = .86, gamma = 6, P1.bar=1, P2.bar=1, betta=.99,
                rho=c(.9,.9, .9, .9), sig.eps=c(.01,.01, .01, .01), eta=2, 
                theta=.05 )
params.0 <- params
params.0$theta <- 0
baseline <- mod.gen(params, check=FALSE)
baseline.0 <- mod.gen(params.0, check=FALSE)

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
             upper = upper, lower=lower, quad=TRUE, n.quad=2,  burn=1000,
             kappa=25, n.sim=10000, eps = 1.2, delta=.05, endog.init=endog.init, 
             c.iter=100, c.tol=1e-07, c.gain=.6,
             k.iter=20, k.tol=1e-07, k.gain=.6,
             n.iter=1, n.tol=1e-05, n.gain=.75, 
             tol=1e-05, iter=5, model='irbc',
             sr=TRUE, adapt.gain=TRUE, adapt.exp=15, image=TRUE,
             exog.names=exog.names, endog.names=endog.names, 
             cont.names=cont.names, fwd.vars=fwd.vars,
             sym.reg=FALSE, n.fwd=4, ys=baseline$ds.sol$ys, mono="m1" )
    # Symmetry regularization doesn't actually improve anything here

## Linear solution
# opt$sym.reg <- TRUE
# opt$l.sym.ave <- list( ave=c(1), sym=list(c(2,3),c(4,5),c(6,7) ) )
# opt$l.pairs <- list( c(1,2) )
# opt$l.pairs.cont <- list( c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12) )
sol.1 <- sol.irbc.iterate( baseline$ds.sol$coeff, opt, params, 
                                baseline$ds.sol$coeff.cont ) #, debug.flag = T )
opt$iter <- 6 # 2
sol.1.0 <- sol.irbc.iterate( sol.1$coeff, opt, params.0, sol.1$coeff.cont )

## Nonlinear solution
opt$iter <- 3 # 2
opt$N <- 2
# opt$n.gain <- 0.025
n.coeff <- idx_count( opt$N, n.exog + n.endog )
idx.coeff <- apply( idx_create(opt$N, n.exog + n.endog), 
                    1, function(x) sum(x) <=1 )
# opt$l.sym.ave <- list( ave=c(1,5,14,27), 
#                        sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), 
#                                  c(10,15), c(16,22), c(17,24), c(18,23),
#                                  c(19,26), c(20,25), c(21,28) ) )
coeff.init <- matrix(0,n.coeff,n.endog)
coeff.init.cont <- matrix(0,n.coeff,n.cont)
coeff.init[ idx.coeff, ] <- sol.1$coeff
coeff.init.cont[ idx.coeff, ] <- sol.1$coeff.cont
    # Set up the new initial guess
sol.2 <- sol.irbc.iterate( coeff.init, opt, params, coeff.init.cont )# , debug.flag = T )
opt$iter <- 1
sol.2.0 <- sol.irbc.iterate( sol.2$coeff, opt, params.0, sol.2$coeff.cont )

## Create simulations
n.sim <- 100000
baseline.sol <- list( params=params, opt=sol.1$opt, 
                      coeff=baseline$ds.sol$coeff, 
                      coeff.cont=baseline$ds.sol$coeff.cont )
sim.baseline <- sim.sol( baseline.sol, n.sim )
sim.exog <- sim.baseline[,1:n.exog]
sim.sol.1 <- sim.sol( sol.1, n.sim, sim.exog )
sim.sol.1.0 <- sim.sol( sol.1.0, n.sim, sim.exog )
sim.sol.2 <- sim.sol( sol.2, n.sim, sim.exog )
sim.sol.2.0 <- sim.sol( sol.2.0, n.sim, sim.exog )

## Analyze errors on **out-of-sample** simulation
extra.args <- list( n.fwd=opt$n.fwd, y1.ss=opt$ys['Y1'] )

# Baseline
pred <- contemp_eqns_irbc_grid( sim.baseline, opt$lags, params, 
                     n.exog, n.endog, n.cont, extra.args, opt$model )
colnames(pred) <- c( endog.names, cont.names )
pred[,fwd.vars] <- 
    euler_hat_grid( baseline.sol$coeff, baseline.sol$coeff.cont, 
                    sim.baseline, opt$lags, params, 
                    n.exog, n.endog, n.cont, opt$n.fwd, params$rho, 
                    params$sig.eps, 0, 1, upper, lower, opt$cheby, 
                    matrix(0,1,1,), TRUE, opt$n.quad, opt$model, opt$mono )
baseline.err <- sim.baseline[,c( endog.names, cont.names )] - pred

# Linear
pred <- contemp_eqns_irbc_grid( sim.sol.1, opt$lags, params, 
                      n.exog, n.endog, n.cont, extra.args, opt$model )
colnames(pred) <- c( endog.names, cont.names )
pred[,fwd.vars] <- 
  euler_hat_grid( sol.1$coeff, sol.1$coeff.cont, 
                  sim.sol.1, opt$lags, params, 
                  n.exog, n.endog, n.cont, opt$n.fwd, params$rho, 
                  params$sig.eps, 0, 1, upper, lower, opt$cheby, 
                  matrix(0,1,1,), TRUE, opt$n.quad, opt$model, opt$mono )
sol.1.err <- sim.sol.1[,c( endog.names, cont.names )] - pred

# Linear, theta=0
pred <- contemp_eqns_irbc_grid( sim.sol.1.0, opt$lags, params.0, 
                                n.exog, n.endog, n.cont, extra.args, opt$model )
colnames(pred) <- c( endog.names, cont.names )
pred[,fwd.vars] <- 
  euler_hat_grid( sol.1.0$coeff, sol.1.0$coeff.cont, 
                  sim.sol.1.0, opt$lags, params, 
                  n.exog, n.endog, n.cont, opt$n.fwd, params$rho, 
                  params$sig.eps, 0, 1, upper, lower, opt$cheby, 
                  matrix(0,1,1,), TRUE, opt$n.quad, opt$model, opt$mono )
sol.1.0.err <- sim.sol.1.0[,c( endog.names, cont.names )] - pred

# Quadratic
pred <- contemp_eqns_irbc_grid( sim.sol.2, opt$lags, params, 
                      n.exog, n.endog, n.cont, extra.args, opt$model )
colnames(pred) <- c( endog.names, cont.names )
pred[,fwd.vars] <- 
  euler_hat_grid( sol.2$coeff, sol.2$coeff.cont, 
                  sim.sol.2, opt$lags, params, 
                  n.exog, n.endog, n.cont, opt$n.fwd, params$rho, 
                  params$sig.eps, 0, 2, upper, lower, opt$cheby, 
                  matrix(0,1,1,), TRUE, opt$n.quad, opt$model, opt$mono )
sol.2.err <- sim.sol.2[,c( endog.names, cont.names )] - pred

# Quadratic, theta=0
pred <- contemp_eqns_irbc_grid( sim.sol.2.0, opt$lags, params, 
                                n.exog, n.endog, n.cont, extra.args, opt$model )
colnames(pred) <- c( endog.names, cont.names )
pred[,fwd.vars] <- 
  euler_hat_grid( sol.2.0$coeff, sol.2.0$coeff.cont, 
                  sim.sol.2.0, opt$lags, params, 
                  n.exog, n.endog, n.cont, opt$n.fwd, params$rho, 
                  params$sig.eps, 0, 2, upper, lower, opt$cheby, 
                  matrix(0,1,1,), TRUE, opt$n.quad, opt$model, opt$mono )
sol.2.0.err <- sim.sol.2.0[,c( endog.names, cont.names )] - pred


# Euler
euler.bias.pc <- rbind( local=apply(baseline.err[,c(fwd.vars,'B22')],2,mean),
                         global.1=apply(sol.1.err[,c(fwd.vars,'B22')],2,mean),
                         global.2=apply(sol.2.err[,c(fwd.vars,'B22')],2,mean),
                         global.1.0=apply(sol.1.0.err[,c(fwd.vars,'B22')],2,mean),
                         global.2.0=apply(sol.2.0.err[,c(fwd.vars,'B22')],2,mean) ) * 100
euler.aad.l10 <- rbind( local=apply(log(abs(baseline.err[,c(fwd.vars,'B22')]), 10),2,mean),
                        global.1=apply(log(abs(sol.1.err[,c(fwd.vars,'B22')]), 10),2,mean),
                        global.2=apply(log(abs(sol.2.err[,c(fwd.vars,'B22')]), 10),2,mean),
                        global.1.0=apply(log(abs(sol.1.0.err[,c(fwd.vars,'B22')]), 10),2,mean),
                        global.2.0=apply(log(abs(sol.2.0.err[,c(fwd.vars,'B22')]), 10),2,mean) )

## Check BS/UIP stats
l.sim <- list( local=sim.baseline, global.1=sim.sol.1, global.2=sim.sol.2,
               global.1.0=sim.sol.1.0, global.2.0=sim.sol.2.0 )
bs.log <- sapply( l.sim, 
    function(sim) cor( diff(sim[,'C1']-sim[,'C2']), diff(sim[,'Q']) ) )
uip.coeff <- sapply( l.sim, 
    function(sim) lm( diff(sim[,'E']) ~ (sim[,'R_1']-sim[,'R_2'])[-n.sim] )$coeff[2] )
names(uip.coeff) <- names(bs.log)


## Wrap as function


## Run for various parameters





