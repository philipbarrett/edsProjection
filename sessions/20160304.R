### A Very simple parameter setting

params <- list( alpha = .75, gamma = 2, P1.bar=1, P2.bar=1, betta=.95,
                rho=c(.25,.25), sig.eps=c(.03,.03), eta=1 )

sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
upper <- c(  3 * sd.x, rep( .1, 2 ) )
lower <- -upper

opt <- list( lags=1, n.exog=2, n.endog=2, n.cont=13, N=1, cheby=FALSE,
             upper = upper, lower=lower, quad=TRUE, n.quad=3,  burn=1000,
             kappa=25, n.sim=10000, eps = .7, delta=.05, endog.init=c(0, 0), 
             c.iter=400, c.tol=1e-07, c.gain=.25,
             k.iter=50, k.tol=2e-04, k.gain=.2,
             n.iter=1, n.tol=1e-05, n.gain=.25, 
             tol=1e-05, iter=40, model='irbc',
             sr=TRUE, adapt.gain=TRUE, adapt.exp=15, image=TRUE,
             sym.reg=TRUE )

coeff.init <- matrix( 0, 5, 2 )
coeff.init[2,] <- c( -.01,.02)
coeff.init[3,] <- c( .02,-.01)
coeff.init[4,] <- c( -.02,.01)
coeff.init[5,] <- c( .01,-.02)

coeff.cont.init <- matrix( 0, 5, opt$n.cont )
c.ss <- params$alpha ^ params$alpha * ( 1 - params$alpha ) ^ ( 1 - params$alpha )
p.ss <- 1 / c.ss
coeff.cont.init[1, ] <- log( c( c.ss, c.ss, 1 / params$betta, 1 / params$betta, 
                                params$alpha, params$alpha, 1 - params$alpha, 1 - params$alpha,
                                p.ss, p.ss, 1, 1, 1 ) )

coeff.cont.init[2, ] <- c( -.1, .5, -.2, .5, rep( .1, 9 ) ) / 10
coeff.cont.init[3, ] <- c(  .5,-.1,  .5,-.2, rep( .1, 9 ) ) / 10
coeff.cont.init[4, ] <- coeff.cont.init[2, ]
coeff.cont.init[5, ] <- coeff.cont.init[3, ]


opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$l.pairs <- list( c(1,2) )
opt$l.pairs.cont <- list( c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12) )
opt$sym.reg <- TRUE


# Rough solution
opt$iter <- 12
sol.irbc.N.1.b <- sol.irbc.iterate( coeff.init, opt, params, coeff.init.cont, debug.flag = F )

# Finer solution
opt$n.gain <- .02
opt$iter <- 15
sol.irbc.N.1.b <- sol.irbc.iterate( sol.irbc.N.1.b$coeff, opt, params, sol.irbc.N.1.b$coeff.cont )
rep.irbc.N.1.b <- report.data( sol.irbc.N.1.b )
report.create( sol.irbc.N.1.b, rep.irbc.N.1.b )

##### THIS WORKS !!!! ###
## Now: Find a nonlinear solution ##
opt$N <- 2
coeff.init <- matrix( 0, 15, 2 )
coeff.init[ c(1,2,4,7,11), ] <- sol.irbc.N.1.b$coeff
coeff.cont.init <- matrix( 0, 15, 13 )
coeff.cont.init[ c(1,2,4,7,11), ] <- sol.irbc.N.1.b$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                        ave=c(1,5,14) )
opt$iter <- 20
# Try a nonlinear solution
sol.irbc.N.1.b.nl <- sol.irbc.iterate( coeff.init, opt, params, coeff.cont.init )
rep.irbc.N.1.b.nl <- report.data( sol.irbc.N.1.b.nl )
report.create( sol.irbc.N.1.b.nl, rep.irbc.N.1.b.nl )


### Now increase gamma ###
params$gamma <- 3
opt <- sol.irbc.N.1.b$opt
opt$n.gain <- .1
opt$iter <- 16
sol.irbc.2.b <- sol.irbc.iterate( sol.irbc.N.1.b$coeff, opt, params, sol.irbc.N.1.b$coeff.cont )
    # This is about as good at the linear solution gets
opt$N <- 2
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
coeff.init[ c(1,2,4,7,11), ] <- sol.irbc.2.b$coeff
coeff.cont.init[ c(1,2,4,7,11), ] <- sol.irbc.2.b$coeff.cont
opt$n.gain <- .03
opt$iter <- 60
sol.irbc.2.b.nl <- sol.irbc.iterate( coeff.init, opt, params, coeff.cont.init )
    # Nonlinear solution
rep.irbc.2.b <- report.data( sol.irbc.2.b )
report.create( sol.irbc.2.b, rep.irbc.2.b )
rep.irbc.2.b.nl <- report.data( sol.irbc.2.b.nl )
report.create( sol.irbc.2.b.nl, rep.irbc.2.b.nl )

params$gamma <- 4
opt <- sol.irbc.2.b$opt
opt$n.gain <- .025
opt$iter <- 50
sol.irbc.3.b <- sol.irbc.iterate( sol.irbc.2.b$coeff, opt, params, sol.irbc.2.b$coeff.cont )
    # This is about as good at the linear solution gets
opt <- sol.irbc.2.b.nl$opt
opt$k.gain <- .1
opt$k.iter <- 30
coeff.init[ c(1,2,4,7,11), ] <- sol.irbc.3.b$coeff
coeff.cont.init[ c(1,2,4,7,11), ] <- sol.irbc.3.b$coeff.cont
sol.irbc.3.b.nl <- sol.irbc.iterate( coeff.init, opt, params, coeff.cont.init )
    # Nonlinear solution
rep.irbc.3.b <- report.data( sol.irbc.3.b )
report.create( sol.irbc.3.b, rep.irbc.3.b )
rep.irbc.3.b.nl <- report.data( sol.irbc.3.b.nl )
report.create( sol.irbc.3.b.nl, rep.irbc.3.b.nl )


### Now increase rho ###
params <- sol.irbc.N.1.b$params
opt <- sol.irbc.N.1.b$opt
opt$n.gain <- .1
opt$iter <- 22
params$rho <- c( .75, .75 )
sol.irbc.4.b <- sol.irbc.iterate( sol.irbc.N.1.b$coeff, opt, params, sol.irbc.N.1.b$coeff.cont )
    # Linear solution
opt$N <- 2
coeff.init <- matrix( 0, 15, 2 )
coeff.init[ c(1,2,4,7,11), ] <- sol.irbc.4.b$coeff
coeff.cont.init <- matrix( 0, 15, 13 )
coeff.cont.init[ c(1,2,4,7,11), ] <- sol.irbc.4.b$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 40
opt$n.gain <- .1
sol.irbc.4.b.nl <- sol.irbc.iterate( coeff.init, opt, params, coeff.cont.init )
    # Try a nonlinear solution
rep.irbc.4.b <- report.data( sol.irbc.4.b )
report.create( sol.irbc.4.b, rep.irbc.4.b )
rep.irbc.4.b.nl <- report.data( sol.irbc.4.b.nl )
report.create( sol.irbc.4.b.nl, rep.irbc.4.b.nl )
    # Report


### Now increase sigma.eps ###
params <- sol.irbc.N.1.b$params
opt <- sol.irbc.N.1.b$opt
opt$n.gain <- .02
opt$iter <- 10
params$sig.eps <- c( .05, .05 )
sol.irbc.5.b <- sol.irbc.iterate( sol.irbc.N.1.b$coeff, opt, params, sol.irbc.N.1.b$coeff.cont )
    # The linear solution.  It is quite poor.
opt$N <- 2
coeff.init <- matrix( 0, 15, 2 )
coeff.init[ c(1,2,4,7,11), ] <- sol.irbc.5.b$coeff
coeff.cont.init <- matrix( 0, 15, 13 )
coeff.cont.init[ c(1,2,4,7,11), ] <- sol.irbc.5.b$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 30
opt$n.gain <- .1
sol.irbc.5.b.nl <- sol.irbc.iterate( coeff.init, opt, params, coeff.cont.init )
    # Try a nonlinear solution
rep.irbc.5.b <- report.data( sol.irbc.5.b )
report.create( sol.irbc.5.b, rep.irbc.5.b )
rep.irbc.5.b.nl <- report.data( sol.irbc.5.b.nl )
report.create( sol.irbc.5.b.nl, rep.irbc.5.b.nl )
    # Report


### Now increase sigma.eps more ###
params$sig.eps <- c(.075,.075)
sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
upper <- c(  3 * sd.x, rep( .4, 2 ) )
lower <- -upper
opt <- sol.irbc.5.b$opt
opt$iter <- 27
opt$n.gain <- .075
opt$upper <- upper
opt$lower <- lower
sol.irbc.6.b <- sol.irbc.iterate( sol.irbc.5.b$coeff, opt, params, sol.irbc.5.b$coeff.cont )
    # The linear solution.
opt <- sol.irbc.5.b.nl$opt
opt$n.gain <- .05
opt$iter <- 30
sol.irbc.6.b.nl <- sol.irbc.iterate( sol.irbc.5.b.nl$coeff, opt, params, sol.irbc.5.b.nl$coeff.cont )
    # The nonlinear solution
rep.irbc.6.b <- report.data( sol.irbc.6.b )
report.create( sol.irbc.6.b, rep.irbc.6.b )
rep.irbc.6.b.nl <- report.data( sol.irbc.6.b.nl )
report.create( sol.irbc.6.b.nl, rep.irbc.6.b.nl )
    # Report


### Now increase alpha ###
params <- sol.irbc.N.1.b$params
opt <- sol.irbc.N.1.b$opt
opt$n.gain <- .05
opt$iter <- 16
params$alpha <- .85
sol.irbc.7.b <- sol.irbc.iterate( sol.irbc.N.1.b$coeff, opt, params, sol.irbc.N.1.b$coeff.cont )
    # Linear solution
opt$N <- 2
coeff.init <- matrix( 0, 15, 2 )
coeff.init[ c(1,2,4,7,11), ] <- sol.irbc.7.b$coeff
coeff.cont.init <- matrix( 0, 15, 13 )
coeff.cont.init[ c(1,2,4,7,11), ] <- sol.irbc.7.b$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 30
opt$n.gain <- .05
opt$n.k <- 30
sol.irbc.7.b.nl <- sol.irbc.iterate( coeff.init, opt, params, coeff.cont.init )
    # The nonlinear solution
rep.irbc.7.b <- report.data( sol.irbc.7.b )
report.create( sol.irbc.7.b, rep.irbc.7.b )
rep.irbc.7.b.nl <- report.data( sol.irbc.7.b.nl )
report.create( sol.irbc.7.b.nl, rep.irbc.7.b.nl )
    # Report
