rm(list=ls())
library(xtable)
library(scales)
library(edsProjection)
setwd("~/code/2016/edsProjection")

##### 0. INITIALIZATION #####

### 0.1 Baseline parameters ###
params <- list( share = .75, gamma = 2, P1.bar=1, P2.bar=1, betta=.95,
                rho=c(.25,.25), sig.eps=c(.01,.01), eta=1 )

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


### 0.2 Initial coefficients ###
coeff.init <- matrix( .02, 5, 2 )
coeff.init[2,] <- c( -.01,.02)
coeff.init[3,] <- c( .02,-.01)
coeff.init[4,] <- c( -.02,.01)
coeff.init[5,] <- c( .01,-.02)

coeff.cont.init <- matrix( 0, 5, opt$n.cont )
c.ss <- params$share ^ params$share * ( 1 - params$share ) ^ ( 1 - params$share )
p.ss <- 1 / c.ss
coeff.cont.init[1, ] <- log( c( c.ss, c.ss, 1 / params$betta, 1 / params$betta, 
                                params$share, params$share, 1 - params$share, 1 - params$share,
                                p.ss, p.ss, 1, 1, 1 ) )

coeff.cont.init[2, ] <- c( -.1, .5, -.2, .5, rep( .1, 9 ) ) / 10
coeff.cont.init[3, ] <- c(  .5,-.1,  .5,-.2, rep( .1, 9 ) ) / 10
coeff.cont.init[4, ] <- coeff.cont.init[2, ]
coeff.cont.init[5, ] <- coeff.cont.init[3, ]


opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$l.pairs <- list( c(1,2) )
opt$l.pairs.cont <- list( c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12) )
opt$sym.reg <- TRUE


##### 1. BASELINE SOLUTION #####

### 1.1 Rough linear solution ###
opt$n.gain <- .25
opt$iter <- 14
sol.base.1.r <- sol.irbc.iterate( coeff.init, opt, params, coeff.cont.init, debug.flag = F )

### 1.2 Finer linear solution ###
opt$n.gain <- .05
opt$iter <- 20
sol.base.1 <- sol.irbc.iterate( sol.base.1.r$coeff, opt, params, sol.base.1.r$coeff.cont )
opt$n.gain <- .01
opt$iter <- 10
sol.base.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.base.1 <- report.data( sol.base.1 )
print( paste0( "err = ", round( max(apply( abs( rep.base.1$err ), 2, mean )) * 100, 4), "pp" ) )

### 1.3 A nonlinear solution ###
opt$N <- 2
coeff.init.2 <- matrix( 0, 15, 2 )
coeff.init.2[ c(1,2,4,7,11), ] <- sol.base.1$coeff
coeff.cont.init.2 <- matrix( 0, 15, 13 )
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.base.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 25
opt$n.gain <- .05
sol.base.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
    # Solve it
rep.base.2 <- report.data( sol.base.2 )
print( paste0( "err = ", round( max(apply( abs( rep.base.2$err ), 2, mean )) * 100, 4), "pp" ) )

save( sol.base.1, sol.base.2, rep.base.1, rep.base.2, file='keyStats.rdata')

##### 2. TRY CHANGING RHO #####

l.rho <- list()

### 2.1 rho = .2 ###
params$rho <- c(.2, .2)
sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
opt$upper <- c(  3 * sd.x, rep( .1, 2 ) )
opt$lower <- -upper
    # Change parameters and bounds
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .05
opt$iter <- 25
opt$N <- 1
sol.rho.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1.r$coeff.cont )
rep.rho.1 <- report.data( sol.rho.1 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 40
opt$n.gain <- .05
sol.rho.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.rho.2 <- report.data( sol.rho.2 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.rho[[1]] <- list( sol.1=sol.rho.1, rep.1=rep.rho.1, sol.2=sol.rho.2, rep.2=rep.rho.2 )
    # Assign to the list

### 2.2 rho = .15 ###
params$rho <- c(.15, .15)
sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
opt$upper <- c(  3 * sd.x, rep( .1, 2 ) )
opt$lower <- -upper
    # Change parameters and bounds
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .04
opt$iter <- 35
opt$N <- 1
sol.rho.1 <- sol.irbc.iterate( l.rho[[1]]$sol.1$coeff, opt, params, l.rho[[1]]$sol.1$coeff.cont )
rep.rho.1 <- report.data( sol.rho.1 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 40
opt$n.gain <- .05
sol.rho.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.rho.2 <- report.data( sol.rho.2 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.rho[[2]] <- list( sol.1=sol.rho.1, rep.1=rep.rho.1, sol.2=sol.rho.2, rep.2=rep.rho.2 )
    # Assign to the list

### 2.3 rho = .1 ###
params$rho <- c(.1, .1)
sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
opt$upper <- c(  3 * sd.x, rep( .1, 2 ) )
opt$lower <- -upper
    # Change parameters and bounds
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .04
opt$iter <- 35
opt$N <- 1
sol.rho.1 <- sol.irbc.iterate( l.rho[[2]]$sol.1$coeff, opt, params, l.rho[[2]]$sol.1$coeff.cont )
rep.rho.1 <- report.data( sol.rho.1 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 40
opt$n.gain <- .05
sol.rho.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.rho.2 <- report.data( sol.rho.2 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.rho[[3]] <- list( sol.1=sol.rho.1, rep.1=rep.rho.1, sol.2=sol.rho.2, rep.2=rep.rho.2 )
    # Assign to the list


### 2.4 rho = .05 ###
params$rho <- c(.05, .05)
sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
opt$upper <- c(  3 * sd.x, rep( .1, 2 ) )
opt$lower <- -upper
    # Change parameters and bounds
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .04
opt$iter <- 35
opt$N <- 1
sol.rho.1 <- sol.irbc.iterate( l.rho[[3]]$sol.1$coeff, opt, params, l.rho[[3]]$sol.1$coeff.cont )
rep.rho.1 <- report.data( sol.rho.1 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 40
opt$n.gain <- .075
sol.rho.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.rho.2 <- report.data( sol.rho.2 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.rho[[4]] <- list( sol.1=sol.rho.1, rep.1=rep.rho.1, sol.2=sol.rho.2, rep.2=rep.rho.2 )
    # Assign to the list

### 2.5 rho = .3 ###
params$rho <- c(.3, .3)
sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
opt$upper <- c(  3 * sd.x, rep( .1, 2 ) )
opt$lower <- -upper
    # Change parameters and bounds
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .05
opt$iter <- 25
opt$N <- 1
sol.rho.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.rho.1 <- report.data( sol.rho.1 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 40
opt$n.gain <- .05
sol.rho.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.rho.2 <- report.data( sol.rho.2 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.rho[[5]] <- list( sol.1=sol.rho.1, rep.1=rep.rho.1, sol.2=sol.rho.2, rep.2=rep.rho.2 )
    # Assign to the list

### 2.6 rho = .35 ###
params$rho <- c(.35, .35)
sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
opt$upper <- c(  3 * sd.x, rep( .1, 2 ) )
opt$lower <- -upper
    # Change parameters and bounds
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .04
opt$iter <- 35
opt$N <- 1
sol.rho.1 <- sol.irbc.iterate( l.rho[[5]]$sol.1$coeff, opt, params, l.rho[[5]]$sol.1$coeff.cont )
rep.rho.1 <- report.data( sol.rho.1 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 40
opt$n.gain <- .075
sol.rho.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.rho.2 <- report.data( sol.rho.2 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.rho[[6]] <- list( sol.1=sol.rho.1, rep.1=rep.rho.1, sol.2=sol.rho.2, rep.2=rep.rho.2 )
    # Assign to the list

### 2.7 rho = .4 ###
params$rho <- c(.4, .4 )
sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
opt$upper <- c(  3 * sd.x, rep( .1, 2 ) )
opt$lower <- -upper
    # Change parameters and bounds
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .03
opt$iter <- 35
opt$N <- 1
sol.rho.1 <- sol.irbc.iterate( l.rho[[6]]$sol.1$coeff, opt, params, l.rho[[6]]$sol.1$coeff.cont )
rep.rho.1 <- report.data( sol.rho.1 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 30
opt$n.gain <- .075
sol.rho.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.rho.2 <- report.data( sol.rho.2 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.rho[[7]] <- list( sol.1=sol.rho.1, rep.1=rep.rho.1, sol.2=sol.rho.2, rep.2=rep.rho.2 )
    # Assign to the list

### 2.8 rho = .45 ###
params$rho <- c(.45, .45 )
sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
opt$upper <- c(  3 * sd.x, rep( .1, 2 ) )
opt$lower <- -upper
    # Change parameters and bounds
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .02
opt$iter <- 28
opt$N <- 1
sol.rho.1 <- sol.irbc.iterate( l.rho[[7]]$sol.1$coeff, opt, params, l.rho[[7]]$sol.1$coeff.cont )
rep.rho.1 <- report.data( sol.rho.1 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 40
opt$n.gain <- .075
sol.rho.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.rho.2 <- report.data( sol.rho.2 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.rho[[8]] <- list( sol.1=sol.rho.1, rep.1=rep.rho.1, sol.2=sol.rho.2, rep.2=rep.rho.2 )
    # Assign to the list

### 2.9 rho = .5 ###
params$rho <- c(.5, .5 )
sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
opt$upper <- c(  3 * sd.x, rep( .1, 2 ) )
opt$lower <- -upper
    # Change parameters and bounds
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .02
opt$iter <- 28
opt$N <- 1
sol.rho.1 <- sol.irbc.iterate( l.rho[[8]]$sol.1$coeff, opt, params, l.rho[[8]]$sol.1$coeff.cont )
rep.rho.1 <- report.data( sol.rho.1 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 40
opt$n.gain <- .075
sol.rho.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.rho.2 <- report.data( sol.rho.2 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.rho[[9]] <- list( sol.1=sol.rho.1, rep.1=rep.rho.1, sol.2=sol.rho.2, rep.2=rep.rho.2 )
    # Assign to the list

### 2.10 rho = .55 ###
params$rho <- c(.55, .55 )
sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
opt$upper <- c(  3 * sd.x, rep( .1, 2 ) )
opt$lower <- -upper
    # Change parameters and bounds
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .05
opt$iter <- 19
opt$N <- 1
sol.rho.1 <- sol.irbc.iterate( l.rho[[9]]$sol.1$coeff, opt, params, l.rho[[9]]$sol.1$coeff.cont )
rep.rho.1 <- report.data( sol.rho.1 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 50
opt$n.gain <- .08
sol.rho.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.rho.2 <- report.data( sol.rho.2 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.rho[[10]] <- list( sol.1=sol.rho.1, rep.1=rep.rho.1, sol.2=sol.rho.2, rep.2=rep.rho.2 )
    # Assign to the list                        

save( l.rho, file='~/Dropbox/outsize/irbc/rho_change.rdata' )

### 2.11 rho = .6 ###
params$rho <- c(.6, .6 )
sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
opt$upper <- c(  3 * sd.x, rep( .1, 2 ) )
opt$lower <- -upper
    # Change parameters and bounds
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .02
opt$iter <- 10
opt$N <- 1
sol.rho.1 <- sol.irbc.iterate( l.rho[[10]]$sol.1$coeff, opt, params, l.rho[[10]]$sol.1$coeff.cont )
rep.rho.1 <- report.data( sol.rho.1 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 40
opt$n.gain <- .08
sol.rho.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.rho.2 <- report.data( sol.rho.2 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.rho[[11]] <- list( sol.1=sol.rho.1, rep.1=rep.rho.1, sol.2=sol.rho.2, rep.2=rep.rho.2 )
    # Assign to the list

### 2.12 rho = .65 ###
params$rho <- c(.65, .65 )
sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
opt$upper <- c(  3 * sd.x, rep( .1, 2 ) )
opt$lower <- -upper
    # Change parameters and bounds
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .02
opt$iter <- 15
opt$N <- 1
sol.rho.1 <- sol.irbc.iterate( l.rho[[11]]$sol.1$coeff, opt, params, l.rho[[11]]$sol.1$coeff.cont )
rep.rho.1 <- report.data( sol.rho.1 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 50
opt$n.gain <- .08
sol.rho.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.rho.2 <- report.data( sol.rho.2 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.rho[[12]] <- list( sol.1=sol.rho.1, rep.1=rep.rho.1, sol.2=sol.rho.2, rep.2=rep.rho.2 )
    # Assign to the list

### 2.13 rho = .7 ###
params$rho <- c(.7, .7 )
sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
opt$upper <- c(  3 * sd.x, rep( .1, 2 ) )
opt$lower <- -upper
    # Change parameters and bounds
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .02
opt$iter <- 15
opt$N <- 1
sol.rho.1 <- sol.irbc.iterate( l.rho[[12]]$sol.1$coeff, opt, params, l.rho[[12]]$sol.1$coeff.cont )
rep.rho.1 <- report.data( sol.rho.1 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 50
opt$n.gain <- .08
sol.rho.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.rho.2 <- report.data( sol.rho.2 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.rho[[13]] <- list( sol.1=sol.rho.1, rep.1=rep.rho.1, sol.2=sol.rho.2, rep.2=rep.rho.2 )
    # Assign to the list

### 2.14 rho = .75 ###
params$rho <- c(.75, .75 )
sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
opt$upper <- c(  3 * sd.x, rep( .1, 2 ) )
opt$lower <- -upper
    # Change parameters and bounds
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .01
opt$iter <- 10
opt$N <- 1
sol.rho.1 <- sol.irbc.iterate( l.rho[[13]]$sol.1$coeff, opt, params, l.rho[[13]]$sol.1$coeff.cont )
rep.rho.1 <- report.data( sol.rho.1 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.rho.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 50
opt$n.gain <- .08
sol.rho.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.rho.2 <- report.data( sol.rho.2 )
print( paste0( "err = ", round( max(apply( abs( rep.rho.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.rho[[14]] <- list( sol.1=sol.rho.1, rep.1=rep.rho.1, sol.2=sol.rho.2, rep.2=rep.rho.2 )
    # Assign to the list
# report.create( sol.rho.2, rep.rho.2 )


save( l.rho, file='~/Dropbox/outsize/irbc/rho_change.rdata' )

### 2.15 rho = .8 ###

## Higher rho fails