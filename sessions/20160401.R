rm(list=ls())
library(xtable)
library(scales)
library(edsProjection)
setwd("~/code/2016/edsProjection")

##### 0. INITIALIZATION #####

### 0.1 Baseline parameters ###
params <- list( share = .75, gamma = 2, P1.bar=1, P2.bar=1, betta=.95,
                rho=c(.5,.5), sig.eps=c(.01,.01), eta=1 )

sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
upper <- c(  3 * sd.x, rep( .5, 2 ) )
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
opt$n.gain <- .02
opt$iter <- 28
sol.base.1 <- sol.irbc.iterate( sol.base.1.r$coeff, opt, params, sol.base.1.r$coeff.cont )
opt$iter <- 40
opt$n.gain <- .075
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

save( sol.base.1, sol.base.2, rep.base.1, rep.base.2, file='~/Dropbox/outsize/irbc/eta_change.rdata')

##### 2. TRY CHANGING ETA #####
l.eta <- list()

#### NEED TO RERUN THIS WITH THE WIDER BOUNDS ON B ####

### 2.1 eta = 1.25 ###
params$eta <- 1.25
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .15
opt$iter <- 25
opt$N <- 1
sol.eta.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.eta.1 <- report.data( sol.eta.1 )
print( paste0( "err = ", round( max(apply( abs( rep.eta.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.eta.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.eta.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 30
opt$n.gain <- .075
sol.eta.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.eta.2 <- report.data( sol.eta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.eta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.eta[[1]] <- list( sol.1=sol.eta.1, rep.1=rep.eta.1, sol.2=sol.eta.2, rep.2=rep.eta.2 )
    # Assign to the list


### 2.2 eta = 1.5 ###
params$eta <- 1.5
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .02
opt$iter <- 5
opt$N <- 1
sol.eta.1 <- sol.irbc.iterate( l.eta[[1]]$sol.1$coeff, opt, params, l.eta[[1]]$sol.1$coeff.cont )
rep.eta.1 <- report.data( sol.eta.1 )
print( paste0( "err = ", round( max(apply( abs( rep.eta.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.eta.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.eta.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 50
opt$n.gain <- .075
sol.eta.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.eta.2 <- report.data( sol.eta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.eta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.eta[[2]] <- list( sol.1=sol.eta.1, rep.1=rep.eta.1, sol.2=sol.eta.2, rep.2=rep.eta.2 )
    # Assign to the list

### 2.2 eta = 1.65 ###
params$eta <- 1.65
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .02
opt$iter <- 4
opt$N <- 1
sol.eta.1 <- sol.irbc.iterate( l.eta[[2]]$sol.1$coeff, opt, params, l.eta[[2]]$sol.1$coeff.cont )
rep.eta.1 <- report.data( sol.eta.1 )
print( paste0( "err = ", round( max(apply( abs( rep.eta.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2 <- l.eta[[2]]$sol.2$coeff
coeff.cont.init.2 <- l.eta[[2]]$sol.2$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 40
opt$n.gain <- .075
opt$k.gain <- .4
opt$tol <- 1e-06
sol.eta.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.eta.2 <- report.data( sol.eta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.eta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.eta[[3]] <- list( sol.1=sol.eta.1, rep.1=rep.eta.1, sol.2=sol.eta.2, rep.2=rep.eta.2 )
    # Assign to the list

### 2.2 eta = 1.75 ###
params$eta <- 1.75
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .01
opt$k.gain <- .05
opt$iter <- 4
opt$N <- 1
# sol.eta.1 <- sol.irbc.iterate( l.eta[[3]]$sol.1$coeff, opt, params, l.eta[[3]]$sol.1$coeff.cont )
# rep.eta.1 <- report.data( sol.eta.1 )
# print( paste0( "err = ", round( max(apply( abs( rep.eta.1$err ), 2, mean )) * 100, 4), "pp" ) )
#     # The linear solution
### Linear solution fails
opt$N <- 2
coeff.init.2 <- l.eta[[3]]$sol.2$coeff
coeff.cont.init.2 <- l.eta[[3]]$sol.2$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 40
opt$n.gain <- .075
opt$k.gain <- .05
opt$tol <- 1e-06
sol.eta.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.eta.2 <- report.data( sol.eta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.eta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.eta[[4]] <- list( sol.2=sol.eta.2, rep.2=rep.eta.2 )
# Assign to the list

save( l.eta, sol.base.1, sol.base.2, rep.base.1, rep.base.2, file='~/Dropbox/outsize/irbc/eta_change.rdata')