rm(list=ls())
library(xtable)
library(scales)
library(edsProjection)
setwd("~/code/2016/edsProjection")

##### 0. INITIALIZATION #####

### 0.1 Baseline parameters ###
params <- list( share = .75, gamma = 2, P1.bar=1, P2.bar=1, betta=.95,
                rho=c(.5,.5), sig.eps=c(.01,.01), eta=1.5 )

sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
upper <- c(  3 * sd.x, rep( .75, 2 ) )
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

### 1.1 Linear solution ###
opt$n.gain <- .25
opt$iter <- 18
sol.base.1 <- sol.irbc.iterate( coeff.init, opt, params, coeff.cont.init, debug.flag = F )
rep.base.1 <- report.data( sol.base.1 )
print( paste0( "err = ", round( max(apply( abs( rep.base.1$err ), 2, mean )) * 100, 4), "pp" ) )

### 1.2 A nonlinear solution ###
opt$N <- 2
coeff.init.2 <- matrix( 0, 15, 2 )
coeff.init.2[ c(1,2,4,7,11), ] <- sol.base.1$coeff
coeff.cont.init.2 <- matrix( 0, 15, 13 )
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.base.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 30
opt$n.gain <- .1
sol.base.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
    # Solve it
rep.base.2 <- report.data( sol.base.2 )
print( paste0( "err = ", round( max(apply( abs( rep.base.2$err ), 2, mean )) * 100, 4), "pp" ) )

save( sol.base.1, sol.base.2, rep.base.1, rep.base.2, file='~/Dropbox/outsize/irbc/gamma_change.rdata')

##### 2. TRY CHANGING GAMMA #####
l.gamma <- list()

### 2.1 gamma = 1.75 ###
params$gamma <- 1.75
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$iter <- 16
opt$N <- 1
sol.gamma.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.gamma.1 <- report.data( sol.gamma.1 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 35
opt$n.gain <- .1
sol.gamma.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.gamma.2 <- report.data( sol.gamma.2 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.gamma[[1]] <- list( sol.1=sol.gamma.1, rep.1=rep.gamma.1, sol.2=sol.gamma.2, rep.2=rep.gamma.2 )
    # Assign to the list

### 2.2 gamma = 1.5 ###
params$gamma <- 1.5
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$iter <- 11
opt$N <- 1
sol.gamma.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.gamma.1 <- report.data( sol.gamma.1 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 33
sol.gamma.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.gamma.2 <- report.data( sol.gamma.2 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.gamma[[2]] <- list( sol.1=sol.gamma.1, rep.1=rep.gamma.1, sol.2=sol.gamma.2, rep.2=rep.gamma.2 )
    # Assign to the list

### 2.3 gamma = 1.25 ###
params$gamma <- 1.25
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$iter <- 21
opt$N <- 1
sol.gamma.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.gamma.1 <- report.data( sol.gamma.1 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$n.gain <- .15
opt$iter <- 35
opt$c.iter <- 50
sol.gamma.2 <- sol.irbc.iterate( l.gamma[[2]]$sol.2$coeff, opt, params, l.gamma[[2]]$sol.2$coeff.cont )
rep.gamma.2 <- report.data( sol.gamma.2 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.gamma[[3]] <- list( sol.1=sol.gamma.1, rep.1=rep.gamma.1, sol.2=sol.gamma.2, rep.2=rep.gamma.2 )
    # Assign to the list

save( l.gamma, sol.base.1, sol.base.2, rep.base.1, rep.base.2, file='~/Dropbox/outsize/irbc/gamma_change.rdata')

### 2.4 gamma = 1 ###
params$gamma <- 1
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$iter <- 21
opt$N <- 1
sol.gamma.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.gamma.1 <- report.data( sol.gamma.1 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 35
opt$c.iter <- 50
sol.gamma.2 <- sol.irbc.iterate( l.gamma[[2]]$sol.2$coeff, opt, params, l.gamma[[2]]$sol.2$coeff.cont )
rep.gamma.2 <- report.data( sol.gamma.2 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.gamma[[4]] <- list( sol.1=sol.gamma.1, rep.1=rep.gamma.1, sol.2=sol.gamma.2, rep.2=rep.gamma.2 )
    # Assign to the list

### 2.5 gamma = 0.75 ###
params$gamma <- 0.75
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$iter <- 30
opt$N <- 1
sol.gamma.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.gamma.1 <- report.data( sol.gamma.1 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 33
opt$c.iter <- 50
sol.gamma.2 <- sol.irbc.iterate( l.gamma[[2]]$sol.2$coeff, opt, params, l.gamma[[2]]$sol.2$coeff.cont )
rep.gamma.2 <- report.data( sol.gamma.2 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.gamma[[5]] <- list( sol.1=sol.gamma.1, rep.1=rep.gamma.1, sol.2=sol.gamma.2, rep.2=rep.gamma.2 )
    # Assign to the list

### 2.6 gamma = 0.5 ###
params$gamma <- 0.5
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$iter <- 30
opt$N <- 1
sol.gamma.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.gamma.1 <- report.data( sol.gamma.1 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 10
opt$c.iter <- 50
sol.gamma.2 <- sol.irbc.iterate( l.gamma[[2]]$sol.2$coeff, opt, params, l.gamma[[2]]$sol.2$coeff.cont )
rep.gamma.2 <- report.data( sol.gamma.2 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.gamma[[6]] <- list( sol.1=sol.gamma.1, rep.1=rep.gamma.1, sol.2=sol.gamma.2, rep.2=rep.gamma.2 )
    # Assign to the list

### 2.7 gamma = 2.25 ###
params$gamma <- 2.25
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$iter <- 4
opt$N <- 1
sol.gamma.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.gamma.1 <- report.data( sol.gamma.1 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 50
opt$n.gain <- .05
sol.gamma.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.gamma.2 <- report.data( sol.gamma.2 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.gamma[[7]] <- list( sol.1=sol.gamma.1, rep.1=rep.gamma.1, sol.2=sol.gamma.2, rep.2=rep.gamma.2 )
    # Assign to the list


### 2.8 gamma = 2.5 ###
params$gamma <- 2.5
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$iter <- 4
opt$N <- 1
sol.gamma.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.gamma.1 <- report.data( sol.gamma.1 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 80
opt$n.gain <- .05
sol.gamma.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.gamma.2 <- report.data( sol.gamma.2 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.gamma[[8]] <- list( sol.1=sol.gamma.1, rep.1=rep.gamma.1, sol.2=sol.gamma.2, rep.2=rep.gamma.2 )
    # Assign to the list

save( l.gamma, sol.base.1, sol.base.2, rep.base.1, rep.base.2, file='~/Dropbox/outsize/irbc/gamma_change.rdata')


### 2.9 gamma = 2.75 ###
params$gamma <- 2.75
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$iter <- 20
opt$N <- 1
sol.gamma.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.gamma.1 <- report.data( sol.gamma.1 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 50
opt$n.gain <- .075
sol.gamma.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.gamma.2 <- report.data( sol.gamma.2 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.gamma[[9]] <- list( sol.1=sol.gamma.1, rep.1=rep.gamma.1, sol.2=sol.gamma.2, rep.2=rep.gamma.2 )
    # Assign to the list


### 2.10 gamma = 3 ###
params$gamma <- 3
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$iter <- 20
opt$N <- 1
sol.gamma.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.gamma.1 <- report.data( sol.gamma.1 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 50
opt$n.gain <- .075
sol.gamma.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.gamma.2 <- report.data( sol.gamma.2 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.gamma[[10]] <- list( sol.1=sol.gamma.1, rep.1=rep.gamma.1, sol.2=sol.gamma.2, rep.2=rep.gamma.2 )
    # Assign to the list

### 2.12 gamma = 3.25 ###
params$gamma <- 3.25
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$iter <- 20
opt$N <- 1
sol.gamma.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.gamma.1 <- report.data( sol.gamma.1 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 50
sol.gamma.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.gamma.2 <- report.data( sol.gamma.2 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.gamma[[11]] <- list( sol.1=sol.gamma.1, rep.1=rep.gamma.1, sol.2=sol.gamma.2, rep.2=rep.gamma.2 )
    # Assign to the list

### 2.12 gamma = 3.5 ###
params$gamma <- 3.5
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$iter <- 20
opt$N <- 1
sol.gamma.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.gamma.1 <- report.data( sol.gamma.1 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.gamma.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 50
sol.gamma.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.gamma.2 <- report.data( sol.gamma.2 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.gamma[[12]] <- list( sol.1=sol.gamma.1, rep.1=rep.gamma.1, sol.2=sol.gamma.2, rep.2=rep.gamma.2 )
    # Assign to the list

save( l.gamma, sol.base.1, sol.base.2, rep.base.1, rep.base.2, file='~/Dropbox/outsize/irbc/gamma_change.rdata')
