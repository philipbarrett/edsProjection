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

save( sol.base.1, sol.base.2, rep.base.1, rep.base.2, file='~/Dropbox/outsize/irbc/betta_change.rdata')

##### 2. TRY CHANGING BETTA #####
l.betta <- list()

### 2.1 betta = .96 ###
params$betta <- .96
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .15
opt$iter <- 25
opt$N <- 1
sol.betta.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.betta.1 <- report.data( sol.betta.1 )
print( paste0( "err = ", round( max(apply( abs( rep.betta.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.betta.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.betta.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 40
opt$n.gain <- .05
sol.betta.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.betta.2 <- report.data( sol.betta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.betta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.betta[[1]] <- list( sol.1=sol.betta.1, rep.1=rep.betta.1, sol.2=sol.betta.2, rep.2=rep.betta.2 )
    # Assign to the list

### 2.2 betta = .97 ###
params$betta <- .97
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .15
opt$iter <- 15
opt$N <- 1
sol.betta.1 <- sol.irbc.iterate( l.betta[[1]]$sol.1$coeff, opt, params, l.betta[[1]]$sol.1$coeff.cont )
rep.betta.1 <- report.data( sol.betta.1 )
print( paste0( "err = ", round( max(apply( abs( rep.betta.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.betta.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.betta.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 70
opt$n.gain <- .05
sol.betta.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.betta.2 <- report.data( sol.betta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.betta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.betta[[2]] <- list( sol.1=sol.betta.1, rep.1=rep.betta.1, sol.2=sol.betta.2, rep.2=rep.betta.2 )
    # Assign to the list

### 2.2 betta = .98 ###
params$betta <- .98
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .075
opt$iter <- 7
opt$N <- 1
sol.betta.1 <- sol.irbc.iterate( l.betta[[2]]$sol.1$coeff, opt, params, l.betta[[2]]$sol.1$coeff.cont )
rep.betta.1 <- report.data( sol.betta.1 )
print( paste0( "err = ", round( max(apply( abs( rep.betta.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.betta.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.betta.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 120
opt$n.gain <- .04
sol.betta.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.betta.2 <- report.data( sol.betta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.betta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.betta[[3]] <- list( sol.1=sol.betta.1, rep.1=rep.betta.1, sol.2=sol.betta.2, rep.2=rep.betta.2 )
    # Assign to the list

### 2.3 betta = .99 ###
params$betta <- .99
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .075
opt$iter <- 11
opt$N <- 1
sol.betta.1 <- sol.irbc.iterate( l.betta[[3]]$sol.1$coeff, opt, params, l.betta[[3]]$sol.1$coeff.cont )
rep.betta.1 <- report.data( sol.betta.1 )
print( paste0( "err = ", round( max(apply( abs( rep.betta.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.betta.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.betta.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 180
opt$n.gain <- .03
sol.betta.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.betta.2 <- report.data( sol.betta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.betta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.betta[[4]] <- list( sol.1=sol.betta.1, rep.1=rep.betta.1, sol.2=sol.betta.2, rep.2=rep.betta.2 )
    # Assign to the list

save( params, opt, sol.base.1, sol.base.2, rep.base.1, rep.base.2, l.betta, file='~/Dropbox/outsize/irbc/betta_change.rdata')

### 2.4 betta = .925 ###
params$betta <- .925
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .1
opt$iter <- 31
opt$N <- 1
sol.betta.1 <- sol.irbc.iterate( sol.base.1$coeff, opt, params, sol.base.1$coeff.cont )
rep.betta.1 <- report.data( sol.betta.1 )
print( paste0( "err = ", round( max(apply( abs( rep.betta.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.betta.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.betta.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 40
opt$n.gain <- .1
sol.betta.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.betta.2 <- report.data( sol.betta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.betta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.betta[[5]] <- list( sol.1=sol.betta.1, rep.1=rep.betta.1, sol.2=sol.betta.2, rep.2=rep.betta.2 )
    # Assign to the list

### 2.5 betta = .9 ###
params$betta <- .9
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .05
opt$iter <- 26
opt$N <- 1
sol.betta.1 <- sol.irbc.iterate( l.betta[[5]]$sol.1$coeff, opt, params, l.betta[[5]]$sol.1$coeff.cont )
rep.betta.1 <- report.data( sol.betta.1 )
print( paste0( "err = ", round( max(apply( abs( rep.betta.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.betta.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.betta.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 65
opt$n.gain <- .1
sol.betta.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.betta.2 <- report.data( sol.betta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.betta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.betta[[6]] <- list( sol.1=sol.betta.1, rep.1=rep.betta.1, sol.2=sol.betta.2, rep.2=rep.betta.2 )
    # Assign to the list

### 2.5 betta = .85 ###
params$betta <- .85
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .075
opt$iter <- 16
opt$N <- 1
sol.betta.1 <- sol.irbc.iterate( l.betta[[6]]$sol.1$coeff, opt, params, l.betta[[6]]$sol.1$coeff.cont )
rep.betta.1 <- report.data( sol.betta.1 )
print( paste0( "err = ", round( max(apply( abs( rep.betta.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.betta.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.betta.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 65
opt$n.gain <- .1
sol.betta.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.betta.2 <- report.data( sol.betta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.betta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.betta[[7]] <- list( sol.1=sol.betta.1, rep.1=rep.betta.1, sol.2=sol.betta.2, rep.2=rep.betta.2 )

### 2.5 betta = .8 ###
params$betta <- .8
    # Change parameters
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$n.gain <- .075
opt$iter <- 14
opt$N <- 1
sol.betta.1 <- sol.irbc.iterate( l.betta[[7]]$sol.1$coeff, opt, params, l.betta[[7]]$sol.1$coeff.cont )
rep.betta.1 <- report.data( sol.betta.1 )
print( paste0( "err = ", round( max(apply( abs( rep.betta.1$err ), 2, mean )) * 100, 4), "pp" ) )
    # The linear solution
opt$N <- 2
coeff.init.2[ c(1,2,4,7,11), ] <- sol.betta.1$coeff
coeff.cont.init.2[ c(1,2,4,7,11), ] <- sol.betta.1$coeff.cont
opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$iter <- 65
opt$n.gain <- .1
sol.betta.2 <- sol.irbc.iterate( coeff.init.2, opt, params, coeff.cont.init.2 )
rep.betta.2 <- report.data( sol.betta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.betta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # The nonlinear solution
l.betta[[8]] <- list( sol.1=sol.betta.1, rep.1=rep.betta.1, sol.2=sol.betta.2, rep.2=rep.betta.2 )

save( params, opt, sol.base.1, sol.base.2, rep.base.1, rep.base.2, l.betta, file='~/Dropbox/outsize/irbc/betta_change.rdata')


