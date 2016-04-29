## Attempt to push up both gamma and eta togather

rm(list=ls())
library(xtable)
library(scales)
library(edsProjection)
setwd("~/code/2016/edsProjection")

load('/home/philip/Dropbox/outsize/irbc/gamma_change.rdata')

base <- l.gamma[[12]]
rm(l.gamma)
    # Retain only the highest-gamma solution

##### 1. CHECK THE BASELINE SOLUTION WORKS #####
opt <- base$sol.2$opt
params <- base$sol.2$params
opt$iter <- 1
sol.gamma.2 <- sol.irbc.iterate( base$sol.2$coeff, opt, params, base$sol.2$coeff.cont )
rep.gamma.2 <- report.data( sol.gamma.2 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # Looks good

l.sol <- list()
l.sol[[1]] <- list( sol.2=sol.gamma.2, rep.2=rep.gamma.2 )

##### 2. INCREASE GAMMA TO 4 #####
params$gamma <- 4
opt$iter <- 50
opt$k.tol <- 5e-05
sol.gamma.2 <- sol.irbc.iterate( base$sol.2$coeff, opt, params, base$sol.2$coeff.cont )
rep.gamma.2 <- report.data( sol.gamma.2 )
print( paste0( "err = ", round( max(apply( abs( rep.gamma.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # Great!  Error is nice and low
l.sol[[2]] <- list( sol.2=sol.gamma.2, rep.2=rep.gamma.2 )


##### 3. INCREASE ETA #####

### 3.0 Change settings for solution to make it work for increased eta ###
opt$adapt.gain <- FALSE
opt$n.gain <- .05 # .02
opt$c.gain <- .5 # .25
opt$c.iter <- 100
opt$c.tol <- 1e-09
opt$k.gain <- .5
opt$k.iter <- 40
opt$k.tol <- 1e-06
opt$tol <- 1e-06

### 3.1 Eta = 2 ###
params$eta <- 2
opt$iter <- 50
sol.eta.2 <- sol.irbc.iterate( l.sol[[2]]$sol.2$coeff, opt, params, l.sol[[2]]$sol.2$coeff.cont )
rep.eta.2 <- report.data( sol.eta.2 )
opt$n.gain <- .02
opt$iter <- 20
sol.eta.2 <- sol.irbc.iterate( sol.eta.2$coeff, opt, params, sol.eta.2$coeff.cont )
rep.eta.2 <- report.data( sol.eta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.eta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # It is harder to get convergence here - but I did it!
l.sol[[3]] <- list( sol.2=sol.eta.2, rep.2=rep.eta.2 )


### 3.1 Eta = 2.5 ###
params$eta <- 2.5
opt$n.gain <- .05
opt$iter <- 40
sol.eta.2 <- sol.irbc.iterate( l.sol[[3]]$sol.2$coeff, opt, params, l.sol[[3]]$sol.2$coeff.cont )
rep.eta.2 <- report.data( sol.eta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.eta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # Ave abs err < 1bp
l.sol[[4]] <- list( sol.2=sol.eta.2, rep.2=rep.eta.2 )

### 3.2 Eta = 3 ###
params$eta <- 3
opt$iter <- 13
sol.eta.2 <- sol.irbc.iterate( l.sol[[4]]$sol.2$coeff, opt, params, l.sol[[4]]$sol.2$coeff.cont )
rep.eta.2 <- report.data( sol.eta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.eta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # Ave abs err < 1bp
l.sol[[5]] <- list( sol.2=sol.eta.2, rep.2=rep.eta.2 )

save( l.sol, file='~/Dropbox/outsize/irbc/eta_high.rdata')

### 3.3 Eta = 3.5 ###
params$eta <- 3.5
opt$iter <- 50
sol.eta.2 <- sol.irbc.iterate( l.sol[[5]]$sol.2$coeff, opt, params, l.sol[[5]]$sol.2$coeff.cont )
rep.eta.2 <- report.data( sol.eta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.eta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # Ave abs err < 1bp
l.sol[[6]] <- list( sol.2=sol.eta.2, rep.2=rep.eta.2 )

### 3.4 Eta = 3.75 ###
params$eta <- 3.75
sol.eta.2 <- sol.irbc.iterate( l.sol[[6]]$sol.2$coeff, opt, params, l.sol[[6]]$sol.2$coeff.cont )
rep.eta.2 <- report.data( sol.eta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.eta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # Ave abs err = .6 bp
l.sol[[7]] <- list( sol.2=sol.eta.2, rep.2=rep.eta.2 )


##### 4. TRY INCREASING GAMMA SOME MORE #####

### 4.1 gamma=4.5 ###
params$gamma <- 4.5
opt$iter <- 80
sol.eta.2 <- sol.irbc.iterate( l.sol[[7]]$sol.2$coeff, opt, params, l.sol[[7]]$sol.2$coeff.cont )
rep.eta.2 <- report.data( sol.eta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.eta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # Ave abs err .7bp
l.sol[[8]] <- list( sol.2=sol.eta.2, rep.2=rep.eta.2 )

### 4.2 gamma=5 ###
params$gamma <- 5
sol.eta.2 <- sol.irbc.iterate( l.sol[[7]]$sol.2$coeff, opt, params, l.sol[[7]]$sol.2$coeff.cont )
rep.eta.2 <- report.data( sol.eta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.eta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # Ave abs err .7bp
l.sol[[9]] <- list( sol.2=sol.eta.2, rep.2=rep.eta.2 )


##### 5. NOW UP ETA A LITTLE #####

### 3.4 Eta = 4 ###
params$eta <- 4
sol.eta.2 <- sol.irbc.iterate( l.sol[[9]]$sol.2$coeff, opt, params, l.sol[[9]]$sol.2$coeff.cont )


rep.eta.2 <- report.data( sol.eta.2 )
print( paste0( "err = ", round( max(apply( abs( rep.eta.2$err ), 2, mean )) * 100, 4), "pp" ) )
    # Ave abs err = .6 bp
l.sol[[7]] <- list( sol.2=sol.eta.2, rep.2=rep.eta.2 )


### Other ideas: Increase rho or sig_eps




# ### TRY A 3rd ORDER SOLUTION
# params$eta <- 3.75
# opt$N <- 3
# opt$iter <- 5
# opt$n.gain <- .01
# idx.3 <- idx_create( 3, 4 )
# coeff.init <- matrix( 0, 35, 2 )
# coeff.init[ which( apply( idx.3, 1, sum ) <= 2 ), ] <- l.sol[[7]]$sol.2$coeff
# coeff.init.cont <- matrix( 0, 35, 13 )
# coeff.init.cont[ which( apply( idx.3, 1, sum ) <= 2 ), ] <- l.sol[[7]]$sol.2$coeff.cont
#     # Initialize the coefficients
# opt$l.sym.ave <- list( sym=list( c(2,5), c(3,8), c(4,10), c(7,9), c(8,13), c(11,21),
#                                  c(12,24), c(13,26), c(14,22), c(15,25), c(16,23),
#                                  c(17,31), c(18,33), c(19,32), c(20,35), c(28,29),
#                                  c(30,34) ),
#                        ave=c(1,6,27) )
# 
# sol.eta.3 <- sol.irbc.iterate( coeff.init, opt, params, coeff.init.cont )