rm(list=ls())
setwd("~/code/2016/edsProjection")
library(edsProjection)
library(xtable)
load('27feb2016.Rdata')
opt <- sol.irbc$opt
params <- sol.irbc$params
opt$model <- 'irbc'

opt$kappa <- 15
opt$c.gain <- .1
opt$n.gain <- .02
opt$c.iter <- 2000
opt$k.tol <- 1e-04
opt$k.iter <- 10
opt$iter <- 10
opt$sym.reg <- TRUE
opt$sr <- TRUE

sd.x <- params$sig.eps / sqrt( 1 - params$rho ^ 2 )
opt$upper[1:2] <- 5 * sd.x
opt$lower <- - opt$upper
    # Fix the range of the solutions

opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                       ave=c(1,5,14) )
opt$l.pairs <- list( c(1,2) )
opt$l.pairs.cont <- list( c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12) )
opt$l.solo.cont <- c(13)
    # Enforce symmetry

## Create the initial guesses
coeff.init <- 0 * sol.irbc$coeff
coeff.init.cont <- 0 * sol.irbc$coeff.cont
for( i in 1:2) 
  coeff.init[,i] <- coeff_reg( sol.irbc$X.cont[,2+i], sol.irbc$X.cont[,c(1,2,7,8)], 
                               opt$N, opt$lower, opt$upper, opt$cheby )
for( i in 1:13) 
  coeff.init.cont[,i] <- coeff_reg( sol.irbc$X.cont[,8+i], sol.irbc$X.cont[,c(1,2,7,8)], 
                                    opt$N, opt$lower, opt$upper, opt$cheby )

# ## Profiling
# Rprof(NULL)
# opt$iter <- 1
# Rprof('03mar.Rprof')
# sol.irbc.2 <- sol.irbc.iterate( coeff.init, opt, params, coeff.init.cont )
# Rprof(NULL)
# summaryRprof('03mar.Rprof')
# stop()

## Solve the model
sol.irbc.2 <- sol.irbc.iterate( coeff.init, opt, params, coeff.init.cont )
rep.irbc.2 <- report.data( sol.irbc.2 )
report.create( sol.irbc.2, rep.irbc.2 )

## Increase the schock variance
params$sig.eps <- c( .02, .02 )
sd.x <- params$sig.eps / sqrt( 1 - params$rho ^ 2 )
opt$upper[1:2] <- 5 * sd.x
opt$lower <- - opt$upper
opt$iter <- 30
sol.irbc.3 <- sol.irbc.iterate( sol.irbc.2$coeff, opt, params, sol.irbc.2$coeff.cont )
opt$iter <- 50
sol.irbc.3 <- sol.irbc.iterate( sol.irbc.3$coeff, opt, params, sol.irbc.3$coeff.cont )
opt$n.gain <- .05
opt$iter <- 10
sol.irbc.3 <- sol.irbc.iterate( sol.irbc.3$coeff, opt, params, sol.irbc.3$coeff.cont )
## STILL NOT QUITE CONVERGED
rep.irbc.3 <- report.data( sol.irbc.3 )
report.create( sol.irbc.3, rep.irbc.3 )
