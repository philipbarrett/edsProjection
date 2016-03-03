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
opt$n.gain <- .1
opt$c.iter <- 2000
opt$k.tol <- 1e-04
opt$k.iter <- 10
opt$iter <- 2
opt$sym.reg <- TRUE
opt$sr <- TRUE

sd.x <- params$sig.eps / sqrt( 1 - params$rho ^ 2 )
opt$upper[1:2] <- 5 * sd.x
opt$lower <- - opt$upper

opt$l.sym.ave <- list( sym=list( c(2,4), c(3,6), c(7,11), c(8,13), c(9,12), c(10,15) ),
                   ave=c(1,5,14) )
opt$l.pairs <- list( c(1,2) )
opt$l.pairs.cont <- list( c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12) )
opt$l.solo.cont <- c(13)
    # Enforce symmetry

params$gamma <- 6.5
sol.irbc.2 <- sol.irbc.iterate( sol.irbc$coeff / 10, opt, params, sol.irbc$coeff.cont / 10 )
rep.irbc.2 <- report.data( sol.irbc.2 )
report.create( sol.irbc.2, rep.irbc.2 )

params$gamma <- 8
sol.irbc.3 <- sol.irbc.iterate( sol.irbc.2$coeff, opt, params, sol.irbc.2$coeff.cont )
rep.irbc.3 <- report.data( sol.irbc.3 )
report.create( sol.irbc.3, rep.irbc.3 )

params$sig.eps <- .1
sd.a <- params$sig.eps / sqrt( 1 - params$rho ^ 2 )
opt$upper
