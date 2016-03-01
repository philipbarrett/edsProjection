rm(list=ls())
library(xtable)
load('27feb2016.Rdata')
opt <- sol.irbc$opt
params <- sol.irbc$params
opt$model <- 'irbc'
opt$kappa <- 15
opt$iter <- 50

params$gamma <- 6.5
sol.irbc.2 <- sol.irbc.iterate( sol.irbc$coeff, opt, params, sol.irbc$coeff.cont )
rep.irbc.2 <- report.data( sol.irbc.2 )
report.create( sol.irbc.2, rep.irbc.2 )

params$gamma <- 8
sol.irbc.3 <- sol.irbc.iterate( sol.irbc.2$coeff, opt, params, sol.irbc.2$coeff.cont )
rep.irbc.3 <- report.data( sol.irbc.3 )
report.create( sol.irbc.3, rep.irbc.3 )
