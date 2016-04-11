
setwd("/home/philip/Dropbox/2016/Research/IRBC puzzles/reports/param_experiments")

plot.stat <- function( st.param, stats, st.var, y.lab, file.name ){
  st.name <- paste0( st.param, '/', file.name, '.pdf' )
  pdf(st.name)
  plot( stats[st.param,], stats[st.var,], lwd=2, type='l', xlab=st.param, ylab = y.lab )
  dev.off()
}

v.st.var <- c( 'err.pp', 'cor.q.c.gth', 'cor.q.c', 'home.bias.coeff.NA', 'asset.diff',
               'abs.tot.assets', 'uip.coeff.r.diff' )
v.y.lab <- c( 'Maximum average equation error (pp)',
              'Correlation of q and c (growth)',
              'Correlation of q and c (level)', 
              'Home bias regression coefficient',
              'Mean net domestic assets (unscaled)',
              'Mean absolute total assets',
              'UIP regression coefficient' )
v.file.name <- c( 'err', 'corQCgth', 'corQC', 'homeBias', 'netDomAsset', 'meanAbsAsset', 'uipReg' )

load('/home/philip/Dropbox/outsize/irbc/eta_change.rdata')
stats <- key.stats( l.eta, sol.base.2, rep.base.2, 'eta' )
sapply( 1:length(v.st.var), function(i) plot.stat( 'eta', stats, v.st.var[i], 
                                                   v.y.lab[i], v.file.name[i] ) )

load('/home/philip/Dropbox/outsize/irbc/gamma_change.rdata')
stats <- key.stats( l.gamma, sol.base.2, rep.base.2, 'gamma' )
sapply( 1:length(v.st.var), function(i) plot.stat( 'gamma', stats, v.st.var[i], 
                                                   v.y.lab[i], v.file.name[i] ) )


load('/home/philip/Dropbox/outsize/irbc/keyStats.rdata')
load('/home/philip/Dropbox/outsize/irbc/rho_change.rdata')
stats <- key.stats( l.rho, sol.base.2, rep.base.2, 'rho' )
sapply( 1:length(v.st.var), function(i) plot.stat( 'rho', stats, v.st.var[i], 
                                                   v.y.lab[i], v.file.name[i] ) )

load('/home/philip/Dropbox/outsize/irbc/betta_change.rdata')
stats <- key.stats( lapply( c(1:3,5:8), function(i) l.betta[[i]] ), sol.base.2, rep.base.2, 'betta' )
sapply( 1:length(v.st.var), function(i) plot.stat( 'betta', stats, v.st.var[i], 
                                                   v.y.lab[i], v.file.name[i] ) )
