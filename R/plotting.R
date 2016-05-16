#########################################################################
# plotting.R
#
# Plotting for the EDS-projection algorithm
# Philip Barrett, Chicago
# 11feb2016
#########################################################################

plot.coeffs <- function( coeff, ... ){
  barplot(t(coeff), beside=T, col=1:ncol(coeff), ... )
}

cloud.density <- function( X, Y, ... ){
  f1 <- kde2d( X, Y, n=250 )
  filled.contour(f1, color.palette =  colorRampPalette(c("white", "black") ), ... )
}

paper.charts <- function( sol, rep.data, 
    loc='/home/philip/Dropbox/2016/Research/IRBC puzzles/Paper/graphs/' ){
  
  tot.assets <- rep.data$endog.sim[,3] -
    rep.data$endog.sim[,4] * exp( rep.data$cont.sim[,13] )
  asset.diff <- rep.data$endog.sim[,3] + 
    rep.data$endog.sim[,4] * exp( rep.data$cont.sim[,13] )
  
  pdf(paste0( loc, 'debt_dist.pdf') )
    plot( density( rep.data$endog.sim[,3] ), col='red', xlab='Assets', ylab='Density', 
          main='', xlim=c(-.5,.5), lwd=2 )
    lines( density( - rep.data$endog.sim[,4] * exp( rep.data$cont.sim[,13] ) ),col='blue', lwd=2 )
    lines( density( tot.assets ),col='black', lwd=2 )
    abline( v=0, lty=2 )
    legend( 'topright', 
            c('Domestic assets','Foreign assets', 'Total assets' ), 
            lty=1, col=c('red','blue', 'black'), bty='n', lwd=2 )
  dev.off()
  
  pdf(paste0( loc, 'debt_dist_diff.pdf') )
    plot( density( asset.diff ), col='black', xlab='Domestic less foreign assets', 
          ylab='Density', main='', lwd=2 )
    abline( v=0, lty=2 )
  dev.off()
  
  qq <- - rep.data$cont.sim[,9] + rep.data$cont.sim[,10] + rep.data$cont.sim[,13]
  
  pdf( paste0( loc, 'cq_corr.pdf') )
    cloud.density( exp( rep.data$cont.sim[,1] ), exp( qq ), xlab='Consumption', 
                   ylab='Real ex rate', xlim=c(.47, .53) )
  dev.off()
  
  pdf( paste0( loc, 'cq_corr_diff.pdf') )
    cloud.density( diff( rep.data$cont.sim[,1] ), diff( qq ), xlim=c(-.04,.04), 
                  xlab='Consumption growth', ylab='Real ex rate growth' )
  dev.off()
  
  # pdf( paste0( loc, 'ce_corr_diff.pdf') )
  #   cloud.density( diff( rep.data$endog.sim[,1] ), diff( rep.data$endog.sim[,13] ), 
  #                  xlab='Consumption growth', ylab='Nominal ex rate growth' )
  # dev.off()

}