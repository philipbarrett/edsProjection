rm(list=ls())

params <- list( share = .75, gamma = 4, P1.bar=1, P2.bar=1, betta=.95,
                rho=c(.5,.5), sig.eps=c(.01,.01), eta=2.5 )
n.burn <- 1e05
n.pds <- 1e05
n.sample <- 100

prs.sim <- prs.sol.sim( n.pds, n.sample, n.burn, params )

q <- exp( prs.sim$cont[,'e.12'] + prs.sim$cont[,'p.2'] - prs.sim$cont[,'p.1'] )
q.check <- ( exp( prs.sim$cont[,'c.2'] - prs.sim$cont[,'c.1'] ) ) ^ ( - params$gamma )
    # These are not the same.  This is a little concerning.

cor( q, exp(prs.sim$cont[,'c.1'] ) )
cor( q.check, exp( prs.sim$cont[,'c.1'] ) )
cor( diff(log(q.check)), diff(prs.sim$cont[,'c.1'] ) )
cor( diff(log(q)), diff(prs.sim$cont[,'c.1'] ) )

cor( 1 / q.check, prs.sim$cont[,'c.2'] )


pdf( '~/Dropbox//2016//Research//IRBC puzzles/Paper/graphs/baselineCharts_20160429/clouds/c_e_CM.pdf')
  f1 <- kde2d(exp( prs.sim$cont[,'c.1'] ), q.check, n=250 )
  filled.contour( f1, color.palette =  colorRampPalette(c("white", "black")),
      xlab='Consumption', ylab='Real ex. rate', ylim=c(.9,1.1))
dev.off()

pdf( '~/Dropbox//2016//Research//IRBC puzzles/Paper/graphs/baselineCharts_20160429/clouds/c_e_CM_diff.pdf')
  f1 <- kde2d( diff( prs.sim$cont[,'c.1'] ), diff(log(q.check)), n=250 )
  filled.contour( f1, color.palette =  colorRampPalette(c("white", "black")),
                  xlab='Consumption growth', ylab='Real ex. rate growth', 
                  xlim=c(-.04,.04) )
dev.off()


load('~/Dropbox/outsize/irbc/eta_high.rdata')

q.sim <- l.sol[[4]]$rep.2$cont.sim[,13] - l.sol[[4]]$rep.2$cont.sim[,9] + l.sol[[4]]$rep.2$cont.sim[,10]

# plot( l.sol[[4]]$rep.2$cont.sim[,1], q.sim, pch=19, col= alpha('black', 0.15) )

# pdf( '~/Dropbox//2016//Research//IRBC puzzles/Paper/graphs/baselineCharts_20160429/clouds/c_e_levels.pdf')
# plot( exp( l.sol[[4]]$rep.2$cont.sim[,1] ), 
#       exp( q.sim ), pch=19, col= alpha('black', 0.15),
#       xlab='Consumption', ylab='Real ex. rate')
# dev.off()

cor(exp( l.sol[[4]]$rep.2$cont.sim[,1] ), exp( q.sim ) )

cor( q.sim[-1], q.sim[-length(q.sim)] )
cor( log( q.check[-1] ), log( q.check[-length(q.check)] ) )
    # RER autocorrel
cor( l.sol[[4]]$rep.2$cont.sim[-1,13], 
     l.sol[[4]]$rep.2$cont.sim[-length(q.sim),13] )
cor( prs.sim$cont[-1,'e.12'], 
     prs.sim$cont[-nrow(prs.sim$cont),'e.12'])
    # NER autocorrel
cor( q.sim, l.sol[[4]]$rep.2$cont.sim[,13] )
cor( q.check, prs.sim$cont[,'e.12'] )
    # RER, NER correl.  (Worrying)

sd(prs.sim$cont[,'c.1'] ) * 100
sd(l.sol[[4]]$rep.2$cont.sim[,1] ) * 100


ql <- quantile( apply( 100 * 100 * abs(l.sol[[4]]$rep.2$err), 1, max ), c( .5, .75, .9, .95, .99 ))

pdf( '~/Dropbox//2016//Research//IRBC puzzles/Paper/graphs/baselineCharts_20160429/err_plot.pdf')
  plot( density( 100^2 * apply( abs(l.sol[[4]]$rep.2$err), 1, max ), n=512*8, from=1e-08 ), 
        xlim=c(0, 5), main='', xlab='Basis points (pp/100)', lwd=2 )
  abline( v=ql, lty=2 )
  text( ql, .4, paste0( 'p=', c(0.5, 0.75, 0.9, 0.95, 0.99 ) ), srt=90, pos=2 )
dev.off()
