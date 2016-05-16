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
  f1 <- kde2d(exp( prs.sim$cont[,'c.1'] ), q.check, n=500 )
  filled.contour( f1, color.palette =  colorRampPalette(c("white", "black")),
      xlab='Consumption', ylab='Real ex. rate')
dev.off()

pdf( '~/Dropbox//2016//Research//IRBC puzzles/Paper/graphs/baselineCharts_20160429/clouds/c_e_CM_diff.pdf')
plot( diff( prs.sim$cont[,'c.1'] ), 
      diff(log(q.check)), pch=19, col= alpha('black', 0.15),
      xlab='Consumption growth', ylab='Real ex. rate growth')
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
