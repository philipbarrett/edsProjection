

load('~/Dropbox/outsize/irbc/eta_high.rdata')

params <- l.sol[[4]]$sol.2$params
opt <- l.sol[[4]]$sol.2$opt
    # Establish parameteres and options
coeff.init <- l.sol[[4]]$sol.2$coeff[ c(1,2,4,7,11), ]
coeff.cont.init <- l.sol[[4]]$sol.2$coeff.cont[ c(1,2,4,7,11), ]
opt$N <- 1

opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5) ), ave=c(1) )
opt$l.pairs <- list( c(1,2) )
opt$l.pairs.cont <- list( c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12) )
opt$sym.reg <- TRUE
opt$iter <- 2

sol.1 <- sol.irbc.iterate( coeff.init, opt, params, coeff.cont.init, debug.flag = F )
rep.1 <- report.data( sol.1 )
print( paste0( "err = ", round( max(apply( abs( rep.1$err ), 2, mean )) * 100, 4), "pp" ) )

ql.1 <- quantile( apply( 100 * 100 * abs(rep.1$err), 1, max ), c( .5, .75, .9, .95, .99 ))
pdf( '~/Dropbox//2016//Research//IRBC puzzles/Paper/graphs/baselineCharts_20160429/err_plot_linear.pdf')
  plot( density( 100^2 * apply( abs(rep.1$err), 1, max ), from=1e-08, n=512*8 ),
        xlim=c(0, 300), main='', xlab='Basis points (pp/100)', lwd=2, col='red' )
  abline( v=ql, lty=2 )
  text( ql.1, .011, paste0( 'p=', c(0.5, 0.75, 0.9, 0.95, 0.99 ) ), srt=90, pos=2 )
dev.off()
