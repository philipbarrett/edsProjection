params <- list( share = .86, gamma = 2, P1.bar=1, P2.bar=1, betta=.95,
                rho=c(.9,.9), sig.eps=c(.025,.025), eta=2 )


nn.1 <- 10
nn.2 <- 20
nn <- nn.1 + nn.2
eta.1.range <- c( 1.001, 1.09)
eta.2.range <- c( 1.1, 12 )
eta.range <- c( eta.1.range, eta.2.range )

v.eta <- c( seq( eta.1.range[1], eta.range[2], length.out=nn.1),
            seq( eta.2.range[1], eta.2.range[2], length.out = nn.2 ) )

all <- matrix( 0, nn, 11 )
    # All the statistics
colnames(all) <- c('af1', 'ds.bias', 'ds.aad', 'ds.mean.max.abs.err', 'nl.bias', 'nl.aad',
                   'nl.mean.max.abs.err', 
                   'ds.bs', 'nl.bs', 'ds.uip', 'nl.uip' )

for( i.eta in 1:nn ){
  
  eta <- v.eta[i.eta]
  
  message('\n****************************')
  message('***  eta = ', eta, ' **************')
  message('****************************')
  
  params$eta <- eta
  all[i.eta,] <- stat.ds(params, "all")
}



# plot( v.eta, all[,'ds.bs'], type='l', lwd=2, col='red' )
# lines( v.eta, all[,'nl.bs'], lwd=2, col='blue' )

# ## And now UIP
# uip <- matrix( 0, nn, 2 )
# colnames(uip) <- c('ds','nl')
# for( i.eta in 1:nn ){
#   
#   eta <- v.eta[i.eta]
#   
#   message('\n****************************')
#   message('***  eta = ', eta, ' **************')
#   message('****************************')
#   
#   params$eta <- eta
#   temp <- stat.ds(params, "uip")
#   uip[i.eta,] <- c( temp$ds[2], temp$nl[2] )
# }
# 
# 
# omit <- c(1,11,16,24,28)
# pdf('~/Dropbox/2016/Research/IRBC puzzles/Paper/graphs/uip_eta.pdf')
#   plot( v.eta[-omit], uip[-omit,'ds'], type='l', lwd=2, col='red', xlab=expression(eta),
#         ylab='UIP regression coefficient')
#   lines( v.eta[-omit], uip[-omit,'nl'], type='l', lwd=2, col='blue' )
#   legend( 'topright', c('Devreux-Sutherland', 'Global'), lwd=2, col=c('red', 'blue'),
#           bty='n')
# dev.off()
# 
# ### NEXT BIT TEMPORARY: NEED TO RERUN FOR bs.cov ###
# v.eta.2 <- seq(2,16, length.out=19)
# omit.2 <- c(5, 15, 18)
# pdf('~/Dropbox/2016/Research/IRBC puzzles/Paper/graphs/bs_eta.pdf')
#   plot( v.eta.2[-omit.2], bs.cov[-omit.2,'nl'], type='l', lwd=2, col='blue', xlab=expression(eta),
#         ylab='RER-consumption correlation')
#   lines( v.eta.2[-omit.2], bs.cov[-omit.2,'ds'], type='l', lwd=2, col='red')
#   legend( 'topright', c('Devreux-Sutherland', 'Global'), lwd=2, col=c('red', 'blue'),
#           bty='n')
# dev.off()

