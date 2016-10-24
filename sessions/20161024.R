params <- list( share = .86, gamma = 2, P1.bar=1, P2.bar=1, betta=.95,
                rho=c(.9,.9, .9, .9), sig.eps=c(.025,.025, .01, .01), eta=2, theta=.05 )

# nn.1 <- 5
# nn.2 <- 10
# nn <- nn.1 + nn.2
# eta.1.range <- c( 1.001, 1.09)
# eta.2.range <- c( 1.1, 5 )
# eta.range <- c( eta.1.range, eta.2.range )
# 
# v.eta <- c( seq( eta.1.range[1], eta.range[2], length.out=nn.1),
#             seq( eta.2.range[1], eta.2.range[2], length.out = nn.2 ) )

nn <- 8
v.sigma.m <- seq(0,.025, length.out=nn)

all <- matrix( 0, nn, 5 )

for( i.eta in 1:nn ){
  
#   eta <- v.eta[i.eta]
  sig.m <- v.sigma.m[i.eta]
  
  message('\n****************************')
#   message('***  eta = ', eta, ' **************')
  message('***  sig.m = ', sig.m, ' **************')
  message('****************************')
  
#   params$eta <- eta
  params$sig.eps <- c( params$sig.eps[1:2], rep( sig.m,2) )
  temp <- mod.gen(params, sim.stats = TRUE)
  print(temp)
  all[i.eta,] <- temp
  
}

colnames(all) <- names(temp)

plot( v.eta, all[,'alpha.tilde'], type='l', lwd=2)
plot( v.eta, all[,'bs.basic'], type='l', lwd=2)

plot( v.sigma.m, all[,'alpha.tilde'], type='l', lwd=2)
plot( v.sigma.m, all[,'bs.basic'], type='l', lwd=2)


