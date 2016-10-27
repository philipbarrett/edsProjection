params <- list( share = .86, gamma = 2, P1.bar=1, P2.bar=1, betta=.95,
                rho=c(.9,.9, .9, .9), sig.eps=c(.01,.01, .0025, .0025), eta=2, 
                theta=.025 )
params.orig <- params

nn.1 <- 5
nn.2 <- 10
nn <- nn.1 + nn.2
eta.1.range <- c( 1.001, 1.09)
eta.2.range <- c( 1.1, 5 )
eta.range <- c( eta.1.range, eta.2.range )

v.eta <- c( seq( eta.1.range[1], eta.range[2], length.out=nn.1),
            seq( eta.2.range[1], eta.2.range[2], length.out = nn.2 ) )

baseline <- mod.gen(params, sim.stats = TRUE)
params.cm <- params
params.cm$sig.eps[3:4] <- 0
cm <- mod.gen(params.cm, sim.stats = TRUE)

baseline.sim <- mod.gen(params, return.sim = TRUE)
cm.sim <- mod.gen(params.cm, return.sim = TRUE)

h <- hexbin(diff(baseline.sim[,'C1']-baseline.sim[,'C2']),diff(baseline.sim[,'Q']), xbins=50)
g <- hexbin(diff(cm.sim[,'C1']-cm.sim[,'C2']),diff(cm.sim[,'Q']), xbins=50)
i <- hexbin(exp(baseline.sim[,'C1'])-exp(baseline.sim[,'C2']),exp(baseline.sim[,'Q']), xbins=50)
j <- hexbin(exp(cm.sim[,'C1'])-exp(cm.sim[,'C2']),exp(cm.sim[,'Q']), xbins=50)
pdf('~/Dropbox//2016/Research/IRBC puzzles/Paper/graphs/bs_dens_florida.pdf')
  plot(h, legend=F, xlab='Consumption growth differential', ylab='Real exchange rate appreciation')
dev.off()
pdf('~/Dropbox//2016/Research/IRBC puzzles/Paper/graphs/bs_dens_florida_cm.pdf')
  plot(g, legend=F, xlab='Consumption growth differential', ylab='Real exchange rate appreciation')
dev.off()
pdf('~/Dropbox//2016/Research/IRBC puzzles/Paper/graphs/bs_dens_florida_levels.pdf')
  plot(i, legend=F, xlab='Consumption growth differential', ylab='Real exchange rate appreciation')
dev.off()
pdf('~/Dropbox//2016/Research/IRBC puzzles/Paper/graphs/bs_dens_florida_cm_levels.pdf')
  plot(j, legend=F, xlab='Consumption growth differential', ylab='Real exchange rate appreciation')
dev.off()
B11 <- baseline.sim[-nrow(baseline.sim),'af1'] * exp( baseline.sim[-1,'rb1'] + baseline.sim[-1,'P1'] )
pdf('~/Dropbox//2016/Research/IRBC puzzles/Paper/graphs/bond_density_florida.pdf')
  plot(density(B11), lwd=2, main='', xlab='Domestic bond holdings')
dev.off()

eta.stats <- matrix( 0, nn, 6 )

for( i.eta in 1:nn ){
  
  eta <- v.eta[i.eta]
  
  message('\n****************************')
  message('***  eta = ', eta, ' **************')
  message('****************************')
  
  params$eta <- eta
  temp <- mod.gen(params, sim.stats = TRUE)
  print(temp)
  eta.stats[i.eta,] <- temp
  
}

colnames(eta.stats) <- names(temp)

pdf('~/Dropbox//2016/Research/IRBC puzzles/Paper/graphs/alpha_tilde_florida.pdf')
  plot( v.eta, params$betta * eta.stats[, 'alpha.tilde'], type='l', lwd=2, 
      xlab=expression(eta), ylab='Mean domestic bond holdings' )
  abline(v=2,lty=2)
dev.off()

pdf('~/Dropbox//2016/Research/IRBC puzzles/Paper/graphs/bs_eta_florida.pdf')
  plot( v.eta, eta.stats[, 'bs.basic'], type='l', lwd=2, 
      xlab=expression(eta), ylab='Backus-Smith correlation' )
  abline(v=2,lty=2)
  abline(h=0)
dev.off()

nn <- 20
v.sigma.m <- seq(0,.02, length.out=nn)
params <- params.orig
sig.stats <- matrix( 0, nn, 6 )

for( i in 1:nn ){

  sig.m <- v.sigma.m[i]
  
  message('\n****************************')
  message('***  sig.m = ', sig.m, ' **************')
  message('****************************')
  
  params$sig.eps <- c( params$sig.eps[1:2], rep( sig.m,2) )
  temp <- mod.gen(params, sim.stats = TRUE)
  print(temp)
  sig.stats[i,] <- temp
}

colnames(sig.stats) <- names(temp)

pdf('~/Dropbox//2016/Research/IRBC puzzles/Paper/graphs/bs_sigma_m_florida.pdf')
  plot( v.sigma.m, sig.stats[, 'bs.basic'], type='l', lwd=2, 
      xlab=expression(sigma), ylab='Backus-Smith correlation' )
  abline(v=.0025,lty=2)
  abline(h=0)
dev.off()

nn <- 19
v.gamma <- seq(.5,5, length.out=nn)
params <- params.orig
gam.stats <- matrix( 0, nn, 6 )

for( i in 1:nn ){
  
  gam <- v.gamma[i]
  
  message('\n****************************')
  message('***  gamma = ', gam, ' **************')
  message('****************************')
  
  params$gamma <- gam
  temp <- mod.gen(params, sim.stats = TRUE)
  print(temp)
  gam.stats[i,] <- temp
}

colnames(gam.stats) <- names(temp)

pdf('~/Dropbox//2016/Research/IRBC puzzles/Paper/graphs/bs_gamma_florida.pdf')
  plot( v.gamma[v.gamma>=1], gam.stats[v.gamma>=1, 'bs.basic'], type='l', lwd=2, 
        xlab=expression(gamma), ylab='B' )
  abline(v=2,lty=2)
  abline(h=0)
dev.off()


pdf('~/Dropbox//2016/Research/IRBC puzzles/Paper/graphs/alpha_tilde_gam_florida.pdf')
  plot( v.gamma[v.gamma>=1], params$betta * gam.stats[v.gamma>=1, 'alpha.tilde'], type='l', lwd=2, 
        xlab=expression(gamma), ylab='Mean domestic bond holdings' )
  abline(v=2,lty=2)
dev.off()

# plot( v.sigma.m, all[,'alpha.tilde'], type='l', lwd=2)
# plot( v.sigma.m, all[,'bs.basic'], type='l', lwd=2)


