## Parameters
rho <- diag(c(.95,.95,.9,.9))
sig <- diag(c(.000016, .000016, .000031, .000031))
base.theta <- .002
params <- list( share = .86, gamma = 1.5, P1.bar=1, P2.bar=1, betta=.99,
                rho=rho, sig.eps=sig, eta=3.5, theta=base.theta, mu=.55, xi=.44 )

v.theta <- .01 * 1.5 ^ (seq(0,-10,length.out=20))

baseline.stats <- sapply( v.theta, function(theta){
  params$theta <- theta ;
  return( mod.gen(params, sim.stats=TRUE )  )
}  )


params$theta <- base.theta

opt.lim <- list( n.gain=c(.25,.1,.25,.1), iter=c(6,10,6,10), n.iter=3, kappa=100, tol=1e-06 )

AA <- mod.eval(params, opt.lim=opt.lim, return.sims.err = TRUE, n.sim = 200000 )

# params.b <- params
# params.b$gamma <- 8
# BB <- mod.eval(params.b, opt.lim=opt.lim, return.sims.err = TRUE, n.sim = 200000 )

par(mfrow=c(1,3))
plot( log(abs(apply(AA$err.abs.sumy[fwd.vars,],2,mean)),10), 
      log(abs(apply(AA$err.ave.sumy[fwd.vars,],2,mean)),10), 
      pch=1:5, cex=2, xlab='Average absolute error', ylab='Average error',
      main='Euler equations')
abline(h=0,lty=2)
# legend( 'topright', 
#         c( 'Devreux-Sutherland', 'Global Linear EDF', 'Global Quadratic EDF',
#            'Global Linear', 'Global Quadratic' ), pch=1:5, bty='n' )

plot( log(abs(AA$err.abs.sumy['B22',]),10), 
      log(abs(AA$err.ave.sumy['B22',]),10),
      pch=1:5, cex=2, xlab='Average absolute error', 
      ylab='Average error', main='Debt law of motion' )
legend( 'bottomright',
        c( 'Devreux-Sutherland', 'Global Linear EDF', 'Global Quadratic EDF',
           'Global Linear', 'Global Quadratic' ), pch=1:5, bty='n' )
abline(h=0,lty=2)

plot( log(abs(apply(AA$err.abs.sumy[-c(1,2,5,6,15),],2,mean)),10), 
      log(abs(apply(AA$err.ave.sumy[-c(1,2,5,6,15),],2,mean)),10), 
      pch=1:5, cex=2, xlab='Average absolute error', ylab='Average error',
      main='Contemporaneous equations')
abline(h=0,lty=2)
# legend( 'topright', 
#         c( 'Devreux-Sutherland', 'Global Linear EDF', 'Global Quadratic EDF',
#            'Global Linear', 'Global Quadratic' ), pch=1:5, bty='n' )
par(mfrow=c(1,1))

par(mfrow=c(1,2))
  plot( v.theta, baseline.stats['bs.basic',], xlab=expression(theta), 
        ylab='Backus-Smith coefficient', lwd=2, type='l')
  abline( h=AA$bs.log['global.2.0'], lwd=2, lty=2 )
  abline( v=base.theta, lty=3 )
  legend('topright', c('Devreux-Sutherland', 'Global quadratic'),
         lwd=2, lty=c(1,2), bty='n' )
  
  plot( v.theta, baseline.stats['b11.sd',], xlab=expression(theta), 
        ylab='Standard deviation of domestic assets', lwd=2, type='l')
  abline( h=sd(AA$l.sim$global.2.0[,'B11']), lwd=2, lty=2 )
  abline( v=base.theta, lty=3 )
  legend('topright', c('Devreux-Sutherland', 'Global quadratic'),
         lwd=2, lty=c(1,2), bty='n' )
par(mfrow=c(1,1))