## Parameters
rho <- diag(c(.95,.95,.9,.9))
sig <- diag(c(.000016, .000016, .000031, .000031))
params <- list( share = .86, gamma = 1.5, P1.bar=1, P2.bar=1, betta=.99,
                rho=rho, sig.eps=sig, eta=3.5, theta=0.001, mu=.55, xi=.44 )

opt.lim <- list( n.gain=c(.25,.1,.25,.1), iter=c(6,6,6,6), n.iter=3, kappa=100, tol=1e-06 )

AA <- mod.eval(params, opt.lim=opt.lim, return.sims.err = TRUE, n.sim = 200000 )

# params.b <- params
# params.b$gamma <- 8
# BB <- mod.eval(params.b, opt.lim=opt.lim, return.sims.err = TRUE, n.sim = 200000 )

par(mfrow=c(1,3))
plot( apply(AA$err.abs.sumy[fwd.vars,],2,mean), apply(AA$err.ave.sumy[fwd.vars,],2,mean), 
      pch=1:5, cex=2, xlab='Average absolute error', ylab='Average error',
      main='Euler equation')
legend( 'topright', 
        c( 'Devreux-Sutherland', 'Global Linear EDF', 'Global Quadratic EDF',
           'Global Linear', 'Global Quadratic' ), pch=1:5, bty='n' )

plot( AA$err.abs.sumy['B22',], AA$err.ave.sumy['B22',], 
      pch=1:5, cex=2, xlab='Average absolute error', 
      ylab='Average error', main='Debt law of motion' )
# legend( 'right', 
#         c( 'Devreux-Sutherland', 'Global Linear EDF', 'Global Quadratic EDF',
#            'Global Linear', 'Global Quadratic' ), pch=1:5, bty='n' )

plot( apply(AA$err.abs.sumy[-c(1,2,5,6,15),],2,mean), 
      apply(AA$err.ave.sumy[-c(1,2,5,6,15),],2,mean), 
      pch=1:5, cex=2, xlab='Average absolute error', ylab='Average error',
      main='Contemporaneous equations')
# legend( 'topright', 
#         c( 'Devreux-Sutherland', 'Global Linear EDF', 'Global Quadratic EDF',
#            'Global Linear', 'Global Quadratic' ), pch=1:5, bty='n' )
par(mfrow=c(1,1))