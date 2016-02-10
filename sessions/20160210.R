#### Code to check that various versions of the NGM model work well


##### 1. Full depreciation case ####
params.full.dep <- list( A=1, alpha=.3, delta=1, gamma=1, betta=.99, 
                         rho=0, sig.eps=.01 )
params.full.dep$A <- ( 1 / params.full.dep$betta - 1 + 
                         params.full.dep$delta ) / params.full.dep$alpha
k.ss <- ( ( 1 / params.full.dep$betta - 1 + params.full.dep$delta ) / 
            ( params.full.dep$A * params.full.dep$alpha )  ) ^ ( 1 / ( params.full.dep$alpha - 1 ) )
# The steady state
sd.x<- sqrt( params.full.dep$sig.eps / ( 1 - params.full.dep$rho ^ 2 ) )
upper <- c(  3 * sd.x, k.ss * 1.4 )
lower <- c( -3 * sd.x, k.ss * 0.6 )
# Create the bounds
opt <- list( model='ngm', lags=1, n.exog=1, n.endog=1, N=1, cheby=FALSE, 
             upper = upper, lower=lower, quad=TRUE, n.quad=5,
             diff.tol=1e-05, n.iter=20, burn=1000, kappa=25, n.sim=10000,
             eps = .3, delta=.02, endog.init=k.ss, gain=1, reg.method=TRUE, 
             sr=FALSE )
# Solution options
a.b.range <- upper - lower
coeff.init <- matrix( c( k.ss, 
                         params.full.dep$alpha * a.b.range[2] / 2, 
                         k.ss * a.b.range[1] / 2 ), 3, 1)
# The linear-approx solution evaluated @ the steady state:
#   k' = alpha * betta * exp(x) * A * k ^ alpha
#     ~= k.ss + k.ss * x + alpha * ( k - k.ss )
#   NB: The range is already centered on k.ss, so can just use 
#       these coeffs

sol <- sol.iterate( coeff.init, opt, params.full.dep )
# The solution
n.irf <- 10
sim.irf <- 1000
irf.init <- irf_create( n.irf+1, sim.irf, opt$N, 0, params.full.dep$rho, 
                        params.full.dep$sig.eps, coeff.init, opt$upper, opt$lower, 
                        c(0,k.ss), opt$n.endog, opt$n.exog, .03, FALSE )
irf <- irf_create( n.irf+1, sim.irf, opt$N, 0, params.full.dep$rho, 
                   params.full.dep$sig.eps, sol, opt$upper, opt$lower, 
                   c(0,k.ss), opt$n.endog, opt$n.exog, .03, FALSE )
plot( 0:n.irf, irf.init[,2], type='l', lwd=2, col='blue' )
lines( 0:n.irf, irf[,2], type='l', lwd=2, col='red' )
lines( 0:n.irf, irf[,1], lwd=2 )


##### 2 Full depreciation case: Bad initial guess ####
coeff.far <- matrix( c( k.ss, 
                        params.full.dep$alpha * a.b.range[2] / 2, 
                        k.ss * a.b.range[1] / 2 ), 3, 1) * c( 1.6, 1.4, .5 )
# The linear-approx solution evaluated @ the steady state:
#   k' = alpha * betta * exp(x) * A * k ^ alpha
#     ~= k.ss + k.ss * x + alpha * ( k - k.ss )
#   NB: The range is already centered on k.ss, so can just use 
#       these coeffs

sol.far <- sol.iterate( coeff.far, opt, params.full.dep )
# The solution
irf.far <- irf_create( n.irf+1, sim.irf, opt$N, 0, params.full.dep$rho, 
                       params.full.dep$sig.eps, sol.far, opt$upper, opt$lower, 
                       c(0,k.ss), opt$n.endog, opt$n.exog, .03, FALSE )
plot( 0:n.irf, irf.init[,2], type='l', lwd=2, col='red', lty=2 )
lines( 0:n.irf, irf.far[,1], lwd=2 )
lines( 0:n.irf, irf.far[,2], type='l', lwd=2, col='red' )


##### 3. Full depreciation case: Nonlinear solution ####
opt$N <- 2
coeff.init.2 <- matrix( c( sol[1:2], 0, sol[3], 0, 0 ), 6,1 )
sol.2 <- sol.iterate( coeff.init.2, opt, params.full.dep )
# The solution
irf.2 <- irf_create( n.irf+1, sim.irf, opt$N, 0, params.full.dep$rho, 
                     params.full.dep$sig.eps, sol.2, opt$upper, opt$lower, 
                     c(0,k.ss), opt$n.endog, opt$n.exog, .03, FALSE )
plot( 0:n.irf, irf.init[,2], type='l', lwd=2, col='red', lty=2 )
lines( 0:n.irf, irf.2[,1], lwd=2 )
lines( 0:n.irf, irf.2[,2], type='l', lwd=2, col='red' )

opt$N <- 3
idx.2 <- idx_create( 2, 2 )
idx.3 <- idx_create( 3, 2 )
init.3 <- apply( idx.3, 1, function(rr) 
  which( apply( idx.2, 1, identical, rr )) )
coeff.init.3 <- matrix( sapply( init.3, function(i) 
  if( length(i)[[1]] == 0 ) 0 else sol.2[i] ), ncol=1 )
sol.3 <- sol.iterate( coeff.init.3, opt, params.full.dep )
# The solution

exact.sol <- function( k, x, params ){
  # The exact solution
  return( params$alpha * params$betta * exp(x) * params$A * k ^ params$alpha )
}
KK <- seq( lower[2], upper[2], length.out=201 )
pol.sol.exact <- sapply(KK, function(k) exact.sol(k,0,params.full.dep) )
plot( KK, pol.sol.exact, type='l', lwd=2 )
abline( v=k.ss, lwd=.5 )
abline( a=0, b=1, lwd=.5 )
add.line.0 <- function( this.coeff, N, this.opt, ... ){
  exog.0 <- matrix( rep( 0, 20 ), ncol=1 )
  endog.prime <- sapply( KK, 
                         function(k) endog_update( 0, k, this.coeff, this.opt$n.exog, 
                                                   this.opt$n.endog, N, this.opt$upper, this.opt$lower, this.opt$cheby ) )
  lines( KK, endog.prime, ... )
}
add.line.0( sol, 1, opt, col='red', lty=2, lwd=2 )
add.line.0( sol.2, 2, opt, col='red', lwd=2 )
add.line.0( sol.3, 3, opt, col='blue', lwd=1 )


##### 4. Full depreciation case: Nonlinear solution, Chebychev polynomial basis ####
opt.c <- opt
opt.c$cheby <- TRUE
sol.3.cheby <- sol.iterate( sol.3, opt.c, params.full.dep )
add.line.0( sol.3, 3, opt.c, col='blue', lwd=2, lty=2 )

### CHEBYCHEV POLYNOMIALS STILL NEED FIXING ###


##### 5. Partial depreciation case ####
params <- list( A=1, alpha=.3, delta=.025, gamma=1, betta=.99, 
                rho=0.95, sig.eps=.01 )
params$A <- ( 1 / params$betta - 1 + params$delta ) / params$alpha
k.ss <- ( ( 1 / params$betta - 1 + params$delta ) / 
            ( params$A * params$alpha )  ) ^ ( 1 / ( params$alpha - 1 ) )
sd.x<- sqrt( params$sig.eps / ( 1 - params$rho ^ 2 ) )
upper <- c(  3 * sd.x, k.ss * 1.4 )
lower <- c( -3 * sd.x, k.ss * 0.6 )
# Create the bounds
opt <- list( model='ngm', lags=1, n.exog=1, n.endog=1, N=1, cheby=FALSE, 
             upper = upper, lower=lower, quad=TRUE, n.quad=5,
             diff.tol=1e-05, n.iter=5, burn=1000, kappa=25, n.sim=10000,
             eps = .3, delta=.02, endog.init=k.ss, gain=.1, reg.method=TRUE, 
             sr=FALSE )
    # Solution options
a.b.range <- upper - lower
coeff.init <- matrix( c( k.ss, 
                         params.full.dep$alpha * a.b.range[2] / 2, 
                         k.ss * a.b.range[1] / 2 ), 3, 1)
    # Initialize with the full-depreciation solution
sol.slow <- sol.iterate( coeff.init, opt, params, FALSE )
opt$gain <- 1
opt$n.iter <- 30
sol.fast <- sol.iterate( sol.slow, opt, params )
    # Set the gain very low initially, then converge fast
n.irf <- 40
irf <- irf_create( n.irf+1, sim.irf, opt$N, 0, params$rho, 
                   params$sig.eps, sol.fast, opt$upper, opt$lower, 
                   c(0,k.ss), opt$n.endog, opt$n.exog, .01, FALSE )
    # Create IRF

opt$N <- 2
coeff.init.2 <- matrix( c( sol.fast[1:2], 0, sol.fast[3], 0, 0 ), 6,1 )
sol.2.fast <- sol.iterate( coeff.init.2, opt, params )
    # The second-order approximate solution
irf.2 <- irf_create( n.irf+1, sim.irf, opt$N, 0, params$rho, 
                     params$sig.eps, sol.2.fast, opt$upper, opt$lower, 
                     c(0,k.ss), opt$n.endog, opt$n.exog, .01, FALSE )
    # The IRF for the second-order approx

plot( 0:n.irf, irf[,1], lwd=2, type='l' )
lines( 0:n.irf, irf[,2], lwd=2, col='red' )
lines( 0:n.irf, irf.2[,2], lwd=2, col='red', lty=2 )

