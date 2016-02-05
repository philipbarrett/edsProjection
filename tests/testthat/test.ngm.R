context('Checking NGM model')

test_that("Integration", {
  
  params <- list( A=1, alpha=1/3, delta=.05, gamma=1, delta=.02, betta=.95 )
  upper <- c( 1, 1.5 )
  lower <- c( -.5, 0 )
  coeff <- matrix( c( 1, .1, .15 ), 3, 1)
  exog <- matrix(0,1,1)
  endog <- matrix(1,2,1)
      # Variables used
  rho <- .8
  sig.eps <- .1
  n.mc <- 100000 # 1000000
      # Number of Monte Carlo shocks.  Works better & slower with bigger numbers.
  n.quad <- 7
      # Number of quadrature nodes
  
  set.seed=(12345)
  exog.innov <- matrix( rnorm( n.mc, 0, sig.eps ), n.mc, 1 )
      # Monte Carlo integration shocks
  err.mc <- err_ngm( exog, endog, exog.innov, params, coeff, 1, 1, rho, 
                     n.mc, 1, upper, lower, FALSE, rep(1,n.mc)/n.mc )
  nodes <- quad_nodes_1d( n.quad, sig.eps, 0 )
  wts <- quad_weights_1d( n.quad )
  err.quad <- err_ngm( exog, endog, nodes, params,
                       coeff, 1, 1, rho, n.quad, 1,
                       upper, lower, FALSE, wts )
  
  expect_true( abs( err.mc - err.quad ) < 2e-04 )
      # Check that the integrals are close(ish).  MCintegration is quite bad!
})

test_that("Derivatives", {
  
  params <- list( A=1, alpha=1/3, delta=.05, gamma=1, delta=.02, betta=.95 )
  upper <- c( 1, 1.5 )
  lower <- c( -.5, 0 )
  N <- 2
      # approximation order
  coeff <- matrix( c( 1, .1, .15, .1, .1, -.05 ), 6, 1)
  exog <- matrix(0,1,1)
  endog <- matrix(1,2,1)
  exog_shk <- .01
      # Variables used
  rho <- .8
  sig.eps <- .1
  n.quad <- 7
      # Number of quadrature nodes

  ##### The integrand #####
  
  ## Ordinary polynomials
  exog <- matrix(0,1,1)
  base <- integrand_ngm( exog, endog, exog_shk, params, coeff, 1, 1, 2, upper, 
                         lower, FALSE)
  deriv.x <- integrand_ngm_D( exog, endog, exog_shk, params, coeff, 1, 1, 2, 
                              upper, lower, FALSE)
      # The exact derivative
  deriv.fd <- 0 * deriv.x
      # Initialize the first-diff approx
  inc <- 1e-06 * diag( 6 )
      # The matrix of increments
  for( i in 1:6 )
    deriv.fd[i] <- 1e06 * ( integrand_ngm( exog, endog, exog_shk, params, 
                                           coeff + inc[i,], 
                                           1, 1, 2, upper, lower, FALSE) - base )
        # Fill in the vector of first derivatives
  expect_equal( deriv.x, deriv.fd, tolerance=1e-04 )
      # Need to have low tolerance because FD is only an approximation
  
  ## Chebychev polynomials
  base <- integrand_ngm( exog, endog, exog_shk, params, coeff, 1, 1, 2, upper, 
                         lower, TRUE )
  deriv.x <- integrand_ngm_D( exog, endog, exog_shk, params, coeff, 1, 1, 2, 
                              upper, lower, TRUE )
      # The exact derivative
  deriv.fd <- 0 * deriv.x
      # Initialize the first-diff approx
  inc <- 1e-06 * diag( 6 )
      # The matrix of increments
  for( i in 1:6 )
    deriv.fd[i] <- 1e06 * ( integrand_ngm( exog, endog, exog_shk, params, 
                                           coeff + inc[i,], 
                                           1, 1, 2, upper, lower, TRUE) - base )
      # Fill in the vector of first derivatives
  expect_equal( deriv.x, deriv.fd, tolerance=1e-05 )
      # Need to have low tolerance because FD is only an approximation

  
  ##### The integral #####
  
  ## Ordinary polynomials
  exog <- matrix(0,1,1)
  nodes <- quad_nodes_1d( n.quad, sig.eps, 0 )
  wts <- quad_weights_1d( n.quad )
  base <- err_ngm( exog, endog, nodes, params, coeff, 1, 1, rho, 
                     n.quad, 2, upper, lower, FALSE, wts )  
      # Case where the sign of the error derivative is important
  deriv.x <- err_ngm_D( exog, endog, nodes, params, coeff, 1, 1, rho, 
                        n.quad, 2, upper, lower, FALSE, wts )  
      # The exact derivative
  deriv.fd <- 0 * deriv.x
      # Initialize the first-diff approx
  inc <- 1e-06 * diag( 6 )
      # The matrix of increments
  for( i in 1:6 )
    deriv.fd[i] <- 1e06 * ( err_ngm( exog, endog, nodes, params, 
                                     coeff + inc[i,], 1, 1, rho, 
                                     n.quad, 2, upper, lower, FALSE, wts ) - base )
      # Fill in the vector of first derivatives
  expect_equal( deriv.x, deriv.fd, tolerance=1e-04 )
      # Need to have low tolerance because FD is only an approximation
  
  exog <- matrix(0.2,1,1)
  base <- err_ngm( exog, endog, nodes, params, coeff, 1, 1, rho, 
                   n.quad, 2, upper, lower, FALSE, wts )  
      # Case where the sign of the error derivative is not important
  deriv.x <- err_ngm_D( exog, endog, nodes, params, coeff, 1, 1, rho, 
                        n.quad, 2, upper, lower, FALSE, wts )  
      # The exact derivative
  deriv.fd <- 0 * deriv.x
      # Initialize the first-diff approx
  inc <- 1e-06 * diag( 6 )
      # The matrix of increments
  for( i in 1:6 )
    deriv.fd[i] <- 1e06 * ( err_ngm( exog, endog, nodes, params, 
                                     coeff + inc[i,], 1, 1, rho, 
                                     n.quad, 2, upper, lower, FALSE, wts ) - base )
      # Fill in the vector of first derivatives
  expect_equal( deriv.x, deriv.fd, tolerance=1e-04 )
      # Need to have low tolerance because FD is only an approximation
  
  ## Chebychev polynomials
  exog <- matrix(0,1,1)
  nodes <- quad_nodes_1d( n.quad, sig.eps, 0 )
  wts <- quad_weights_1d( n.quad )
  base <- err_ngm( exog, endog, nodes, params, coeff, 1, 1, rho, 
                   n.quad, 2, upper, lower, TRUE, wts )  
      # Case where the sign of the error derivative is important
  deriv.x <- err_ngm_D( exog, endog, nodes, params, coeff, 1, 1, rho, 
                        n.quad, 2, upper, lower, TRUE, wts )  
      # The exact derivative
  deriv.fd <- 0 * deriv.x
      # Initialize the first-diff approx
  inc <- 1e-06 * diag( 6 )
      # The matrix of increments
  for( i in 1:6 )
    deriv.fd[i] <- 1e06 * ( err_ngm( exog, endog, nodes, params, 
                                     coeff + inc[i,], 1, 1, rho, 
                                     n.quad, 2, upper, lower, TRUE, wts ) - base )
      # Fill in the vector of first derivatives
  expect_equal( deriv.x, deriv.fd, tolerance=1e-05 )
      # Need to have low tolerance because FD is only an approximation
  
  exog <- matrix(0.6,1,1)
  base <- err_ngm( exog, endog, nodes, params, coeff, 1, 1, rho, 
                   n.quad, 2, upper, lower, TRUE, wts )  
      # Case where the sign of the error derivative is not important
  deriv.x <- err_ngm_D( exog, endog, nodes, params, coeff, 1, 1, rho, 
                        n.quad, 2, upper, lower, TRUE, wts )  
      # The exact derivative
  deriv.fd <- 0 * deriv.x
      # Initialize the first-diff approx
  inc <- 1e-06 * diag( 6 )
      # The matrix of increments
  for( i in 1:6 )
    deriv.fd[i] <- 1e06 * ( err_ngm( exog, endog, nodes, params, 
                                     coeff + inc[i,], 1, 1, rho, 
                                     n.quad, 2, upper, lower, TRUE, wts ) - base )
      # Fill in the vector of first derivatives
  expect_equal( deriv.x, deriv.fd, tolerance=1e-04 )
      # Need to have low tolerance because FD is only an approximation
  
} )

test_that("Checking full-depreciation solution", {
  
  params.full.dep <- list( A=1, alpha=.3, delta=1, gamma=1, betta=.99, 
                           rho=.0, sig.eps=.01 )
  k.ss <- ( ( 1 / params.full.dep$betta - 1 + params.full.dep$delta ) / 
              ( params.full.dep$A * params.full.dep$alpha )  ) ^ ( 1 / ( params.full.dep$alpha - 1 ) )
      # The steady state
  sd.x<- sqrt( params.full.dep$sig.eps / ( 1 - params.full.dep$rho ^ 2 ) )
  upper <- c( 3 * sd.x, k.ss * 1.4 )
  lower <- c( -3 * sd.x, k.ss * 0.6 )
      # Create the bounds
  opt <- list( model='ngm', lags=1, n.exog=1, n.endog=1, N=1, cheby=FALSE, 
               upper = upper, lower=lower, quad=TRUE, n.quad=7,
               err.tol=1e-05, n.iter=10, burn=0, kappa=1, n.sim=10000,
               eps = .3, delta=.02, endog.init=k.ss )
      # Solution options
  a.b.range <- upper - lower
  coeff.init <- matrix( c( k.ss, 
                           params.full.dep$alpha * a.b.range[2], 
                           k.ss * a.b.range[1] ), 3, 1)
      # The linear-approx solution evaluated @ the steady state:
      #   k' = alpha * betta * exp(x) * A * k ^ alpha
      #     ~= k.ss + k.ss * x + alpha * ( k - k.ss )
      #   NB: The range is already centered on k.ss, so can just use 
      #       these coeffs
  
  k.ss.prime <- endog_update( exog = c(0), k.ss, coeff.init, 1, 1, 1, upper, lower, FALSE )
  expect_true( k.ss == k.ss.prime )
  
  k.ss.euler <- 
    integrand_ngm( matrix( 0, 1, 1), matrix( k.ss, nrow=2, ncol=1), 0, 
                   params.full.dep, coeff.init, 1, 1, 1, upper, lower, FALSE )
  expect_equal( k.ss.euler * params.full.dep$betta, 1 )
  
  nodes <- quad_nodes_1d( opt$n.quad, params.full.dep$sig.eps, 0 )
  wts <- quad_weights_1d( opt$n.quad )
  
  single.integral <- 
    err_ngm( matrix( 0, 1, 1), matrix( k.ss, nrow=2, ncol=1), nodes, 
             params.full.dep, coeff.init, opt$n.exog, opt$n.endog, 
             params.full.dep$rho, opt$n.quad, opt$N, upper, lower, FALSE, wts )
  
  mult.pts <- 
    eval_err( coeff.init, matrix( c( 0, k.ss, 0, k.ss ), nrow=1 ), 'ngm', opt$lags, 
              params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
              params.full.dep$sig.eps,  0, opt$N, upper, lower, FALSE, 
              matrix(0,1,1,), TRUE, opt$n.quad )
  
  expect_equal( single.integral, mult.pts )
  
  exog.sim <- ar1_sim( opt$n.sim * opt$kappa + opt$burn, 
                       params.full.dep$rho, params.full.dep$sig.eps )
  endog.sim <- endog_sim( opt$n.sim, exog.sim, coeff.init, opt$N, opt$upper,
                          opt$lower, k.ss, FALSE, opt$kappa, opt$burn, opt$lags )
  indices <- p_eps_cheap_const_idx( endog.sim[,c(1,2,4)], opt$eps, opt$delta )
  sum(indices)
  X <- endog.sim[indices==1,]
      # State reduction
  grid.init <- eval_err( coeff.init, X, 'ngm', opt$lags, 
                         params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
                         params.full.dep$sig.eps,  0, opt$N, upper, lower, FALSE, 
                         matrix(0,1,1,), TRUE, opt$n.quad )
  grid.init.grad <- eval_err_D( coeff.init, X, 'ngm', opt$lags, 
                         params.full.dep, opt$n.exog, opt$n.endog, params.full.dep$rho, 
                         params.full.dep$sig.eps,  0, opt$N, upper, lower, FALSE, 
                         matrix(0,1,1,), TRUE, opt$n.quad )

  
#   sol <- sol.iterate( coeff.init, opt, params.full.dep, TRUE )
  
})