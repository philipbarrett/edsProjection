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
  nodes <- quad_nodes_1d( n.quad, sig.eps, 0 )
  wts <- quad_weights_1d( n.quad )
  base <- err_ngm( exog, endog, nodes, params, coeff, 1, 1, rho, 
                     n.quad, 2, upper, lower, FALSE, wts, TRUE )  
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
                   n.quad, 2, upper, lower, FALSE, wts, TRUE )  
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
  nodes <- quad_nodes_1d( n.quad, sig.eps, 0 )
  wts <- quad_weights_1d( n.quad )
  base <- err_ngm( exog, endog, nodes, params, coeff, 1, 1, rho, 
                   n.quad, 2, upper, lower, TRUE, wts, TRUE )  
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
                   n.quad, 2, upper, lower, TRUE, wts, TRUE )  
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
  expect_equal( deriv.x, deriv.fd, tolerance=1e-05 )
      # Need to have low tolerance because FD is only an approximation
  
} )