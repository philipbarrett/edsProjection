context('Checking Simulation code')

test_that("AR(1), Linear", {
# Linear case with a single forcing variable

  
  skip("These tests conductied in test.poly.R")
  
  n <- 20
  rho <- .6
  sig.eps <- .01
  N <- 1
  k.mean <- 0
      # Some parameters
  
  exog.sim <- ar1_sim( n, rho, sig.eps )
      # The exogenous simulation
  upper <- c( 1.01*max(exog.sim), k.mean + 4 * sig.eps / sqrt( 1 - rho ^ 2 ) )
  lower <- c( 1.01*min(exog.sim), k.mean - 4 * sig.eps / sqrt( 1 - rho ^ 2 ) )
      # Bounds
  range <- upper - lower
      # The range
  coeffs <- matrix( c( k.mean, 0, 1 ) * c( 0, range[2] / 2, range[1] / 2 ) )
      # Coefficients
  first.pt <- endog_update( exog.sim[1], 0, coeffs, 1, 1, N, upper, lower, FALSE )
      # The first point in the simulaiton
  expect_equal( first.pt, exog.sim[1] )
      # Becuase the coefficient should just be one on the exogenous variable
  
  
  endog.sim <- endog_sim( n, exog.sim, coeffs, N, upper, lower, k.mean, FALSE, 1, 0, 1 )
  # The simulation
  
  
} )