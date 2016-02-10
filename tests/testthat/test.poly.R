context('Checking polynomials')

test_that("Polynomial evaluation", {
# Checks the multi-dimensional polynomial evaluation code
  
  # A very simple case
  set.seed(1234)
  n <- 20
  rho <- .6
  sig.eps <- .01
  N <- 1
  k.mean <- 0
  cheby <- FALSE
  K <- 2
      # Some parameters
  
  exog.sim <- ar1_sim( n, rho, sig.eps )
      # A series of forceing variables. use this to generate reasonable
      # upper/lower bounds
  upper <- c( 1.01 * max(abs(exog.sim)), k.mean + 4 * sig.eps / sqrt( 1 - rho ^ 2 ) )
  lower <- c( - 1.01 * max(abs(exog.sim)), k.mean - 4 * sig.eps / sqrt( 1 - rho ^ 2 ) )
      # Bounds
  range <- upper - lower
      # The range
  coeffs <- matrix( c( k.mean, 0, 1 ) * c( 1, range[2] / 2, range[1] / 2 ) )
      # So k_t = a_t
  X <- matrix( c( exog.sim, rep( k.mean, n ) ), ncol=2 )
      # Create a matrix of test shocks
  
  zero.poly <- poly_eval( coeffs[,1], 0*X, N, lower, upper, cheby )
      # Simple test case where all elements are zero
  expect_equal( zero.poly, matrix(0, n, 1) )
  
  
  X.limits <- X_limits( X, lower, upper, nrow(X), ncol(X) )
      # X rescaled to match limits
  X.limits.test <- ( X - matrix( ( upper + lower ) / 2, nrow(X), ncol(X), 
                  byrow=T ) ) / matrix( range / 2, nrow(X), ncol(X), byrow=T )
  expect_equal( X.limits, X.limits.test )
      # Check the rescaling first
  
  basis <- basis_cube( X.limits, N, K, cheby )
      # The polynomial basis
  for( i in 1:n ){
    test.mat <- matrix( c( 1, 1, X.limits[i,] ), 2, 2, byrow=T )
    expect_equal( test.mat, basis[, , i] )
  }   # Now guarantee that the basis function works
  
  y <- poly_eval( coeffs[,1], X, N, lower, upper, cheby )
      # Try polynomial evaluation
  expect_equal( exog.sim, y )
      # This should be true
  endog.sim <- endog_sim( n, exog.sim, coeffs, N, upper, lower, 0,
                            FALSE, 1, 0, TRUE )
      # The endogenous simulation
  expect_equal( c(exog.sim), endog.sim[,1] )
      # The endogenous simulation should generate the just the exogenous
      # variables
  
    ## LESSON: RANGES FOR RESCALING NEED TO BE SYMMETRIC (ABOUT WHAT!?)
  
  #### PART II: ALLOW FOR AUTOCORRELATION ####
  k.mean <- 1
      # Also allowing for level shift
  upper <- c( 1.01 * max(abs(exog.sim)), k.mean + 4 * sig.eps / sqrt( 1 - rho ^ 2 ) )
  lower <- c( - 1.01 * max(abs(exog.sim)), k.mean - 4 * sig.eps / sqrt( 1 - rho ^ 2 ) )
      # Bounds
  range <- upper - lower
      # The range
  coeffs <- matrix( c( k.mean, .5, 1 ) * c( 1, range[2] / 2, range[1] / 2 ) )
      # Change the oefficients to allow for persistence
  n <- 20000
      # Increae the sample
  exog.sim <- ar1_sim( n, rho, sig.eps )
  endog.sim <- endog_sim( n, exog.sim, coeffs, N, upper, lower, 0,
                          FALSE, 1, 0, TRUE )
  coeff.reg <- coeff_reg( endog.sim[,2], endog.sim[,c(1,4)], N, lower, upper, cheby )
  expect_equal( coeffs, coeff.reg )
      # Recovers the coefficients perfectly
  
} )

test_that("Regression", {
# Test the regression coefficient construction
  
  set.seed(1234)
  n.elems <- 1000000
  K <- 3
  upper <- rep( 1, K )
  lower <- - upper
  X <- matrix( runif( K * n.elems, lower[1], upper[1] ), ncol=K )
      # Some high(ish)-dimensional data.  Should not get rescaled.
  err <- rnorm( n.elems )
      # The error
  N <- 1
  coeff <- 1:K
      # The coefficients
  y <- X %*% coeff + err
  coeff.reg <- coeff_reg( y, X, N, lower, upper, FALSE )
  expect_true( max( abs( coeff.reg - c( 0, K:1 ) ) ) < 4e-3 )
  
  
  N <- 2
  n.terms <- idx_count( N, K )
  coeff <- matrix( 0, ncol=1, nrow=n.terms )
  coeff[ 2, 1 ] <- coeff[ 4, 1 ] <- coeff[ 7, 1 ] <- 1
      # Linear terms only
  y <- poly_eval( coeff, X, N, lower, upper, FALSE ) + err
  coeff.reg <- coeff_reg( y, X, N, lower, upper, FALSE )
  expect_true( max( abs( coeff.reg - coeff ) ) < 4e-3 )
  
  set.seed(1234)
  coeff <- matrix( runif( n.terms, -1, 1 ), ncol=1, nrow=n.terms )
      # Try with some random polynomial coefficients
  y <- poly_eval( coeff, X, N, lower, upper, FALSE ) + err
  coeff.reg <- coeff_reg( y, X, N, lower, upper, FALSE )
  expect_true( max( abs( coeff.reg - coeff ) ) < 4e-3 )
  
  upper <- rep( 10, K )
  lower <- rep( 2, K )
  n.elems <- 1000000
  X <- matrix( runif( K * n.elems, lower[1], upper[1] ), ncol=K )
      # Some new data.  Should get rescaled.
  err <- rnorm( n.elems )
  coeff <- matrix( runif( n.terms, -1, 1 ), ncol=1, nrow=n.terms )
      # Reset the coefficients
  y <- poly_eval( coeff, X, N, lower, upper, FALSE ) + err
  coeff.reg <- coeff_reg( y, X, N, lower, upper, FALSE )
  expect_true( max( abs( coeff.reg - coeff ) ) < 1e-2 )
      # Generous tolerance because regression coefficients just have some error
  
} )