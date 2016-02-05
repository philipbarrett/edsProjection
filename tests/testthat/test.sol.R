context('Checking the generic error evaluation code in sol.cpp & related')


test_that("Set reduction and error evaluation: NGM model", {
  
  params <- list( A=1, alpha=1/3, delta=.05, gamma=1, delta=.02, betta=.95 )
  upper <- c( 1, 1.5 )
  lower <- c( -.5, 0 )
  coeff <- matrix( c( 1, .1, .15 ), 3, 1)
  # Variables used
  rho <- .8
  sig.eps <- .1
  n.quad <- 7
  # Number of quadrature nodes
  n.sim.out <- 1e04
  # The number of points extracted from the simulation
  burn <- 1000
  kappa <- 100
  endog.init <- c(1)
  exog.sim <- ar1_sim( n.sim.out * kappa + burn, rho, sig.eps )
  # The simulated exogenous state
  endog.sim <- endog_sim( n.sim.out, exog.sim, coeff, 1, upper, lower, 
                          endog.init, FALSE, kappa, burn, TRUE )
  idces <- p_eps_cheap_const_idx( endog.sim[,c(1:2,4)], .3, .01 )
  # Need to include the lagged state in the evaluation set for the
  # equilibrium condition
  X <- endog.sim[ idces == 1, ]
  # Set reduction
  nodes <- quad_nodes_1d( n.quad, sig.eps, 0 )
  wts <- quad_weights_1d( n.quad )
  # Quadrature weights and nodes
  sum.err <- 0
  for( i in 1:nrow(X) ){
    exog.1 <- matrix( X[i,1], 1, 1 )
    endog.1 <- matrix( X[i,c(2,4)], 2, 1 )
    # The first elements of the states, as matrices 
    err.1 <- err_ngm( exog.1, endog.1, nodes, params, coeff, 1, 1, rho, 
                      n.quad, 1, upper, lower, FALSE, wts )
    # The first line of the error
    err.all.1 <- eval_err( coeff, matrix( X[i,], nrow=1) , 'ngm', 1, params, 1, 1, rho, sig.eps, 0, 1,
                           upper, lower, FALSE, matrix(0,1,1), TRUE, n.quad )
    # The same, just using the matrix of all terms but restricted to the
    # first line
    expect_equal( err.1, err.all.1 )
    sum.err <- sum.err + err.1 / nrow(X)
  }
  
  err.all <- eval_err( coeff, X, 'ngm', 1, params, 1, 1, rho, sig.eps, 0, 1,
                       upper, lower, FALSE, matrix(0,1,1), TRUE, n.quad )
  # Usage
  expect_equal( err.all, sum.err )
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
      # Variables used
  rho <- .8
  sig.eps <- .1
  n.quad <- 7
      # Number of quadrature nodes
  n.sim.out <- 1e04
      # The number of points extracted from the simulation
  burn <- 1000
  kappa <- 100
  endog.init <- c(1)
  exog.sim <- ar1_sim( n.sim.out * kappa + burn, rho, sig.eps )
      # The simulated exogenous state
  endog.sim <- endog_sim( n.sim.out, exog.sim, coeff, N, upper, lower, 
                          endog.init, FALSE, kappa, burn, TRUE )
  idces <- p_eps_cheap_const_idx( endog.sim[,c(1:2,4)], .3, .01 )
      # Need to include the lagged state in the evaluation set for the
      # equilibrium condition
  X <- endog.sim[ idces == 1, ]
      # The restricted simulation
  base <- eval_err( coeff, X, 'ngm', 1, params, 1, 1, rho, sig.eps, 0, N,
                    upper, lower, FALSE, matrix(0,1,1), TRUE, n.quad )
  deriv.x <- eval_err_D( coeff, X, 'ngm', 1, params, 1, 1, rho, sig.eps, 0, N,
                       upper, lower, FALSE, matrix(0,1,1), TRUE, n.quad )
      # The exact derivative
  deriv.fd <- 0 * deriv.x
      # Initialize the first-diff approx
  inc <- 1e-06 * diag( 6 )
      # The matrix of increments
  for( i in 1:6 )
    deriv.fd[i] <- 1e06 * 
                    ( eval_err( coeff + inc[i,], X, 'ngm', 1, params, 1, 1, 
                                  rho, sig.eps, 0, N, upper, lower, FALSE, 
                                  matrix(0,1,1), TRUE, n.quad ) - base )
      # Fill in the vector of first derivatives
  expect_equal( deriv.x, deriv.fd, tolerance=1e-04 )
      # Need to have low tolerance because FD is only an approximation
  
  
})

test_that("Minimizing coefficients", {
  
  ## THIS TEST SKIPPED - JUST USAGE ##
  
#   params <- list( A=1, alpha=1/3, delta=.05, gamma=1, delta=.02, betta=.95 )
#   upper <- c( 1, 1.5 )
#   lower <- c( -.5, 0 )
#   N <- 2
#       # approximation order
#   coeff <- matrix( c( 1, .1, .15, .1, .1, -.05 ), 6, 1)
#   exog <- matrix(0,1,1)
#   endog <- matrix(1,2,1)
#       # Variables used
#   rho <- .8
#   sig.eps <- .1
#   cheby <- TRUE
#   n.quad <- 7
#       # Number of quadrature nodes
#   n.sim.out <- 1e04
#       # The number of points extracted from the simulation
#   burn <- 1000
#   kappa <- 100
#   endog.init <- c(1)
#   exog.sim <- ar1_sim( n.sim.out * kappa + burn, rho, sig.eps )
#       # The simulated exogenous state
#   endog.sim <- endog_sim( n.sim.out, exog.sim, coeff, N, upper, lower, 
#                           endog.init, cheby, kappa, burn, TRUE )
#   idces <- p_eps_cheap_const_idx( endog.sim[,c(1:2,4)], .3, .01 )
#       # Need to include the lagged state in the evaluation set for the
#       # equilibrium condition
#   X <- endog.sim[ idces == 1, ]
#       # The restricted simulation
#   err.min(  coeff, X, 'ngm', 1, params, 1, 1, rho, sig.eps, 0, N,
#             upper, lower, cheby, matrix(0,1,1), TRUE, n.quad )
#         # Just a usage example
  
})