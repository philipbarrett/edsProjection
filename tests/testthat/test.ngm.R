context('Checking integration for NGM model')

test_that("Quadrature vs. MC integration", {
  
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

test_that("Set reduction and error evaluation", {
  
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
  idces <- p_eps_cheap_const_idx( endog.sim[,1:2], .3, .01 )
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