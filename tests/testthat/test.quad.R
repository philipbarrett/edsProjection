
context('Checking multi-dimensional quadrature rules')

test_that("One-dimensional quadrature", {
  
  ## 1. For N(0,1) ##
  nds <- quad_nodes_1d( 10, 1, 0 )
  wts <- quad_weights_1d( 10 )
  mu <- sum( nds * wts )
  sig <- sqrt( sum( nds ^ 2 * wts ) - mu ^ 2 )
  expect_equal(mu, 0 )
  expect_equal(sig, 1 )

  ## 2. For N(1,1) ##
  nds <- quad_nodes_1d( 10, 1, 1 )
  wts <- quad_weights_1d( 10 )
  mu <- sum( nds * wts )
  sig <- sqrt( sum( nds ^ 2 * wts ) - mu ^ 2 )
  expect_equal(mu, 1 )
  expect_equal(sig, 1 )
  
  ## 3. For N(0,3) ##
  nds <- quad_nodes_1d( 10, 3, 0 )
  wts <- quad_weights_1d( 10 )
  mu <- sum( nds * wts )
  sig <- sqrt( sum( nds ^ 2 * wts ) - mu ^ 2 )
  expect_equal(mu, 0 )
  expect_equal(sig, 3)
  
  ## 3. For N(2,3) ##
  nds <- quad_nodes_1d( 10, 3, 2 )
  wts <- quad_weights_1d( 10 )
  mu <- sum( nds * wts )
  sig <- sqrt( sum( nds ^ 2 * wts ) - mu ^ 2 )
  expect_equal(mu, 2 )
  expect_equal(sig, 3 )
  
})

test_that("Multi-dimensional quadrature", {
  
  ## 1. Spherical errors, zero mean ##
  quad <- quad_nodes_weights( 5, 2, c(1,1), c(0,0) )
  mu <- t( quad$nodes ) %*% quad$weights
  sig <- sqrt( t( quad$nodes ^ 2 ) %*% quad$weights - mu ^ 2 )
  cov <- t( apply( quad$nodes, 1, prod ) ) %*% quad$weights - prod( mu )
  expect_true( max( abs( mu - c( 0, 0 ) ) ) < 1e-16 )
  expect_true( max( abs( sig - c( 1, 1 ) ) ) < 1e-09 )
  expect_true( max( abs( cov - 0 ) ) < 1e-16 )

  ## 2. Spherical errors, non-zero mean ##
  quad <- quad_nodes_weights( 5, 2, c(1,1), c(2,2) )
  mu <- t( quad$nodes ) %*% quad$weights
  sig <- sqrt( t( quad$nodes ^ 2 ) %*% quad$weights - mu ^ 2 )
  cov <- t( apply( quad$nodes, 1, prod ) ) %*% quad$weights - prod( mu )
  expect_true( max( abs( mu - c( 2, 2 ) ) ) < 1e-08 )
  expect_true( max( abs( sig - c( 1, 1 ) ) ) < 1e-08 )
  expect_true( max( abs( cov - 0 ) ) < 1e-08 )
  
  ## 3. Independent, non-spherical errors, zero mean ##
  quad <- quad_nodes_weights( 5, 2, c(2,.5), c(0,0) )
  mu <- t( quad$nodes ) %*% quad$weights
  sig <- sqrt( t( quad$nodes ^ 2 ) %*% quad$weights - mu ^ 2 )
  cov <- t( apply( quad$nodes, 1, prod ) ) %*% quad$weights - prod( mu )
  expect_true( max( abs( mu - c( 0, 0 ) ) ) < 1e-08 )
  expect_true( max( abs( sig - c( 2, .5 ) ) ) < 1e-08 )
  expect_true( max( abs( cov - 0 ) ) < 1e-08 )
  
  ## 4. Independent, non-spherical errors, nonzero mean ##
  quad <- quad_nodes_weights( 5, 2, c(2,.5), c(1,-1) )
  mu <- t( quad$nodes ) %*% quad$weights
  sig <- sqrt( t( quad$nodes ^ 2 ) %*% quad$weights - mu ^ 2 )
  cov <- t( apply( quad$nodes, 1, prod ) ) %*% quad$weights - prod( mu )
  expect_true( max( abs( mu - c( 1, -1 ) ) ) < 1e-08 )
  expect_true( max( abs( sig - c( 2, .5 ) ) ) < 1e-08 )
  expect_true( max( abs( cov - 0 ) ) < 1e-08 )
  
})