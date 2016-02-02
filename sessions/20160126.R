exog <- cbind( ar1_sim( 100, .8, .5 ), ar1_sim( 100, .5, .8 ) )
idx_create( 1, 4 )
range(exog)
upper <- c( rep( 4, 4 ))
lower <- c( rep( -4, 4 ))
coeff <- matrix( c( .5, 0, .5 * 4, 0, .1, 0, 0, 4, 1, 0 ), 5, 2 )
oo <- endog_sim( 100, exog, coeff, 1, upper, lower, c(0,0), FALSE )
qq <- endog_sim( 100, exog, coeff, 1, upper, lower, c(0,0), FALSE, lag=TRUE )
plot( c(1,100), range(oo), type='n')
for( i in 1:4 ) lines( 1:100, oo[,i], col=i )
plot( c(1,100), range(qq), type='n')
for( i in 1:4 ) lines( 1:100, qq[,i], col=i )
for( i in 1:4 ) lines( 1:100, qq[,4+i], col=i, lty=2 )

par(mfrow=c(2,2))
pp <- irf_create( 20, 1000, 1, 0, c(.8,.5), c(.5, .8), coeff, upper, lower, 
                  c( 0, 0, 1, 1 ), 2, 2, .01 )
qq <- irf_create( 20, 1000, 1, 1, c(.8,.5), c(.5, .8), coeff, upper, lower, 
                  c( 0, 0, 1, 1 ), 2, 2, .01 )
for( i in 1:4 ){
  plot( 1:20, pp[, i], type='l' )
  lines( 1:20, qq[, i], col=2 )
}
par(mfrow=c(1,1))

params <- list( A=1, alpha=1/3, delta=.05, gamma=1, delta=.02, betta=.95 )
integrand_ngm( matrix(0,1,1), matrix(1,2,1), c(.01), params, 
               matrix( c( 1, .1, .15 ), 3, 1), 1, 1, 1, c(1, 1.5), c( -.5, 0 ) )

exog.innov <- matrix( rnorm( 10000 ), 10000, 1 )
integrand_ngm( matrix(0,1,1), matrix(1,2,1), exog.innov[1], params, 
               matrix( c( 1, .1, .15 ), 3, 1), 1, 1, 1, c(1, 1.5), c( -.5, 0 ) )

1 - params$betta * integrand_ngm( matrix(0,1,1), matrix(1,2,1), .1 * exog.innov[1], params, 
               matrix( c( 1, .1, .15 ), 3, 1), 1, 1, 1, c(1, 1.5), c( -.5, 0 ) )
err_ngm( matrix(0,1,1), matrix(1,2,1), exog.innov, params,
         matrix( c( 1, .1, .15 ), 3, 1), 1, 1, .8, .1, 1, 1, 
         c(1, 1.5), c( -.5, 0 ), FALSE, rep(1,1) )
    # Check that the individual evaluation is ok

err.mc <- err_ngm( matrix(0,1,1), matrix(1,2,1), exog.innov, params,
                    matrix( c( 1, .1, .15 ), 3, 1), 1, 1, .8, .1, 10000, 1, 
                    c(1, 1.5), c( -.5, 0 ), FALSE, rep(1,10000)/10000 )
v.err.quad <- rep(NULL,10)
for( i in 1:10 ){
  nodes <- quad_nodes_1d( i )
  wts <- quad_weights_1d( i )
  v.err.quad[i] <- err_ngm( matrix(0,1,1), matrix(1,2,1), nodes, params,
                            matrix( c( 1, .1, .15 ), 3, 1), 1, 1, .8, .1, i, 1, 
                            c(1, 1.5), c( -.5, 0 ), FALSE, wts )
}
err.quad <- err_ngm( matrix(0,1,1), matrix(1,2,1), nodes, params,
                      matrix( c( 1, .1, .15 ), 3, 1), 1, 1, .8, .1, i, 1, 
                      c(1, 1.5), c( -.5, 0 ), FALSE, wts )

##  NB: Gauss-Hermite integration is for int e^(-x^2)f(x) on R  ##
##      So I need to rescale for the constant of integration    ##

sum( sqrt(2) * nodes * wts )
sum( ( sqrt(2) * nodes ) ^ 2 * wts )


## NEED TO PROPERLY DESIGN THE MC/QUADRATURE INTEGRATION ##
## TO DO: ADD CHECKS OF THE EVAL_ERR FUNCTION ##