exog <- cbind( ar1_sim( 100, .8, .5 ), ar1_sim( 100, .5, .8 ) )
idx_create( 1, 4 )
range(exog)
upper <- c( rep( 4, 4 ))
lower <- c( rep( -4, 4 ))
coeff <- matrix( c( .5, 0, .5 * 4, 0, .1, 0, 0, 4, 1, 0 ), 5, 2 )
oo <- endog_sim( 100, exog, coeff, 1, upper, lower, c(0,0), TRUE )
plot( c(1,100), range(oo), type='n')
for( i in 1:4 ) lines( 1:100, oo[,i], col=i )

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
err_ngm_mc( matrix(0,1,1), matrix(1,2,1), exog.innov, params,
                matrix( c( 1, .1, .15 ), 3, 1), 1, 1, .8, .1, 10000, 1, 
                c(1, 1.5), c( -.5, 0 ) )

