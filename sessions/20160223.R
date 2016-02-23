params <- list( alpha = .85, gamma = 5, P1.bar=1, P2.bar=1, betta=.99,
                rho=c(.95,.95), sig.eps=c(.01,.01), eta=1 )
    # Parameters

c.ss <- params$alpha ^ params$alpha * ( 1 - params$alpha ) ^ ( 1 - params$alpha )
p.ss <- 1 / c.ss

exog <- matrix( 0, 2, 2 )
endog <- exog
cont <- log( c( c.ss, c.ss, 1 / params$betta, 1 / params$betta, 
                params$alpha, params$alpha, 1 - params$alpha, 1 - params$alpha, 
                p.ss, p.ss, 1, 1, 1 ) )
# exog <- matrix( X[1,c(1,2,5,6)], 2, 2, byrow=T )
# endog <- matrix( X[1,c(3,4,7,8)], 2, 2, byrow=T )
# cont <- cont.sim[1,]

contemp_eqns_irbc( exog, endog, cont, params )

X.22 <- X.11 <- log( seq( .01, 1, length.out=101 ) )
x.11.out <- function(x.11)
  return( contemp_eqns_irbc( exog, endog, c( cont[1:4], x.11, cont[5:13]), params )[5] )
x.22.out <- function(x.22)
  return( contemp_eqns_irbc( exog, endog, c( cont[1:5], x.22, cont[6:13]), params )[6] )

x.11.out(log(params$alpha)) - log( params$alpha )
x.22.out(log(params$alpha)) - log( params$alpha )


plot( X.11, sapply(X.11, x.11.out), type='l' )
abline( 0, 1 )
abline(v=log(params$alpha), lty=2)
abline(h=log(params$alpha), lty=2)

# for( i in 1:20 ){
#   exog[1,1] <- exog[1,1] - .02
#   lines( X.11, sapply(X.11, x.11.out), col=i+1 )
# }

plot( X.22, sapply(X.22, x.22.out), type='l' )
abline( 0, 1 )
abline(v=log(params$alpha), lty=2)
abline(h=log(params$alpha), lty=2)


