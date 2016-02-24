alpha <- .85
c.ss <- alpha ^ alpha * ( 1 - alpha ) ^ ( 1 - alpha )
c.1 <- log(c.ss)
c.2 <- log(c.ss)
A.1 <- 0
A.2 <- 0

x12.fn <- function( x.11 )
  return( ( c.1 - alpha * x.11 ) / ( 1 - alpha ) )

x21.fn <- function( x.22 )
  return( ( c.2 - alpha * x.22 ) / ( 1 - alpha ) )

x11.fn <- function( x.21 )
  return( log( exp(A.1) - exp(x.21) ) )

x22.fn <- function( x.12 )
  return( log( exp(A.2) - exp(x.12) ) )

x11.all <- function( x.11 )
  return( x11.fn( x21.fn( x22.fn( x12.fn( x.11 ) ) ) ) )

x11.all( log( alpha ) )
log(alpha)

X11 <- seq( .8, 1.2, length.out = 101 ) * log(alpha)
plot( X11, sapply( X11, x11.all ), type='l' )
abline( 0, 1 )
abline( v=log(alpha), lty=2 )

x11.conv <- log( alpha ) * .9
for( i in 2:100 )
  x11.conv <- c( x11.conv, x11.all( x11.conv[i-1] ) )

A.1 <- -.02
X11 <- seq( .5, 1.5, length.out = 101 ) * log(alpha)
plot( X11, sapply( X11, x11.all ), type='l' )
c.1 <- log( c.ss ) * .98
lines( X11, sapply( X11, x11.all ), col=2 )
abline( 0, 1 )
abline( v=log(alpha), lty=2 )
