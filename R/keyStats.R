#########################################################################
# keyStats.R
#
# Code to extract the key statistics from a list of solutions
# Philip Barrett, Chicago
# 25mar2016
#########################################################################

key.stats.sol <- function( rep.data ){
# Extracts the key statistics of a given model
  
  err <-  max(apply( abs( rep.data$err ), 2, mean ) ) * 100
      # The error
  
  cons <- rep.data$cont.sim[,1]
  r.cons <- rep.data$r.cont[,1]
  cons.gth <- r.cons - cons
      # Consumption and consumption growth
  q <- rep.data$cont.sim[,13] - rep.data$cont.sim[,9] + rep.data$cont.sim[,10]
  r.q <- rep.data$r.cont[,13] - rep.data$r.cont[,9] + rep.data$r.cont[,10]
  q.gth <- r.q - q
      # Real exchange rate and RER growth
  cov.q.c <- cov( q, cons )
  cor.q.c <- cor( q, cons )
  cov.q.c.gth <- cov( q.gth, cons.gth )
  cor.q.c.gth <- cor( q.gth, cons.gth )
      # Statistics for q & c
  
  tot.assets <- rep.data$endog.sim[,3] * exp( - rep.data$cont.sim[,3] ) -
                        rep.data$endog.sim[,4] * exp( rep.data$cont.sim[,13] - rep.data$cont.sim[,4] )
  dom.assets <- rep.data$endog.sim[,3] * exp( - rep.data$cont.sim[,3] )
  for.assets <- rep.data$endog.sim[,4] * exp( - rep.data$cont.sim[,4] )
  fit <- lm( dom.assets ~ tot.assets )
      # The home bias regression
  home.bias.coeff <- fit$coefficients[1]
  names(home.bias.coeff) <- NA
  dom.mean <- mean( dom.assets )
  for.mean <- mean( for.assets )
  abs.tot.assets <- mean( abs( tot.assets ) )
  asset.diff <- mean( rep.data$endog.sim[,3] * exp( - rep.data$cont.sim[,3] ) +
    rep.data$endog.sim[,4] * exp( rep.data$cont.sim[,13] - rep.data$cont.sim[,4] ) )
      # The statistics for home bias.  Last one is difference in assets (sum as use neg've for assets)
  
  e.e <- rep.data$e.cont[,13] - rep.data$cont.sim[,13]
  r.diff <- rep.data$cont.sim[,3] - rep.data$cont.sim[,4]
  uip <- lm( e.e~r.diff )
      # The uip regression
  uip.coeff <- uip$coefficients['r.diff']
      # The key statistic
  
  return( c( cov.q.c=cov.q.c, cor.q.c=cor.q.c, cov.q.c.gth=cov.q.c.gth, cor.q.c.gth=cor.q.c.gth, 
             home.bias.coeff=home.bias.coeff, dom.mean=dom.mean, for.mean=for.mean, 
             uip.coeff=uip.coeff, err.pp=err, asset.diff=asset.diff, abs.tot.assets=abs.tot.assets ) )
}

key.stats <- function( l.sol, sol.base, base.rep, st.param ){
  
  data <- sapply( l.sol, function(x) c( x$sol.2$params[[st.param]][1], key.stats.sol( x$rep.2 ) ) )
  data <- cbind( c( st.param=sol.base$params[[st.param]][1], key.stats.sol( base.rep ) ), data )
  rownames(data)[1] <- st.param
  data <- data[, order(data[1,])]
  return( data )
}




