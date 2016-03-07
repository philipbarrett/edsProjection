rm(list=ls())
setwd("~/code/2016/edsProjection")
library(edsProjection)
library(xtable)
load('02mar2016.rdata')
opt <- sol.irbc.2$opt
params <- sol.irbc.2$params

opt$N <- 1
coeff.init <- sol.irbc.2$coeff[ c(1,2,4,7,11), ]
coeff.init.cont <- sol.irbc.2$coeff.cont[ c(1,2,4,7,11), ]
coeff.init.cont[2:5, 13] <- c( rep( mean( abs( coeff.init.cont[ 2:3, 13] ) ), 2 ),
                               rep( mean( abs( coeff.init.cont[ 4:5, 13] ) ), 2 ) ) * 
                            sign( coeff.init.cont[2:5, 13] )
opt$l.sym.ave <- list( sym=list(c(2,3),c(4,5), ave=1 ) )
opt$sym.reg <- TRUE
opt$l.sym.ave <- list( sym=list( c(2,3), c(4,5) ), ave=c(1) )


n.draws <- 250
a.draw <- .01 * rnorm( n.draws )
m.a <- cbind( c( a.draw, rep(0, n.draws) ), c( rep(0, n.draws), a.draw ) )
    # The simulation of points for a
m.b.upper <- matrix( rnorm( 2 * n.draws ), n.draws, 2 )
    # Simulated points for B

endog.sim.l <- cbind( m.a, rbind(m.b.upper, m.b.upper[,2:1]) )
endog.sim <- cbind(endog.sim.l, endog.sim.l )
    # The simulation (assuming B & a are the same in consecutive periods - 
    # inconsequential).  BUT: The simulation *is* symmetric.  Every observation
    # corresponds to one with the a's and b's reversed


#### WITH SYMMETRY IN B ####
cont.sim <- cont_sim( endog.sim, coeff.init.cont, opt$N, opt$n.endog, 
                      opt$n.exog, opt$upper, opt$lower, opt$cheby )
    # Let's first see if cont.sim provides evaluations from which we can recover coefficients

print( c.coeff.lm <- 
        cbind( lm( cont.sim[,1] ~ endog.sim[,1:4] )$coeff,
               lm( cont.sim[,2] ~ endog.sim[,1:4] )$coeff ) )
print( c.coeff <- 
         cbind( coeff_reg( cont.sim[,1], endog.sim[,1:4], opt$N, opt$lower, opt$upper, opt$cheby ),
                coeff_reg( cont.sim[,2], endog.sim[,1:4], opt$N, opt$lower, opt$upper, opt$cheby ) ) )
print( c.coeff[,1] - c.coeff[c(1,3,2,5,4),2] )
print( c.coeff[,1:2] - coeff.init.cont[,1:2] )
    # So consumption is fine

print( p.coeff.lm <- 
         cbind(lm( cont.sim[,9] ~ endog.sim[,1:4] )$coeff,
               lm( cont.sim[,10] ~ endog.sim[,1:4] )$coeff ) )
    # Inconclusive.  Could be symmetry, could be inverse.

cont.hat <- contemp_eqns_irbc_grid( cbind( endog.sim, cont.sim ), 
                                    opt$lags, params, opt$n.exog, opt$n.endog, opt$n.cont )
    # Now create the equation errors
print( c.coeff.hat.lm <- 
         cbind( lm( cont.hat[,1] ~ endog.sim[,1:4] )$coeff,
                lm( cont.hat[,2] ~ endog.sim[,1:4] )$coeff ) )
print( p.coeff.hat.lm <- 
         cbind(lm( cont.hat[,9] ~ endog.sim[,1:4] )$coeff,
               lm( cont.hat[,10] ~ endog.sim[,1:4] )$coeff ) )
    # This looks correct, actually

#### WITH ASYMMETRY IN B ####
endog.sim.l.2 <- cbind( m.a, matrix( rnorm( 4 * n.draws ), ncol=2 ) )
endog.sim.2 <- cbind(endog.sim.l.2, endog.sim.l.2 )
    # Now try with variation in B across the draws
cont.sim.2 <- cont_sim( endog.sim.2, coeff.init.cont, opt$N, opt$n.endog, 
                        opt$n.exog, opt$upper, opt$lower, opt$cheby )

print( c.coeff.lm.2 <- 
         cbind( lm( cont.sim.2[,1] ~ endog.sim.2[,1:4] )$coeff,
                lm( cont.sim.2[,2] ~ endog.sim.2[,1:4] )$coeff ) )
print( p.coeff.lm.2 <- 
         cbind(lm( cont.sim.2[,9] ~ endog.sim.2[,1:4] )$coeff,
               lm( cont.sim.2[,10] ~ endog.sim.2[,1:4] )$coeff ) )
    # Looks good

print( c.coeff.2 <- 
         cbind( coeff_reg( cont.sim.2[,1], endog.sim.2[,1:4], opt$N, opt$lower, opt$upper, opt$cheby ),
                coeff_reg( cont.sim.2[,2], endog.sim.2[,1:4], opt$N, opt$lower, opt$upper, opt$cheby ) ) )
print( p.coeff.2 <- 
     cbind( coeff_reg( cont.sim.2[,9], endog.sim.2[,1:4], opt$N, opt$lower, opt$upper, opt$cheby ),
            coeff_reg( cont.sim.2[,10], endog.sim.2[,1:4], opt$N, opt$lower, opt$upper, opt$cheby ) ) )
    # Same for the coeff_reg versions
coeff.init.cont[,9:10]

cont.hat.2 <- contemp_eqns_irbc_grid( cbind( endog.sim.2, cont.sim.2 ), 
                                      opt$lags, params, opt$n.exog, opt$n.endog, opt$n.cont )


print( c.coeff.hat.lm.2 <- 
         cbind( lm( cont.hat.2[,1] ~ endog.sim.2[,1:4] )$coeff,
                lm( cont.hat.2[,2] ~ endog.sim.2[,1:4] )$coeff ) )
print( p.coeff.hat.lm.2 <- 
         cbind(lm( cont.hat.2[,9] ~ endog.sim.2[,1:4] )$coeff,
               lm( cont.hat.2[,10] ~ endog.sim.2[,1:4] )$coeff ) )
    # Looks good again

print( c.coeff.hat.2 <- 
         cbind( coeff_reg( cont.hat.2[,1], endog.sim.2[,1:4], opt$N, opt$lower, opt$upper, opt$cheby ),
                coeff_reg( cont.hat.2[,2], endog.sim.2[,1:4], opt$N, opt$lower, opt$upper, opt$cheby ) ) )
print( p.coeff.hat.2 <- 
         cbind( coeff_reg( cont.hat.2[,9], endog.sim.2[,1:4], opt$N, opt$lower, opt$upper, opt$cheby ),
                coeff_reg( cont.hat.2[,10], endog.sim.2[,1:4], opt$N, opt$lower, opt$upper, opt$cheby ) ) )
    # And... same for the coeff_reg versions


#### NOW: TRY WITH A GENUINE SIMULAITON FROM THIS SET OF COEFFICIENTS ####
kappa <- 25
exog.sim <- cbind( ar1_sim( n.draws * kappa, params$rho[1], params$sig.eps[1] ),
                    ar1_sim( n.draws * kappa, params$rho[1], params$sig.eps[1] ) )
    # Simulate the exogenous state



opt$iter <- 1
opt$n.gain <- 1
opt$sym.reg <- FALSE
opt$sr <- FALSE
sol.irbc.N.1 <- sol.irbc.iterate( coeff.init, opt, params, coeff.init.cont )
