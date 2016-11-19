#########################################################################
# simSol.R
#
# Code to create a standardized simulation from a solution object
# Philip Barrett, Washington DC
# 09nov2016
#########################################################################

sim.sol <- function(sol, n.sim=100000, sim.exog=NULL ){
# Creates a standardized simulation
  
  opt <- sol$opt
  params <- sol$params
      # Copy the solution objects
  
  if(is.null(sim.exog)) 
    sim.exog <- var1_sim( n.sim, params$rho, params$sig.eps )
    # sim.exog <- sapply( 1:n.exog, function(i) ar1_sim( n.sim, 
    #                              params$rho[i], params$sig.eps[i] ) )
      # Create the initial solution if required
  endog.sim <- endog_sim( n.sim, sim.exog, sol$coeff, opt$N, opt$upper, 
                          opt$lower, c(sol$coeff[1,]), opt$cheby, 1, 0, TRUE )
      # The endogenous variables
  nn <- c( opt$exog.names, opt$endog.names )
  colnames(endog.sim) <- c( nn, paste0( nn, '(-1)') )
      # Name the variables
  cont.sim <- cont_sim( endog.sim, sol$coeff.cont, opt$N, opt$n.endog, opt$n.exog,
                        opt$upper, opt$lower, opt$cheby )
      # The contemporaneous variables
  colnames( cont.sim ) <- opt$cont.names
      # Names
  return( cbind( endog.sim, cont.sim ) )
}

sim.err <- function( sim, sol, extra.args=NULL ){
# Computes the errors on a simulation
  opt <- sol$opt
  params <- sol$params
      # Copy option and parameters
  pred <- contemp_eqns_irbc_grid( sim, opt$lags, params, opt$n.exog, 
                                  opt$n.endog, opt$n.cont, extra.args, opt$model )
  colnames(pred) <- c( opt$endog.names, opt$cont.names )
  pred[,fwd.vars] <- 
    euler_hat_grid( sol$coeff, sol$coeff.cont, sim, opt$lags, params, 
                    opt$n.exog, opt$n.endog, opt$n.cont, opt$n.fwd, params$rho, 
                    params$sig.eps, 0, opt$N, opt$upper, opt$lower, opt$cheby, 
                    matrix(0,1,1), TRUE, opt$n.quad, opt$model, opt$mono )
  err <- sim[,c( endog.names, cont.names )] - pred
  return(err)
}


