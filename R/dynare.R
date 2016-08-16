#########################################################################
# dynare.R
#
# Code to utilize dynare to produce approximate local solutions via the 
# Devreux-Sutherland method
# Philip Barrett, Chicago
# 11aug2016
#########################################################################

params.convert <- function( st ){
# Converts the parameters to the dynare names
  if( st == 'share' ) return('alph')
  if( st == 'gamma' ) return('RH')
  if( st == 'P1.bar' ) return('p1bar')
  if( st == 'P2.bar' ) return('p2bar')
  if( st == 'betta' ) return('BT')
  if( st == 'rho' ) return('rho')
  if( st == 'sig.eps' ) return('sigeps')
  if( st == 'eta' ) return('eta')
  
}

mod.create <- function( params ){
# Creates the .mod file for a list of parameters (plus some cleaning up)
  
  if( file.exists('matlab/adamsBarrett.mod') ){
    tstamp <- file.info('R/dynare.R')$mtime
        # Timestamp of the old file
    new.name <- paste0( 'matlab/old/', format(tstamp, '%Y-%m-%d-%H%M%S'), '.mod' )
    file.rename( 'matlab/adamsBarrett.mod', new.name )
        # Move the old file
  }
  
  write( paste( paste(rep('%', 40 ), collapse=''), 
                '\n% AUTO-GENERATED CODE FROM DYNARE.R \n% CREATED ',
                gsub( ':', '', gsub( ' ', '-', Sys.time() ) ), '\n',
                paste( rep('%', 40 ), collapse=''), '\n',
                sep=' ' ), file='matlab/adamsBarrett.mod', append=FALSE)
  write( readLines('matlab/head.txt',warn=FALSE ), 
         file='matlab/adamsBarrett.mod', append=TRUE )
      # Top part of the code
  par.names <- names(params)
      # Parameter names
  for( i in 1:length( params ) ){
    st.par <- params.convert( par.names[i] )
    if( length(params[[i]]) == 1){
      write( paste( st.par, '=', params[[i]], ';', sep=' ' ), 
             file='matlab/adamsBarrett.mod', append=TRUE )
    }else{
      for( j in 1:length(params[[i]]) ){
        write( paste( st.par, j, ' = ', params[[i]][j], ' ;', sep='' ), 
               file='matlab/adamsBarrett.mod', append=TRUE )
      }
    }
  }
      # Write parameter values
  write( readLines('matlab/tail.txt', warn=FALSE ), 
         file='matlab/adamsBarrett.mod', append=TRUE )
      # Bottom part of the code
}

mod.run <- function(){
# Runs the DS model
  system( 'matlab -nojvm < matlab/launcher.m' )
}

mod.read <- function(){
# Reads the saved parameteres from the DS-style solution
  return( readMat('matlab/Rvars.mat') )
}

mod.gen <- function(params, nsim=1e6, burn=1e4, cheby=FALSE, check=FALSE ){
# Generates the coefficient matrices and vector of equation errors
  
  ### 1. Run and extract the model ###
  N <- 1
      # The order of the model
  mod.create(params)
      # Create the model
  mod.run()
      # Run the model
  mod <- mod.read()
      # The model save file as a list
  sim <- stoch_simDS( mod, params$sig.eps, params$betta, nsim, burn, 0 )
      # The simulation (Change the hard-coded zero if the index of y1 is other than 1)
  mod$varo <- c( sub("\\s+$", "", (unlist(mod$varo))), 'af1' )
  colnames(sim) <- mod$varo
  mod$ys <- c( mod$ys, mod$alpha.tilde )
  names(mod$ys) <- mod$varo
      # Assign names
  exog.order <- c('A1','A2')
  endog.order <- c( 'NFA', 'Z1', 'Z2' )
  cont.order <- c( 'C1', 'C2', 'rb1', 'rb2', 'X11', 'X22', 'X12', 'X21', 
                   'P1', 'P2', 'P11', 'P22', 'P12', 'P21', 'E', 'Q', 'af1',
                   'Y1', 'Y2', 'cd', 'cg' )
      # Variable names for the DS-style solution
  sim.exog <- sim[, exog.order]
  sim.endog <- sim[, endog.order]
  sim.cont <- sim[, cont.order]
      # Separate out the exogenous and endogenous states and the static variables
  n.exog <- length( exog.order )
  n.endog <- length( endog.order )
  n.cont <- length( cont.order )
      # Numbers of types of variables
  
  ### 2. Run regressions to generate coefficients ###
  sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
      # The standard deviation of the exogenous state
  upper <- c(  3 * sd.x, mod$ys[endog.order] + .25 )
  lower <- c( -3 * sd.x, mod$ys[endog.order] - .25 )
      # The bounds of the states
  endog.reg <- sim.endog[-nsim,]
  exog.reg <- sim.exog[-1,]
  X <- cbind( exog.reg, endog.reg )
      # The X-variables for the regressions
  n.X <- 1 + length(endog.order) + length(exog.order) 
      # The number of X variables is the number of states plus one (the constant
      # term)
  endog.coeff <- matrix( 0, n.X, n.endog )
  cont.coeff <- matrix( 0, n.X, n.cont )
      # The coefficient matrices  
  for(i in 1:n.endog) endog.coeff[,i] <- coeff_reg( sim.endog[,i][-1], X, N,
                                                    lower, upper, cheby )
  for(i in 1:n.cont) cont.coeff[,i] <- coeff_reg( sim.cont[,i][-1], X, N, 
                                                  lower, upper, cheby )
      # Populate the coefficient matrices

  if(check){
    sim.endog.alt <- endog_sim( nsim, sim.exog, endog.coeff, N, upper, lower, 
                                   sim.endog[1,], cheby, 1, 0, TRUE )
    browser()
        # Do checks in debug - looks fine to me! :)
  }
  
  ds.sol <- list( endog.coeff=endog.coeff, cont.coeff=cont.coeff, 
                  upper=upper, lower=lower )
      # Details of the Devreux-Sutherland solution
  
  ### 3. Evaluate equation errors ###
  X <- cbind( sim.exog[-1,], sim.endog[-1,], sim.exog[-nsim,], sim.endog[-nsim,], 
              sim.cont[-1,] )
      # Redefine X for evaluating the errors
  extra.args <- list(n.fwd=4, y1.ss=mod$ys['Y1'])
      # The number of forward-looking equations
  cont.nfa <- contemp_eqns_irbc_grid( X, 1, params, n.exog, n.endog, 
                                      n.cont, extra.args, "ds" )
      # The controls + NFA predictions
  colnames(cont.nfa) <- c( cont.order, 'NFA' )
  
  browser()
  
  ## Euler equation errors here
  
  
  return( list( ds.sol=ds.sol, err=NA, alt.states=NA ) )
}


# X <- cbind( sim.endog[-1,], sim.endog[-nsim,], sim.cont[-1,], sim.cont[-nsim,] )
# # The matrix of points on which the error will be assessed
# cont.nfa <- contemp_eqns_irbc_grid( X, 1, params)






