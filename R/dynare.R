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

mod.gen <- function(params, nsim=1e6, burn=1e4, cheby=FALSE, check=FALSE, n.nodes=3 ){
# Generates the coefficient matrices and vector of equation errors
  
  ### 0. Hard-coding ###
  N <- 1
      # The order of the model
  n.lag <- 1
      # Number of lags
  n.fwd <- 4
      # Number of Euler equations
  
  ### 1. Run and extract the model ###
  message('Creating & solving dynare model')
      # Screen updating
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
  message("Mean home asset holding = ", round( mod$alpha.tilde, 4 ) )
  
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
  message('Generating coefficient rules for error evaluation')
      # Screen updating
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
  message('Evaluating errors from the dynare solution')
      # Screen updating
  X <- cbind( sim.exog[-1,], sim.endog[-1,], sim.exog[-nsim,], sim.endog[-nsim,], 
              sim.cont[-1,] )
      # Redefine X for evaluating the errors
  extra.args <- list(n.fwd=4, y1.ss=mod$ys['Y1'])
      # The number of forward-looking equations
  cont.nfa <- contemp_eqns_irbc_grid( X, 1, params, n.exog, n.endog, 
                                      n.cont, extra.args, "ds" )
      # The controls + NFA predictions
  colnames(cont.nfa) <- c( cont.order, 'NFA' )
      # Rename the columns
  cont.pred <- cont.nfa[,1:n.cont]
      # The control predictors (Q and af1 to be added after the Euler eqs)
  cont.err <- abs( sim.cont[-1,] - cont.pred )[,(!cont.order %in% c('Q', 'af1') )]
      # The errors on the static equations
  nfa.err <- abs( cont.nfa[, 1+n.cont] - sim.endog[-1,'NFA'] )
      # The NFA error
  
  euler.pred <- euler_hat_grid( endog.coeff, cont.coeff, X, n.lag, params,
                                n.exog, n.endog, n.cont, n.fwd, params$rho, 
                                params$sig.eps, 0, N, upper, lower, cheby,
                                matrix(0,1,1), TRUE, n.nodes, "ds" )
      # Generate predictions from the Euler equation errors
  euler.orig <- cbind( sim.endog[-1,'Z1'], sim.cont[-1,'af1'], 
                       sim.endog[-1,'Z2'], sim.cont[-1,'Q'] )
  euler.err <- abs( euler.pred - euler.orig )
      # The Euler equation errors
  err <- apply( cbind( cont.err, nfa.err, euler.err ), 1, max )
      # Extract the distribution of maximum equation errors
  
  ### 4. Now for the alternative-state representation ###
  # States are now A1, A2, NFA, and af1
  endog.order.alt <- c( 'NFA', 'af1' )
  cont.order.alt <- c( 'C1', 'C2', 'rb1', 'rb2', 'X11', 'X22', 'X12', 'X21', 
                   'P1', 'P2', 'P11', 'P22', 'P12', 'P21', 'Z1', 'Z2', 'E', 'Q', 
                   'Y1', 'Y2', 'cd', 'cg' )
      # Variable names for the DSalternative-state respresentation
  sim.endog.alt <- sim[, endog.order.alt]
  sim.cont.alt <- sim[, cont.order.alt]
      # Separate out the exogenous and endogenous states and the static variables
  n.endog <- length( endog.order.alt )
  n.cont <- length( cont.order.alt )
      # Numbers of types of variables
  message('Generating coefficient rules for alternative-state representation')
      # Screen updating
  upper <- c(  3 * sd.x, mod$ys[endog.order.alt] + .5 )
  lower <- c( -3 * sd.x, mod$ys[endog.order.alt] - .5 )
      # The bounds of the states
  endog.reg <- sim.endog.alt[-nsim,]
  exog.reg <- sim.exog[-1,]
  X <- cbind( exog.reg, endog.reg )
      # The X-variables for the regressions
  n.X <- 1 + length(endog.order.alt) + length(exog.order) 
      # The number of X variables is the number of states plus one (the constant
      # term)
  endog.coeff <- matrix( 0, n.X, n.endog )
  cont.coeff <- matrix( 0, n.X, n.cont )
      # The coefficient matrices  
  
  for(i in 1:n.endog) endog.coeff[,i] <- coeff_reg( sim.endog.alt[,i][-1], X, N,
                                                    lower, upper, cheby )
  for(i in 1:n.cont) cont.coeff[,i] <- coeff_reg( sim.cont.alt[,i][-1], X, N, 
                                                  lower, upper, cheby )
      # Populate the coefficient matrices
  alt.sol <- list( endog.coeff=endog.coeff, cont.coeff=cont.coeff, 
                  upper=upper, lower=lower )
      # Details of the alternative solution

  return( list( mod=mod, ds.sol=ds.sol, err=err, alt.sol=alt.sol ) )
}


# X <- cbind( sim.endog[-1,], sim.endog[-nsim,], sim.cont[-1,], sim.cont[-nsim,] )
# # The matrix of points on which the error will be assessed
# cont.nfa <- contemp_eqns_irbc_grid( X, 1, params)






