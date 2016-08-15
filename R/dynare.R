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

mod.gen <- function(params, lower, upper, nsim=1e6, burn=1e4 ){
# Generates the coefficient matrices 
  mod.create(params)
      # Create the model
  mod.run()
      # Run the model
  mod <- mod.read()
      # The model save file as a list
  sim <- stoch_simDS( mod, params$sig.eps, params$betta, nsim, burn, 0 )
      # The simulation (Change the hard-coded zero if the index of y1 is other than 1)
  nn <- c( sub("\\s+$", "", (unlist(mod$varo))), 'af1' )
  colnames(sim) <- nn
      # Name the columns of the simulation
  endog.order <- c( 'NFA', 'Z1', 'Z2' )
  cont.order <- c( 'C1', 'C2', 'rb1', 'rb2', 'X11', 'X22', 'X12', 'X21', 
                   'P1', 'P2', 'P11', 'P22', 'P12', 'P21', 'E', 'Q', 'af1',
                   'Y1', 'Y2', 'cd', 'cg' )
      # Variable names for the DS-style solution
  sim.exog <- sim[, c('A1','A2')]
  sim.endog <- sim[, endog.order]
  sim.cont <- sim[, cont.order]
      # Separate out the exogenous and endogenous states and the static variables
  
  
  browser()
}










