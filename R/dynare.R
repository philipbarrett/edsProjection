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
  if( st == 'theta' ) return('theta')
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
# Reads the saved parameters from the DS-style solution
  return( readMat('matlab/Rvars.mat') )
}

mod.gen <- function(params, nsim=1e6, burn=1e4, cheby=FALSE, check=TRUE, 
                    n.nodes=3, err.deets=FALSE, sim.stats=FALSE, return.sim=FALSE ){
# Generates the coefficient matrices and vector of equation errors
  
  ### 0. Hard-coding ###
  N <- 1
      # The order of the model
  n.lag <- 1
      # Number of lags
  n.fwd <- 4
      # Number of Euler equations
  n.sample <- 20000
      # The number of rows of the simulation to sample for the coefficient rules
      # & such
  sr <- FALSE
  eps <- 1
  delta <- .025
      # State reduction controls
  
  ### 1. Run and extract the model ###
  message('Creating & solving dynare model')
      # Screen updating
  mod.create(params)
      # Create the model
  mod.run()
      # Run the model
  mod <- mod.read()
      # The model save file as a list

  message("****** Mean home asset holding = ", round( mod$alpha.tilde, 4 ), " ******" )
      # Screen updating
  
  y1.idx <- which( unlist( mod$varo ) == "Y1 " ) - 1
      # The location of Y1
  sim <- stoch_simDS( mod, params$sig.eps, params$betta, nsim, burn, y1.idx )
      # The simulation
    
  varo <- c( sub("\\s+$", "", (unlist(mod$varo))), 'af1' )
  ys <- c( mod$ys, mod$alpha.tilde )
  names(ys) <- varo
  colnames(sim) <- varo
      # Assign names
  B11 <- exp(ys[y1.idx+1]) *# params$betta * #exp( -params$theta * sim[,'C1'] ) *
    sim[,'af1'] * exp( sim[,'R_1'] + sim[,'P1'] )
  B22 <- - ( sim[,'NFA'] - #params$betta * #exp( -params$theta * sim[,'C1'] ) * 
             sim[,'af1'] ) * exp(ys[y1.idx+1] ) * exp( sim[,'R_2'] + sim[,'P1']  - sim[,'E'])
      # Create nominal debt series:  NB THESE ARE END OF PERIOD
      # B22 is negative because it is held by households in country 2
  sim <- cbind( sim, B11, B22 )
      # Append B11, B22 to simulation
  varo <- c( varo, 'B11', 'B22' )
  ys <- c( ys, B11=mean(B11), B22=mean(B22) )
      # Assign names
  
  exog.order <- c('A1','A2', 'P11', 'P22')
  endog.order <- c( 'B11', 'B22' )
  cont.order <- c( 'C1', 'C2', 'R_1', 'R_2', #'rb1', 'rb2', 
                   'X11', 'X22', 'X12', 'X21', 
                   'P1', 'P2', 
                   'P12', 'P21', 'E', 'Q' ) #,
                   # 'Y1', 'Y2', 'cd', 'cg' ) # , 'NFA', 'af1' )
      # Variable names for the DS-style solution
  fwd.order <- c( 'B11', 'E', 'R_1', 'R_2')
      # The forward-looking equation variables
  sim.exog <- sim[, exog.order]
  sim.endog <- sim[, endog.order]
  sim.cont <- sim[, cont.order]
      # Separate out the exogenous and endogenous states and the static variables
  n.exog <- length( exog.order )
  n.endog <- length( endog.order )
  n.cont <- length( cont.order )
      # Numbers of types of variables
  
  ## 1a. Return some stats from the simulation
#   browser()
  if( sim.stats ){
    bs.basic <- cor(diff(sim.cont[,'C1'] - sim.cont[,'C2']), diff(sim.cont[,'Q']) )
    bs.soph <- cor(diff(sim.cont[,'C1'] - sim.cont[,'C2']-sim.cont[,'Q'] / params$gamma), 
                   diff(sim.cont[,'E']) )
    bs.level <- cor(exp(sim.cont[,'C1']) - exp(sim.cont[,'C2']), exp(sim.cont[,'Q']) )
    sim.R1 <- sim.cont[,'R_1']
    sim.R2 <- sim.cont[,'R_2']
    sim.rdiff <- sim.R1 - sim.R2
    sim.e.app <- diff(sim.cont[,'E'])
    uip.coeffs <- lm(sim.e.app~sim.rdiff[-length(sim.rdiff)])$coeff
    sim.uip.err <- sim.e.app - sim.rdiff[-length(sim.rdiff)]
    uip.nfa <- lm( sim.uip.err ~ sim.endog[,'NFA'][-nrow(sim.endog)])$coeff 
        # The UIP coefficients
    out <- c( mod$alpha.tilde, bs.basic, bs.soph, bs.level, uip.coeffs[2], uip.nfa[2] )
    names(out) <- c('alpha.tilde', 'bs.basic', 'bs.soph', 'bs.level', 'uip.coeff', 'uip.nfa' )
    return( out )
  }

  if( return.sim ){
    return( cbind( sim.exog, sim.endog, sim.cont ) )
  }
  
  ### 2. Run regressions to generate coefficients ###
  message('Generating coefficient rules for error evaluation')
      # Screen updating
  if(sr){
    # browser()
    idces <- p_eps_cheap_const_idx( cbind( sim.exog[-1,], sim.endog[-nsim,]),
                                    eps, delta )
        # Need to include the lagged state in the evaluation set for the
        # equilibrium condition
    endog.reg <- sim.endog[-nsim,][idces,]
    exog.reg <- sim.exog[-1,][idces,]
  }else{
    endog.reg <- sim.endog[-nsim,]
    exog.reg <- sim.exog[-1,]
  }
  X <- cbind( exog.reg, endog.reg )
      # The X-variables for the regressions
  
  var.x <- solve( diag(n.exog) - params$rho %*% params$rho, params$sig.eps )
  sd.x <- sqrt(diag(var.x))
      # The standard deviation of the exogenous state
  
  var.check <- max(abs(var(sim.exog) - var.x))
  message( "Simulated and analytical exogenous variance differ by at most ", 
           round(var.check,8) )
  
  upper <- c(   3 * sd.x, ys[endog.order] + 3 * apply(sim.endog,2,sd) )
  lower <- c( - 3 * sd.x, ys[endog.order] - 3 * apply(sim.endog,2,sd) )
      # The bounds of the states
  n.X <- 1 + length(endog.order) + length(exog.order) 
      # The number of X variables is the number of states plus one (the constant
      # term)
  coeff <- matrix( 0, n.X, n.endog )
  coeff.cont <- matrix( 0, n.X, n.cont )
  colnames(coeff) <- endog.order
  colnames(coeff.cont) <- cont.order
      # The coefficient matrices  
  for(i in 1:n.endog) coeff[,i] <- coeff_reg( sim.endog[,i][-1], X, N,
                                                    lower, upper, cheby )
  for(i in 1:n.cont) coeff.cont[,i] <- coeff_reg( sim.cont[,i][-1], X, N, 
                                                  lower, upper, cheby )
      # Populate the coefficient matrices
# browser()
  if(check){
    message('Verifying dynare solution conversion')
    sim.endog.alt <- endog_sim( nsim-1, sim.exog[-1,], coeff, N, upper, lower, 
                                sim.endog[1,], cheby, 1, 0, TRUE )
        # Do checks in debug - looks fine to me! :)
    pc.err <- abs( sim.endog[-1,] / sim.endog.alt[,n.exog+1:n.endog] - 1 ) * 100
        # The percentage error
    message( "  % absolute error on endogenous state rules:" )
    print( t( apply( pc.err, 2, function(x) quantile( x, c(.25, .5, .75, .9, .95, .99 ) ) ) ) )
    # browser()
  }
  ds.sol <- list( coeff=coeff, coeff.cont=coeff.cont, 
                  upper=upper, lower=lower, ys=ys )
      # Details of the Devreux-Sutherland solution
  
  if( !err.deets ){
    l.err <- NA
    return( list( mod=mod, ds.sol=ds.sol, l.err=l.err ) )
  }
  
  ### 3. Evaluate equation errors ###
  # browser()
  message('Evaluating errors from the dynare solution')
      # Screen updating
  X <- cbind( sim.exog[-1,], sim.endog[-1,], sim.exog[-nsim,], sim.endog[-nsim,], 
              sim.cont[-1,] )[sample(nrow(X),size=n.sample,replace=TRUE),]
      # Redefine X for evaluating the errors.  Timing is good.
  
  extra.args <- list(n.fwd=4, y1.ss=ys['Y1'])
      # The number of forward-looking equations
  pred <- contemp_eqns_irbc_grid( X, 1, params, n.exog, n.endog, 
                                      n.cont, extra.args, "irbc" )
      # The contemporaneous predictors
  colnames(pred) <- c( endog.order, cont.order )
      # Rename the columns
  pred[,fwd.order] <- euler_hat_grid( coeff, coeff.cont, X, n.lag, params,
                                      n.exog, n.endog, n.cont, n.fwd, params$rho, 
                                      params$sig.eps, 0, N, upper, lower, cheby,
                                      matrix(0,1,1), TRUE, n.nodes, "irbc" )
      # Generate predictions from the Euler equation errors
  err <- pred - X[,c(n.exog+1:n.endog, 
                          (1+n.lag)*(n.exog+n.endog)+1:n.cont)]
  
  bias <- apply(err, 2, mean)
  aad <- apply(abs(err), 2, mean)
  mean.max.abs.err <- mean(apply(abs(err),1,max))
  max.abs.err.dist <- density( apply( abs( err ), 1, max ) )
  log.max.abs.err.dist <- density( log( apply( abs( err ), 1, max ), 10 ) )
        # The distribution of maximum equation errors
  if( err.deets ){
    l.err <- list( bias=bias, aad=aad, max.abs.err.dist=max.abs.err.dist,
                   mean.max.abs.err=mean.max.abs.err,
                   log.max.abs.err.dist=log.max.abs.err.dist, pred=pred, err=err, 
                   X=X[,c(n.exog+1:n.endog,(1+n.lag)*(n.exog+n.endog)+1:n.cont)] )
  }else{
    l.err <- list( bias=bias, aad=aad, mean.max.abs.err=mean.max.abs.err,
                   max.abs.err.dist=max.abs.err.dist,
                   log.max.abs.err.dist=log.max.abs.err.dist )
  }
    

#                  max.abs.err.dist=max.abs.err.dist,
#                  err=err)
        # The list of different error measures
  
#   ### 4. Now for the alternative-state representation ###
#   # States are now A1, A2, NFA, and af1
#   endog.order.alt <- c( 'NFA', 'af1' )
#   cont.order.alt <- c( 'C1', 'C2', 'rb1', 'rb2', 'X11', 'X22', 'X12', 'X21', 
#                    'P1', 'P2', 'P11', 'P22', 'P12', 'P21', 'Z1', 'Z2', 'E', 'Q', 
#                    'Y1', 'Y2', 'cd', 'cg' )
#       # Variable names for the DSalternative-state respresentation
#   sim.endog.alt <- sim[, endog.order.alt]
#   sim.cont.alt <- sim[, cont.order.alt]
#       # Separate out the exogenous and endogenous states and the static variables
#   n.endog <- length( endog.order.alt )
#   n.cont <- length( cont.order.alt )
#       # Numbers of types of variables
#   message('Generating coefficient rules for alternative-state representation')
#       # Screen updating
#   upper <- c(  3 * sd.x, ys[endog.order.alt] + .5 )
#   lower <- c( -3 * sd.x, ys[endog.order.alt] - .5 )
#       # The bounds of the states
#   idx.sample <- sample(nsim-1,size=n.sample,replace=TRUE)
#   endog.reg <- sim.endog.alt[-nsim,][idx.sample,]
#   cont.reg <- sim.cont.alt[-1,][idx.sample,]
#   exog.reg <- sim.exog[-1,][idx.sample,]
#   X <- cbind( exog.reg, endog.reg )
#       # The X-variables for the regressions
#   n.X <- 1 + length(endog.order.alt) + length(exog.order) 
#       # The number of X variables is the number of states plus one (the constant
#       # term)
#   coeff <- matrix( 0, n.X, n.endog )
#   coeff.cont <- matrix( 0, n.X, n.cont )
#       # The coefficient matrices  
#   browser()
#   for(i in 1:n.endog) coeff[,i] <- coeff_reg( endog.reg[,i], X, N,
#                                                     lower, upper, cheby )
#   for(i in 1:n.cont) coeff.cont[,i] <- coeff_reg( cont.reg[,i], X, N, 
#                                                   lower, upper, cheby )
#       # Populate the coefficient matrices
#   alt.sol <- list( coeff=coeff, coeff.cont=coeff.cont, 
#                   upper=upper, lower=lower )
#       # Details of the alternative solution
  alt.sol <- NA

  return( list( mod=mod, ds.sol=ds.sol, l.err=l.err ) )
}


# X <- cbind( sim.endog[-1,], sim.endog[-nsim,], sim.cont[-1,], sim.cont[-nsim,] )
# # The matrix of points on which the error will be assessed
# cont.nfa <- contemp_eqns_irbc_grid( X, 1, params)






