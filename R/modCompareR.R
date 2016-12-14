mod.eval <- function( params, opt.lim=NULL, n.sim=100000, return.sims.err=FALSE ){
# Evaluates the model for the DS, linear and quadratic solutions and with theta=0
  
  #### 1. SET UP ####
  params.0 <- params
  params.0$theta <- 0
      # Establish the parameters
  
  #### 2. BASELINE SOLUTION ####
  message("**** Creating DS solutions ****")
  message("** theta > 0 **")
  baseline <- mod.gen(params, check=FALSE ) #, err.deets = TRUE )
  message("** theta = 0 **")
  baseline.0 <- mod.gen(params.0, check=FALSE ) #, err.deets = TRUE )
      # The DS solutions
  
  #### 3. SET UP OPTIONS ####
  upper <- baseline$ds.sol$upper
  lower <- baseline$ds.sol$lower
  n.endog <- ncol(baseline$ds.sol$coeff)
  n.exog <- nrow(baseline$ds.sol$coeff) - n.endog - 1
  n.cont <- ncol(baseline$ds.sol$coeff.cont)
  endog.init <- tail(baseline$ds.sol$ys, n.endog)
      # Extract from baseline
  
  exog.names <- c('shk1','shk2', 'P11', 'P22')
  endog.names <- c( 'B11', 'B22' )
  cont.names <- c( 'C1', 'C2', 'R_1', 'R_2', 'X11', 'X22', 'X12', 'X21', 
                   'P1', 'P2', 'P12', 'P21', 'E', 'Q', 'A1', 'A2', 'PN1', 'PN2', 
                   'PT1', 'PT2', 'CN1', 'CN2', 'CT1', 'CT2' )
  fwd.vars <- c('B11', 'E', 'R_1', 'R_2')
      # Model fundamentals
  
  opt <- list( lags=1, n.exog=n.exog, n.endog=n.endog, n.cont=n.cont, N=1, cheby=FALSE,
               upper = upper, lower=lower, quad=TRUE, n.quad=3,  burn=1000,
               kappa=25, n.sim=10000, eps = 1.0, delta=.02, endog.init=endog.init, 
               c.iter=100, c.tol=1e-07, c.gain=.8,
               k.iter=20, k.tol=1e-07, k.gain=.7,
               n.iter=2, n.tol=1e-05, n.gain=.5, 
               tol=1e-05, iter=6, model='irbc', cont.coeff.gain=.05,
               sr=TRUE, adapt.gain=TRUE, adapt.exp=15, image=FALSE,
               exog.names=exog.names, endog.names=endog.names, 
               cont.names=cont.names, fwd.vars=fwd.vars, mono="m1",
               sym.reg=FALSE, n.fwd=length(fwd.vars), ys=baseline$ds.sol$ys )
      # Set default global solution options
  
  baseline.sol <- list( params=params, opt=opt, 
                        coeff=baseline$ds.sol$coeff, 
                        coeff.cont=baseline$ds.sol$coeff.cont )
      # Create a solution object for the baseline
  baseline.sol.0 <- list( params=params, opt=opt, 
                          coeff=baseline.0$ds.sol$coeff, 
                          coeff.cont=baseline.0$ds.sol$coeff.cont )
  
  #### 4. LINEAR SOLUTIONS ####
  message("**** Creating linear solutions ****")
  message("\n** theta > 0 **")
  sol.idx <- 1
  if(!is.null(opt.lim)){
    for( nn in names(opt.lim) )
      opt[[nn]] <- opt.lim[[nn]][min(length(opt.lim[[nn]]),sol.idx)]
    # Paste the limited options over the top of the defaults
  }
  sol.1 <- sol.irbc.iterate( baseline$ds.sol$coeff, opt, params, 
                             baseline$ds.sol$coeff.cont )
  message("\n** theta = 0 **")
  sol.idx <- 2
  if(!is.null(opt.lim)){
    for( nn in names(opt.lim) )
      opt[[nn]] <- opt.lim[[nn]][min(length(opt.lim[[nn]]),sol.idx)]
    # Paste the limited options over the top of the defaults
  }
  sol.1.0 <- sol.irbc.iterate( baseline.sol$coeff, opt, params.0, baseline.sol$coeff.cont )
      # The linear solutions
  
  #### 4. QUADRATIC SOLUTIONS ####
  message("**** Creating quadratic solutions ****")
  message("\n** theta > 0 **")
  opt$N <- 2
  sol.idx <- 3
  if(!is.null(opt.lim)){
    for( nn in names(opt.lim) )
      opt[[nn]] <- opt.lim[[nn]][min(length(opt.lim[[nn]]),sol.idx)]
    # Paste the limited options over the top of the defaults
  }
  n.coeff <- idx_count( opt$N, n.exog + n.endog )
  idx.coeff <- apply( idx_create(opt$N, n.exog + n.endog), 
                      1, function(x) sum(x) <=1 )
  coeff.init <- matrix(0,n.coeff,n.endog)
  coeff.init.cont <- matrix(0,n.coeff,n.cont)
  coeff.init[ idx.coeff, ] <- baseline.sol$coeff # sol.1$coeff # 
  coeff.init.cont[ idx.coeff, ] <-  baseline.sol$coeff.cont # sol.1$coeff.cont # 
      # Set up the new initial guess
  sol.idx <- 4
  if(!is.null(opt.lim)){
    for( nn in names(opt.lim) )
      opt[[nn]] <- opt.lim[[nn]][min(length(opt.lim[[nn]]),sol.idx)]
    # Paste the limited options over the top of the defaults
  }
  sol.2 <- sol.irbc.iterate( coeff.init, opt, params, coeff.init.cont )
      # Quadratic solution
  message("\n** theta = 0 **")
  # coeff.init[ idx.coeff, ] <- sol.1.0$coeff
  # coeff.init.cont[ idx.coeff, ] <- sol.1.0$coeff.cont
      # Set up initial guess again
  sol.2.0 <- sol.irbc.iterate( coeff.init, opt, params.0, coeff.init.cont )
      # The solution
  
  #### 5. SIMULATIONS ####
  message( '**** Computing simulations ... ****')
  message( '      ... Devreux-Sutherland, theta > 0 ... ')
  sim.baseline <- sim.sol( baseline.sol, n.sim )
  message( '      ... Devreux-Sutherland, theta = 0 ... ')
  sim.baseline.0 <- sim.sol( baseline.sol.0, n.sim )
  sim.exog <- sim.baseline[,1:n.exog]
  message( '      ... linear, theta > 0 ... ')
  sim.sol.1 <- sim.sol( sol.1, n.sim, sim.exog )
  message( '      ... linear, theta = 0 ... ')
  sim.sol.1.0 <- sim.sol( sol.1.0, n.sim, sim.exog )
  message( '      ... quadratic, theta > 0 ... ')
  sim.sol.2 <- sim.sol( sol.2, n.sim, sim.exog )
  message( '      ... quadratic, theta = 0. ')
  sim.sol.2.0 <- sim.sol( sol.2.0, n.sim, sim.exog )
  
  #### 6. MEASURING THE ERRORS ####
  message( '**** Measuring errors ... ****')
  extra.args <- list( n.fwd=opt$n.fwd, y1.ss=opt$ys['Y1'] )
  message( '      ... Devreux-Sutherland, theta > 0 ... ')
  baseline.err <- sim.err( sim.baseline, baseline.sol, extra.args, params = params.0 )
  message( '      ... linear, theta > 0 ... ')
  err.sol.1 <- sim.err( sim.sol.1, sol.1, extra.args, params = params.0 )
  message( '      ... linear, theta = 0 ... ')
  err.sol.1.0 <- sim.err( sim.sol.1.0, sol.1.0, extra.args, params = params.0 )
  message( '      ... quadratic, theta > 0 ... ')
  err.sol.2 <- sim.err( sim.sol.2, sol.2, extra.args, params = params.0 )
  message( '      ... quadratic, theta = 0. ')
  err.sol.2.0 <- sim.err( sim.sol.2.0, sol.2.0, extra.args, params = params.0 )
      # Should always have params.0 because that is the true model to compare to.
  
  #### 7. FORMATTING THE OUTPUT ####
  baseline.0.sumy <- summary(sim.baseline.0)
  l.sol <- list( local=baseline.sol, global.1=sol.1, 
                 global.2=sol.2, global.1.0=sol.1.0, global.2.0=sol.2.0 )
  l.sim <- list( local=sim.baseline, global.1=sim.sol.1, 
                 global.2=sim.sol.2, global.1.0=sim.sol.1.0, global.2.0=sim.sol.2.0 )
  l.err <- list( local=baseline.err, global.1=err.sol.1, 
                 global.2=err.sol.2, global.1.0=err.sol.1.0, global.2.0=err.sol.2.0 )
  bs.log <- sapply( l.sim, 
                  function(sim) cor( diff(sim[,'C1']-sim[,'C2']), diff(sim[,'Q']) ) )
  bs.level <- sapply( l.sim, 
                  function(sim) cor( exp(sim[,'C1'])-exp(sim[,'C2']), exp(sim[,'Q']) ) )
  uip.coeff <- sapply( l.sim, 
                  function(sim) lm( diff(sim[,'E']) ~ (sim[,'R_1']-sim[,'R_2'])[-n.sim] )$coeff[2] )
  names(uip.coeff) <- names(bs.log)
  err.ave.sumy <- sapply( l.err, function(x) apply(x, 2, mean) )
  err.abs.sumy <- sapply( l.err, function(x) apply(abs(x), 2, mean) )
  alpha.tilde <- sapply( l.sol, function(sol) mean(sol$coeff[1,]) )
  
  if(!return.sims.err){
    return( list( bs.log=bs.log, bs.level=bs.level, uip.coeff=uip.coeff,
                  err.ave.sumy=err.ave.sumy, err.abs.sumy=err.abs.sumy,
                  alpha.tilde=alpha.tilde, baseline.0.sumy=baseline.0.sumy ) )
  }else{
    return( list( bs.log=bs.log, bs.level=bs.level, uip.coeff=uip.coeff,
                  err.ave.sumy=err.ave.sumy, err.abs.sumy=err.abs.sumy,
                  alpha.tilde=alpha.tilde, l.sol=l.sol,
                  l.sim=l.sim, l.err=l.err, baseline.0.sumy=baseline.0.sumy ) )
  }
  
  
  
}

