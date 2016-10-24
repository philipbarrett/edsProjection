



stat.ds <- function( params, stat="err", opt=NULL, bo.plot = FALSE, bo.nl=TRUE ){
# Computes a given statistic from a sample for both the DS & NL model solutions.
# Stat may be: err, bs for Bakus-Smith correlation, hb for home bias, or uip 
# for UIP test.
  
  ## 1. THE DS SOLUTION ##
  l.coeffs <- mod.gen(params, err.deets=TRUE, n.nodes=2 )
      # The dynare solution
  exog.names <- c('A1','A2')
  endog.names <- c( 'NFA', 'af1' )
  cont.names <- c( 'C1', 'C2', 'rb1', 'rb2', 'X11', 'X22', 'X12', 'X21', 
                   'P1', 'P2', 'P11', 'P22', 'P12', 'P21', 'E', 'Q', 
                   'Y1', 'Y2', 'cd', 'cg' )
  if( is.null(opt)){
    opt <- list( lags=1, n.exog=2, n.endog=4, n.fwd=4, n.cont=20, N=2, cheby=FALSE,
#                  lags=1, n.exog=2, n.endog=2, n.fwd=4, n.cont=20, N=2, cheby=FALSE,
                 upper = l.coeffs$ds.sol$upper, lower=l.coeffs$ds.sol$lower, quad=TRUE, 
                 #                n.quad=4,  burn=1000, kappa=25, n.sim=20000, eps = 1, delta=.01, 
                 n.quad=2,  burn=1000, kappa=40, n.sim=20000, eps = .8, delta=.025, # eps=.5
                 endog.init=l.coeffs$mod$ys[c('NFA','af1')], 
                 fwd.vars=c('rb1','rb2','af1','Q'),
                 exog.names=exog.names, endog.names=endog.names, cont.names=cont.names,
                 c.iter=20, c.tol=1e-07, c.gain=.8,
                 k.iter=20, k.tol=1e-06, k.gain=.05, # k.gain=.05,
                 n.iter=4, n.tol=1e-05, n.gain=.1, 
                 #                tol=1e-05, iter=1, model='ds',
                 tol=1e-05, iter=3, model='ds',
                 sr=TRUE, adapt.gain=FALSE, adapt.exp=15, image=FALSE,
                 sym.reg=FALSE, ys=l.coeffs$mod$ys )
    # Solution options
  }

  
  if( bo.nl ){

    ## 2. THE NL SOLUTION ##
    message("*** Creating the nonlinear solution ***")

        # Variable names
    
#     opt$N <- 1
#     idx.1 <- idx_create( opt$N, opt$n.endog+opt$n.exog )
#         # The indices for the parameters of the DS solution
#     coeff.init <- matrix( 0, nrow(idx.1), opt$n.endog )
#     coeff.init[ which( apply( idx.1, 1, sum ) < 2 ), ] <- l.coeffs$ds.sol$coeff
#     coeff.init.cont <- matrix( 0,  nrow(idx.1), opt$n.cont )
#     coeff.init.cont[ which( apply( idx.1, 1, sum ) < 2 ), ] <- l.coeffs$ds.sol$coeff.cont
#         # The initial coefficients
#     sol.1 <- sol.irbc.iterate( coeff.init, opt, params, 
#                                coeff.init.cont, FALSE, TRUE )
#         # Linear solution
  
    opt$N <- 2

    idx.2 <- idx_create( opt$N, opt$n.endog+opt$n.exog )
        # The indices for the parameters of the DS solution
    coeff.init <- matrix( 0, nrow(idx.2), opt$n.endog )
#     coeff.init[ which( apply( idx.2, 1, sum ) < 2 ), ] <- sol.1$coeff
    coeff.init[ which( apply( idx.2, 1, sum ) < 2 ), ] <- l.coeffs$ds.sol$coeff
    coeff.init.cont <- matrix( 0,  nrow(idx.2), opt$n.cont )
    coeff.init.cont[ which( apply( idx.2, 1, sum ) < 2 ),  ] <- l.coeffs$ds.sol$coeff.cont
#     coeff.init.cont[ which( apply( idx.2, 1, sum ) < 2 ),  ] <- sol.1$coeff.cont
  
    bo.nl <- tryCatch( 
      sol.2 <- sol.irbc.iterate( coeff.init, opt, params, 
                                 coeff.init.cont, FALSE, TRUE ),
                  error=function(cond) return(FALSE) )
    if( !is.logical(bo.nl) ) bo.nl <- TRUE
        # Create the solution
    if( !bo.nl )
      warning("NONLINEAR SOLUTION FAILED")
  #   browser()
  }
  
  ## 3. OPTION 1: COMPUTE THE ERRORS ##
  if( stat %in% c( "err", "all" ) ){
    n.sample <- nrow(l.coeffs$l.err$err)
#     n.sample <- 20000
        # Number of sample points
    message("*** Evaluating errors ***")
    if( bo.nl ){
      nl.check <- sol.irbc.check(sol.2)
          # Resimulate
      x.sample <- sample(nrow(nl.check$endog.sim)-1,size=n.sample,replace=TRUE)
      nfa.nl.idx <- opt$n.exog + 1
          # Sample 
      nl.bias <- mean( nl.check$err[x.sample,], na.rm = TRUE )
      nl.aad <- mean( abs( nl.check$err[x.sample,] ), na.rm = TRUE )
      nl.mean.max.abs.err <- mean( apply( abs( nl.check$err[x.sample,] ), 1, max) , na.rm = TRUE )
          # The key statistics for the errors
    }else{
      nl.bias <- nl.aad <- nl.mean.max.abs.err <- NA
    }
    
    if( bo.plot ){
      plt.data.ds <- data.frame( NFA=l.coeffs$l.err$X[,'NFA'], 
                                 NFA.err=log(abs(l.coeffs$l.err$err[,'NFA']), 10 ) )
      plt.data.nl <- data.frame( NFA=nl.check$endog.sim[x.sample,nfa.nl.idx], 
                                 NFA.err=log(abs(nl.check$err[x.sample,'NFA']), 10 ) )
      plt.data <- data.frame( rbind( plt.data.ds, plt.data.nl) , 
                              type=c(rep("DS", n.sample), rep("nonlinear", n.sample)))
      ggplot(plt.data, aes(NFA,NFA.err, colour=type) ) + geom_point(alpha=.2) +
        geom_smooth(aes(NFA,NFA.err, colour=type) ) + #, se=FALSE ) +
        facet_wrap(~type, scales="free_x") +
        labs(x = "NFA", y = "log 10 NFA error")
    }
    
    if( stat == "err" ){
      return( list( ds=list( bias=l.coeffs$l.err$bias, aad=l.coeffs$l.err$aad,
                             mean.max.abs.err=l.coeffs$l.err$mean.max.abs.err ),
                    nl=list( bias=nl.bias, aad=nl.aad,
                             mean.max.abs.err=nl.mean.max.abs.err ) ) )
    }
  }
  
  ## 4. OPTION 2: COMPUTE THE REAL EXCHANGE RATE CORRELATION OR THE UIP REGRESSION ##
  if( stat %in% c( "bs", "uip" , "all") ){
    
    ## 4.1 SIMULATE ##
    n.sim <- 80000
        # Number of simulation points
    exog.sim <- sapply( 1:opt$n.exog, function(i) ar1_sim( n.sim, params$rho[i], 
                                                           params$sig.eps[i] ) )
        # Create a fresh exogenous simulation
    if( bo.nl ){
      nl.sim <- endog_sim( n.sim, exog.sim, sol.2$coeff, 2, opt$upper,
                           opt$lower, opt$endog.init, FALSE, 1, 0, TRUE )
      nl.cont.sim <- cont_sim( nl.sim, sol.2$coeff.cont, 2, opt$n.endog, opt$n.exog, 
                               sol.2$opt$upper, sol.2$opt$lower, FALSE )
      colnames(nl.cont.sim) <- cont.names
        # The nonlinear simulation
    }
    ds.sim <- endog_sim( n.sim, exog.sim, l.coeffs$ds.sol$coeff, 1, opt$upper,
                          opt$lower, opt$endog.init, FALSE, 1, 0, TRUE )
    ds.cont.sim <- cont_sim( ds.sim, l.coeffs$ds.sol$coeff.cont, 1, opt$n.endog, opt$n.exog, 
                              opt$upper, opt$lower, FALSE )
    colnames(ds.cont.sim) <- cont.names
        # The DS simulation
    
    browser()
    
    ## 4.2 BACKUS-SMITH ##
    if( stat %in% c( "bs", "all") ){
      
      message(" * Creating Backus-Smith stats *")
      bs.nl <- if(bo.nl) cor( diff(nl.cont.sim[,'E']), 
                    diff(nl.cont.sim[,'C1'] - nl.cont.sim[,'C2'] -
                      nl.cont.sim[,'Q'] / params$gamma ) ) else NA
      bs.ds <- cor( diff(ds.cont.sim[,'E']), 
                    diff(ds.cont.sim[,'C1'] - ds.cont.sim[,'C2'] -
                      ds.cont.sim[,'Q'] / params$gamma ) )
          # The correlations
      
      nfa.qt.nl <- if(bo.nl) quantile( nl.sim[,3], c(.25, .75 ), na.rm=T ) else NA
      nfa.qt.ds <- quantile( ds.sim[,3], c(.25, .75 ), na.rm=T )
          # The upper and lower quartile boundaries
      
      nfa.lower.idx.nl <- if(bo.nl) nl.sim[,3] < nfa.qt.nl[1] else NA
      nfa.upper.idx.nl <- if(bo.nl) nl.sim[,3] > nfa.qt.nl[2] else NA
      nfa.lower.idx.ds <- ds.sim[,3] < nfa.qt.ds[1]
      nfa.upper.idx.ds <- ds.sim[,3] > nfa.qt.ds[2]
          # Indices for upper and lower quartiles
      
      bs.upper.nl <- if(bo.nl) cor( diff(nl.cont.sim[nfa.upper.idx.nl,'E']), 
                          nl.cont.sim[nfa.upper.idx.nl,'C1'][-1] - 
                            nl.cont.sim[nfa.upper.idx.nl,'C2'][-1] +
                            diff(nl.cont.sim[nfa.upper.idx.nl,'Q']) 
                          / params$gamma ) else NA
      bs.upper.ds <- cor( diff(ds.cont.sim[nfa.upper.idx.ds,'E']), 
                          ds.cont.sim[nfa.upper.idx.ds,'C1'][-1] - 
                            ds.cont.sim[nfa.upper.idx.ds,'C2'][-1] +
                            diff(ds.cont.sim[nfa.upper.idx.ds,'Q']) / params$gamma )
      bs.lower.nl <- if(bo.nl) cor( diff(nl.cont.sim[nfa.lower.idx.nl,'E']), 
                          nl.cont.sim[nfa.lower.idx.nl,'C1'][-1] - 
                            nl.cont.sim[nfa.lower.idx.nl,'C2'][-1] +
                            diff(nl.cont.sim[nfa.lower.idx.nl,'Q']) 
                          / params$gamma ) else NA
      bs.lower.ds <- cor( diff(ds.cont.sim[nfa.lower.idx.ds,'E']), 
                          ds.cont.sim[nfa.lower.idx.ds,'C1'][-1] - 
                            ds.cont.sim[nfa.lower.idx.ds,'C2'][-1] +
                            diff(ds.cont.sim[nfa.lower.idx.ds,'Q']) / params$gamma )    
      ## Correlations for top and bottom quartiles
      
      if( stat == "bs" ){
        return( c( nl=bs.nl, ds=bs.ds, bs.upper.nl=bs.upper.nl, bs.upper.ds=bs.upper.ds,
                   bs.lower.nl=bs.lower.nl, bs.lower.ds=bs.lower.ds) )
      }
    }
    
#     browser()
    
    if( stat %in% c( "uip", "all")){
      ## 4.2 UIP REGRESSION ##
      
      message( " * Creating UIP stats * ")
            
      if(bo.nl){
        nl.R1 <- nl.cont.sim[,'rb1'][-1] + diff( nl.cont.sim[,'P1'] )
        nl.R2 <- nl.cont.sim[,'rb2'][-1] + diff( nl.cont.sim[,'P1'] ) - diff( nl.cont.sim[,'E'] )
            ## NEED TO ADD THE R BACK IN HERE SOMEHOW
        nl.rdiff <- nl.R1 - nl.R2
        nl.e.app <- diff(nl.cont.sim[,'E'])
        nl.coeffs <- tryCatch( lm(nl.e.app~nl.rdiff)$coeff , 
                               error=function(cond) return(c(NA) ) )
        nl.uip.err <- nl.e.app - nl.rdiff 
        nl.uip.nfa <- tryCatch( lm( nl.uip.err ~ nl.sim[,3][-nrow(nl.sim)])$coeff , 
                                error=function(cond) return(c(NA) ) )
      }else{
        nl.coeffs <- NA
      }
      ds.R1 <- ds.cont.sim[,'rb1'][-1] + diff( ds.cont.sim[,'P1'] )
      ds.R2 <- ds.cont.sim[,'rb2'][-1] + diff( ds.cont.sim[,'P1'] ) - diff( ds.cont.sim[,'E'] )
      ds.rdiff <- ds.R1 - ds.R2
      ds.e.app <- diff(ds.cont.sim[,'E'])
      ds.coeffs <- tryCatch( lm(ds.e.app~ds.rdiff)$coeff,
                            error=function(cond) return(c(NA) ) )
          # The UIP coefficients
      ds.uip.err <- ds.e.app - ds.rdiff 
      ds.uip.nfa <- tryCatch( lm( ds.uip.err ~ ds.sim[,3][-nrow(ds.sim)])$coeff , 
                              error=function(cond) return(c(NA) ) )
      if( stat == "uip" ){
        return( list( nl=nl.coeffs, ds=ds.coeffs ) )
      }
    }
    
    
    
    return( c( af1=l.coeffs$mod$alpha.tilde,
               ds.bias=l.coeffs$l.err$bias, ds.aad=l.coeffs$l.err$aad,
               ds.mean.max.abs.err=l.coeffs$l.err$mean.max.abs.err,
               nl.bias=nl.bias, nl.aad=nl.aad, nl.mean.max.abs.err=nl.mean.max.abs.err,
               ds.bs=bs.ds, nl.bs=bs.nl, 
               ds.uip=ds.coeffs[1], nl.uip=nl.coeffs[1] ) )

  }
  
  
  
  
}
                     
