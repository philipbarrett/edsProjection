#########################################################################
# reports.R
#
# Code to produce summary reports of the model
# Philip Barrett, Chicago
# 29feb2016
#########################################################################

## TO DO:
## 1. Create list with: endog sim, cont sim, error, expected conts, 
##     lead conts [DONE]
## 1a. Return the realized shocks for the exogenous variables.
## 2. Generic charting function [DONE]
## 3. Apply to a set of charts.  Save output.
## 4. Compute the correlations
## 5. Create latex tables of correlations
## 6. Create latex table of parameters
## 7. Create wrapper tex file and compile into a single latex report.
## 8. Wrap for a single solution

report.data <- function( sol ){
# Computes the list of data used in the report
  opt <- sol$opt
  params <- sol$params
      # Local copy of the options and parameters
  check <- if( sol$opt$model == 'irbc') sol.irbc.check( sol ) else sol.check( sol, sol$opt, sol$params )
      # Create the check list
  e.cont <- e_cont( sol$coeff.cont, check$endog.sim[,1:4], opt$n.exog, 
                    opt$n.endog, opt$n.cont, params$rho, 
                    params$sig.eps, 0, opt$N, opt$upper, opt$lower, 
                    opt$cheby, matrix(0,1,1), TRUE, 5 )
      # The expected continuations
  r.cont <- real_cont( sol$coeff.cont, check$endog.sim[,1:4], opt$n.exog, 
                        opt$n.endog, opt$n.cont, params$rho, 
                        params$sig.eps, opt$N, opt$upper, opt$lower, opt$cheby )
  out <- check
  out$r.cont <- r.cont$r.cont
  out$r.exog <- r.cont$r.exog
  out$e.cont <- e.cont
      # Set up the output
  return(out)
}

###### CHARTS #######
cht.generic <- function( X, Y, xlab, ylab, file, loc ){
  pdf( paste0( loc, '/', file, '.pdf' ) )
  plot( X, Y, xlab=xlab, ylab=ylab )
  dev.off()
}
  
cht.rep <- function( rep.data, loc=NULL ){
# Plots the report data
  
  if( is.null( loc ) ) loc <- '/home/philip/Dropbox/2016/Research/IRBC puzzles'
  
  ## Create real exchange rate ##
  rep.data$cont.sim <- cbind( rep.data$cont.sim, rep.data$cont.sim[,13] - 
                                rep.data$cont.sim[,9] + rep.data$cont.sim[,10] )
  rep.data$r.cont <- cbind( rep.data$r.cont, rep.data$r.cont[,13] - 
                            rep.data$r.cont[,9] + rep.data$r.cont[,10],
                            rep.data$cont.sim[,4] + rep.data$e.cont[,13] 
                              - rep.data$cont.sim[,13] )
      # Create real exchange rates and add the foreign one mult by expected
      # appreciation
  
  v.y.gth <- c( rep( T, 5 ), F, F, rep( T, 6 ) )
      # Indicator for growth for Y
  v.x.cont <- c( rep( F, 7 ), rep( T, 6 ) )
      # Whether the x index is a control or not
  v.x.idx <- c( rep( 1, 11 ), 13, 14 )
      # The indices
  v.x.lab <- c( rep( 'Income growth', 7 ), rep( 'Consumption growth', 4), 
                'Nom exchange rate growth', 'Real exchange rate growth' )
      # The x labels
  v.y.lab <- c( 'Domestic inflation', 'Foreign inflation', 'Consumption growth',
                'Nom exchange rate growth', 'Real exchange rate growth',
                'Log dom int rate', 'Log: (for int rate) x (ex rate gth)', 
                'Nom ex rate gth', 'Real ex rate gth', 'Dom inflation',
                'For inflation', 'Dom inflation', 'For inflation' )
      # The y labels
  v.file <- paste0( 'chart', 1:13 )
      # The vector of file names
  for( i in 1:13 ){
    if( v.x.cont ){
      X <- rep.data$cont.sim[, v.x.idx[i]]
    }else{
      X <- rep.data$cont.sim[, v.x.idx[i]]
    }
      
  }
    
  
}


