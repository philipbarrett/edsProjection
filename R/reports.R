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
  pdf( paste0( loc, file, '.pdf' ) )
  plot( X, Y, xlab=xlab, ylab=ylab )
  dev.off()
  return( c( cov(X,Y), cor(X,Y), lm(Y~X)$coefficients['X'] ) )
}
  
report.corr <- function( rep.data, loc=NULL ){
# Plots the report data
  
  if( is.null( loc ) ) loc <- '/home/philip/Dropbox/2016/Research/IRBC puzzles/'
  loc.chts <- paste0( loc, 'charts/' )
  
  ## Create real exchange rate ##
  rep.data$cont.sim <- cbind( rep.data$cont.sim, rep.data$cont.sim[,13] - 
                                rep.data$cont.sim[,9] + rep.data$cont.sim[,10] )
  rep.data$r.cont <- cbind( rep.data$r.cont, 
                        rep.data$r.cont[,13] - rep.data$r.cont[,9] + rep.data$r.cont[,10],
                        rep.data$cont.sim[,4] + rep.data$e.cont[,13] - rep.data$cont.sim[,13] )
      # Create real exchange rates and add the foreign one mult by expected
      # appreciation.  Also total assets
  
  v.y.gth <- c( rep( T, 5 ), F, F, rep( T, 6 ) )
      # Indicator for growth for Y
  v.x.cont <- c( rep( F, 7 ), rep( T, 6 ) )
      # Whether the x index is a control or not
  v.x.idx <- c( rep( 1, 11 ), 13, 14 )
  v.y.idx <- c( 9, 10, 1, 13, 14, 4, 15, 13, 14, 9, 10, 9, 10 )
      # The indices      
  v.x.lab <- c( rep( 'Income growth', 7 ), rep( 'Consumption growth', 4), 
                'Nom exchange rate growth', 'Real exchange rate growth' )
      # The x labels
  v.y.lab <- c( 'Domestic inflation', 'Foreign inflation', 'Consumption growth',
                'Nom ex rate growth', 'Real ex rate growth',
                'Log dom int rate', 'Log: (for int) * (NER gth)', 
                'Nom ex rate gth', 'Real ex rate gth', 'Dom inflation',
                'For inflation', 'Dom inflation', 'For inflation' )
      # The y labels
  v.file <- paste0( 'chart', 1:13 )
      # The vector of file names
  out.df <- data.frame( x = v.x.lab, y = v.y.lab, cov=NA, cor=NA, reg=NA )
      # Empty dataframe
  for( i in 1:13 ){
    if( v.x.cont[i] ){
      X <- rep.data$r.cont[, v.x.idx[i]] - rep.data$cont.sim[, v.x.idx[i]]
    }else{
      X <- rep.data$r.exog[, v.x.idx[i]] - rep.data$endog.sim[, v.x.idx[i]]
    }
    if( v.y.gth[i] ){
      Y <- rep.data$r.cont[, v.y.idx[i]] - rep.data$cont.sim[, v.y.idx[i]]
    }else{
      Y <- rep.data$r.cont[, v.y.idx[i]]
    }
    out.df[i,3:5] <- cht.generic( X, Y, v.x.lab[i], v.y.lab[i], v.file[i], loc.chts )
  }
  
  write( print(xtable(out.df, digits=c( NA, NA,NA, 6, 2, 2)), include.rownames=F), file = paste0( loc, 'reports/correl.tex' ) )
  
  # Add here the chart for home/total assets + a 45 degree line
  pdf(paste0( loc, 'charts/assets.pdf') )
  plot( rep.data$endog.sim[,3] - rep.data$endog.sim[,4] * exp( rep.data$cont.sim[,13] ),
        rep.data$endog.sim[,3], xlab='Total assets', ylab='Domestically held assets' )
  abline( 0, .5, lty=2, col='blue' )
  dev.off()
  
  # Add the log10 errors here
  pdf(paste0( loc, 'charts/err.pdf') )
  plot( 1:ncol(rep.data$err), log( apply( abs( rep.data$err ), 2, mean ), 10 ),
        xlab='Equation number', ylab='Log(10) mean absolute error' )
  dev.off()
  
  return( out.df )  
}

###### CREATE THE REPORT ######
report.chart.latex <- function( cht.file, out.file ){
# Creaetes the latex code for including charts
  write( '\\begin{figure}[phtb]\n\\centering', file=out.file, append=T )
  write( paste0('\\includegraphics[width=3.5in]{', cht.file, '}'), file=out.file, append=T )
  write( '\\end{figure}', file=out.file, append=T )
}


report.create <- function( sol, rep.data=NULL, loc=NULL ){
# Create the report
  
  if( is.null( rep.data ) ) rep.data <- report.data( sol )
  rep.corr <- report.corr( rep.data )
      # The data for the charts
  st.time <- gsub( ':', '', gsub( ' ', '-', Sys.time() ) )
      # The timestamp
  if( is.null( loc ) ) loc <- '/home/philip/Dropbox/2016/Research/IRBC puzzles/'
  out.file <- paste0( loc, 'reports/', st.time, '.tex' )
      # The output file
  
  ## Header
  write( '\\documentclass[12pt]{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage{graphicx,ctable,booktabs}', file=out.file )
  write( paste0( '\\title{Model solution Report: ', st.time, '}\n' ), file=out.file, append=T )  
  write( '\\begin{document}\n', file=out.file, append=T )  
  write( '\\maketitle\n', file=out.file, append=T )  
  
  
  ## Parameters table
  write( '\\begin{table}[htb]\n\\centering\n\\begin{tabular}{cccccccc}', file=out.file, append=T )
  write( '$\\alpha$ & $\\beta$ & $\\gamma$ & $\\eta$ & $\\rho$ & $\\sigma_\\epsilon$ & $\\bar{p}_1$ & $\\bar{p}_2$ \\\\', 
         file=out.file, append=T )
  write( '\\hline', file=out.file, append=T )
  write( paste0( sol$params$alpha , ' & ', sol$params$betta , ' & ', sol$params$gamma , ' & ',
                 sol$params$eta , ' & ', sol$params$rho[1] , ' & ', sol$params$sig.eps[1] , ' & ',
                 sol$params$P1.bar , ' & ', sol$params$P2.bar ), 
         file=out.file, append=T )
  write( '\\end{tabular} \n\\caption{Model parameters} \n\\end{table} \n\n', file=out.file, append=T )

  ## Correlation tables
  write( '\\input{correl.tex}', file=out.file, append=T )
  
  ## Variance table
  v.sd <- sapply( c(9, 10, 13), function(i) sd( rep.data$r.cont[, i] - rep.data$cont.sim[, i] ) )
  v.sd[4] <- sd( rep.data$r.cont[,13] - rep.data$r.cont[,9] + rep.data$r.cont[,10] - 
                   ( rep.data$cont.sim[,13] - rep.data$cont.sim[,9] + rep.data$cont.sim[,10] ) )
  v.sd <- round( v.sd, 4 )
      # Vector of control standard deviations
  write( '\\begin{table}[htb]\n\\centering\n\\begin{tabular}{ccccc}', file=out.file, append=T )
  write( ' & Domestic inflation & Foreign inflation & NER growth  & RER growth \\\\', 
         file=out.file, append=T )
  write( '\\hline', file=out.file, append=T )
  write( paste0( 'Sample std dev & ', v.sd[1], ' & ', v.sd[2], ' & ', v.sd[3], ' & ', 
                 v.sd[4], ' \\\\ '), 
         file=out.file, append=T )
  write( '\\end{tabular} \n\\caption{Model parameters} \n\\end{table} \n\n', file=out.file, append=T )
  
  ## Euler eq decomp
  write( '\n***Euler equation decomposition to go here***\n', file=out.file, append=T)
  
  ## Charts
  report.chart.latex('../charts/err.pdf', out.file)
  for( i in 1:13) report.chart.latex( paste0('../charts/chart', i, '.pdf'), out.file)
  report.chart.latex('../charts/assets.pdf', out.file)
  
  ## Footer
  write( '\\end{document}', file=out.file, append=T )
  
  setwd( paste0( loc, 'reports/' ) )
  system( paste0( 'pdflatex ', st.time, '.tex' ) )
}


