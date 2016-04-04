#########################################################################
# reports.R
#
# Code to produce summary reports of the model
# Philip Barrett, Chicago
# 29feb2016
#########################################################################

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
  plot( X, Y, xlab=xlab, ylab=ylab, pch=19, col= alpha('black', 0.15) )
  dev.off()
  return( c( cov(X,Y), cor(X,Y), lm(Y~X)$coefficients[-1] ) )
}
  
report.corr <- function( rep.data, loc=NULL ){
# Plots the report data
  
  if( is.null( loc ) ) loc <- '~/Dropbox/2016/Research/IRBC puzzles/'
  loc.chts <- paste0( loc, 'charts/' )
  
  ## Create real exchange rate ##
  rep.data$cont.sim <- cbind( rep.data$cont.sim, 
                        rep.data$cont.sim[,13] - rep.data$cont.sim[,9] + rep.data$cont.sim[,10],
                        rep.data$cont.sim[,4] + rep.data$e.cont[,13] - rep.data$cont.sim[,13])
  rep.data$r.cont <- cbind( rep.data$r.cont, 
                        rep.data$r.cont[,13] - rep.data$r.cont[,9] + rep.data$r.cont[,10] )
      # Create real exchange rates and add the foreign one mult by expected
      # appreciation.
  
  v.x.gth <- c( rep( T, 13 ), rep( F, 3 ) )
  v.y.gth <- c( rep( T, 5 ), F, F, rep( T, 6 ), rep(F,3) )
      # Indicator for growth for X & Y
  v.x.cont <- c( rep( F, 7 ), rep( T, 6 ), F, T, T )
      # Whether the x index is a control or not
  v.x.idx <- c( rep( 1, 11 ), 13, 14, rep(1,3) )
  v.y.idx <- c( 9, 10, 1, 13, 14, 3, 15, 13, 14, 9, 10, 9, 10, 1, 13, 14 )
      # The indices
  v.x.lab <- c( rep( 'Technology growth', 7 ), rep( 'Consumption growth', 4), 
                'Nom exchange rate growth', 'Real exchange rate growth',
                'Technology', rep( 'Consumption', 2 ) )
      # The x labels
  v.y.lab <- c( 'Domestic inflation', 'Foreign inflation', 'Consumption growth',
                'Nom ex rate growth', 'Real ex rate growth',
                'Log dom int rate', 'Log: (for int) * (NER gth)', 
                'Nom ex rate gth', 'Real ex rate gth', 'Dom inflation',
                'For inflation', 'Dom inflation', 'For inflation',
                'Consumption', 'Nom ex rate', 'Real ex rate' )
      # The y labels
  v.file <- paste0( 'chart', 1:16 )
      # The vector of file names
  out.df <- data.frame( x = v.x.lab, y = v.y.lab, cov=NA, cor=NA, reg=NA )
      # Empty dataframe
  for( i in 1:16 ){
    if( v.x.gth[i] ){
      if( v.x.cont[i] ){
        X <- rep.data$r.cont[, v.x.idx[i]] - rep.data$cont.sim[, v.x.idx[i]]
      }else{
        X <- rep.data$r.exog[, v.x.idx[i]] - rep.data$endog.sim[, v.x.idx[i]]
      }
    }else{
      if( v.x.cont[i] ){
        X <- rep.data$cont.sim[, v.x.idx[i]]
      }else{
        X <- rep.data$endog.sim[, v.x.idx[i]]
      }
    }
    
    if( v.y.gth[i] ){
      Y <- rep.data$r.cont[, v.y.idx[i]] - rep.data$cont.sim[, v.y.idx[i]]
    }else{
      Y <- rep.data$cont.sim[, v.y.idx[i]]
    }
    out.df[i,3:5] <- cht.generic( X, Y, v.x.lab[i], v.y.lab[i], v.file[i], loc.chts )
  }
  
  write( print(xtable(out.df, digits=c( NA, NA,NA, 6, 2, 2 ), 
                      caption='Main cross-correlations of the model'), 
               include.rownames=F), file = paste0( loc, 'reports/correl.tex' ) )
  
  # UIP Regression here
#   q <- rep.data$r.cont[,14]
  e.e <- rep.data$e.cont[,13] - rep.data$cont.sim[,13]
  r.diff <- rep.data$cont.sim[,3] - rep.data$cont.sim[,4]
#   write( print(xtable(summary( lm( q~r.diff ) ), digits=3, caption='Real UIP regression') ), 
#          file = paste0( loc, 'reports/UIP.tex' ) )
  write( print(xtable(summary( lm( e.e~r.diff ) ), digits=3, caption='Nominal UIP regression')), 
         file = paste0( loc, 'reports/UIP_nom.tex' ) )
  pdf(paste0( loc, 'charts/uip.pdf') )
  plot( r.diff, e.e, xlab='Domestic bond premium', 
        ylab='Expected exchange rate growth', pch=19, col= alpha('black', 0.15) )
  abline( 0, 1, lty=2, col='blue' )
  dev.off()
  
  # Chart for home/total assets + a 45 degree line
  pdf(paste0( loc, 'charts/assets.pdf') )
  tot.assets <- rep.data$endog.sim[,3] -
                    rep.data$endog.sim[,4] * exp( rep.data$cont.sim[,13] )
  dom.assets <- rep.data$endog.sim[,3] # * exp( - rep.data$cont.sim[,3] )
  plot( tot.assets, dom.assets, xlab='Total assets', 
        ylab='Domestically held assets', pch=19, col= alpha('black', 0.15) )
  abline( 0, .5, lty=2, col='blue' )
  fit <- lm( dom.assets ~ tot.assets )
  abline( fit, lty=2, col='red' )
  legend( 'topleft', c('22.5 degree line', 'Regression line'),
          col=c('blue', 'red'), lty=2, bty='n' )
  dev.off()
  write( print( xtable( summary( fit ), digits=3, 
                        caption = 'Regression of domestic on total assets' ) ), 
         file = paste0( loc, 'reports/assets.tex' ) )
  
  asset.ratio <- rep.data$endog.sim[,3] / ( rep.data$endog.sim[,3] - rep.data$endog.sim[,4] * exp( rep.data$cont.sim[,13] ) )
  asset.diff <- rep.data$endog.sim[,3] + rep.data$endog.sim[,4] * exp( rep.data$cont.sim[,13] )
  write( print( xtable( summaryfunction( asset.ratio[abs(asset.ratio)<1000] ), digits=3, 
                      caption = 'Ratio of home to total assets' ) ), 
       file = paste0( loc, 'reports/assetRatio.tex' ) )
      # Drop the outliers for the asset ratio
  write( print( xtable( summaryfunction( asset.diff ), digits=3, 
                      caption = 'Difference between home and foreign assets' ) ), 
       file = paste0( loc, 'reports/assetDiff.tex' ) )
  write( print( xtable( summaryfunction( tot.assets ), digits=3, 
                      caption = 'Total assets' ) ), 
       file = paste0( loc, 'reports/assetTot.tex' ) )

  # Asset densities
  pdf(paste0( loc, 'charts/debt_dist.pdf') )
  plot( density( rep.data$endog.sim[,3] ), col='red', xlab='Domestic assets', ylab='Density' )
  lines( density( rep.data$endog.sim[,4] ),col='blue' )
  abline( v=mean(rep.data$endog.sim[,3]), lty=2, col='red' )
  abline( v=mean(rep.data$endog.sim[,4]), lty=2, col='blue' )
  legend( 'topright', 
          c('Country 1','Country 2','Country 1 mean','Country 2 mean'), 
          lty=c(1,1,2,2), col=c('red','blue','red','blue') )
  dev.off()
  
  # The log10 errors
  pdf(paste0( loc, 'charts/err.pdf') )
  plot( 1:ncol(rep.data$err), log( apply( abs( rep.data$err ), 2, mean ), 10 ),
        xlab='Equation number', ylab='Log(10) mean absolute error' )
  dev.off()

  # The relationship between the real exchange rate and consumptio n differentials
  pdf(paste0( loc, 'charts/cons_diff.pdf') )
  cons.diff <- ( rep.data$r.cont[, 3] - rep.data$cont.sim[, 3] ) - 
                    ( rep.data$r.cont[, 4] - rep.data$cont.sim[, 4] )
  q.gth <- rep.data$r.cont[, 14] - rep.data$cont.sim[, 14]
  plot( cons.diff, q.gth, xlab='International consumption differential growth', 
        ylab='Real exchange rate growth', pch=19, col= alpha('black', 0.15) )
  fit.cons.diff <- lm( q.gth ~ cons.diff )
  abline( fit.cons.diff, lty=2, col='red' )
  dev.off()
  write( print( xtable( summary( fit.cons.diff ), digits=3, 
          caption = 'Regression of real exchange rate on international consumption differential' ) ), 
          file = paste0( loc, 'reports/cons_diff.tex' ) )
  
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
  if( is.null( loc ) ) loc <- '~/Dropbox/2016/Research/IRBC puzzles/'
  out.file <- paste0( loc, 'reports/', st.time, '.tex' )
      # The output file
  
  ## Header
  write( '\\documentclass[12pt]{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage{graphicx,ctable,booktabs}', file=out.file )
  write( paste0( '\\title{Model solution Report: ', st.time, '}\n' ), file=out.file, append=T )  
  write( '\\begin{document}\n', file=out.file, append=T )  
  write( '\\maketitle\n', file=out.file, append=T )  
  
  
  ## Parameters table
  write( '\\begin{table}[htb]\n\\centering\n\\begin{tabular}{cccccccccc}', file=out.file, append=T )
  write( 'Share & $\\hat\\alpha$ & $\\beta$ & $\\gamma$ & $\\eta$ & $\\rho$ & $\\sigma_\\epsilon$ & $\\bar{p}_1$ & $\\bar{p}_2$ & $N$ \\\\', 
         file=out.file, append=T )
  write( '\\hline', file=out.file, append=T )
  write( paste0( sol$params$share, ' & ', round( sol$params$alphahat, 3) , ' & ', sol$params$betta , ' & ', 
                 sol$params$gamma , ' & ',sol$params$eta , ' & ', sol$params$rho[1] , ' & ', 
                 sol$params$sig.eps[1] , ' & ', sol$params$P1.bar , ' & ', sol$params$P2.bar, ' & ',
                 sol$opt$N ), 
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
  write( '\\end{tabular} \n\\caption{Model standard deviations} \n\\end{table} \n\n', file=out.file, append=T )
  
  ## Consumption differential
#   write( '\\input{UIP.tex}', file=out.file, append=T )
  report.chart.latex('../charts/cons_diff.pdf', out.file)
  write( '\\input{cons_diff.tex}', file=out.file, append=T )
  write( '\n\\clearpage\n', file=out.file, append=T )

  ## Law of one price
  loop <- log( sol$params$P1.bar ) - rep.data$cont.sim[,11 ] + 
    log( sol$params$P2.bar ) - rep.data$cont.sim[,12 ]
  write( print(xtable(summaryfunction( loop ), digits=6, 
                      caption='Check on the law of one price (should all be zero)'), 
               include.rownames=F), file = paste0( loc, 'reports/loop.tex' ) )
  write( '\\input{loop.tex}', file=out.file, append=T )

  ## UIP
  #   write( '\\input{UIP.tex}', file=out.file, append=T )
  report.chart.latex('../charts/uip.pdf', out.file)
  write( '\\input{UIP_nom.tex}', file=out.file, append=T )


#   write( '\n\\clearpage\n', file=out.file, append=T )

  ## Charts
  report.chart.latex('../charts/err.pdf', out.file)
  report.chart.latex('../charts/assets.pdf', out.file)
  write( '\\input{assets.tex}', file=out.file, append=T )
  write( '\\input{assetRatio.tex}', file=out.file, append=T )
  write( '\\input{assetDiff.tex}', file=out.file, append=T )
  write( '\\input{assetTot.tex}', file=out.file, append=T )
  report.chart.latex('../charts/debt_dist.pdf', out.file)
  write( '\n\\clearpage\n', file=out.file, append=T )
  for( i in 1:16) report.chart.latex( paste0('../charts/chart', i, '.pdf'), out.file)
    
  ## Footer
  write( '\\end{document}', file=out.file, append=T )
  
  setwd( paste0( loc, 'reports/' ) )
  system( paste0( 'pdflatex ', st.time, '.tex' ) )
}

summaryfunction= function (x){
  if( is.numeric(x)!=TRUE) {stop("Supplied X is not numeric")}
  mysummary = data.frame(
    "Min." =as.numeric( min(x)),
    "1st Qu." = quantile(x)[2],
    "Median" = median(x),
    "Mean" = mean(x),
    "3rd Qu." = quantile(x)[4],
    "Max." = max(x),
    row.names=""
    
  )
  names(mysummary) = c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")
  return( mysummary )
}
