#########################################################################
# plotting.R
#
# Plotting for the EDS-projection algorithm
# Philip Barrett, Chicago
# 11feb2016
#########################################################################

plot.coeffs <- function( coeff, ... ){
  barplot(t(coeff), beside=T, col=1:ncol(coeff), ... )
}