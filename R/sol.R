#########################################################################
# sol.R
#
# Generic solution algorithm for the EDS-projection algorithm
# Philip Barrett, Chicago
# 02feb2016
#########################################################################

err.min <- function( coeffs.init, X, model, lags, params, n.exog, 
                     n.endog, rho, sig.eps, n.integ, N, upper, lower, cheby,
                     exog.innov.mc, quad, n.nodes ){
## Minimizes the error on the model equations given a grid X
  
  browser()
  
  n.terms <- idx_count( N, n.exog + n.endog )
      # Number of terms in the polynomial approximation
  f <- function(x) eval_err( matrix( x, nrow=n.terms, ncol=n.endog ), 
                             X, model, lags, params, n.exog, n.endog,
                             rho, sig.eps, n.integ, N, upper, lower, cheby,
                             exog.innov.mc, quad, n.nodes )
      # The error evaluated on the restricted grid X
  f.grad <- function(x) eval_err_D( matrix( x, nrow=n.terms, ncol=n.endog ), 
                                    X, model, lags, params, n.exog, n.endog,
                                    rho, sig.eps, n.integ, N, upper, lower, 
                                    cheby,exog.innov.mc, quad, n.nodes )
      # The gradient of the error
  
  opts <- list("algorithm"="NLOPT_LD_LBFGS",
               "xtol_rel"=1.0e-8, print_level=1 )
  
  sol <- nloptr( x0=coeffs.init, eval_f = f, eval_grad_f = f.grad, opts=opts )
  
  
  sol <- ipoptr( x0=coeffs.init, eval_f = f, eval_grad_f = f.grad )
  
  browser()
  
  return(sol)
  
}
# 
# 
# eval_err( arma::mat  )
# 
# eval_err( arma::mat coeffs, arma::mat X, std::string model, 
#           int lags, List params, int n_exog, int n_endog,
#           arma::rowvec rho, arma::rowvec sig_eps, int n_integ,
#           int N, arma::rowvec upper, arma::rowvec lower, bool cheby,
#           arma::mat exog_innov_mc, bool quad=true, int n_nodes=0 )