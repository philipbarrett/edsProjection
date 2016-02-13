/***********************************************************************************
 * ngmCont2.cpp
 * 
 * Code to compute the errors on the 2-country neoclassical growth model with controls
 * 
 * 13feb2016
 * Philip Barrett, Chicago
 * 
 ***********************************************************************************/

#include "ngmCont2.hpp"
#include "sim.hpp"

// [[Rcpp::export]]
arma::rowvec integrand_ngm_cont_2( 
               arma::mat exog, arma::mat endog, arma::rowvec cont, 
               arma::rowvec exog_lead, List params, arma::mat coeffs, 
               arma::mat coeffs_cont, int n_exog, int n_endog, int n_cont, 
               int N, arma::rowvec upper, arma::rowvec lower, bool cheby=false ){

  double A = params["A"] ;
  double alpha = params["alpha"] ;
  double delta = params["delta"] ;
  double gamma = params["gamma"] ;
  
  rowvec cont_lead = endog_update( exog_lead, endog.row(0), coeffs_cont, n_exog, 
                                    n_endog, N, upper, lower, cheby ) ;
      // Create the next-period control.  Remember, the current-period state is
      // an end-of-period variable, so the control depends directly on its lag.

  rowvec integrand(2) ;
  integrand(0) = pow( cont_lead(0) / cont(0), - gamma ) * 
          ( 1 - delta + exp( exog_lead( 0 ) ) * alpha * A * pow( endog( 0, 0 ), alpha - 1 ) ) ;
  integrand(1) = pow( cont_lead(0) / cont(0), - gamma ) * 
          ( 1 - delta + exp( exog_lead( 1 ) ) * alpha * A * pow( endog( 0, 1 ), alpha - 1 ) ) ;
      // Two countries have same consumption but diefferent capital stocks and
      // technologies
  return integrand ;
}



