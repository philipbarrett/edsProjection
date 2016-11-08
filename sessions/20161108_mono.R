########## INFORMAL TEST OF MONOMIAL INTEGRATION ##########

# Create coefficients from dynare
params <- list( share = .86, gamma = 2, P1.bar=1, P2.bar=1, betta=.95,
                rho=c(.9,.9, .9, .9), sig.eps=c(.01,.01, .0025, .0025), eta=2, 
                theta=.025 )
baseline <- mod.gen(params, check=FALSE)

# Some definitions
upper <- baseline$ds.sol$upper
lower <- baseline$ds.sol$lower
n.endog <- ncol(baseline$ds.sol$coeff)
n.exog <- nrow(baseline$ds.sol$coeff) - n.endog - 1
n.cont <- ncol(baseline$ds.sol$coeff.cont)
n.fwd <- n.exog
endog.init <- tail(baseline$ds.sol$ys, n.endog)
coeff <- baseline$ds.sol$coeff
coeff.cont <- baseline$ds.sol$coeff.cont
    # From Dynare
exog.names <- c('A1','A2', 'P11', 'P22')
endog.names <- c( 'B11', 'B22' )
cont.names <- c( 'C1', 'C2', 'R_1', 'R_2', 'X11', 'X22', 'X12', 'X21', 
                 'P1', 'P2', 'P12', 'P21', 'E', 'Q' )
fwd.vars <- c('B11', 'E', 'R_1', 'R_2')
    # For Euler equations evaluation
rho <- params$rho
sig.eps <- params$sig.eps
    # From parameters
n.sim <- 10000
kappa <- 1
burn <- 0
lags <- 1
N <- 1
cheby <- FALSE
endog.init <- baseline$ds.sol$ys[c('B11','B22')]
    # Fresh

# Create the simulation
exog.sim <- sapply( 1:n.exog, function(i) ar1_sim( n.sim * kappa + burn, 
                                                   rho[i], sig.eps[i] ) )
colnames(exog.sim) <- exog.names
    # Exogenous simulation
endog.sim <- endog_sim( n.sim, exog.sim, coeff, N, upper, lower, 
                        endog.init, cheby, kappa, burn, (lags>0) )
colnames(endog.sim) <- c( c( exog.names, endog.names ), 
                          paste0(c( exog.names, endog.names), '(-1)' ) )
    # Endogenous simulation
cont.sim <- cont_sim( endog.sim, coeff.cont, N, n.endog, n.exog, upper, lower, cheby )
colnames(cont.sim) <- cont.names
    # The control simulation
X.cont <- cbind( endog.sim, cont.sim )
    # The full simulation grid

# Evaluation the Euler equation errors
system.time(
  k.hat.quad.1 <- euler_hat_grid( coeff, coeff.cont, X.cont, lags, params, 
                                 n.exog, n.endog, n.cont, n.fwd, rho, sig.eps, 
                                 0, N, upper, lower, cheby, matrix(0,1,1,), 
                                 TRUE, 1, "irbc" ) )
    # .4 sec
system.time(
  k.hat.quad.2 <- euler_hat_grid( coeff, coeff.cont, X.cont, lags, params, 
                                n.exog, n.endog, n.cont, n.fwd, rho, sig.eps, 
                                0, N, upper, lower, cheby, matrix(0,1,1,), 
                                TRUE, 2, "irbc" ) )
    # 5.6 sec
system.time(
  k.hat.quad.3 <- euler_hat_grid( coeff, coeff.cont, X.cont, lags, params, 
                                  n.exog, n.endog, n.cont, n.fwd, rho, sig.eps, 
                                  0, N, upper, lower, cheby, matrix(0,1,1,), 
                                  TRUE, 3, "irbc" ) )
    # 32 sec
system.time(
  k.hat.quad.m1 <- euler_hat_grid( coeff, coeff.cont, X.cont, lags, params, 
                                  n.exog, n.endog, n.cont, n.fwd, rho, sig.eps, 
                                  0, N, upper, lower, cheby, matrix(0,1,1,), 
                                  TRUE, 2, "irbc", "m1" ) )
    # 5.7 sec

