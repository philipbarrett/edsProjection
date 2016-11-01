params <- list( share = .86, gamma = 2, P1.bar=1, P2.bar=1, betta=.95,
                rho=c(.9,.9, .9, .9), sig.eps=c(.01,.01, .0025, .0025), eta=2, 
                theta=.025 )

### 0. Hard-coding ###
N <- 1
n.lag <- 1
n.fwd <- 4
n.sample <- 4000
nsim <- 1e6
burn <- 1e4

### 1. Run and extract the model ###
mod.create(params)
mod.run()
mod <- mod.read()
y1.idx <- which( unlist( mod$varo ) == "Y1 " ) - 1
    # The location of Y1
sim <- stoch_simDS( mod, params$sig.eps, params$betta, nsim, burn, y1.idx )
    # The simulation
varo <- c( sub("\\s+$", "", (unlist(mod$varo))), 'af1' )
ys <- c( mod$ys, mod$alpha.tilde )
names(ys) <- varo
colnames(sim) <- varo
    # Assign names
B11 <- exp(mod$ys[y1.idx+1]) * params$betta * exp( -params$theta * sim[,'C1'] ) *
          sim[,'af1'] * exp( sim[,'R_1'] + sim[,'P1'] )
B22 <- ( sim[,'NFA'] - params$betta * exp( -params$theta * sim[,'C1'] ) * 
           sim[,'af1'] ) * exp(mod$ys[y1.idx+1] ) * exp( sim[,'R_2'] + sim[,'P1']  - sim[,'E'])
    # Create nominal debt series:  NB THESE ARE END OF PERIOD
sim <- cbind( sim, B11, B22 )
    # Append B11, B22 to simulation
varo <- c( varo, 'B11', 'B22' )
ys <- c( ys, B11=mean(B11), B22=mean(B22) )
    # Assign names

exog.order <- c('A1','A2', 'P11', 'P22')
endog.order <- c( 'B11', 'B22' )
cont.order <- c( 'C1', 'C2', 'rb1', 'rb2', 'X11', 'X22', 'X12', 'X21', 
                 'P1', 'P2', 'R_1', 'R_2',
                 'P12', 'P21', 'E', 'Q',
                 'Y1', 'Y2', 'cd', 'cg', 'NFA', 'af1' )
    # Variable names for the DS-style solution
fwd.order <- c( 'rb1', 'rb2', 'B11', 'B22')
    # The forward-looking equation variables
sim.exog <- sim[, exog.order]
sim.endog <- sim[, endog.order]
sim.cont <- sim[, cont.order]
    # Separate out the exogenous and endogenous states and the static variables
n.exog <- length( exog.order )
n.endog <- length( endog.order )
n.cont <- length( cont.order )
    # Numbers of types of variables

sd.x <- params$sig.eps / sqrt( ( 1 - params$rho ^ 2 ) )
    # The standard deviation of the exogenous state
upper <- c(   3 * sd.x, ys[endog.order] + 3 * apply(sim.endog,2,sd) )
lower <- c( - 3 * sd.x, ys[endog.order] - 3 * apply(sim.endog,2,sd) )
    # The bounds of the states
endog.reg <- sim.endog[-nsim,]
exog.reg <- sim.exog[-1,]
X <- cbind( exog.reg, endog.reg )
    # The X-variables for the regressions
n.X <- 1 + length(endog.order) + length(exog.order) 

coeff <- matrix( 0, n.X, n.endog )
coeff.cont <- matrix( 0, n.X, n.cont )
    # The coefficient matrices  
cheby <- FALSE
X.reg <- coeff_reg_X( X, N, lower, upper, cheby )
for(i in 1:n.endog) coeff[,i] <- coeff_reg( sim.endog[,i][-1], X, N,
                                            lower, upper, cheby )
    # Populate the coefficient matrices (Endogenous states only)

sim.endog.alt <- endog_sim( nsim-1, sim.exog[-1,], coeff, N, upper, lower, 
                            sim.endog[1,], cheby, 1, 0, TRUE )
    # Alternative endogenous simulation




