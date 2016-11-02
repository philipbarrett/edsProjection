params <- list( share = .86, gamma = 2, P1.bar=1, P2.bar=1, betta=.95,
                rho=c(.9,.9, .9, .9), sig.eps=c(.01,.01, .0025, .0025), eta=2, 
                theta=.025 )

baseline <- mod.gen(params, check=FALSE)
