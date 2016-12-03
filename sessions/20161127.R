## Parameters
rho <- diag(rep(.99,4))
sig <- diag(c(.01,.01,.005,.005)^2)
params <- list( share = .9, gamma = 2, P1.bar=1, P2.bar=1, betta=.99,
                rho=rho, sig.eps=sig, eta=3, theta=0.01 )

v.rho <- seq(.95,.99, length.out=5)
ll <- list()
for( i in 1:length(v.rho) ){
  params$rho <- diag(c(rep(v.rho[i],4)))
  ll[[i]] <- mod.eval(params)
  ll[[i]]$params <- params
}

v.gamma <- seq(2, 8, length.out=7)
kk <- list()
params$rho <- diag(.98,4)
for( i in 1:length(v.gamma) ){
  params$gamma <- v.gamma[i]
  kk[[i]] <- mod.eval(params)
  kk[[i]]$params <- params
}

v.share <- seq(.8, .95, length.out=7)
mm <- list()
params$gamma <- 2
for( i in 1:length(v.share) ){
  params$share <- v.share[i]
  mm[[i]] <- mod.eval(params, opt.lim = list(n.gain=0.25, iter=5))
  mm[[i]]$params <- params
}

v.eta <- seq(2, 10, length.out=9)
nn <- list()
params$share <- .9
for( i in 1:length(v.eta) ){
  params$eta <- v.eta[i]
  nn[[i]] <- mod.eval(params, opt.lim = list(n.gain=0.25, iter=5))
}

v.beta <- seq(.9, .99, length.out=10)
oo <- list()
params$eta <- 3
for( i in 1:length(v.betta) ){
  params$betta <- v.betta[i]
  oo[[i]] <- mod.eval(params)
}