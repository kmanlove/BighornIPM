TuljarFun2 <- function(alpha, gamma, sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost){
  Tuljar1.out <- TuljarFun1(alpha, gamma)
  
  L.healthy <- UpdateLeslieFun(current.state = "healthy", sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
  L.sp <- UpdateLeslieFun(current.state = "spillover", sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
  L.infected <- UpdateLeslieFun(current.state = "infected", sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
  
  b <- L.healthy * Tuljar1.out$w.naught[1] + L.sp * Tuljar1.out$w.naught[2] + L.infected * Tuljar1.out$w.naught[3]
  
  b.eigen <- eigen(b)
  lambda.naught <- Re(b.eigen$value[1])
  u.naught <- Re(b.eigen$vector[ , 1])
  v.naught <- Re(Conj(ginv(b.eigen$vector))[1, ])
  check <- t(v.naught) %*% u.naught
  
  return(list(L.healthy = L.healthy, L.sp = L.sp, L.infected = L.infected, b = b, lambda.naught = lambda.naught, u.naught = u.naught, v.naught = v.naught))
}
