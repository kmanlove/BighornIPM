TuljarFun6 <- function(n.states, n.environs, Tuljar1.out, Tuljar2.out, Tuljar3.out){
  #  Tuljar1.out <- TuljarFun1(alpha, gamma)
  #  Tuljar2.out <- TuljarFun2(alpha, gamma, sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
  #  Tuljar3.out <- TuljarFun3(Tuljar2.out)
  I.s <- diag(1, n.states)
  I.sk <- diag(1, n.states * n.environs)
  
  out <- ginv(I.sk - kronecker(Tuljar1.out$P.mat, Tuljar3.out$D))
  
  return(list(I.s = I.s, I.sk = I.sk, out = out))
}