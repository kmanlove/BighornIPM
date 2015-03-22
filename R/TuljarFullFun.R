TuljarFullFun <- function(n.states, n.environs, alpha, gamma, sex.ratio, samples.to.draw, 
                          tot.chains, joint.posterior.coda, posterior.names, intro.cost)
{
  Tuljar1.out <- TuljarFun1(alpha, gamma)
  Tuljar2.out <- TuljarFun2(alpha, gamma, sex.ratio, samples.to.draw, tot.chains, 
                            joint.posterior.coda, posterior.names, intro.cost)
  Tuljar3.out <- TuljarFun3(Tuljar2.out)
  Tuljar4.out <- TuljarFun4(Tuljar1.out, Tuljar2.out, Tuljar3.out)
  Tuljar5.out <- TuljarFun5(Tuljar1.out, Tuljar2.out, Tuljar3.out, Tuljar4.out)
  Tuljar6.out <- TuljarFun6(n.states, n.environs, Tuljar1.out, Tuljar2.out, Tuljar3.out)
  
  Q.naught <- Tuljar2.out$u.naught %*% t(Tuljar2.out$v.naught)
  W2 <- t(Tuljar2.out$v.naught) %*% Tuljar5.out$g1 %*% ((kronecker(Tuljar1.out$P.mat, Tuljar6.out$I.s) %*% (Tuljar6.out$out)) - (kronecker(Tuljar1.out$P.mat, Q.naught))) %*% Tuljar5.out$g2 %*% Tuljar2.out$u.naught
  
  return(list(W1 = Tuljar4.out$W1, W2 = W2))
}