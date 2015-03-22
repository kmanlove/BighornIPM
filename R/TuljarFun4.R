TuljarFun4 <- function(Tuljar1.out, Tuljar2.out, Tuljar3.out){
  
  v.naught.kron <- t(kronecker(Tuljar2.out$v.naught, Tuljar2.out$v.naught))
  u.naught.kron <- kronecker(Tuljar2.out$u.naught, Tuljar2.out$u.naught)
  
  C.he.kron <- kronecker(Tuljar3.out$Cx.he, Tuljar3.out$Cx.he)
  C.sp.kron <- kronecker(Tuljar3.out$Cx.sp, Tuljar3.out$Cx.sp)
  C.inf.kron <- kronecker(Tuljar3.out$Cx.inf, Tuljar3.out$Cx.inf)
  
  Sigma.naught <- Tuljar1.out$w.naught[1] * C.he.kron + Tuljar1.out$w.naught[2] * C.sp.kron + Tuljar1.out$w.naught[3] * C.inf.kron
  
  W1 <- (1 / 2 * Tuljar2.out$lambda.naught ^ 2) * v.naught.kron %*% Sigma.naught %*% u.naught.kron 
  
  return(list(v.naught.kron = v.naught.kron, u.naught.kron = u.naught.kron, C.he.kron = C.he.kron, C.sp.kron = C.sp.kron, C.inf.kron = C.inf.kron, Sigma.naught = Sigma.naught, W1 = W1))
}
