TuljarFun1 <- function(alpha, gamma){
  # transitions out of healthy into spillover and infected
  H.trans <- c((1 - alpha), alpha, 0)
  Sp.trans <- c((1 - alpha) * gamma, alpha, (1 - alpha) * (1 - gamma))
  I.trans <- c((1 - alpha) * gamma, alpha, (1 - alpha) * (1 - gamma))
  P.mat <- cbind(H.trans, Sp.trans, I.trans)
  w.naught <- Re(eigen(P.mat)$vectors)[, 1]
  return(list(w.naught = w.naught, P.mat = P.mat))
}