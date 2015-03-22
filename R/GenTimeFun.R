GenTimeFun <- function(current.state, sex.ratio, samples.to.draw, tot.chains, 
                       joint.posterior.coda, posterior.names, intro.cost)
{
  leslie <- UpdateLeslieFun(current.state, sex.ratio, samples.to.draw, tot.chains, 
                            joint.posterior.coda, posterior.names, intro.cost)
  Tmat <- rbind(rep(0, dim(leslie)[1]), leslie[2:dim(leslie)[1], ])
  Fmat <- rbind(leslie[1, ], matrix(0, nrow = dim(Tmat)[1] - 1, ncol = dim(Tmat)[1]))
  I <- diag(rep(1, dim(Tmat)[1]))
  N <- ginv(I - Tmat)
  R <- Fmat %*% N
  Ro <- eigen(R)$values[1]
  lambda <- Re(eigen(leslie)$values[1])[1]
  gen.time <- log(Ro) / log(lambda)[1]
  return(list(Ro = Ro, lambda = lambda, gen.time = gen.time))
}