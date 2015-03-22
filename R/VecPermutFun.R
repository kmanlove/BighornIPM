VecPermutFun <- function(n.ages, n.classes, reps, alpha, gamma, intro.cost, sex.ratio, 
                         samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)
{
  #  intro.out.elast <- fade.out.elast <- rep(NA, reps)
  # 19 age classes; 2 environmental states. 
  fade.out.elast <- intro.elast <- rep(NA, reps)
  
  p = n.classes
  s = n.ages
  
  for(i in 1:reps){
    #     # P is the vec-permutation matrix (sensu Hunter and Caswell 2005 Ecological Modelling)
    #   P <- matrix(NA, nrow = 19 * 2, ncol = 19 * 2)
    #   e.full <- array(NA, dim = c(38, 38, 38))
    #   for(k in 1:2){
    #   for(j in 1:19){
    #       e.kj <- matrix(0,  nrow = 2, ncol = 19)
    #       e.kj[k, j] <- 1
    #       e.full[, , (k - 1) * 19 + j ] <- kronecker(e.kj, t(e.kj))
    #     }}
    #   P <- apply(e.full, c(1, 2), sum)
    
    P = matrix(0, s * p, s * p)
    for(l in 1:s){
      for(j in 1:p){
        E = matrix(0, s, p)
        E[l, j] = 1
        E
        P = P + kronecker(E, t(E))
      }
    }
    
    # B is block diagonal, with 2 19x19 blocks for the 2 environmental states.
    healthy.leslie <- UpdateLeslieFun(current.state = "healthy", sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
    spillover.leslie <- UpdateLeslieFun(current.state = "spillover", sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
    endemic.leslie <- UpdateLeslieFun(current.state = "infected", sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
    leslie.list <- list(healthy.leslie, spillover.leslie, endemic.leslie)
    B <- bdiag(leslie.list)
    
    # M is block diagonal with 19 2x2 blocks for the 19 demographic states
    small.M <- cbind(c(1 - alpha, alpha, 0), c(gamma, 0, 1 - gamma), c(gamma, 0, 1 - gamma)) # set up to be columns into rows, so that columns of M sum to 1. 
    M <- bdiag(list(small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M))
    
    # A is population projection matrix with environmental stochasticity
    A <-  B %*% t(P) %*% M %*% P
    A_eigens <- eigen(A) # complex????
    S_a <- sensitivity(A)
    
    # Sensitivity (environmental transitions)
    # assume disease status is updated before demography.
    S_m <- P %*% t(B) %*% S_a %*% t(P)
    #    round(S_m, 2)[1:10, 1:10]
    E_m <- (1 / Re(A_eigens$value)[1]) * M * S_m # regular * because it's a Hadamard production
    #    E_m <- (1 / .92) * M * S_m # regular * because it's a Hadamard production
    #    round(E_m, 2)[1:10, 1:10]
    # compare elasticities of fade-out to elasticity of reintroduction
    #    even.indices <- seq(1:19) * 2
    #    odd.indices <- seq(1:19) * 2 - 1
    upper.left.indices <- seq(1 : 19) * 3 - 2
    # elasticites are sum (off-diagonal elasts in 2x2 blocks) - sum(main-diag elasts)
    fade.out.elast[i] <- sum(E_m[upper.left.indices, upper.left.indices + 1] + E_m[upper.left.indices, upper.left.indices + 2]) - sum(E_m[upper.left.indices + 2, upper.left.indices + 1] + E_m[upper.left.indices + 2, upper.left.indices + 2])
    #    intro.indices <- seq(1:19) * 2 - 1
    intro.elast[i] <- sum(E_m[upper.left.indices, upper.left.indices]) - sum(E_m[upper.left.indices + 1, upper.left.indices])
  }
  return(list(fade.out.elast = fade.out.elast, intro.elast = intro.elast))
}
