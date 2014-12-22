#-- simulation source functions --#

# function to update environmental state
update.status.fun <- function(alpha, gamma, current.state){
  current.state.new <- rep(NA, 1)
  if(current.state == "healthy"){
    gets.infected <- rbinom(1, 1, alpha)
#    current.state.new[1] <- ifelse(gets.infected == 0, "healthy", "spillover")
    current.state.new[1] <- ifelse(gets.infected == 0, "healthy", "infected")
  }
  else if(current.state == "infected"){
#    else if(current.state == "spillover"){
#      current.state.new[1] <- "infected"
#  } else if(current.state == "infected"){
    fade.out <- rbinom(1, 1, gamma)
    current.state.new[1] <- ifelse(fade.out == 1, "healthy", "infected")
  }
  return(list(current.state.new = current.state.new))
}

# function to update Leslie matrix parameters
update.leslie.fun <- function(current.state, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names){
  chain <- ceiling(runif(1, 0, tot.chains))
  post.draw <- as.data.frame(t(joint.posterior.coda[[chain]][sample(x = samples.to.draw, size = 1), ]))
  names(post.draw) <- posterior.names
  
  post.draw$beta.repro.1.2 <- exp(post.draw$beta.repro.1.2) / (1 + exp(post.draw$beta.repro.1.2))
  post.draw$beta.repro.1.3 <- exp(post.draw$beta.repro.1.3) / (1 + exp(post.draw$beta.repro.1.3))
  post.draw$beta.repro.1.4 <- exp(post.draw$beta.repro.1.4) / (1 + exp(post.draw$beta.repro.1.4))
  post.draw$beta.repro.1.5 <- exp(post.draw$beta.repro.1.5) / (1 + exp(post.draw$beta.repro.1.5))
  post.draw$beta.wean.1.2 <- exp(post.draw$beta.wean.1.2) / (1 + exp(post.draw$beta.wean.1.2))
  post.draw$beta.wean.1.3 <- exp(post.draw$beta.wean.1.3) / (1 + exp(post.draw$beta.wean.1.3))
  post.draw$beta.wean.1.4 <- exp(post.draw$beta.wean.1.4) / (1 + exp(post.draw$beta.wean.1.4))
  post.draw$beta.wean.1.5 <- exp(post.draw$beta.wean.1.5) / (1 + exp(post.draw$beta.wean.1.5))
#  post.draw$beta.overwinter.1 <- exp(post.draw$beta.overwinter.1) / (1 + exp(post.draw$beta.overwinter.1))
#  post.draw$beta.adsurv.1.1 <- exp(post.draw$beta.adsurv.1.1) / (1 + exp(post.draw$beta.adsurv.1.1))
  post.draw$beta.adsurv.1.2 <- exp(post.draw$beta.adsurv.1.2) / (1 + exp(post.draw$beta.adsurv.1.2))
  post.draw$beta.adsurv.1.3 <- exp(post.draw$beta.adsurv.1.3) / (1 + exp(post.draw$beta.adsurv.1.3))
  post.draw$beta.adsurv.1.4 <- exp(post.draw$beta.adsurv.1.4) / (1 + exp(post.draw$beta.adsurv.1.4))
  post.draw$beta.adsurv.1.5 <- exp(post.draw$beta.adsurv.1.5) / (1 + exp(post.draw$beta.adsurv.1.5))
  
  post.draw$beta.repro.2.2 <- exp(post.draw$beta.repro.2.2) / (1 + exp(post.draw$beta.repro.2.2))
  post.draw$beta.repro.2.3 <- exp(post.draw$beta.repro.2.3) / (1 + exp(post.draw$beta.repro.2.3))
  post.draw$beta.repro.2.4 <- exp(post.draw$beta.repro.2.4) / (1 + exp(post.draw$beta.repro.2.4))
  post.draw$beta.repro.2.5 <- exp(post.draw$beta.repro.2.5) / (1 + exp(post.draw$beta.repro.2.5))
  post.draw$beta.wean.2.2 <- exp(post.draw$beta.wean.2.2) / (1 + exp(post.draw$beta.wean.2.2))
  post.draw$beta.wean.2.3 <- exp(post.draw$beta.wean.2.3) / (1 + exp(post.draw$beta.wean.2.3))
  post.draw$beta.wean.2.4 <- exp(post.draw$beta.wean.2.4) / (1 + exp(post.draw$beta.wean.2.4))
  post.draw$beta.wean.2.5 <- exp(post.draw$beta.wean.2.5) / (1 + exp(post.draw$beta.wean.2.5))
#  post.draw$beta.overwinter.2 <- exp(post.draw$beta.overwinter.2) / (1 + exp(post.draw$beta.overwinter.2))
#  post.draw$beta.adsurv.2.1 <- exp(post.draw$beta.adsurv.2.1) / (1 + exp(post.draw$beta.adsurv.2.1))
  post.draw$beta.adsurv.2.2 <- exp(post.draw$beta.adsurv.2.2) / (1 + exp(post.draw$beta.adsurv.2.2))
  post.draw$beta.adsurv.2.3 <- exp(post.draw$beta.adsurv.2.3) / (1 + exp(post.draw$beta.adsurv.2.3))
  post.draw$beta.adsurv.2.4 <- exp(post.draw$beta.adsurv.2.4) / (1 + exp(post.draw$beta.adsurv.2.4))
  post.draw$beta.adsurv.2.5 <- exp(post.draw$beta.adsurv.2.5) / (1 + exp(post.draw$beta.adsurv.2.5))
  
#  leslie.out <- rep(NA, 6, 6)
  if(current.state == "healthy"){
#     repros <- c(0, rep((post.draw$beta.repro.1.2 * post.draw$beta.wean.1.2 * post.draw$beta.overwinter.1), 3), rep((post.draw$beta.repro.1.3 * post.draw$beta.wean.1.3 * post.draw$beta.overwinter.1), 6), rep((post.draw$beta.repro.1.4 * post.draw$beta.wean.1.4 * post.draw$beta.overwinter.1), 5), rep((post.draw$beta.repro.1.5 * post.draw$beta.wean.1.5 * post.draw$beta.overwinter.1), 4))    
#     survs <- c(1, rep(post.draw$beta.adsurv.1.2, 2), rep(post.draw$beta.adsurv.1.3, 6), rep(post.draw$beta.adsurv.1.4, 5), rep(post.draw$beta.adsurv.1.5, 3), 0)    
#     leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs))))
    repros <- c(0, rep((post.draw$beta.repro.1.2 * post.draw$beta.wean.1.2 ), 3), rep((post.draw$beta.repro.1.3 * post.draw$beta.wean.1.3 ), 6), rep((post.draw$beta.repro.1.4 * post.draw$beta.wean.1.4 ), 5), rep((post.draw$beta.repro.1.5 * post.draw$beta.wean.1.5 ), 4))    
    survs <- c(1, rep(post.draw$beta.adsurv.1.2, 2), rep(post.draw$beta.adsurv.1.3, 6), rep(post.draw$beta.adsurv.1.4, 5), rep(post.draw$beta.adsurv.1.5, 3), 0)    
    leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs))))
  }
  else {
#    repros <- c(0, rep((post.draw$beta.repro.2.2 * post.draw$beta.wean.2.2 * post.draw$beta.overwinter.2), 1), rep((post.draw$beta.repro.2.3 * post.draw$beta.wean.2.3 * post.draw$beta.overwinter.2), 6), rep((post.draw$beta.repro.2.4 * post.draw$beta.wean.2.4 * post.draw$beta.overwinter.2), 6), rep((post.draw$beta.repro.2.5 * post.draw$beta.wean.2.5 * post.draw$beta.overwinter.2), 4))    
#    survs <- c(1, rep(post.draw$beta.adsurv.2.2, 1), rep(post.draw$beta.adsurv.2.3, 6), rep(post.draw$beta.adsurv.2.4, 6), rep(post.draw$beta.adsurv.2.5, 3), 0)    
#    leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs))))
    repros <- c(0, rep((post.draw$beta.repro.2.2 * post.draw$beta.wean.2.2 ), 3), rep((post.draw$beta.repro.2.3 * post.draw$beta.wean.2.3 ), 6), rep((post.draw$beta.repro.2.4 * post.draw$beta.wean.2.4 ), 5), rep((post.draw$beta.repro.2.5 * post.draw$beta.wean.2.5 ), 4))    
    survs <- c(1, rep(post.draw$beta.adsurv.2.2, 2), rep(post.draw$beta.adsurv.2.3, 6), rep(post.draw$beta.adsurv.2.4, 5), rep(post.draw$beta.adsurv.2.5, 3), 0)    
    leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs))))
  }
  return(leslie)
}

# project.fun <- function(timesteps, ages.init, alpha, gamma, he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr){
#   N <- matrix(NA, nrow = length(ages.init), ncol = timesteps)
#   N[, 1] <- ages.init
#   tot.pop.size <- log.lambda.s <- rep(NA, length = timesteps)
#   disease.status <- rep(NA, timesteps)
#   disease.status[1:11] <- c(rep("healthy", 10), "spillover")
#   for(i in 1:10){
#     new.leslie <- update.leslie.fun(current.state = disease.status[i], he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
#     N[, i + 1] <- round(t(N[ , i]) %*% new.leslie) 
#     tot.pop.size[i] <- sum(N[ , i])
#     if(i == 1){
#       log.lambda.s[i] <- NA
#     } else {
#       log.lambda.s[i] <- ifelse(tot.pop.size[i] == 0, NA, log(tot.pop.size[i] / tot.pop.size[i - 1]))
#     }
#   }
#   for(i in 11:(timesteps - 1)){
#     new.leslie <- update.leslie.fun(current.state = disease.status[i], he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
#     N[, i + 1] <- round(t(N[ , i]) %*% new.leslie)
#     disease.status[i + 1] <- update.status.fun(alpha, gamma, current.state = disease.status[i])$current.state.new[1]
#     tot.pop.size[i] <- sum(N[ , i])
#     if(i == 1){
#       log.lambda.s[i] <- NA
#     } else {
#       log.lambda.s[i] <- ifelse(tot.pop.size[i] == 0, NA, log(tot.pop.size[i] / tot.pop.size[i - 1]))
#     }
#   }
#   out.list <- list(N = N, disease.status = disease.status, tot.pop.size = tot.pop.size, log.lambda.s = log.lambda.s)
#   return(out.list)
# }

project.fun <- function(timesteps, ages.init, alpha, gamma, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names){
  N <- matrix(NA, nrow = length(ages.init), ncol = timesteps)
  N[, 1] <- ages.init
  tot.pop.size <- log.lambda.s <- rep(NA, length = timesteps)
  disease.status <- rep("healthy", timesteps)
  for(i in 1:10){
    new.leslie <- update.leslie.fun(current.state = disease.status[i], samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)
    N[, i + 1] <- t(N[ , i]) %*% new.leslie  
    tot.pop.size[i] <- sum(N[ , i])
    if(i == 1){
      log.lambda.s[i] <- NA
    } else {
      log.lambda.s[i] <- log(tot.pop.size[i] / tot.pop.size[i - 1])
    }
  }
  for(i in 11:(timesteps - 1)){
    new.leslie <- update.leslie.fun(current.state = disease.status[i], samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)
    N[, i + 1] <- t(N[ , i]) %*% new.leslie  
    disease.status[i + 1] <- update.status.fun(alpha, gamma, current.state = disease.status[i])$current.state.new[1]
    tot.pop.size[i] <- sum(N[ , i])
    if(i == 1){
      log.lambda.s[i] <- NA
    } else {
      log.lambda.s[i] <- log(tot.pop.size[i] / tot.pop.size[i - 1])
    }
  }
  out.list <- list(N = N, disease.status = disease.status, tot.pop.size = tot.pop.size, log.lambda.s = log.lambda.s)
  return(out.list)
}

# projection function in absence of disease
# healthy.project.fun <- function(timesteps, ages.init, alpha, gamma, he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr){
healthy.project.fun <- function(timesteps, ages.init, alpha, gamma, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names){
  N <- matrix(NA, nrow = length(ages.init), ncol = timesteps)
  N[, 1] <- ages.init
  tot.pop.size <- log.lambda.s <- rep(NA, length = timesteps)
  disease.status <- rep("healthy", timesteps)
  for(i in 1:(timesteps - 1)){
    new.leslie <- update.leslie.fun(current.state = disease.status[i], samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)
    N[, i + 1] <- t(N[ , i]) %*% new.leslie  
    tot.pop.size[i] <- sum(N[ , i])
    if(i == 1){
      log.lambda.s[i] <- NA
    } else {
      log.lambda.s[i] <- log(tot.pop.size[i] / tot.pop.size[i - 1])
    }
  }
  out.list <- list(N = N, disease.status = disease.status, tot.pop.size = tot.pop.size, log.lambda.s = log.lambda.s)
  return(out.list)
}

popgrowth.sim.fun <- function(timesteps, reps, alpha.range, gamma.range, alpha.steps, gamma.steps){
  alphas <- 1 / seq(min(alpha.range), max(alpha.range), length.out = alpha.steps)  
  gammas <- 1 / seq(min(gamma.range), max(gamma.range), length.out = gamma.steps)
  alpha.gamma.frame <- expand.grid(alphas, gammas)
  popsize.ij <- loglambda.ij <- vector("list", dim(alpha.gamma.frame)[1])
  mean.lnlambda <- rep(NA, dim(alpha.gamma.frame)[1])
  for(i in 1:length(popsize.ij)){
    popsize.ij[[i]] <- loglambda.ij[[i]] <- matrix(NA, nrow = timesteps, ncol = reps)
    for(j in 1:reps){
      popgrowthsim.ij <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = alpha.gamma.frame[i, 1], gamma = alpha.gamma.frame[i, 2], he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
      popsize.ij[[i]][, j] <- popgrowthsim.ij$tot.pop.size
      loglambda.ij[[i]][, j] <- popgrowthsim.ij$log.lambda.s
    }
    mean.lnlambda[i] <- mean(na.omit(unlist(as.vector(loglambda.ij[[i]][-c(11, 12), ]))))
  }
  outlist <- list(alpha.gamma.frame = alpha.gamma.frame, popsize.ij = popsize.ij, loglambda.ij = loglambda.ij, mean.lnlambda = mean.lnlambda)
  return(outlist)
}

vec.permut.fun <- function(reps, alpha, gamma){
#  intro.out.elast <- fade.out.elast <- rep(NA, reps)
  # 19 age classes; 2 environmental states. 
  fade.out.elast <- intro.elast <- rep(NA, reps)
  for(i in 1:reps){
    # P is the vec-permutation matrix (sensu Hunter and Caswell 2005 Ecology)
  P <- matrix(NA, nrow = 19 * 2, ncol = 19 * 2)
  e.full <- array(NA, dim = c(38, 38, 38))
  for(k in 1:2){
  for(j in 1:19){
      e.kj <- matrix(0,  nrow = 2, ncol = 19)
      e.kj[k, j] <- 1
      e.full[, , (k - 1) * 19 + j ] <- kronecker(e.kj, t(e.kj))
    }}
  P <- apply(e.full, c(1, 2), sum)
    
    # B is block diagonal, with 3 19x19 blocks for the 3 environmental states.
    healthy.leslie <- update.leslie.fun(current.state = "healthy", samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)
    endemic.leslie <-  update.leslie.fun(current.state = "infected", samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)
    leslie.list <- list(healthy.leslie, endemic.leslie)
  B <- bdiag(leslie.list)
    
    # M is block diagonal with 20 2x2 blocks for the 20 demographic states
    small.M <- rbind(c(1 - alpha, alpha), c(gamma, 1 - gamma))
    M <- bdiag(list(small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M))
    
    # A is population projection matrix with environmental stochasticity
    A <- t(P) %*% M %*% P %*% B
    A_eigens <- eigen(A) # complex????
    S_a <- sensitivity(A)
    
    # Sensitivity (environmental transitions)
    # assume disease status is updated before demography.
    S_m <- P %*% t(B) %*% S_a %*% t(P)
    round(S_m, 2)[1:10, 1:10]
    E_m <- (1 / A_eigens$value[1]) * M * S_m # regular * because it's a Hadamard production
#    E_m <- (1 / .92) * M * S_m # regular * because it's a Hadamard production
    round(E_m, 2)[1:10, 1:10]
    # compare elasticities of fade-out to elasticity of reintroduction
    fade.indices <- seq(1:19) * 2
    fade.out.elast[i] <- sum(E_m[fade.indices, fade.indices])
    intro.indices <- seq(1:19) * 2 - 1
    intro.elast[i] <- sum(E_m[intro.indices, intro.indices])
  }
  return(list(fade.out.elast = fade.out.elast, intro.elast = intro.elast))
}


addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}