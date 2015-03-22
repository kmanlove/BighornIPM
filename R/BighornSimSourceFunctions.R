#-- simulation source functions --#

# # function to update environmental state
# update.status.fun <- function(alpha, gamma, current.state){
#   current.state.new <- rep(NA, 1)
#   if(current.state == "healthy"){
#     gets.infected <- rbinom(1, 1, alpha)
#     current.state.new[1] <- ifelse(gets.infected == 0, "healthy", "spillover")
#   }
# #  else if(current.state == "infected"){
#     else if(current.state == "spillover"){
#       gets.infected <- rbinom(1, 1, alpha)
#       fade.out <- rbinom(1, 1, gamma)
#       current.state.new[1] <- ifelse(gets.infected == 1, "spillover", ifelse((gets.infected == 0 & fade.out == 1), "healthy", "infected"))
#   } else if(current.state == "infected"){
#     gets.infected <- rbinom(1, 1, alpha)
#     fade.out <- rbinom(1, 1, gamma)
#     current.state.new[1] <- ifelse(gets.infected == 1, "spillover", ifelse((gets.infected == 0 & fade.out == 1), "healthy", "infected"))
#   }
#   return(list(current.state.new = current.state.new))
# }
# 
# # function to update Leslie matrix parameters
# update.leslie.fun <- function(current.state, sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost){
#   chain <- ceiling(runif(1, 0, tot.chains))
#   post.draw <- as.data.frame(t(joint.posterior.coda[[chain]][sample(x = samples.to.draw, size = 1), ]))
#   names(post.draw) <- posterior.names
# 
#   post.draw$beta.wean.1.2 <- exp(post.draw$beta.wean.1.2) / (1 + exp(post.draw$beta.wean.1.2))
#   post.draw$beta.wean.1.3 <- exp(post.draw$beta.wean.1.3) / (1 + exp(post.draw$beta.wean.1.3))
#   post.draw$beta.wean.1.4 <- exp(post.draw$beta.wean.1.4) / (1 + exp(post.draw$beta.wean.1.4))
#   post.draw$beta.wean.1.5 <- exp(post.draw$beta.wean.1.5) / (1 + exp(post.draw$beta.wean.1.5))
#   post.draw$beta.wean.1.6 <- exp(post.draw$beta.wean.1.6) / (1 + exp(post.draw$beta.wean.1.6))
#   post.draw$beta.adsurv.1.2 <- exp(post.draw$beta.adsurv.1.2) / (1 + exp(post.draw$beta.adsurv.1.2))
#   post.draw$beta.adsurv.1.3 <- exp(post.draw$beta.adsurv.1.3) / (1 + exp(post.draw$beta.adsurv.1.3))
#   post.draw$beta.adsurv.1.4 <- exp(post.draw$beta.adsurv.1.4) / (1 + exp(post.draw$beta.adsurv.1.4))
#   post.draw$beta.adsurv.1.5 <- exp(post.draw$beta.adsurv.1.5) / (1 + exp(post.draw$beta.adsurv.1.5))
#   post.draw$beta.adsurv.1.6 <- exp(post.draw$beta.adsurv.1.6) / (1 + exp(post.draw$beta.adsurv.1.6))
# 
#   post.draw$beta.wean.2.2 <- exp(post.draw$beta.wean.2.2) / (1 + exp(post.draw$beta.wean.2.2))
#   post.draw$beta.wean.2.3 <- exp(post.draw$beta.wean.2.3) / (1 + exp(post.draw$beta.wean.2.3))
#   post.draw$beta.wean.2.4 <- exp(post.draw$beta.wean.2.4) / (1 + exp(post.draw$beta.wean.2.4))
#   post.draw$beta.wean.2.5 <- exp(post.draw$beta.wean.2.5) / (1 + exp(post.draw$beta.wean.2.5))
#   post.draw$beta.adsurv.2.2 <- exp(post.draw$beta.adsurv.2.2) / (1 + exp(post.draw$beta.adsurv.2.2))
#   post.draw$beta.adsurv.2.3 <- exp(post.draw$beta.adsurv.2.3) / (1 + exp(post.draw$beta.adsurv.2.3))
#   post.draw$beta.adsurv.2.4 <- exp(post.draw$beta.adsurv.2.4) / (1 + exp(post.draw$beta.adsurv.2.4))
#   post.draw$beta.adsurv.2.5 <- exp(post.draw$beta.adsurv.2.5) / (1 + exp(post.draw$beta.adsurv.2.5))
#   
#   if(current.state == "healthy"){
#     repros <- c(0, rep((post.draw$beta.wean.1.2 * sex.ratio), 1), rep((post.draw$beta.wean.1.3  * sex.ratio), 5), rep((post.draw$beta.wean.1.4  * sex.ratio), 6), rep((post.draw$beta.wean.1.5  * sex.ratio), 6))    
#     survs <- c(1, rep(post.draw$beta.adsurv.1.2, 1), rep(post.draw$beta.adsurv.1.3, 5), rep(post.draw$beta.adsurv.1.4, 6), rep(post.draw$beta.adsurv.1.5, 5))    
#     leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, dim(diag(survs))[1])))
#   }
#   else if(current.state == "infected"){
#     repros <- c(0, rep((post.draw$beta.wean.2.2 * sex.ratio), 1), rep((post.draw$beta.wean.2.3  * sex.ratio), 5), rep((post.draw$beta.wean.2.4  * sex.ratio), 6), rep((post.draw$beta.wean.2.5 * sex.ratio), 6))    
#     survs <- c(1, rep(post.draw$beta.adsurv.2.2, 1), rep(post.draw$beta.adsurv.2.3, 5), rep(post.draw$beta.adsurv.2.4, 6), rep(post.draw$beta.adsurv.2.5, 5))    
#     leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs))))
#   } 
#   else if(current.state == "spillover"){
#     repros <- c(0, rep((post.draw$beta.wean.2.2 * sex.ratio), 1), rep((post.draw$beta.wean.2.3  * sex.ratio), 5), rep((post.draw$beta.wean.2.4  * sex.ratio), 6), rep((post.draw$beta.wean.2.5 * sex.ratio), 6))    
#     survs <- c(1, rep(post.draw$beta.adsurv.1.2, 1), rep(post.draw$beta.adsurv.1.3, 5), rep(post.draw$beta.adsurv.1.4, 6), rep(post.draw$beta.adsurv.1.5, 5)) * (1 - intro.cost)
#     leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs)))) 
#   }
#   return(leslie)
# }
# 
# # Johnson age-structure version of update.leslie.fun
# johnson.update.leslie.fun <- function(current.state, sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost){
#   chain <- ceiling(runif(1, 0, tot.chains))
#   post.draw <- as.data.frame(t(joint.posterior.coda[[chain]][sample(x = samples.to.draw, size = 1), ]))
#   names(post.draw) <- posterior.names
#   
#   post.draw$beta.wean.1.2 <- exp(post.draw$beta.wean.1.2) / (1 + exp(post.draw$beta.wean.1.2))
#   post.draw$beta.wean.1.3 <- exp(post.draw$beta.wean.1.3) / (1 + exp(post.draw$beta.wean.1.3))
#   post.draw$beta.wean.1.4 <- exp(post.draw$beta.wean.1.4) / (1 + exp(post.draw$beta.wean.1.4))
# #  post.draw$beta.wean.1.5 <- exp(post.draw$beta.wean.1.5) / (1 + exp(post.draw$beta.wean.1.5))
# #  post.draw$beta.wean.1.6 <- exp(post.draw$beta.wean.1.6) / (1 + exp(post.draw$beta.wean.1.6))
#   post.draw$beta.adsurv.1.2 <- exp(post.draw$beta.adsurv.1.2) / (1 + exp(post.draw$beta.adsurv.1.2))
#   post.draw$beta.adsurv.1.3 <- exp(post.draw$beta.adsurv.1.3) / (1 + exp(post.draw$beta.adsurv.1.3))
#   post.draw$beta.adsurv.1.4 <- exp(post.draw$beta.adsurv.1.4) / (1 + exp(post.draw$beta.adsurv.1.4))
# #  post.draw$beta.adsurv.1.5 <- exp(post.draw$beta.adsurv.1.5) / (1 + exp(post.draw$beta.adsurv.1.5))
# #  post.draw$beta.adsurv.1.6 <- exp(post.draw$beta.adsurv.1.6) / (1 + exp(post.draw$beta.adsurv.1.6))
#   
#   post.draw$beta.wean.2.2 <- exp(post.draw$beta.wean.2.2) / (1 + exp(post.draw$beta.wean.2.2))
#   post.draw$beta.wean.2.3 <- exp(post.draw$beta.wean.2.3) / (1 + exp(post.draw$beta.wean.2.3))
#   post.draw$beta.wean.2.4 <- exp(post.draw$beta.wean.2.4) / (1 + exp(post.draw$beta.wean.2.4))
# #  post.draw$beta.wean.2.5 <- exp(post.draw$beta.wean.2.5) / (1 + exp(post.draw$beta.wean.2.5))
#   post.draw$beta.adsurv.2.2 <- exp(post.draw$beta.adsurv.2.2) / (1 + exp(post.draw$beta.adsurv.2.2))
#   post.draw$beta.adsurv.2.3 <- exp(post.draw$beta.adsurv.2.3) / (1 + exp(post.draw$beta.adsurv.2.3))
#   post.draw$beta.adsurv.2.4 <- exp(post.draw$beta.adsurv.2.4) / (1 + exp(post.draw$beta.adsurv.2.4))
# #  post.draw$beta.adsurv.2.5 <- exp(post.draw$beta.adsurv.2.5) / (1 + exp(post.draw$beta.adsurv.2.5))
#   
#   if(current.state == "healthy"){
# #    repros <- c(0, rep((post.draw$beta.wean.1.2 * sex.ratio), 1), rep((post.draw$beta.wean.1.3  * sex.ratio), 5), rep((post.draw$beta.wean.1.4  * sex.ratio), 6), rep((post.draw$beta.wean.1.5  * sex.ratio), 6))    
# #    survs <- c(1, rep(post.draw$beta.adsurv.1.2, 1), rep(post.draw$beta.adsurv.1.3, 5), rep(post.draw$beta.adsurv.1.4, 6), rep(post.draw$beta.adsurv.1.5, 5))    
# #    leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, dim(diag(survs))[1])))
#     repros <- c(0, rep((post.draw$beta.wean.1.2 * sex.ratio), 1), rep((post.draw$beta.wean.1.3  * sex.ratio), 17))  
#     survs <- c(1, rep(post.draw$beta.adsurv.1.2, 1), rep(post.draw$beta.adsurv.1.3, 16))
#     leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, dim(diag(survs))[1])))
#   }
#   else if(current.state == "infected"){
# #    repros <- c(0, rep((post.draw$beta.wean.2.2 * sex.ratio), 1), rep((post.draw$beta.wean.2.3  * sex.ratio), 5), rep((post.draw$beta.wean.2.4  * sex.ratio), 6), rep((post.draw$beta.wean.2.5 * sex.ratio), 6))    
# #    survs <- c(1, rep(post.draw$beta.adsurv.2.2, 1), rep(post.draw$beta.adsurv.2.3, 5), rep(post.draw$beta.adsurv.2.4, 6), rep(post.draw$beta.adsurv.2.5, 5))    
# #    leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs))))
#     repros <- c(0, rep((post.draw$beta.wean.2.2 * sex.ratio), 1), rep((post.draw$beta.wean.2.3  * sex.ratio), 17))    
#     survs <- c(1, rep(post.draw$beta.adsurv.2.2, 1), rep(post.draw$beta.adsurv.2.3, 16))    
#     leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs))))
#   } 
#   else if(current.state == "spillover"){
# #    repros <- c(0, rep((post.draw$beta.wean.2.2 * sex.ratio), 1), rep((post.draw$beta.wean.2.3  * sex.ratio), 5), rep((post.draw$beta.wean.2.4  * sex.ratio), 6), rep((post.draw$beta.wean.2.5 * sex.ratio), 6))    
# #    survs <- c(1, rep(post.draw$beta.adsurv.1.2, 1), rep(post.draw$beta.adsurv.1.3, 5), rep(post.draw$beta.adsurv.1.4, 6), rep(post.draw$beta.adsurv.1.5, 5)) * (1 - intro.cost)
# #    leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs)))) 
#     repros <- c(0, rep((post.draw$beta.wean.2.2 * sex.ratio), 1), rep((post.draw$beta.wean.2.3  * sex.ratio), 17))    
#     survs <- c(1, rep(post.draw$beta.adsurv.1.2, 1), rep(post.draw$beta.adsurv.1.3, 16)) * (1 - intro.cost)
#     leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs)))) 
#   }
#   return(leslie)
# }

project.fun <- function(johnson, timesteps, sex.ratio, ages.init, alpha, gamma, intro.cost, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, fixed.start.time){
  N <- matrix(NA, nrow = length(ages.init), ncol = timesteps)
  N[, 1] <- ages.init
  tot.pop.size <- log.lambda.s <- rep(NA, length = timesteps)
  disease.status <- rep(NA, timesteps)
  disease.status[1:10] <- "healthy"
  if(johnson == F){
  for(i in 1:10){
    new.leslie <- UpdateLeslieFun(current.state = disease.status[i], sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
    N[, i + 1] <- t(N[ , i]) %*% new.leslie  
    tot.pop.size[i] <- sum(N[ , i])
    if(i == 1){
      log.lambda.s[i] <- NA
    } else {
      log.lambda.s[i] <- log(tot.pop.size[i] / tot.pop.size[i - 1])
    }
  }
  if(fixed.start.time == T){
    for(i in c(11)){
      disease.status[i] <- "spillover"
      new.leslie <- UpdateLeslieFun(current.state = disease.status[i], sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
      N[, i + 1] <- t(N[ , i]) %*% new.leslie  
      disease.status[i + 1] <- UpdateStatusFun(alpha, gamma, current.state = disease.status[i])$current.state.new[1]
      tot.pop.size[i] <- sum(N[ , i])
      if(i == 1){
        log.lambda.s[i] <- NA
      } else {
        log.lambda.s[i] <- log(tot.pop.size[i] / tot.pop.size[i - 1])
      }
    }
    for(i in 12:(timesteps - 1)){
      new.leslie <- UpdateLeslieFun(current.state = disease.status[i], sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
      N[, i + 1] <- t(N[ , i]) %*% new.leslie  
      disease.status[i + 1] <- UpdateStatusFun(alpha, gamma, current.state = disease.status[i])$current.state.new[1]
      tot.pop.size[i] <- sum(N[ , i])
      if(i == 1){
        log.lambda.s[i] <- NA
      } else {
        log.lambda.s[i] <- log(tot.pop.size[i] / tot.pop.size[i - 1])
      }
    }
  } else{
  for(i in 11:(timesteps - 1)){
    disease.status[11] <- UpdateStatusFun(alpha, gamma, current.state = disease.status[10]$current.state.new[1])
    new.leslie <- UpdateLeslieFun(current.state = disease.status[i], sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
    N[, i + 1] <- t(N[ , i]) %*% new.leslie  
    disease.status[i + 1] <- UpdateStatusFun(alpha, gamma, current.state = disease.status[i])$current.state.new[1]
    tot.pop.size[i] <- sum(N[ , i])
    if(i == 1){
      log.lambda.s[i] <- NA
    } else {
      log.lambda.s[i] <- log(tot.pop.size[i] / tot.pop.size[i - 1])
    }
  }
  }
  }
  else{
    for(i in 1:10){
      new.leslie <- JohnsonUpdateLeslieFun(current.state = disease.status[i], sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
      N[, i + 1] <- t(N[ , i]) %*% new.leslie  
      tot.pop.size[i] <- sum(N[ , i])
      if(i == 1){
        log.lambda.s[i] <- NA
      } else {
        log.lambda.s[i] <- log(tot.pop.size[i] / tot.pop.size[i - 1])
      }
    }
    if(fixed.start.time == T){
      for(i in c(11)){
        disease.status[i] <- "spillover"
        new.leslie <- JohnsonUpdateLeslieFun(current.state = disease.status[i], sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
        N[, i + 1] <- t(N[ , i]) %*% new.leslie  
        disease.status[i + 1] <- johnson.update.status.fun(alpha, gamma, current.state = disease.status[i])$current.state.new[1]
        tot.pop.size[i] <- sum(N[ , i])
        if(i == 1){
          log.lambda.s[i] <- NA
        } else {
          log.lambda.s[i] <- log(tot.pop.size[i] / tot.pop.size[i - 1])
        }
      }
      for(i in 12:(timesteps - 1)){
        new.leslie <- JohnsonUpdateLeslieFun(current.state = disease.status[i], sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
        N[, i + 1] <- t(N[ , i]) %*% new.leslie  
        disease.status[i + 1] <- johnson.update.status.fun(alpha, gamma, current.state = disease.status[i])$current.state.new[1]
        tot.pop.size[i] <- sum(N[ , i])
        if(i == 1){
          log.lambda.s[i] <- NA
        } else {
          log.lambda.s[i] <- log(tot.pop.size[i] / tot.pop.size[i - 1])
        }
      }
    } else{
      for(i in 11:(timesteps - 1)){
        disease.status[11] <- johnson.update.status.fun(alpha, gamma, current.state = disease.status[10]$current.state.new[1])
        new.leslie <- JohnsonUpdateLeslieFun(current.state = disease.status[i], sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
        N[, i + 1] <- t(N[ , i]) %*% new.leslie  
        disease.status[i + 1] <- johnson.update.status.fun(alpha, gamma, current.state = disease.status[i])$current.state.new[1]
        tot.pop.size[i] <- sum(N[ , i])
        if(i == 1){
          log.lambda.s[i] <- NA
        } else {
          log.lambda.s[i] <- log(tot.pop.size[i] / tot.pop.size[i - 1])
        }
      }
    }
  }
  out.list <- list(N = N, disease.status = disease.status, tot.pop.size = tot.pop.size, log.lambda.s = log.lambda.s)
  return(out.list)
}

# projection function in absence of disease
# healthy.project.fun <- function(timesteps, ages.init, alpha, gamma, he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr){
healthy.project.fun <- function(timesteps, sex.ratio, ages.init, alpha, gamma, intro.cost, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names){
  N <- matrix(NA, nrow = length(ages.init), ncol = timesteps)
  N[, 1] <- ages.init
  tot.pop.size <- log.lambda.s <- rep(NA, length = timesteps)
  disease.status <- rep("healthy", timesteps)
  for(i in 1:(timesteps - 1)){
    new.leslie <- UpdateLeslieFun(current.state = disease.status[i], sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
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

# popgrowth.sim.fun <- function(timesteps, reps, alpha.range, gamma.range, alpha.steps, gamma.steps, sex.ratio, intro.cost){
#   alphas <- 1 / seq(min(alpha.range), max(alpha.range), length.out = alpha.steps)  
#   gammas <- 1 / seq(min(gamma.range), max(gamma.range), length.out = gamma.steps)
#   alpha.gamma.frame <- expand.grid(alphas, gammas)
#   popsize.ij <- loglambda.ij <- vector("list", dim(alpha.gamma.frame)[1])
#   mean.lnlambda <- rep(NA, dim(alpha.gamma.frame)[1])
#   for(i in 1:length(popsize.ij)){
#     popsize.ij[[i]] <- loglambda.ij[[i]] <- matrix(NA, nrow = timesteps, ncol = reps)
#     for(j in 1:reps){
#       popgrowthsim.ij <- project.fun(timesteps = timesteps, sex.ratio, ages.init = ages.init, alpha = alpha.gamma.frame[i, 1], gamma = alpha.gamma.frame[i, 2], he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
#       popsize.ij[[i]][, j] <- popgrowthsim.ij$tot.pop.size
#       loglambda.ij[[i]][, j] <- popgrowthsim.ij$log.lambda.s
#     }
#     mean.lnlambda[i] <- mean(na.omit(unlist(as.vector(loglambda.ij[[i]][-c(11, 12), ]))))
#   }
#   outlist <- list(alpha.gamma.frame = alpha.gamma.frame, popsize.ij = popsize.ij, loglambda.ij = loglambda.ij, mean.lnlambda = mean.lnlambda)
#   return(outlist)
# }

vec.permut.fun <- function(n.ages, n.classes, reps, alpha, gamma, intro.cost, sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names){
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

# function to calculate generation time (following Caswell 2001 pg 126-128 + 112 (for definition of N))
GenTimeFun <- function(current.state, sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost){
  leslie <- UpdateLeslieFun(current.state, sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
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

# Johnson age-structure vsn of GenTimeFun
Jo.GenTimeFun <- function(current.state, sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost){
  leslie <- JohnsonUpdateLeslieFun(current.state, sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
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

# functions from Tuljapurkar and Haridas Ecol Letters 2006 for interannual variation vs. temp autocorr
# 1) build environ state-transition matrix and extract eigenvectors
TuljarFun1 <- function(alpha, gamma){
  # transitions out of healthy into spillover and infected
  H.trans <- c((1 - alpha), alpha, 0)
  Sp.trans <- c((1 - alpha) * gamma, alpha, (1 - alpha) * (1 - gamma))
  I.trans <- c((1 - alpha) * gamma, alpha, (1 - alpha) * (1 - gamma))
  P.mat <- cbind(H.trans, Sp.trans, I.trans)
  w.naught <- Re(eigen(P.mat)$vectors)[, 1]
  return(list(w.naught = w.naught, P.mat = P.mat))
}

# TuljarFun1(alpha, gamma)
# 2) compute average projection matrix b
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

# TuljarFun2(alpha, gamma, current.state, sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)

# compute D, Cx
TuljarFun3 <- function(Tuljar2.out){
  D <- Tuljar2.out$b / Tuljar2.out$lambda.naught - Tuljar2.out$u.naught %*% t(Tuljar2.out$v.naught)
  
  Cx.he <- Tuljar2.out$L.healthy - Tuljar2.out$b
  Cx.sp <- Tuljar2.out$L.sp - Tuljar2.out$b
  Cx.inf <- Tuljar2.out$L.infected - Tuljar2.out$b
  
  return(list(D = D, Cx.he = Cx.he, Cx.sp = Cx.sp, Cx.inf = Cx.inf))
}
# TuljarFun3(alpha, gamma, current.state, sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)

# 
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

# TuljarFun4(alpha, gamma, current.state, sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
TuljarFun5 <- function(Tuljar1.out, Tuljar2.out, Tuljar3.out, Tuljar4.out){
  g1 <- (1 / Tuljar2.out$lambda.naught) * cbind(Tuljar3.out$Cx.he, Tuljar3.out$Cx.sp, Tuljar3.out$Cx.inf)
  g2 <- (1 / Tuljar2.out$lambda.naught) * rbind(Tuljar1.out$w.naught[1] * t(Tuljar3.out$Cx.he), Tuljar1.out$w.naught[2] * t(Tuljar3.out$Cx.sp), Tuljar1.out$w.naught[3] * t(Tuljar3.out$Cx.inf))

  return(list(g1 = g1, g2 = g2))
}

TuljarFun6 <- function(n.states, n.environs, Tuljar1.out, Tuljar2.out, Tuljar3.out){
#  Tuljar1.out <- TuljarFun1(alpha, gamma)
#  Tuljar2.out <- TuljarFun2(alpha, gamma, sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
#  Tuljar3.out <- TuljarFun3(Tuljar2.out)
  I.s <- diag(1, n.states)
  I.sk <- diag(1, n.states * n.environs)
  
  out <- ginv(I.sk - kronecker(Tuljar1.out$P.mat, Tuljar3.out$D))
  
  return(list(I.s = I.s, I.sk = I.sk, out = out))
}

TuljarFunFull <- function(n.states, n.environs, alpha, gamma, sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost){
  Tuljar1.out <- TuljarFun1(alpha, gamma)
  Tuljar2.out <- TuljarFun2(alpha, gamma, sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
  Tuljar3.out <- TuljarFun3(Tuljar2.out)
  Tuljar4.out <- TuljarFun4(Tuljar1.out, Tuljar2.out, Tuljar3.out)
  Tuljar5.out <- TuljarFun5(Tuljar1.out, Tuljar2.out, Tuljar3.out, Tuljar4.out)
  Tuljar6.out <- TuljarFun6(n.states, n.environs, Tuljar1.out, Tuljar2.out, Tuljar3.out)
  
  Q.naught <- Tuljar2.out$u.naught %*% t(Tuljar2.out$v.naught)
  W2 <- t(Tuljar2.out$v.naught) %*% Tuljar5.out$g1 %*% ((kronecker(Tuljar1.out$P.mat, Tuljar6.out$I.s) %*% (Tuljar6.out$out)) - (kronecker(Tuljar1.out$P.mat, Q.naught))) %*% Tuljar5.out$g2 %*% Tuljar2.out$u.naught
  
  return(list(W1 = Tuljar4.out$W1, W2 = W2))
}

#TuljarFunFull(n.states, n.environs, alpha, gamma, current.state, sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost)
pupilplot <- function (wf, cp = NULL, col = topo.colors(256), addContours = FALSE, 
                       cscale = TRUE, ...) 
{
  if (cscale) {
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2]) * par("csi") * 2.54
    layout(matrix(c(2, 1), ncol = 2), widths = c(1, lcm(w)))
    par(las = 1)
    mar <- mar.orig
    mar[4] <- mar[2]
    mar[2] <- 1
    par(mar = mar) 
    thelist <- list(...)  
    findz <- which(names(thelist) == 'zlim')  
    if (length(findz) > 0 ) {   
      zlim <- thelist$zlim  
    }else{  
      zlim <- range(wf, finite = TRUE) #the original line  
    } 
    # end of my hack  
    levels <- seq(zlim[1], zlim[2], length = length(col))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1], col = col,  density = NA)
    axis(4)
    box()
    mar <- mar.orig
    mar[4] <- 0
    par(mar = mar)
  }
  if (is.null(cp)) {
    axis1 <- 1:nrow(wf)
    axis2 <- 1:ncol(wf)
  }
  else {
    axis1 <- ((1:nrow(wf)) - cp$xc)/cp$rx
    axis2 <- ((1:ncol(wf)) - cp$yc)/cp$ry
  }
  image(axis1, axis2, wf, col = col, asp = 1, xlab = "Years to Reintroduction", ylab = "Years to Fade-out",  ...)
  if (addContours) 
    contour(axis1, axis2, wf, add = TRUE)
}