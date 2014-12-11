#-- simulation source functions --#

# function to update environmental state
update.status.fun <- function(alpha, gamma, current.state){
  current.state.new <- rep(NA, 1)
  if(current.state == "healthy"){
    gets.infected <- rbinom(1, 1, alpha)
    current.state.new[1] <- ifelse(gets.infected == 0, "healthy", "spillover")
  }
  else if(current.state == "spillover"){
    current.state.new[1] <- "infected"
  } else if(current.state == "infected"){
    fade.out <- rbinom(1, 1, gamma)
    current.state.new[1] <- ifelse(fade.out == 1, "healthy", "infected")
  }
  return(list(current.state.new = current.state.new))
}

# function to update Leslie matrix parameters
# need to expand Leslie to be 18x18... also, need individuals to age....
update.leslie.fun <- function(current.state, he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr){
  leslie.out <- rep(NA, 6, 6)
  if(current.state == "healthy"){
    repros <- c(0, rep(he.repro.post[sample(1:3000, 1), 2], 2), rep(he.repro.post[sample(1:3000, 1), 3], 4), rep(he.repro.post[sample(1:3000, 1), 4], 6), rep(he.repro.post[sample(1:3000, 1), 5], 5))    
    survs <- c(he.recr[sample(1:length(he.recr), 1)], rep(he.surv.post[sample(1:3000, 1), 2], 2), rep(he.surv.post[sample(1:3000, 1), 3], 4), rep(he.surv.post[sample(1:3000, 1), 4], 6), rep(he.surv.post[sample(1:3000, 1), 5], 4))
    leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs))))
  }
  else {
    repros <- c(0, rep(inf.repro.post[sample(1:3000, 1), 2], 2), rep(inf.repro.post[sample(1:3000, 1), 3], 4), rep(inf.repro.post[sample(1:3000, 1), 4], 6), rep(inf.repro.post[sample(1:3000, 1), 5], 5))    
    survs <- c(pn.recr[sample(1:length(pn.recr), 1)], rep(inf.surv.post[sample(1:3000, 1), 2], 2), rep(inf.surv.post[sample(1:3000, 1), 3], 4), rep(inf.surv.post[sample(1:3000, 1), 4], 6), rep(inf.surv.post[sample(1:3000, 1), 5], 4))
    leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs))))
  }
  return(leslie)
}


project.fun <- function(timesteps, ages.init, alpha, gamma, he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr){
  N <- matrix(NA, nrow = length(ages.init), ncol = timesteps)
  N[, 1] <- ages.init
  tot.pop.size <- log.lambda.s <- rep(NA, length = timesteps)
  disease.status <- rep(NA, timesteps)
  disease.status[1:11] <- c(rep("healthy", 10), "spillover")
  for(i in 1:10){
    new.leslie <- update.leslie.fun(current.state = disease.status[i], he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
    N[, i + 1] <- round(t(N[ , i]) %*% new.leslie) 
    tot.pop.size[i] <- sum(N[ , i])
    if(i == 1){
      log.lambda.s[i] <- NA
    } else {
      log.lambda.s[i] <- ifelse(tot.pop.size[i] == 0, NA, log(tot.pop.size[i] / tot.pop.size[i - 1]))
    }
  }
  for(i in 11:(timesteps - 1)){
    new.leslie <- update.leslie.fun(current.state = disease.status[i], he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
    N[, i + 1] <- round(t(N[ , i]) %*% new.leslie)
    disease.status[i + 1] <- update.status.fun(alpha, gamma, current.state = disease.status[i])$current.state.new[1]
    tot.pop.size[i] <- sum(N[ , i])
    if(i == 1){
      log.lambda.s[i] <- NA
    } else {
      log.lambda.s[i] <- ifelse(tot.pop.size[i] == 0, NA, log(tot.pop.size[i] / tot.pop.size[i - 1]))
    }
  }
  out.list <- list(N = N, disease.status = disease.status, tot.pop.size = tot.pop.size, log.lambda.s = log.lambda.s)
  return(out.list)
}

# projection function in absence of disease
healthy.project.fun <- function(timesteps, ages.init, alpha, gamma, he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr){
  N <- matrix(NA, nrow = length(ages.init), ncol = timesteps)
  N[, 1] <- ages.init
  tot.pop.size <- log.lambda.s <- rep(NA, length = timesteps)
  disease.status <- rep("healthy", timesteps)
  for(i in 1:(timesteps - 1)){
    new.leslie <- update.leslie.fun(current.state = disease.status[i], he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
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