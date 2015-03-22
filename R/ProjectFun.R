ProjectFun <- function(johnson, timesteps, sex.ratio, ages.init, alpha, gamma, 
                       intro.cost, samples.to.draw, tot.chains, joint.posterior.coda, 
                       posterior.names, fixed.start.time)
{
  # This function projections a leslie matrix with a specified Markov environmental state-transition
  # structure forward through time.
  #
  # Args:
  # 
  # johnson = indicator for whether to use age-structure specified in Johnson et al.,
  #   2010, Ecological Applications
  # timesteps = integer, how many years to run the simulation over
  # sex.ratio = proportion of lambs that are male
  # ages.init = vector containing initial age structure
  # alpha = pathogen introduction rate
  # gamma = pathogen fade-out rate
  # intro.cost = increase in overall mortality rate in introduction years (over healthy)
  # samples.to.draw = steps in Markov chain(s) from which samples should be drawn
  # tot.chains = number of chains in joint.posterior.coda
  # joint.posterior.coda = posterior from IPM (as coda object)
  # posterior.names = vector of names associated with columns in joint.posterior.coda
  # fixed.start.time = indicator: does first introduction happen in a particular year
  #
  # Returns:
  # N = N, 
  # disease.status = disease.status, 
  # tot.pop.size = tot.pop.size, 
  # log.lambda.s = log.lambda.s
  
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
  out.list <- list(N = N, disease.status = disease.status, tot.pop.size = tot.pop.size, 
                   log.lambda.s = log.lambda.s)
  return(out.list)
}
