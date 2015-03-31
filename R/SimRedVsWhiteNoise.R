SimRedVsWhiteNoise <- function(timesteps,
                              reps,
                              alpha.range,
                              gamma.range,
                              alpha.steps,
                              gamma.steps,
                              samples.to.use.in,
                              pop.size, 
                              ages.init.in)
{
  ages.init <- pop.size * ages.init.in
  intro.cost.in <- .3
  
  popsize.gam.2.alph.05 <- matrix(NA, ncol = timesteps, nrow = reps)
  popsize.gam.2.alph.2 <- matrix(NA, ncol = timesteps, nrow = reps)
  popsize.gam.2.alph1 <- matrix(NA, ncol = timesteps, nrow = reps)
  
  unstr.popsize.gam.2.alph.05 <- matrix(NA, ncol = timesteps, nrow = reps)
  unstr.popsize.gam.2.alph.2 <- matrix(NA, ncol = timesteps, nrow = reps)
  unstr.popsize.gam.2.alph1 <- matrix(NA, ncol = timesteps, nrow = reps)
  
  unstr.intros.popsize.gam.2.alph.05 <- matrix(NA, ncol = timesteps, nrow = reps)
  unstr.intros.popsize.gam.2.alph.2 <- matrix(NA, ncol = timesteps, nrow = reps)
  unstr.intros.popsize.gam.2.alph1 <- matrix(NA, ncol = timesteps, nrow = reps)
  
  loglambda.gam.2.alph.05 <- matrix(NA, ncol = timesteps, nrow = reps)
  loglambda.gam.2.alph.2 <- matrix(NA, ncol = timesteps, nrow = reps)
  loglambda.gam.2.alph1 <- matrix(NA, ncol = timesteps, nrow = reps)
  
  disstat.gam.2.alph.05 <- matrix(NA, ncol = timesteps, nrow = reps)
  disstat.gam.2.alph.2 <- matrix(NA, ncol = timesteps, nrow = reps)
  disstat.gam.2.alph1 <- matrix(NA, ncol = timesteps, nrow = reps)
  
  disstat.intros.gam.2.alph.05 <- matrix(NA, ncol = timesteps, nrow = reps)
  disstat.intros.gam.2.alph.2 <- matrix(NA, ncol = timesteps, nrow = reps)
  disstat.intros.gam.2.alph1 <- matrix(NA, ncol = timesteps, nrow = reps)
  
  for(i in 1:reps){
    project.gam.2.alph.05 <- ProjectFun(johnson = F, timesteps = timesteps, intro.cost = intro.cost.in, sex.ratio = .5, ages.init = ages.init, alpha = .05, gamma = .05, samples.to.draw = samples.to.use.in, tot.chains = 3, joint.posterior.coda = ipm11.coda, posterior.names = posterior.names, fixed.start.time = T)
    project.gam.2.alph.2 <- ProjectFun(johnson = F, timesteps = timesteps, intro.cost = intro.cost.in, sex.ratio = .5, ages.init = ages.init, alpha = .2, gamma = .05, samples.to.draw = samples.to.use.in, tot.chains = 3, joint.posterior.coda = ipm11.coda, posterior.names = posterior.names, fixed.start.time = T)
    project.gam.2.alph1 <- ProjectFun(johnson = F, timesteps = timesteps, intro.cost = intro.cost.in, sex.ratio = .5, ages.init = ages.init, alpha = .9, gamma = .05, samples.to.draw = samples.to.use.in, tot.chains = 3, joint.posterior.coda = ipm11.coda, posterior.names = posterior.names, fixed.start.time = T)
    
    popsize.gam.2.alph.05[i, ] <- project.gam.2.alph.05$tot.pop.size
    popsize.gam.2.alph.2[i, ] <- project.gam.2.alph.2$tot.pop.size
    popsize.gam.2.alph1[i, ] <- project.gam.2.alph1$tot.pop.size
    
    #     loglambda.gam.2.alph.05[i, ] <- project.gam.2.alph.05$log.lambda.s
    #     loglambda.gam.2.alph.2[i, ] <- project.gam.2.alph.2$log.lambda.s
    #     loglambda.gam.2.alph1[i, ] <- project.gam.2.alph1$log.lambda.s
    
    # extract vector of disease statuses for each sim
    disstat.gam.2.alph.05[i, ] <- sample(project.gam.2.alph.05$disease.status, length(project.gam.2.alph.05$disease.status), replace = F)
    disstat.gam.2.alph.2[i, ] <- sample(project.gam.2.alph.2$disease.status, length(project.gam.2.alph.2$disease.status), replace = F)
    disstat.gam.2.alph1[i, ] <- sample(project.gam.2.alph1$disease.status, length(project.gam.2.alph1$disease.status), replace = F)
    
    unstr.project.gam.2.alph.05 <- UnstrProjectFun(johnson = F, timesteps = timesteps, 
                                                   intro.cost = intro.cost.in, sex.ratio = .5, ages.init = ages.init, 
                                                   alpha = .05, gamma = .05, samples.to.draw = samples.to.use.in, 
                                                   tot.chains = 3, joint.posterior.coda = ipm11.coda, 
                                                   posterior.names = posterior.names, fixed.start.time = T, states = disstat.gam.2.alph.05[i, ])
    
    unstr.project.gam.2.alph.2 <- UnstrProjectFun(johnson = F, timesteps = timesteps, 
                                                  intro.cost = intro.cost.in, sex.ratio = .5, ages.init = ages.init, 
                                                  alpha = .2, gamma = .05, samples.to.draw = samples.to.use.in, 
                                                  tot.chains = 3, joint.posterior.coda = ipm11.coda, 
                                                  posterior.names = posterior.names, fixed.start.time = T, states = disstat.gam.2.alph.2[i, ])
    
    unstr.project.gam.2.alph1 <- UnstrProjectFun(johnson = F, timesteps = timesteps, 
                                                 intro.cost = intro.cost.in, sex.ratio = .5, ages.init = ages.init, 
                                                 alpha = .05, gamma = .05, samples.to.draw = samples.to.use.in, 
                                                 tot.chains = 3, joint.posterior.coda = ipm11.coda, 
                                                 posterior.names = posterior.names, fixed.start.time = T, states = disstat.gam.2.alph1[i, ])
    
    unstr.popsize.gam.2.alph.05[i, ] <- unstr.project.gam.2.alph.05$tot.pop.size
    unstr.popsize.gam.2.alph.2[i, ] <- unstr.project.gam.2.alph.2$tot.pop.size
    unstr.popsize.gam.2.alph1[i, ] <- unstr.project.gam.2.alph1$tot.pop.size
    
    disstat.intros.gam.2.alph.05[i, ] <- ifelse(sample(project.gam.2.alph.05$disease.status, length(project.gam.2.alph.05$disease.status), replace = F) == "spillover", "spillover", "healthy")
    disstat.intros.gam.2.alph.2[i, ] <- ifelse(sample(project.gam.2.alph.2$disease.status, length(project.gam.2.alph.2$disease.status), replace = F) == "spillover", "spillover", "healthy")
    disstat.intros.gam.2.alph1[i, ] <- ifelse(sample(project.gam.2.alph1$disease.status, length(project.gam.2.alph1$disease.status), replace = F) == "spillover", "spillover", "healthy")
    
    unstr.intros.project.gam.2.alph.05 <- UnstrProjectFun(johnson = F, timesteps = timesteps, 
                                                          intro.cost = intro.cost.in, sex.ratio = .5, ages.init = ages.init, 
                                                          alpha = .05, gamma = .05, samples.to.draw = samples.to.use.in, 
                                                          tot.chains = 3, joint.posterior.coda = ipm11.coda, 
                                                          posterior.names = posterior.names, fixed.start.time = T, states = disstat.intros.gam.2.alph.05[i, ])
    
    unstr.intros.project.gam.2.alph.2 <- UnstrProjectFun(johnson = F, timesteps = timesteps, 
                                                         intro.cost = intro.cost.in, sex.ratio = .5, ages.init = ages.init, 
                                                         alpha = .2, gamma = .05, samples.to.draw = samples.to.use.in, 
                                                         tot.chains = 3, joint.posterior.coda = ipm11.coda, 
                                                         posterior.names = posterior.names, fixed.start.time = T, states = disstat.intros.gam.2.alph.2[i, ])
    
    unstr.intros.project.gam.2.alph1 <- UnstrProjectFun(johnson = F, timesteps = timesteps, 
                                                        intro.cost = intro.cost.in, sex.ratio = .5, ages.init = ages.init, 
                                                        alpha = .05, gamma = .05, samples.to.draw = samples.to.use.in, 
                                                        tot.chains = 3, joint.posterior.coda = ipm11.coda, 
                                                        posterior.names = posterior.names, fixed.start.time = T, states = disstat.intros.gam.2.alph1[i, ])
    
    unstr.intros.popsize.gam.2.alph.05[i, ] <- unstr.intros.project.gam.2.alph.05$tot.pop.size
    unstr.intros.popsize.gam.2.alph.2[i, ] <- unstr.intros.project.gam.2.alph.2$tot.pop.size
    unstr.intros.popsize.gam.2.alph1[i, ] <- unstr.intros.project.gam.2.alph1$tot.pop.size
    
    print(i)
  }  
  
  list.out <- list(popsize.gam.2.alph.05 = popsize.gam.2.alph.05,
                   popsize.gam.2.alph.2 = popsize.gam.2.alph.2,
                   popsize.gam.2.alph1 = popsize.gam.2.alph1,
                   unstr.popsize.gam.2.alph.05 = unstr.popsize.gam.2.alph.05,
                   unstr.popsize.gam.2.alph.2 = unstr.popsize.gam.2.alph.2,
                   unstr.popsize.gam.2.alph1 = unstr.popsize.gam.2.alph1,
                   unstr.intros.popsize.gam.2.alph.05 = unstr.intros.popsize.gam.2.alph.05,
                   unstr.intros.popsize.gam.2.alph.2 = unstr.intros.popsize.gam.2.alph.2,
                   unstr.intros.popsize.gam.2.alph1 = unstr.intros.popsize.gam.2.alph1
  )
  
  return(list.out)
}
