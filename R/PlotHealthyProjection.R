PlotHealthyProjection <- function(timesteps, ages.init, alpha, gamma, samples.to.draw,
                                  tot.chains, joint.posterior.coda, posterior.names, reps)
{

  popsize.he <- log.lambda.s.he <- matrix(NA, ncol = timesteps, nrow = reps)
  
  for(i in 1:reps){
    he.project <- HealthyProjectFun(timesteps, sex.ratio = .5, ages.init, alpha = 0, gamma = 0, samples.to.draw = seq(500:1000), tot.chains = 3, joint.posterior.coda = ipm11.coda, posterior.names, intro.cost = 0)
    popsize.he[i, ] <- he.project$tot.pop.size 
    log.lambda.s.he[i, ] <- he.project$log.lambda.s
  }  
    
  he.quants <- which(popsize.he[, timesteps - 1] %in% as.numeric(quantile(popsize.he[, timesteps - 1], c(0.025, .25, .75, 0.975), type = 3)))
  he.med <- which(popsize.he[, timesteps - 1] %in% as.numeric(quantile(popsize.he[, timesteps - 1], .5, type = 3)))
  
  par(oma = c(1, 1, 1, 1), mar = c(4, 5, 2, 2))
  layout(matrix(c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3), nrow = 2, byrow = T))
  plot(popsize.he[1, -c(1)] ~ seq(2, timesteps), type = "l", ylim = c(0, 1200), xlab = "year", ylab = "population size")
  for(i in 2:reps){
    lines(popsize.he[i, -c(1)] ~ seq(2, timesteps), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in 1:3){
    lines(popsize.he[he.quants[j], -c(1, 2)] ~ seq(3, timesteps), type = "l", col = "black", lwd = 2)  
  }
  lines(popsize.he[he.med, -c(1, 2)] ~ seq(3, timesteps), type = "l", col = "red", lwd = 2)
  
  plot(log.lambda.s.he[1, -c(1, 2)] ~ seq(3, timesteps), type = "l", ylim = c(-.5, .5), xlab = "year", ylab = expression(paste("log(", lambda, "s)", sep = "")))
  for(i in 2:reps){
    lines(log.lambda.s.he[i, -c(1, 2)] ~ seq(3, timesteps), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  abline( h = 0, lty = 2, col = "red", lwd = 2)
  boxplot(as.vector(log.lambda.s.he[, -c(1:10)]), col = "grey80", ylim = c(-.5, .5), ylab = expression(paste("log(", lambda, "s)", sep = "")))
  abline(h = 0, lty = 2, col = "red", lwd = 2)
}