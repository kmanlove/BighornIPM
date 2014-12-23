#---------------------------------------------#
#-- Code for recruitment vs. adult survival --#
#---------------------------------------------#
# require(rjags)
# require(runjags)
require(popbio)
require(Matrix)
source("~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Code/BighornIPM_GIT/BighornSimSourceFunctions.R")
ipm11.coda <- dget("~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Posteriors/IPM/22Dec2014")


#--------------------------------------------------#
#-- Figure 1. IPM posterior estimates -------------#
#--------------------------------------------------#
coda.summary.obj.11 <- summary(ipm11.coda)
row.names(coda.summary.obj.11[[2]])

beta.posts.adsurv <- coda.summary.obj.11[[2]][1:18, ]
beta.posts.repro <- coda.summary.obj.11[[2]][19:36, ]
beta.posts.wean <- coda.summary.obj.11[[2]][37:54, ]

# reminder of age-structure for IPM:
age.class.ind <- c(1, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6) 

# restructure posteriors so that age-classes are represented as many times as they contain years
he.adsurv <- beta.posts.adsurv[c(1, 4, 7, 10, 13, 16), ]
inf.adsurv <- beta.posts.adsurv[c(1, 4, 7, 10, 13, 16) + 1, ]
he.repro <- beta.posts.repro[c(1, 4, 7, 10, 13, 16), ]
inf.repro <- beta.posts.repro[c(1, 4, 7, 10, 13, 16) + 1, ]
he.wean <- beta.posts.wean[c(1, 4, 7, 10, 13, 16), ]
inf.wean <- beta.posts.wean[c(1, 4, 7, 10, 13, 16) + 1, ]
he.adsurv.2 <- he.adsurv[age.class.ind, ]
inf.adsurv.2 <- inf.adsurv[age.class.ind, ]
he.repro.2 <- he.repro[age.class.ind, ]
inf.repro.2 <- inf.repro[age.class.ind, ]
he.wean.2 <- he.wean[age.class.ind, ]
inf.wean.2 <- inf.wean[age.class.ind, ]

# plot betas out
plot.cols <- c("white", "black", "red")
#svg("~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Plots/FigsV1/Posteriors_22Dec2014.svg", width = 3, height = 2, pointsize = 8)
par(mfrow = c(1, 3), mex = 1, oma = c(2, 2, 0, 0), mar = c(4, 5, 1, 1))
plot(-1, -1, ylim = c(0, 1), xlim = c(2, 19), ylab = "Probability of survival", xlab = "Age (years)")
for(i in 2:19){
  segments(x0 = i, x1 = i, y0 = (exp(he.adsurv.2[i, 1])) / (1 + exp(he.adsurv.2[i, 1])), y1 = exp(he.adsurv.2[i, 5]) / (1 + exp(he.adsurv.2[i, 5])), col = "black", lty = 2, lwd = 2)
  segments(x0 = i - 0.25, x1 = i + 0.25, y0 = (exp(he.adsurv.2[i, 3])) / (1 + exp(he.adsurv.2[i, 3])), y1 = exp(he.adsurv.2[i, 3]) / (1 + exp(he.adsurv.2[i, 3])), col = "black", lty = 2, lwd = 2)
  segments(x0 = i + .25, x1 = i + .25, y0 = (exp(inf.adsurv.2[i, 1])) / (1 + exp(inf.adsurv.2[i, 1])), y1 = exp(inf.adsurv.2[i, 5]) / (1 + exp(inf.adsurv.2[i, 5])), col = "red", lty = 1, lwd = 2)
  segments(x0 = i + .25 - 0.25, x1 = i + .25 + 0.25, y0 = (exp(inf.adsurv.2[i, 3])) / (1 + exp(inf.adsurv.2[i, 3])), y1 = exp(inf.adsurv.2[i, 3]) / (1 + exp(inf.adsurv.2[i, 3])), col = "red", lty = 1, lwd = 2)
}
leg.text <- c("Healthy", "Diseased")
legend("bottomleft", leg.text, col = c("black", "red"), lty = c(2, 1), bty = "n", lwd = c(2, 2))

plot(-1, -1, ylim = c(0, 1), xlim = c(2, 19), ylab = "Probability of reproducing", xlab = "Age (years)")
for(i in 2:19){
  segments(x0 = i, x1 = i, y0 = (exp(he.repro.2[i, 1])) / (1 + exp(he.repro.2[i, 1])), y1 = exp(he.repro.2[i, 5]) / (1 + exp(he.repro.2[i, 5])), col = "black", lty = 2, lwd = 2)
  segments(x0 = i - 0.25, x1 = i + 0.25, y0 = (exp(he.repro.2[i, 3])) / (1 + exp(he.repro.2[i, 3])), y1 = exp(he.repro.2[i, 3]) / (1 + exp(he.repro.2[i, 3])), col = "black", lty = 2, lwd = 2)
  segments(x0 = i + .25, x1 = i + .25, y0 = (exp(inf.repro.2[i, 1])) / (1 + exp(inf.repro.2[i, 1])), y1 = exp(inf.repro.2[i, 5]) / (1 + exp(inf.repro.2[i, 5])), col = "red", lty = 1, lwd = 2)
  segments(x0 = i + .25 - 0.25, x1 = i + .25 + 0.25, y0 = (exp(inf.repro.2[i, 3])) / (1 + exp(inf.repro.2[i, 3])), y1 = exp(inf.repro.2[i, 3]) / (1 + exp(inf.repro.2[i, 3])), col = "red", lty = 1, lwd = 2)
}

plot(-1, -1, ylim = c(0, 1), xlim = c(2, 19), ylab = "Probability of weaning", xlab = "Age (years)")
for(i in 2:19){
  segments(x0 = i, x1 = i, y0 = (exp(he.wean.2[i, 1])) / (1 + exp(he.wean.2[i, 1])), y1 = exp(he.wean.2[i, 5]) / (1 + exp(he.wean.2[i, 5])), col = "black", lty = 2, lwd = 2)
  segments(x0 = i - 0.25, x1 = i + 0.25, y0 = (exp(he.wean.2[i, 3])) / (1 + exp(he.wean.2[i, 3])), y1 = exp(he.wean.2[i, 3]) / (1 + exp(he.wean.2[i, 3])), col = "black", lty = 2, lwd = 2)
  segments(x0 = i + .25, x1 = i + .25, y0 = (exp(inf.wean.2[i, 1])) / (1 + exp(inf.wean.2[i, 1])), y1 = exp(inf.wean.2[i, 5]) / (1 + exp(inf.wean.2[i, 5])), col = "red", lty = 1, lwd = 2)
  segments(x0 = i + .25 - 0.25, x1 = i + .25 + 0.25, y0 = (exp(inf.wean.2[i, 3])) / (1 + exp(inf.wean.2[i, 3])), y1 = exp(inf.wean.2[i, 3]) / (1 + exp(inf.wean.2[i, 3])), col = "red", lty = 1, lwd = 2)
}
#dev.off()

#---------------------------------------------#
#-- build population projection structures ---#
#---------------------------------------------#
timesteps <- 60
reps <- 100
ages.init <- c(100, 20, rep(10, 8), rep(5, 5), rep(3, 4))
alpha <- .1
gamma <- 1
samples.to.draw <- seq(500:1000)
tot.chains <- 3
joint.posterior.coda <- ipm11.coda

posterior.names <- c("beta.adsurv.1.1", "beta.adsurv.2.1", "beta.adsurv.3.1", "beta.adsurv.1.2", 
                     "beta.adsurv.2.2", "beta.adsurv.3.2", "beta.adsurv.1.3", "beta.adsurv.2.3",
                     "beta.adsurv.3.3", "beta.adsurv.1.4", "beta.adsurv.2.4", "beta.adsurv.3.4",  
                     "beta.adsurv.1.5", "beta.adsurv.2.5", "beta.adsurv.3.5", "beta.adsurv.1.6", "beta.adsurv.2.6", "beta.adsurv.3.6",
                     "beta.repro.1.1", "beta.repro.2.1",
                     "beta.repro.3.1", "beta.repro.1.2", "beta.repro.2.2",  "beta.repro.3.2", 
                     "beta.repro.1.3", "beta.repro.2.3", "beta.repro.3.3", "beta.repro.1.4",  
                     "beta.repro.2.4", "beta.repro.3.4", "beta.repro.1.5", "beta.repro.2.5",  
                     "beta.repro.3.5", "beta.repro.1.6", "beta.repro.2.6",  
                     "beta.repro.3.6", "beta.wean.1.1", "beta.wean.2.1", "beta.wean.3.1",  
                     "beta.wean.1.2",  "beta.wean.2.2", "beta.wean.3.2", "beta.wean.1.3",  
                     "beta.wean.2.3", "beta.wean.3.3", "beta.wean.1.4", "beta.wean.2.4", 
                     "beta.wean.3.4", "beta.wean.1.5", "beta.wean.2.5", "beta.wean.3.5", "beta.wean.1.6", "beta.wean.2.6", "beta.wean.3.6")

healthy.test <- healthy.project.fun(timesteps, sex.ratio = .6, ages.init, alpha, gamma, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)
  
reps <- 100
popsize.he <- log.lambda.s.he <- matrix(NA, ncol = timesteps, nrow = reps)

for(i in 1:reps){
  he.project <- healthy.project.fun(timesteps, sex.ratio = .6, ages.init, alpha, gamma, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)
  popsize.he[i, ] <- he.project$tot.pop.size 
  log.lambda.s.he[i, ] <- he.project$log.lambda.s
}  

#write.csv(popsize.he, "~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Simulations/Healthy/HealthyPopsize_30Sept2014.csv", row.names = F)
#popsize.he <- as.matrix(read.csv("~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Simulations/Healthy/HealthyPopsize_30Sept2014.csv"))

he.quants <- which(popsize.he[, 40] %in% as.numeric(quantile(popsize.he[, 40], c(0.025, 0.975), type = 3)))
he.med <- which(popsize.he[, 40] %in% as.numeric(quantile(popsize.he[, 40], .5, type = 3)))

#par(mfrow = c(2, 1))
layout(matrix(c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3), nrow = 2, byrow = T))
plot(popsize.he[1, -c(1)] ~ seq(2, timesteps), type = "l", ylim = c(0, 3000), xlab = "year", ylab = "population size")
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

#--------------------------------------------------------------------------#
#-- Exploration of lambda and age structure in each environmental state ---#
#--------------------------------------------------------------------------#

healthy.leslie.list <- infected.leslie.list <- vector("list", length = 1000)
healthy.eigenval1 <- infected.eigenval1 <- rep(NA, 1000)
juvsurv.elast.he <- juvsurv.elast.inf <- adsurv.elast.he <- adsurv.elast.inf <- fecund.elast.he <- fecund.elast.inf <- rep(NA, 1000)
age.struct.he <- age.struct.inf <- matrix(NA, nrow = 19, ncol = 1000)

for(i in 1:1000){
    healthy.leslie.list[[i]] <- update.leslie.fun(current.state = "healthy", samples.to.draw = 500:1000, tot.chains = 3, joint.posterior.coda, posterior.names = posterior.names)
    infected.leslie.list[[i]] <- update.leslie.fun(current.state = "infected", samples.to.draw = 500:1000, tot.chains = 3, joint.posterior.coda, posterior.names = posterior.names)
    healthy.eigenval1[i] <- Re(eigen(healthy.leslie.list[[i]])$values[1]) # strip off only real part of eigenvalue
    infected.eigenval1[i] <- Re(eigen(infected.leslie.list[[i]])$values[1])
    he.eigen.rescale <- sum(Re(eigen(healthy.leslie.list[[i]])$vectors[, 1]))
    inf.eigen.rescale <- sum(Re(eigen(infected.leslie.list[[i]])$vectors[, 1]))
    age.struct.he[, i] <- Re(eigen(healthy.leslie.list[[i]])$vectors[, 1]) / he.eigen.rescale
    age.struct.inf[, i] <- Re(eigen(infected.leslie.list[[i]])$vectors[, 1]) / inf.eigen.rescale
    
    he.sens <- sensitivity(healthy.leslie.list[[i]])
    he.elast <- (1 / Re(eigen(healthy.leslie.list[[i]])$values[1])) * he.sens * healthy.leslie.list[[i]]
    fecund.elast.he[i] <- sum(he.elast[1, ])
    juvsurv.elast.he[i] <- sum(he.elast[2, 1], he.elast[3, 2])
    adsurv.elast.he[i] <- sum(he.elast[4, 3], he.elast[5, 4], he.elast[6, 5], he.elast[7, 6], he.elast[8, 7], he.elast[9, 8], he.elast[10, 9], he.elast[11, 10], he.elast[12, 11], he.elast[13, 12], he.elast[14, 13], he.elast[15, 14], he.elast[16, 15], he.elast[17, 16], he.elast[18, 17])
    
    inf.sens <- sensitivity(infected.leslie.list[[i]])
    inf.elast <- (1 / Re(eigen(infected.leslie.list[[i]])$values[1])) * inf.sens * infected.leslie.list[[i]]
    fecund.elast.inf[i] <- sum(inf.elast[1, ])
    juvsurv.elast.inf[i] <- sum(inf.elast[2, 1], inf.elast[3, 2])
    adsurv.elast.inf[i] <- sum(inf.elast[4, 3], inf.elast[5, 4], inf.elast[6, 5], inf.elast[7, 6], inf.elast[8, 7], inf.elast[9, 8], inf.elast[10, 9], inf.elast[11, 10], inf.elast[12, 11], inf.elast[13, 12], inf.elast[14, 13], inf.elast[15, 14], inf.elast[16, 15], inf.elast[17, 16], inf.elast[18, 17])
}
  
  # 95% intervals on lambda 
  healthy.cred<- quantile(healthy.eigenval1, c(0.25, 0.5, 0.75))
  spillover.cred<- quantile(spillover.eigenval1, c(0.25, 0.5, 0.75))
  infected.cred<- quantile(infected.eigenval1, c(0.25, 0.5, 0.75))
  
  # 95% intervals for fecundity, juvenile survival, and adult survival in each environment (healthy and infected)
  fecund.he <- quantile(fecund.elast.he, c(0.05, 0.5, 0.975))
  fecund.inf <- quantile(fecund.elast.inf, c(.05, 0.5, 0.975))
  juvsurv.he <- quantile(juvsurv.elast.he, c(0.05, 0.5, 0.975))
  juvsurv.inf <- quantile(juvsurv.elast.inf, c(0.05, 0.5, 0.975))
  adsurv.he <- quantile(adsurv.elast.he, c(0.05, 0.5, 0.975))
  adsurv.inf <- quantile(adsurv.elast.inf, c(0.05, 0.5, 0.975))
  
  par(mfrow = c(1, 1), mar = c(4, 8, 2, 2), las = 1, cex.lab = 1.0)
  plot(c(1, 1) ~ c(fecund.he[1], fecund.he[3]), lty = 1, xlim = c(0, 1), ylim = c(0, 7), type = "l", xlab = expression(paste("Elasticity of ", lambda, " to rate", sep = "")), ylab = "", yaxt = "n", lwd = 2)
  lines(c(2, 2) ~ c(fecund.inf[1], fecund.inf[3]), lty = 2, col = "red", lwd = 2)
  lines(c(3, 3) ~ c(juvsurv.he[1], juvsurv.he[3]), lty = 1, col = "black", lwd = 2)
  lines(c(4, 4) ~ c(juvsurv.inf[1], juvsurv.inf[3]), lty = 2, col = "red", lwd = 2)
  lines(c(5, 5) ~ c(adsurv.he[1], adsurv.he[3]), lty = 1, col = "black", lwd = 2)
  lines(c(6, 6) ~ c(adsurv.inf[1], adsurv.inf[3]), lty = 2, col = "red", lwd = 2)
  lines(c(0.75, 1.25) ~ c(fecund.he[2], fecund.he[2]), lty = 1, col = "black", lwd = 2)
  lines(c(1.75, 2.25) ~ c(fecund.inf[2], fecund.inf[2]), lty = 2, col = "red", lwd = 2)
  lines(c(2.75, 3.25) ~ c(juvsurv.he[2], juvsurv.he[2]), lty = 1, col = "black", lwd = 2)
  lines(c(3.75, 4.25) ~ c(juvsurv.inf[2], juvsurv.inf[2]), lty = 2, col = "red", lwd = 2)
  lines(c(4.75, 5.25) ~ c(adsurv.he[2], adsurv.he[2]), lty = 1, col = "black", lwd = 2)
  lines(c(5.75, 6.25) ~ c(adsurv.inf[2], adsurv.inf[2]), lty = 2, col = "red", lwd = 2)
  axis(side = 2, at = c(1:6), cex.axis = 1.0, labels = c("Fecundity (he)", "Fecundity (pers)", "Juvenile survival (he)", "Juvenile surv (pers)", "Adult survival (he)", "Adult survival (pers)"))
  
  k <- hist(c(healthy.eigenval1, infected.eigenval1), breaks = 15, plot = F)
  # par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.5)
  hist(healthy.eigenval1, xlim = c(min(min(healthy.eigenval1), min(infected.eigenval1)), max(max(healthy.eigenval1), max(infected.eigenval1))), breaks = k$breaks, xlab = expression(lambda), main = "Healthy", col = "grey80")
  abline(v = 1, col = "red", lwd = 3)
  abline(v = 1, col = "red", lwd = 3)
  hist(infected.eigenval1, xlim = c(min(min(healthy.eigenval1), min(infected.eigenval1)), max(max(healthy.eigenval1), max(infected.eigenval1))), breaks = k$breaks, ylab = "", xlab = expression(lambda), main = "Persistence", col = "grey80")
  abline(v = 1, col = "red", lwd = 3)

  # environ-state-specific elasticities
  he.sens <- sensitivity(healthy.leslie.list[[1]])
  he.elast <- (1 / Re(eigen(healthy.leslie.list[[1]])$values[1])) * he.sens * healthy.leslie.list[[1]]
  round(he.sens, 2)
  round(he.elast, 2)
  fecund.elast.he <- sum(he.elast[1, ])
  juvsurv.elast.he <- sum(he.elast[2, 1], he.elast[3, 2])
  adsurv.elast.he <- sum(he.elast[4, 3], he.elast[5, 4], he.elast[6, 5], he.elast[7, 6], he.elast[8, 7], he.elast[9, 8], he.elast[10, 9], he.elast[11, 10], he.elast[12, 11], he.elast[13, 12], he.elast[14, 13], he.elast[15, 14], he.elast[16, 15], he.elast[17, 16], he.elast[18, 17])
  
  inf.sens <- sensitivity(infected.leslie.list[[1]])
  inf.elast <- (1 / Re(eigen(infected.leslie.list[[1]])$values[1])) * inf.sens * infected.leslie.list[[1]]
  round(inf.sens, 2)
  round(inf.elast, 2)
  fecund.elast.inf <- sum(inf.elast[1, ])
  juvsurv.elast.inf <- sum(inf.elast[2, 1], inf.elast[3, 2])
  #sum(inf.elast[5, 4], inf.elast[6, 5], inf.elast[7, 6], inf.elast[8, 7], inf.elast[9, 8], inf.elast[10, 9])
  adsurv.elast.inf <- sum(inf.elast[4, 3], inf.elast[5, 4], inf.elast[6, 5], inf.elast[7, 6], inf.elast[8, 7], inf.elast[9, 8], inf.elast[10, 9], inf.elast[11, 10], inf.elast[12, 11], inf.elast[13, 12], inf.elast[14, 13], inf.elast[15, 14], inf.elast[16, 15], inf.elast[17, 16], inf.elast[18, 17])
  
  # intervals on age-structure in healthy and infected environments 
  age.bounds.he <- age.bounds.inf <- matrix(NA, nrow = 18, ncol = 3)
  for(i in 1:18){
    age.bounds.he[i, ] <- quantile(abs(age.struct.he[i, ]), c(0.025, 0.5, 0.975))
    age.bounds.inf[i, ] <- quantile(abs(age.struct.inf[i, ]), c(0.025, 0.5, 0.975))
  }
  
  par(mfrow = c(1, 2), cex.axis = .8, cex.lab = 1, oma = c(2, 0, 0, 0), mar = c(2, 6, 1, 1), las = 0)
  plot(density(healthy.eigenval1), main = "", ylim = c(0, 12), xlim = c(0.7, 1.12), lwd = 2, lty = 2, )
  lines(density(infected.eigenval1), col = "red", lwd = 2, lty = 1)
  mtext(side = 1, line = 3, outer = F, cex = 1, expression(paste(lambda)))
  plot(age.bounds.he[1, c(1, 3)] ~ c(1, 1), lty = 2, type = "l", xlim = c(0.5, 18.5), ylim = c(0, .2), lwd = 2, xlab = "age (years)", ylab = "Proportion of \n population")
  segments(x0 = 1 + .3, x1 = 1 + .3, y0 = age.bounds.inf[1, 1], y1 = age.bounds.inf[1, 2], lwd = 2, col = "red")
  for(i in 2: dim(age.bounds.he)[1]){
    segments(x0 = i, x1 = i, y0 = age.bounds.he[i, 1], y1 = age.bounds.he[i, 2], lwd = 2, lty = 2)
    segments(x0 = i + .3, x1 = i + .3, y0 = age.bounds.inf[i, 1], y1 = age.bounds.inf[i, 2], lwd = 2, col = "red")
  }
  leg.text <- c("healthy", "infected")
  legend("topright", leg.text, col = c("black", "red"), lwd = c(2, 2), bty = "n")
  mtext(side = 1, line = 3, outer = F, cex = 1, expression("age"))


#--------------------------------------------------------#
#-- Sensitivies and Elasticities using Vec-permutation --#
#--------------------------------------------------------#

vec.permut.test0 <- vec.permut.fun(reps = 100, alpha = 0.9, gamma = 0.1)
vec.permut.test1 <- vec.permut.fun(reps = 100, alpha = 0.2, gamma = 0.1)
vec.permut.test2 <- vec.permut.fun(reps = 100, alpha = 0.1, gamma = 0.1)
vec.permut.test3 <- vec.permut.fun(reps = 100, alpha = 0.05, gamma = 0.1)

vec.permut.fixedalpha.test0 <- vec.permut.fun(reps = 100, alpha = 0.1, gamma = 0.9)
vec.permut.fixedalpha.test1 <- vec.permut.fun(reps = 100, alpha = 0.1, gamma = 0.2)
vec.permut.fixedalpha.test2 <- vec.permut.fun(reps = 100, alpha = 0.1, gamma = 0.1)
vec.permut.fixedalpha.test3 <- vec.permut.fun(reps = 100, alpha = 0.1, gamma = 0.05)

fade.out.elasts <- quantile(vec.permut.test$fade.out.elast, c(0.025, 0.975))
intro.elasts <- quantile(vec.permut.test$intro.elast, c(0.025, 0.975))

par(mfrow = c(4, 4), cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.4)
hist(vec.permut.test3$fade.out.elast, ylab = "frequency", xlab = "fade-out elasticity", col = "grey80", main = "E(time to reintro) = 20yrs")
hist(vec.permut.test2$fade.out.elast, ylab = "frequency", xlab = "fade-out elasticity", col = "grey80", main = "E(time to reintro) = 10yrs")
hist(vec.permut.test1$fade.out.elast, ylab = "frequency", xlab = "fade-out elasticity", col = "grey80", main = "E(time to reintro) = 5yrs")
hist(vec.permut.test0$fade.out.elast, ylab = "frequency", xlab = "fade-out elasticity", col = "grey80", main = "E(time to reintro) = 1yrs")
hist(vec.permut.test3$intro.elast, main = "", ylab = "frequency", xlab = "introduction elasticity", col = "grey80")
hist(vec.permut.test2$intro.elast, main = "", ylab = "frequency", xlab = "introduction elasticity", col = "grey80")
hist(vec.permut.test1$intro.elast, main = "", ylab = "frequency", xlab = "introduction elasticity", col = "grey80")
hist(vec.permut.test0$intro.elast, main = "", ylab = "frequency", xlab = "introduction elasticity", col = "grey80")

hist(vec.permut.fixedalpha.test3$fade.out.elast, ylab = "frequency", xlab = "fade-out elasticity", col = "grey80", main = "E(time to fadeout) = 20yrs")
hist(vec.permut.fixedalpha.test2$fade.out.elast, ylab = "frequency", xlab = "fade-out elasticity", col = "grey80", main = "E(time to fadeout) = 10yrs")
hist(vec.permut.fixedalpha.test1$fade.out.elast, ylab = "frequency", xlab = "fade-out elasticity", col = "grey80", main = "E(time to fadeout) = 5yrs")
hist(vec.permut.fixedalpha.test0$fade.out.elast, ylab = "frequency", xlab = "fade-out elasticity", col = "grey80", main = "E(time to fadeout) = 1yrs")
hist(vec.permut.fixedalpha.test3$intro.elast, main = "", ylab = "frequency", xlab = "introduction elasticity", col = "grey80")
hist(vec.permut.fixedalpha.test2$intro.elast, main = "", ylab = "frequency", xlab = "introduction elasticity", col = "grey80")
hist(vec.permut.fixedalpha.test1$intro.elast, main = "", ylab = "frequency", xlab = "introduction elasticity", col = "grey80")
hist(vec.permut.fixedalpha.test0$intro.elast, main = "", ylab = "frequency", xlab = "introduction elasticity", col = "grey80")

# Sensitivity (demographic transitions)
S_b <- S_a %*% t(P) %*% t(M) %*% P
round(S_b, 2)[1:20, 1:20]
E_b <- (1 / .92) * B * S_b # regular * because it's a Hadamard production
round(E_b, 2)



#-------------------------------------------------------#
#-- costs of altering alpha/gamma on 30-year pop size --#
#-------------------------------------------------------#
#-- better to think about gammas as expected time to fade-out.... --#
#-- F(x) = 1/beta * (exp (-x / beta))
#-- where beta = expected survival of system
#-- so, gammas of interest are 1, 1/2, 1/3, 1/4

timesteps <- 40
reps <- 100
yearstofadeout.seq <- seq(1, 15)
gamma.seq <- 1 / yearstofadeout.seq
popsize.30.gamma.mat <- matrix(NA, reps, length(gamma.seq))
popsize.30.gamma.alpha.33.mat <- matrix(NA, reps, length(gamma.seq))

for(i in 1:reps){
  for(j in 1:length(gamma.seq)){
    popsize.30.gamma.mat[i, j] <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = .1, gamma = gamma.seq[j], samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)$tot.pop.size[timesteps - 1]
    popsize.30.gamma.alpha.33.mat[i, j] <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = .33, gamma = gamma.seq[j], samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)$tot.pop.size[timesteps - 1]
  }
}

popsize.30.quants <- matrix(NA, length(gamma.seq), 5)
popsize.30.alpha.33.quants <- matrix(NA, length(gamma.seq), 5)
for(j in 1:length(gamma.seq)){
  popsize.30.quants[j, ] <- quantile(popsize.30.gamma.mat[, j], c(0.025, 0.25, 0.5, 0.75, 0.975))
  popsize.30.alpha.33.quants[j, ] <- quantile(popsize.30.gamma.alpha.33.mat[, j], c(0.025, 0.25, 0.5, 0.75, 0.975))
}

par(mfrow = c(1, 1))
plot(x = 0, y = 0, xlim = c(0, 16), ylim = c(0, 500), xlab = "Expected years to fade-out (alpha = 0.1)", ylab = "Simulated Pop Size After 30 years")
for(i in 1:length(yearstofadeout.seq)){
  segments(x0 = yearstofadeout.seq[i], x1 = yearstofadeout.seq[i], y0 = popsize.30.quants[i, 1], y1 = popsize.30.quants[i, 5])
  segments(x0 = yearstofadeout.seq[i] - 0.1, x1 = yearstofadeout.seq[i] + 0.1, y0 = popsize.30.quants[i, 3], y1 = popsize.30.quants[i, 3])
  segments(x0 = yearstofadeout.seq[i] + .2, x1 = yearstofadeout.seq[i] + .2, y0 = popsize.30.alpha.33.quants[i, 1], y1 = popsize.30.alpha.33.quants[i, 5], col = "red", lwd = 2)
  segments(x0 = yearstofadeout.seq[i] + .1, x1 = yearstofadeout.seq[i] + .3, y0 = popsize.30.alpha.33.quants[i, 3], y1 = popsize.30.alpha.33.quants[i, 3], col = "red", lwd = 2)
}

#-- Same thing, but over levels of alpha, with gamma fixed at e(fade-out) = 3 and e(fade-out) = 10
timesteps <- 40
reps <- 100
yearstoreintro<- seq(1, 15)
alpha.seq <- 1 / yearstoreintro
popsize.30.alpha.gamma.33.mat <- matrix(NA, reps, length(alpha.seq))
popsize.30.alpha.gamma.1.mat <- matrix(NA, reps, length(alpha.seq))

for(i in 1:reps){
  for(j in 1:length(alpha.seq)){
    popsize.30.alpha.gamma.33.mat[i, j] <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = alpha.seq[j], gamma = .33, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)$tot.pop.size[timesteps - 1]
    popsize.30.alpha.gamma.1.mat[i, j] <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = alpha.seq[j], gamma = .1, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)$tot.pop.size[timesteps - 1]
  }
}

popsize.alpha.gamma.33.30.quants <- matrix(NA, length(alpha.seq), 5)
popsize.alpha.gamma.1.30.quants <- matrix(NA, length(alpha.seq), 5)
for(j in 1:length(gamma.seq)){
  popsize.alpha.gamma.33.30.quants[j, ] <- quantile(popsize.30.alpha.gamma.33.mat[, j], c(0.025, 0.25, 0.5, 0.75, 0.975))
  popsize.alpha.gamma.1.30.quants[j, ] <- quantile(popsize.30.alpha.gamma.1.mat[, j], c(0.025, 0.25, 0.5, 0.75, 0.975))
}

par(mfrow = c(1, 1))
plot(x = 0, y = 0, xlim = c(0, 16), ylim = c(0, 500), xlab = "Expected years to reintroduction", ylab = "Simulated Pop Size After 30 years")
for(i in 1:length(yearstoreintro)){
  segments(x0 = yearstoreintro[i], x1 = yearstoreintro[i], y0 = popsize.alpha.gamma.33.30.quants[i, 1], y1 = popsize.alpha.gamma.33.30.quants[i, 5])
  segments(x0 = yearstoreintro[i] - 0.1, x1 = yearstoreintro[i] + 0.1, y0 = popsize.alpha.gamma.33.30.quants[i, 3], y1 = popsize.alpha.gamma.33.30.quants[i, 3])
  segments(x0 = yearstoreintro[i] + .2, x1 = yearstoreintro[i] + .2, y0 = popsize.alpha.gamma.1.30.quants[i, 1], y1 = popsize.alpha.gamma.1.30.quants[i, 5], col = "red")
  segments(x0 = yearstoreintro[i] + .1, x1 = yearstoreintro[i] + .3, y0 = popsize.alpha.gamma.1.30.quants[i, 3], y1 = popsize.alpha.gamma.1.30.quants[i, 3], col = "red")
}

par(mfrow = c(1, 2), oma = c(0, 0, 0, 0), mar = c(4, 4, 1, 1))
plot(x = 0, y = 0, xlim = c(1, 16), ylim = c(0, 500), xlab = "Expected years to fade-out", ylab = "Pop Size After 30 years")
lines(popsize.30.quants[, 3] ~ yearstoreintro, col = "black", lty = 2)
lines(popsize.30.alpha.33.quants[, 3] ~ yearstoreintro, col = "red")
for(i in 1:length(yearstofadeout.seq)){
  segments(x0 = yearstofadeout.seq[i], x1 = yearstofadeout.seq[i], y0 = popsize.30.quants[i, 1], y1 = popsize.30.quants[i, 5], lty = 2)
  segments(x0 = yearstofadeout.seq[i] - 0.1, x1 = yearstofadeout.seq[i] + 0.1, y0 = popsize.30.quants[i, 3], y1 = popsize.30.quants[i, 3])
  segments(x0 = yearstofadeout.seq[i] + .3, x1 = yearstofadeout.seq[i] + .3, y0 = popsize.30.alpha.33.quants[i, 1], y1 = popsize.30.alpha.33.quants[i, 5], col = "red", lwd = 1)
  segments(x0 = yearstofadeout.seq[i] + .2, x1 = yearstofadeout.seq[i] + .4, y0 = popsize.30.alpha.33.quants[i, 3], y1 = popsize.30.alpha.33.quants[i, 3], col = "red", lwd = 1)
}
leg.text <- c("alpha = 0.33", "alpha = 0.1")
legend("topright", leg.text, col = c("black", "red"), lwd = c(1, 1), lty = c(2, 1), bty = "n", cex = .8)

plot(x = 0, y = 0, xlim = c(1, 16), ylim = c(0, 500), xlab = "Expected years to reintroduction", ylab = "Pop Size After 30 years")
lines(popsize.alpha.gamma.33.30.quants[, 3] ~ yearstoreintro, col = "black", lty = 2)
lines(popsize.alpha.gamma.1.30.quants[, 3] ~ yearstoreintro, col = "red")
for(i in 1:length(yearstoreintro)){
  segments(x0 = yearstoreintro[i], x1 = yearstoreintro[i], y0 = popsize.alpha.gamma.33.30.quants[i, 1], y1 = popsize.alpha.gamma.33.30.quants[i, 5], lty = 2)
  segments(x0 = yearstoreintro[i] - 0.1, x1 = yearstoreintro[i] + 0.1, y0 = popsize.alpha.gamma.33.30.quants[i, 3], y1 = popsize.alpha.gamma.33.30.quants[i, 3])
  segments(x0 = yearstoreintro[i] + .3, x1 = yearstoreintro[i] + .3, y0 = popsize.alpha.gamma.1.30.quants[i, 1], y1 = popsize.alpha.gamma.1.30.quants[i, 5], col = "red", lwd = 1)
  segments(x0 = yearstoreintro[i] + .2, x1 = yearstoreintro[i] + .4, y0 = popsize.alpha.gamma.1.30.quants[i, 3], y1 = popsize.alpha.gamma.1.30.quants[i, 3], col = "red", lwd = 1)
}
leg.text2 <- c("gamma = 0.33", "gamma = 0.1")
legend("topright", leg.text2, col = c("black", "red"), lwd = c(1, 1), lty = c(2, 1), bty = "n", cex = .8)



#------------------------------------------------------------------#
#-- project population sizes in loop over gammas of .1, .5, .9 ----#
#------------------------------------------------------------------#
timesteps <- 30
reps <- 100
alpha.range <- c(1, 100)
gamma.range <- c(1, 100)
alpha.steps <- 10
gamma.steps <- 10

timesteps <- 60
reps <- 100
popsize.30.gamma.05 <- matrix(NA, ncol = timesteps, nrow = reps)
popsize.30.gamma.1 <- matrix(NA, ncol = timesteps, nrow = reps)
popsize.30.gamma.2 <- matrix(NA, ncol = timesteps, nrow = reps)
popsize.30.gamma.5 <- matrix(NA, ncol = timesteps, nrow = reps)
popsize.30.gamma.9 <- matrix(NA, ncol = timesteps, nrow = reps)
popsize.30.gamma1 <- matrix(NA, ncol = timesteps, nrow = reps)

loglambda.30.gamma.05 <- matrix(NA, ncol = timesteps, nrow = reps)
loglambda.30.gamma.1 <- matrix(NA, ncol = timesteps, nrow = reps)
loglambda.30.gamma.2 <- matrix(NA, ncol = timesteps, nrow = reps)
loglambda.30.gamma.5 <- matrix(NA, ncol = timesteps, nrow = reps)
loglambda.30.gamma.9 <- matrix(NA, ncol = timesteps, nrow = reps)
loglambda.30.gamma1 <- matrix(NA, ncol = timesteps, nrow = reps)

for(i in 1:reps){
  project.30.gamma.05 <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = .1, gamma = .05, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)
  project.30.gamma.1 <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = .1, gamma = .1,  samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)
  project.30.gamma.2 <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = .1, gamma = .2,  samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)
  project.30.gamma.5 <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = .1, gamma = .5,  samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)
  project.30.gamma.9 <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = .1, gamma = .9,  samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)
  project.30.gamma1 <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = .1, gamma = 1,  samples.to.draw, tot.chains, joint.posterior.coda, posterior.names)
  
  popsize.30.gamma.05[i, ] <- project.30.gamma.05$tot.pop.size
  popsize.30.gamma.1[i, ] <-  project.30.gamma.1$tot.pop.size
  popsize.30.gamma.2[i, ] <-  project.30.gamma.2$tot.pop.size
  popsize.30.gamma.5[i, ] <-  project.30.gamma.5$tot.pop.size
  popsize.30.gamma.9[i, ] <-  project.30.gamma.9$tot.pop.size
  popsize.30.gamma1[i, ]  <-  project.30.gamma1$tot.pop.size
  
  loglambda.30.gamma.05[i, ] <- project.30.gamma.05$log.lambda.s
  loglambda.30.gamma.1[i, ] <-  project.30.gamma.1$log.lambda.s
  loglambda.30.gamma.2[i, ] <-  project.30.gamma.2$log.lambda.s
  loglambda.30.gamma.5[i, ] <-  project.30.gamma.5$log.lambda.s
  loglambda.30.gamma.9[i, ] <-  project.30.gamma.9$log.lambda.s
  loglambda.30.gamma1[i, ]  <-  project.30.gamma1$log.lambda.s
}  


# extract 2.5th and 97.5th quantiles --#
gamma.05.quants <- which(popsize.30.gamma.05[, timesteps - 1] %in% as.numeric(quantile(popsize.30.gamma.05[, timesteps - 1], c(0.025, 0.25, 0.5, 0.75, 0.975), type = 3, na.rm = T)))
gamma.1.quants <- which(popsize.30.gamma.1[, timesteps - 1] %in% as.numeric(quantile(popsize.30.gamma.1[, timesteps - 1], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))
gamma.2.quants <- which(popsize.30.gamma.2[, timesteps - 1] %in% as.numeric(quantile(popsize.30.gamma.2[, timesteps - 1], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))
gamma.5.quants <- which(popsize.30.gamma.5[, timesteps - 1] %in% as.numeric(quantile(popsize.30.gamma.5[, timesteps - 1], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))
gamma.9.quants <- which(popsize.30.gamma.9[, timesteps - 1] %in% as.numeric(quantile(popsize.30.gamma.9[, timesteps - 1], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))
gamma1.quants <- which(popsize.30.gamma1[, timesteps - 1] %in% as.numeric(quantile(popsize.30.gamma1[, timesteps - 1], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))

gamma.05.meds <- which(popsize.30.gamma.05[, timesteps - 1] %in% as.numeric(quantile(popsize.30.gamma.05[, timesteps - 1], .5, type = 3, na.rm = T)))
gamma.1.meds <- which(popsize.30.gamma.1[, timesteps - 1] %in% as.numeric(quantile(popsize.30.gamma.1[, timesteps - 1], .5, type = 3, na.rm = T)))
gamma.2.meds <- which(popsize.30.gamma.2[, timesteps - 1] %in% as.numeric(quantile(popsize.30.gamma.2[, timesteps - 1], .5, type = 3, na.rm = T)))
gamma.5.meds <- which(popsize.30.gamma.5[, timesteps - 1] %in% as.numeric(quantile(popsize.30.gamma.5[, timesteps - 1], .5, type = 3, na.rm = T)))
gamma.9.meds <- which(popsize.30.gamma.9[, timesteps - 1] %in% as.numeric(quantile(popsize.30.gamma.9[, timesteps - 1], .5, type = 3, na.rm = T)))
gamma1.meds <- which(popsize.30.gamma1[, timesteps - 1] %in% as.numeric(quantile(popsize.30.gamma1[, timesteps - 1], .5, type = 3, na.rm = T)))

plot.reps <- 100
par(mfrow = c(1, 5), oma = c(2, 1, 1, 1), mar = c(3, 4, 2, 0))
layout(matrix(c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 3, 5, byrow = T))
par(mfrow = c(1, 5), oma = c(2, 1, 1, 1), mar = c(3, 4, 2, 0))
plot(popsize.30.gamma.05[1, 11:59] ~ seq(11:59), type = "l", ylim = c(0, 800), xlab = "year", ylab = "population size", main = "20 years", bty = "n")
for(i in 2:plot.reps){
  lines(popsize.30.gamma.05[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
for(j in c(1, 2, 3, 4)){
  lines(popsize.30.gamma.05[gamma.05.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
}
lines(popsize.30.gamma.05[gamma.05.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2)  

plot(popsize.30.gamma.1[1, 11:59] ~ seq(11:59), type = "l", ylim = c(0, 800), xlab = "year", ylab = "", main = "10 years", bty = "n")
for(i in 2:plot.reps){
  lines(popsize.30.gamma.1[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
for(j in 1:4){
  lines(popsize.30.gamma.1[gamma.1.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
}
lines(popsize.30.gamma.1[gamma.1.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2)  

plot(popsize.30.gamma.2[1, 11:59] ~ seq(11:59), type = "l", ylim = c(0, 800), xlab = "year", ylab = "", main = "5 years", bty = "n")
for(i in 2:plot.reps){
  lines(popsize.30.gamma.2[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
for(j in 1:4){
  lines(popsize.30.gamma.2[gamma.2.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
}
lines(popsize.30.gamma.2[gamma.2.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2) 

plot(popsize.30.gamma1[1, 11:59] ~ seq(11:59), type = "l", ylim = c(0, 800), xlab = "year", ylab = "", main = "1 year", bty = "n")
for(i in 2:plot.reps){
  lines(popsize.30.gamma1[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
for(j in 1:4){
  lines(popsize.30.gamma1[gamma1.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
}
lines(popsize.30.gamma1[gamma1.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2) 

plot(popsize.he[1, 11:59] ~ seq(11:59), type = "l", ylim = c(0, 800), xlab = "year", ylab = "", main = "No disease", bty = "n")
for(i in 2:plot.reps){
  lines(popsize.he[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
for(j in 1:3){
  lines(popsize.he[he.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
}
lines(popsize.he[he.med, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2)  

#--------------------#
# log-lambda row
#--------------------#
# plot(loglambda.30.gamma.05[1, 11:59] ~ seq(11:59), type = "l", ylim = c(-.3, .3), xlab = "year", ylab = "population size", main = "", bty = "n")
# for(i in 2:plot.reps){
#   lines(loglambda.30.gamma.05[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
# }
# abline(h = 0, col = "red", lwd = 2)
# 
# plot(loglambda.30.gamma.1[1, 11:59] ~ seq(11:59), type = "l",  ylim = c(-.3, .3), xlab = "year", ylab = "", main = "", bty = "n")
# for(i in 2:plot.reps){
#   lines(loglambda.30.gamma.1[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
# }
# abline(h = 0, col = "red", lwd = 2)
# 
# plot(loglambda.30.gamma.2[1, 11:59] ~ seq(11:59), type = "l",  ylim = c(-.3, .3), xlab = "year", ylab = "", main = "", bty = "n")
# for(i in 2:plot.reps){
#   lines(loglambda.30.gamma.2[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
# }
# abline(h = 0, col = "red", lwd = 2)
# 
# plot(loglambda.30.gamma1[1, 11:59] ~ seq(11:59), type = "l",  ylim = c(-.3, .3), xlab = "year", ylab = "", main = "", bty = "n")
# for(i in 2:plot.reps){
#   lines(loglambda.30.gamma1[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
# }
# abline(h = 0, col = "red", lwd = 2)
# 
# plot( log.lambda.s.he[1, 11:59] ~ seq(11:59), type = "l", ylim = c(-.3, .3), xlab = "year", ylab = "", main = "", bty = "n")
# for(i in 2:plot.reps){
#   lines( log.lambda.s.he[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
# }
# abline(h = 0, col = "red", lwd = 2)

par(mfrow = c(1, 6))
boxplot(na.omit(as.vector(loglambda.30.gamma.05)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)
boxplot(na.omit(as.vector(loglambda.30.gamma.1)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)
boxplot(na.omit(as.vector(loglambda.30.gamma.2)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)
boxplot(na.omit(as.vector(loglambda.30.gamma.5)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)
boxplot(na.omit(as.vector(loglambda.30.gamma.9)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)
boxplot(na.omit(as.vector(loglambda.30.gamma1)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)



#-----------------------------------------------#
#-- OLD ----------------------------------------#
#-----------------------------------------------#

#     }
#   return(list(current.state.new = current.state.new))
# }
# 
# # function to update Leslie matrix parameters
#   #  need to expand Leslie to be 18x18... also, need individuals to age....
# update.leslie.fun <- function(current.state, he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr){
#   leslie.out <- rep(NA, 6, 6)
#   if(current.state == "healthy"){
#     repros <- c(0, rep(he.repro.post[sample(1:3000, 1), 2], 2), rep(he.repro.post[sample(1:3000, 1), 3], 4), rep(he.repro.post[sample(1:3000, 1), 4], 6), rep(he.repro.post[sample(1:3000, 1), 5], 5))    
# #    survs <- c(he.repro.post[sample(1:3000, 1), 1], rep(he.surv.post[sample(1:3000, 1), 2], 2), rep(he.surv.post[sample(1:3000, 1), 3], 4), rep(he.surv.post[sample(1:3000, 1), 4], 6), rep(he.surv.post[sample(1:3000, 1), 5], 4))
#     survs <- c(he.recr[sample(1:length(he.recr), 1)], rep(he.surv.post[sample(1:3000, 1), 2], 2), rep(he.surv.post[sample(1:3000, 1), 3], 4), rep(he.surv.post[sample(1:3000, 1), 4], 6), rep(he.surv.post[sample(1:3000, 1), 5], 4))
#     leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs))))
#   }
#   else {
#     repros <- c(0, rep(inf.repro.post[sample(1:3000, 1), 2], 2), rep(inf.repro.post[sample(1:3000, 1), 3], 4), rep(inf.repro.post[sample(1:3000, 1), 4], 6), rep(inf.repro.post[sample(1:3000, 1), 5], 5))    
# #    survs <- c(inf.repro.post[sample(1:3000, 1), 1], rep(inf.surv.post[sample(1:3000, 1), 2], 2), rep(inf.surv.post[sample(1:3000, 1), 3], 4), rep(inf.surv.post[sample(1:3000, 1), 4], 6), rep(inf.surv.post[sample(1:3000, 1), 5], 4))
#     survs <- c(pn.recr[sample(1:length(pn.recr), 1)], rep(inf.surv.post[sample(1:3000, 1), 2], 2), rep(inf.surv.post[sample(1:3000, 1), 3], 4), rep(inf.surv.post[sample(1:3000, 1), 4], 6), rep(inf.surv.post[sample(1:3000, 1), 5], 4))
#     leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs))))
#   }
#   return(leslie)
# }
# 
# 
# project.fun <- function(timesteps, ages.init, alpha, gamma, he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr){
#   N <- matrix(NA, nrow = length(ages.init), ncol = timesteps)
#   N[, 1] <- ages.init
#   tot.pop.size <- log.lambda.s <- rep(NA, length = timesteps)
#   disease.status <- rep(NA, timesteps)
# #  disease.status[1] <- "healthy"
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
# #  tot.pop.size <- apply(N, 2, sum)
# #  out.list <- list(N = N, disease.status = disease.status, tot.pop.size = tot.pop.size)
#   out.list <- list(N = N, disease.status = disease.status, tot.pop.size = tot.pop.size, log.lambda.s = log.lambda.s)
#   return(out.list)
# }
# 
# project.fun.out <- project.fun(timesteps = 20, ages.init = ages.init, alpha = .01, gamma = 1, he.repro = he.repro.post, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
# 
# #--------------------------------------------------------#
# #-- Model check: population trajectory with no disease --#
# #--------------------------------------------------------#
# he.project.fun <- function(timesteps, ages.init, alpha, gamma, he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr){
#   N <- matrix(NA, nrow = length(ages.init), ncol = timesteps)
#   N[, 1] <- ages.init
#   tot.pop.size <- log.lambda.s <- rep(NA, length = timesteps)
#   disease.status <- rep("healthy", timesteps)
#   #  disease.status[1] <- "healthy"
# #  disease.status[1:11] <- c(rep("healthy", 10), "spillover")
#   for(i in 1:(timesteps - 1)){
#     new.leslie <- update.leslie.fun(current.state = disease.status[i], he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
#     N[, i + 1] <- t(N[ , i]) %*% new.leslie  
#     tot.pop.size[i] <- sum(N[ , i])
#     if(i == 1){
#       log.lambda.s[i] <- NA
#     } else {
#       log.lambda.s[i] <- log(tot.pop.size[i] / tot.pop.size[i - 1])
#     }
#   }
# #   for(i in 11:(timesteps - 1)){
# #     new.leslie <- update.leslie.fun(current.state = disease.status[i], he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
# #     N[, i + 1] <- t(N[ , i]) %*% new.leslie
# #     disease.status[i + 1] <- update.status.fun(alpha, gamma, current.state = disease.status[i])$current.state.new[1]
# #   }
# #  tot.pop.size <- apply(N, 2, sum)
#   out.list <- list(N = N, disease.status = disease.status, tot.pop.size = tot.pop.size, log.lambda.s = log.lambda.s)
#   return(out.list)
# # }
# 
# sp.repro.post <- pn.repro
# sp.surv.post <- he.surv * runif(1, .3, 1)
# inf.surv.post <- pn.surv
# he.surv.post <- he.surv
# 
# he.project.fun.out <- he.project.fun(timesteps = 20, ages.init = ages.init, alpha = .01, gamma = 1, he.repro = he.repro.post, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
# 
# timesteps <- 60
# reps <- 100
# popsize.he <- log.lambda.s.he <- matrix(NA, ncol = timesteps, nrow = reps)
# 
# for(i in 1:reps){
#   he.project <- he.project.fun(timesteps = timesteps, ages.init = ages.init, alpha = .01, gamma = .1, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
#   popsize.he[i, ] <- he.project$tot.pop.size 
#   log.lambda.s.he[i, ] <- he.project$log.lambda.s
# }  
# 
# #write.csv(popsize.he, "~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Simulations/Healthy/HealthyPopsize_30Sept2014.csv", row.names = F)
# #popsize.he <- as.matrix(read.csv("~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Simulations/Healthy/HealthyPopsize_30Sept2014.csv"))
# 
# he.quants <- which(popsize.he[, 40] %in% as.numeric(quantile(popsize.he[, 40], c(0.025, 0.975), type = 3)))
# he.med <- which(popsize.he[, 40] %in% as.numeric(quantile(popsize.he[, 40], .5, type = 3)))
# 
# #par(mfrow = c(2, 1))
# layout(matrix(c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3), nrow = 2, byrow = T))
# plot(popsize.he[1, -c(1)] ~ seq(2, timesteps), type = "l", ylim = c(0, 3000), xlab = "year", ylab = "population size")
# for(i in 2:reps){
#   lines(popsize.he[i, -c(1)] ~ seq(2, timesteps), type = "l", col = rgb(.35, .35, .35, alpha = .25))
# }
# for(j in 1:3){
#   lines(popsize.he[he.quants[j], -c(1, 2)] ~ seq(3, timesteps), type = "l", col = "black", lwd = 2)  
# }
# lines(popsize.he[he.med, -c(1, 2)] ~ seq(3, timesteps), type = "l", col = "red", lwd = 2)
# 
# plot(log.lambda.s.he[1, -c(1, 2)] ~ seq(3, timesteps), type = "l", ylim = c(-.5, .5), xlab = "year", ylab = expression(paste("log(", lambda, "s)", sep = "")))
# for(i in 2:reps){
#   lines(log.lambda.s.he[i, -c(1, 2)] ~ seq(3, timesteps), type = "l", col = rgb(.35, .35, .35, alpha = .25))
# }
# abline( h = 0, lty = 2, col = "red", lwd = 2)
# boxplot(as.vector(log.lambda.s.he[, -c(1, 2)]), col = "grey80", ylim = c(-.5, .5), ylab = expression(paste("log(", lambda, "s)", sep = "")))
# abline(h = 0, lty = 2, col = "red", lwd = 2)

#-----------------------------------------------------#
#-- Environmental model 1: Movi presence definition --#
#-----------------------------------------------------#

#------------------------------------------------------------------#
#-- project population sizes in loop over gammas of .1, .5, .9 ----#
#------------------------------------------------------------------#
# popgrowth.sim.fun <- function(timesteps, reps, alpha.range, gamma.range, alpha.steps, gamma.steps){
#   alphas <- 1 / seq(min(alpha.range), max(alpha.range), length.out = alpha.steps)  
#   gammas <- 1 / seq(min(gamma.range), max(gamma.range), length.out = gamma.steps)
#   alpha.gamma.frame <- expand.grid(alphas, gammas)
#   popsize.ij <- loglambda.ij <- vector("list", dim(alpha.gamma.frame)[1])
#   mean.lnlambda <- rep(NA, dim(alpha.gamma.frame)[1])
#   for(i in 1:length(popsize.ij)){
#     popsize.ij[[i]] <- loglambda.ij[[i]] <- matrix(NA, nrow = timesteps, ncol = reps)
#     for(j in 1:reps){
#       popgrowthsim.ij <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = alpha.gamma.frame[i, 1], gamma = alpha.gamma.frame[i, 2], he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
#       popsize.ij[[i]][, j] <- popgrowthsim.ij$tot.pop.size
#       loglambda.ij[[i]][, j] <- popgrowthsim.ij$log.lambda.s
#     }
#     mean.lnlambda[i] <- mean(na.omit(unlist(as.vector(loglambda.ij[[i]][-c(11, 12), ]))))
#   }
#   outlist <- list(alpha.gamma.frame = alpha.gamma.frame, popsize.ij = popsize.ij, loglambda.ij = loglambda.ij, mean.lnlambda = mean.lnlambda)
#   return(outlist)
# }

timesteps <- 30
reps <- 100
alpha.range <- c(1, 100)
gamma.range <- c(1, 100)
alpha.steps <- 10
gamma.steps <- 10
popgrowth.test <- popgrowth.sim.fun(timesteps, reps, alpha.range, gamma.range, alpha.steps, gamma.steps)
popgrowth.test <- outlist
x.vals <- 1 / popgrowth.test$alpha.gamma.frame[, 1]
y.vals <- 1 / popgrowth.test$alpha.gamma.frame[, 2]

require(grDevices)
color.ramp <- colorRampPalette(c("red", "grey20", "blue"))
cols <- color.ramp(100)
color.val <- cols[round(100 * (popgrowth.test$mean.lnlambda - min(popgrowth.test$mean.lnlambda)) / (abs(min(popgrowth.test$mean.lnlambda )) - (min(popgrowth.test$mean.lnlambda)))) ]
#color.val2 <- color.val / max(abs(color.val))
pt.size <- 10 * (popgrowth.test$mean.lnlambda - min(popgrowth.test$mean.lnlambda)) / (abs(min(popgrowth.test$mean.lnlambda )) - (min(popgrowth.test$mean.lnlambda)))
#color.in <- paste("grey", round(100 * color.val), sep = "")

par(mfrow = c(1, 1))
plot(y.vals  ~ x.vals, cex = pt.size, xlab = expression(paste("expected waiting time until reintroduction (", alpha, ")", sep = "")), ylab = expression(paste("waiting time until expected fade-out (", gamma, ")", sep = "")), col = color.val, pch = 16)
axis(side = 1, at = as.numeric(levels(factor(popgrowth.test$alpha.gamma.frame[, 1]))), labels = round(1 / as.numeric(levels(factor(popgrowth.test$alpha.gamma.frame[, 1]))), 2))
axis(side = 2, at = as.numeric(levels(factor(popgrowth.test$alpha.gamma.frame[, 2]))), labels = round(1 / as.numeric(levels(factor(popgrowth.test$alpha.gamma.frame[, 2]))), 2))

timesteps <- 60
reps <- 1000
popsize.30.gamma.05 <- matrix(NA, ncol = timesteps, nrow = reps)
popsize.30.gamma.1 <- matrix(NA, ncol = timesteps, nrow = reps)
popsize.30.gamma.2 <- matrix(NA, ncol = timesteps, nrow = reps)
popsize.30.gamma.5 <- matrix(NA, ncol = timesteps, nrow = reps)
popsize.30.gamma.9 <- matrix(NA, ncol = timesteps, nrow = reps)
popsize.30.gamma1 <- matrix(NA, ncol = timesteps, nrow = reps)

loglambda.30.gamma.05 <- matrix(NA, ncol = timesteps, nrow = reps)
loglambda.30.gamma.1 <- matrix(NA, ncol = timesteps, nrow = reps)
loglambda.30.gamma.2 <- matrix(NA, ncol = timesteps, nrow = reps)
loglambda.30.gamma.5 <- matrix(NA, ncol = timesteps, nrow = reps)
loglambda.30.gamma.9 <- matrix(NA, ncol = timesteps, nrow = reps)
loglambda.30.gamma1 <- matrix(NA, ncol = timesteps, nrow = reps)

for(i in 1:reps){
  project.30.gamma.05 <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = 0, gamma = .05, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
  project.30.gamma.1 <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = 0, gamma = .1, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
  project.30.gamma.2 <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = 0, gamma = .2, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
  project.30.gamma.5 <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = 0, gamma = .5, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
  project.30.gamma.9 <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = 0, gamma = .9, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
  project.30.gamma1 <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = 0, gamma = 1, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
  
  popsize.30.gamma.05[i, ] <- project.30.gamma.05$tot.pop.size
  popsize.30.gamma.1[i, ] <-  project.30.gamma.1$tot.pop.size
  popsize.30.gamma.2[i, ] <-  project.30.gamma.2$tot.pop.size
  popsize.30.gamma.5[i, ] <-  project.30.gamma.5$tot.pop.size
  popsize.30.gamma.9[i, ] <-  project.30.gamma.9$tot.pop.size
  popsize.30.gamma1[i, ]  <-  project.30.gamma1$tot.pop.size
  
  loglambda.30.gamma.05[i, ] <- project.30.gamma.05$log.lambda.s
  loglambda.30.gamma.1[i, ] <-  project.30.gamma.1$log.lambda.s
  loglambda.30.gamma.2[i, ] <-  project.30.gamma.2$log.lambda.s
  loglambda.30.gamma.5[i, ] <-  project.30.gamma.5$log.lambda.s
  loglambda.30.gamma.9[i, ] <-  project.30.gamma.9$log.lambda.s
  loglambda.30.gamma1[i, ]  <-  project.30.gamma1$log.lambda.s
}  

#subset(popsize.30.gamma.1, (popsize.30.gamma.1[, 30]) %in% c(floor(quantile(popsize.30.gamma.1[, 30], 0.025)), ceiling(quantile(popsize.30.gamma.1[, 30], 0.025))))
#                                           , 0.975))))

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

# extract 2.5th and 97.5th quantiles --#
gamma.05.quants <- which(popsize.30.gamma.05[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.05[, timesteps], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))
gamma.1.quants <- which(popsize.30.gamma.1[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.1[, timesteps], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))
gamma.2.quants <- which(popsize.30.gamma.2[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.2[, timesteps], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))
gamma.5.quants <- which(popsize.30.gamma.5[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.5[, timesteps], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))
gamma.9.quants <- which(popsize.30.gamma.9[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.9[, timesteps], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))
gamma1.quants <- which(popsize.30.gamma1[, timesteps] %in% as.numeric(quantile(popsize.30.gamma1[, timesteps], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))

gamma.05.meds <- which(popsize.30.gamma.05[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.05[, timesteps], .5, type = 3, na.rm = T)))
gamma.1.meds <- which(popsize.30.gamma.1[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.1[, timesteps], .5, type = 3, na.rm = T)))
gamma.2.meds <- which(popsize.30.gamma.2[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.2[, timesteps], .5, type = 3, na.rm = T)))
gamma.5.meds <- which(popsize.30.gamma.5[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.5[, timesteps], .5, type = 3, na.rm = T)))
gamma.9.meds <- which(popsize.30.gamma.9[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.9[, timesteps], .5, type = 3, na.rm = T)))
gamma1.meds <- which(popsize.30.gamma1[, timesteps] %in% as.numeric(quantile(popsize.30.gamma1[, timesteps], .5, type = 3, na.rm = T)))

plot.reps <- 100
par(mfrow = c(2, 5))
plot(popsize.30.gamma.05[1, 11:59] ~ seq(11:59), type = "l", ylim = c(0, 1500), xlab = "year", ylab = "population size", main = "20 years", bty = "n")
for(i in 2:plot.reps){
  lines(popsize.30.gamma.05[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
# for(j in 1:4){
#   lines(popsize.30.gamma.05[gamma.05.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(popsize.30.gamma.05[gamma.05.meds, ] ~ seq(11:59), type = "l", col = "red", lwd = 2)  

plot(popsize.30.gamma.1[1, 11:59] ~ seq(11:59), type = "l", ylim = c(0, 1500), xlab = "year", ylab = "", main = "10 years", bty = "n")
for(i in 2:plot.reps){
  lines(popsize.30.gamma.1[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
# for(j in 1:4){
#   lines(popsize.30.gamma.1[gamma.1.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(popsize.30.gamma.1[gamma.1.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2)  

plot(popsize.30.gamma.2[1, 11:59] ~ seq(11:59), type = "l", ylim = c(0, 1500), xlab = "year", ylab = "", main = "5 years", bty = "n")
for(i in 2:plot.reps){
  lines(popsize.30.gamma.2[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
# for(j in 1:4){
#   lines(popsize.30.gamma.2[gamma.2.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(popsize.30.gamma.2[gamma.2.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2) 

plot(popsize.30.gamma1[1, 11:59] ~ seq(11:59), type = "l", ylim = c(0, 1500), xlab = "year", ylab = "", main = "1 year", bty = "n")
for(i in 2:plot.reps){
  lines(popsize.30.gamma1[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
# for(j in 1:4){
#   lines(popsize.30.gamma1[gamma1.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(popsize.30.gamma1[gamma1.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2) 

plot(popsize.he[1, 11:59] ~ seq(11:59), type = "l", ylim = c(0, 1000), xlab = "year", ylab = "", main = "No disease", bty = "n")
for(i in 2:plot.reps){
  lines(popsize.he[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
# for(j in 1:3){
#   lines(popsize.he[he.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(popsize.he[he.med, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2)  

# log-lambda row
plot(loglambda.30.gamma.05[1, 11:59] ~ seq(11:59), type = "l", ylim = c(-1, 1), xlab = "year", ylab = "population size", main = "20 years", bty = "n")
for(i in 2:plot.reps){
  lines(loglambda.30.gamma.05[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
abline(h = 0, col = "red", lwd = 2)
# for(j in 1:4){
#   lines(loglambda.30.gamma.05[gamma.05.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(loglambda.30.gamma.05[gamma.05.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2)  

plot(loglambda.30.gamma.1[1, 11:59] ~ seq(11:59), type = "l", ylim = c(-1, 1), xlab = "year", ylab = "", main = "10 years", bty = "n")
for(i in 2:plot.reps){
  lines(loglambda.30.gamma.1[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
abline(h = 0, col = "red", lwd = 2)
# for(j in 1:4){
#   lines(loglambda.30.gamma.1[gamma.1.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(loglambda.30.gamma.1[gamma.1.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2)  

plot(loglambda.30.gamma.2[1, 11:59] ~ seq(11:59), type = "l", ylim = c(-1, 1), xlab = "year", ylab = "", main = "5 years", bty = "n")
for(i in 2:plot.reps){
  lines(loglambda.30.gamma.2[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
abline(h = 0, col = "red", lwd = 2)
# for(j in 1:4){
#   lines(loglambda.30.gamma.2[gamma.2.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(loglambda.30.gamma.2[gamma.2.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2) 

plot(loglambda.30.gamma1[1, 11:59] ~ seq(11:59), type = "l", ylim = c(-1, 1), xlab = "year", ylab = "", main = "1 year", bty = "n")
for(i in 2:plot.reps){
  lines(loglambda.30.gamma1[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
abline(h = 0, col = "red", lwd = 2)
# for(j in 1:4){
#   lines(loglambda.30.gamma1[gamma1.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(loglambda.30.gamma1[gamma1.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2) 

plot( log.lambda.s.he[1, 11:59] ~ seq(11:59), type = "l", ylim = c(-1, 1), xlab = "year", ylab = "", main = "No disease", bty = "n")
for(i in 2:plot.reps){
  lines( log.lambda.s.he[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
abline(h = 0, col = "red", lwd = 2)
# for(j in 1:3){
#   lines(loglambda.he[he.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(loglambda.he[he.med, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2)  

par(mfrow = c(1, 6))
boxplot(na.omit(as.vector(loglambda.30.gamma.05)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)
boxplot(na.omit(as.vector(loglambda.30.gamma.1)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)
boxplot(na.omit(as.vector(loglambda.30.gamma.2)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)
boxplot(na.omit(as.vector(loglambda.30.gamma.5)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)
boxplot(na.omit(as.vector(loglambda.30.gamma.9)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)
boxplot(na.omit(as.vector(loglambda.30.gamma1)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)

#-------------------------------------------------------#
#-- costs of altering alpha/gamma on 30-year pop size --#
#-------------------------------------------------------#
#-- better to think about gammas as expected time to fade-out.... --#
#-- F(x) = 1/beta * (exp (-x / beta))
#-- where beta = expected survival of system
#-- so, gammas of interest are 1, 1/2, 1/3, 1/4

timesteps <- 40
reps <- 1000
yearstofadeout.seq <- seq(1, 15)
gamma.seq <- 1 / yearstofadeout.seq
popsize.30.gamma.mat <- matrix(NA, reps, length(gamma.seq))
popsize.30.gamma.alpha.33.mat <- matrix(NA, reps, length(gamma.seq))

for(i in 1:reps){
  for(j in 1:length(gamma.seq)){
    popsize.30.gamma.mat[i, j] <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = .1, gamma = gamma.seq[j], he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)$tot.pop.size[timesteps]
    popsize.30.gamma.alpha.33.mat[i, j] <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = .33, gamma = gamma.seq[j], he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)$tot.pop.size[timesteps]
  }
}

popsize.30.quants <- matrix(NA, length(gamma.seq), 5)
popsize.30.alpha.33.quants <- matrix(NA, length(gamma.seq), 5)
for(j in 1:length(gamma.seq)){
  popsize.30.quants[j, ] <- quantile(popsize.30.gamma.mat[, j], c(0.025, 0.25, 0.5, 0.75, 0.975))
  popsize.30.alpha.33.quants[j, ] <- quantile(popsize.30.gamma.alpha.33.mat[, j], c(0.025, 0.25, 0.5, 0.75, 0.975))
}

par(mfrow = c(1, 1))
plot(x = 0, y = 0, xlim = c(0, 16), ylim = c(0, 200), xlab = "Expected years to fade-out (alpha = 0.1)", ylab = "Simulated Pop Size After 30 years")
for(i in 1:length(yearstofadeout.seq)){
  segments(x0 = yearstofadeout.seq[i], x1 = yearstofadeout.seq[i], y0 = popsize.30.quants[i, 1], y1 = popsize.30.quants[i, 5])
  segments(x0 = yearstofadeout.seq[i] - 0.1, x1 = yearstofadeout.seq[i] + 0.1, y0 = popsize.30.quants[i, 3], y1 = popsize.30.quants[i, 3])
  segments(x0 = yearstofadeout.seq[i] + .2, x1 = yearstofadeout.seq[i] + .2, y0 = popsize.30.alpha.33.quants[i, 1], y1 = popsize.30.alpha.33.quants[i, 5], col = "red")
  segments(x0 = yearstofadeout.seq[i] + .1, x1 = yearstofadeout.seq[i] + .3, y0 = popsize.30.alpha.33.quants[i, 3], y1 = popsize.30.alpha.33.quants[i, 3], col = "red")
}

#-- Same thing, but over levels of alpha, with gamma fixed at e(fade-out) = 3 and e(fade-out) = 10
timesteps <- 40
reps <- 1000
yearstoreintro<- seq(1, 15)
alpha.seq <- 1 / yearstoreintro
popsize.30.alpha.gamma.33.mat <- matrix(NA, reps, length(alpha.seq))
popsize.30.alpha.gamma.1.mat <- matrix(NA, reps, length(alpha.seq))

for(i in 1:reps){
  for(j in 1:length(alpha.seq)){
    popsize.30.alpha.gamma.33.mat[i, j] <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = alpha.seq[j], gamma = .33, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)$tot.pop.size[timesteps]
    popsize.30.alpha.gamma.1.mat[i, j] <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = alpha.seq[j], gamma = .1, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)$tot.pop.size[timesteps]
  }
}

popsize.alpha.gamma.33.30.quants <- matrix(NA, length(alpha.seq), 5)
popsize.alpha.gamma.1.30.quants <- matrix(NA, length(alpha.seq), 5)
for(j in 1:length(gamma.seq)){
  popsize.alpha.gamma.33.30.quants[j, ] <- quantile(popsize.30.alpha.gamma.33.mat[, j], c(0.025, 0.25, 0.5, 0.75, 0.975))
  popsize.alpha.gamma.1.30.quants[j, ] <- quantile(popsize.30.alpha.gamma.1.mat[, j], c(0.025, 0.25, 0.5, 0.75, 0.975))
}

par(mfrow = c(1, 1))
plot(x = 0, y = 0, xlim = c(0, 16), ylim = c(0, 100), xlab = "Expected years to reintroduction", ylab = "Simulated Pop Size After 30 years")
for(i in 1:length(yearstoreintro)){
  segments(x0 = yearstoreintro[i], x1 = yearstoreintro[i], y0 = popsize.alpha.gamma.33.30.quants[i, 1], y1 = popsize.alpha.gamma.33.30.quants[i, 5])
  segments(x0 = yearstoreintro[i] - 0.1, x1 = yearstoreintro[i] + 0.1, y0 = popsize.alpha.gamma.33.30.quants[i, 3], y1 = popsize.alpha.gamma.33.30.quants[i, 3])
  segments(x0 = yearstoreintro[i] + .2, x1 = yearstoreintro[i] + .2, y0 = popsize.alpha.gamma.1.30.quants[i, 1], y1 = popsize.alpha.gamma.1.30.quants[i, 5], col = "red")
  segments(x0 = yearstoreintro[i] + .1, x1 = yearstoreintro[i] + .3, y0 = popsize.alpha.gamma.1.30.quants[i, 3], y1 = popsize.alpha.gamma.1.30.quants[i, 3], col = "red")
}

# #--------------------------------------------------------------------------#
# #-- Exploration of lambda and age structure in each environmental state ---#
# #--------------------------------------------------------------------------#
# healthy.leslie.list <- spillover.leslie.list <- infected.leslie.list <- vector("list", length = 1000)
# healthy.eigenval1 <- spillover.eigenval1 <- infected.eigenval1 <- rep(NA, 1000)
# juvsurv.elast.he <- juvsurv.elast.inf <- adsurv.elast.he <- adsurv.elast.inf <- fecund.elast.he <- fecund.elast.inf <- rep(NA, 1000)
# age.struct.he <- age.struct.inf <- matrix(NA, nrow = 18, ncol = 1000)
# 
# for(i in 1:1000){
#   function(current.state, samples.to.draw = 1, tot.chains = 3, joint.posterior.coda = ipmcoda.11, posterior.names = posterior.names){
#     
#   healthy.leslie.list[[i]] <- update.leslie.fun(current.state = "healthy", he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
#   spillover.leslie.list[[i]] <- update.leslie.fun(current.state = "spillover", he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
#   infected.leslie.list[[i]] <- update.leslie.fun(current.state = "infected", he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
#   healthy.eigenval1[i] <- Re(eigen(healthy.leslie.list[[i]])$values[1]) # strip off only real part of eigenvalue
#   spillover.eigenval1[i] <- Re(eigen(spillover.leslie.list[[i]])$values[1])
#   infected.eigenval1[i] <- Re(eigen(infected.leslie.list[[i]])$values[1])
#   he.eigen.rescale <- sum(Re(eigen(healthy.leslie.list[[i]])$vectors[, 1]))
#   inf.eigen.rescale <- sum(Re(eigen(infected.leslie.list[[i]])$vectors[, 1]))
#   age.struct.he[, i] <- Re(eigen(healthy.leslie.list[[i]])$vectors[, 1]) / he.eigen.rescale
#   age.struct.inf[, i] <- Re(eigen(infected.leslie.list[[i]])$vectors[, 1]) / inf.eigen.rescale
#   
#   he.sens <- sensitivity(healthy.leslie.list[[i]])
#   he.elast <- (1 / Re(eigen(healthy.leslie.list[[i]])$values[1])) * he.sens * healthy.leslie.list[[i]]
# #  round(he.sens, 2)
# #  round(he.elast, 2)
#   fecund.elast.he[i] <- sum(he.elast[1, ])
#   juvsurv.elast.he[i] <- sum(he.elast[2, 1], he.elast[3, 2])
#   adsurv.elast.he[i] <- sum(he.elast[4, 3], he.elast[5, 4], he.elast[6, 5], he.elast[7, 6], he.elast[8, 7], he.elast[9, 8], he.elast[10, 9], he.elast[11, 10], he.elast[12, 11], he.elast[13, 12], he.elast[14, 13], he.elast[15, 14], he.elast[16, 15], he.elast[17, 16], he.elast[18, 17])
#   
#   inf.sens <- sensitivity(infected.leslie.list[[i]])
#   inf.elast <- (1 / Re(eigen(infected.leslie.list[[i]])$values[1])) * inf.sens * infected.leslie.list[[i]]
# #  round(inf.sens, 2)
# #  round(inf.elast, 2)
#   fecund.elast.inf[i] <- sum(inf.elast[1, ])
#   juvsurv.elast.inf[i] <- sum(inf.elast[2, 1], inf.elast[3, 2])
#   #sum(inf.elast[5, 4], inf.elast[6, 5], inf.elast[7, 6], inf.elast[8, 7], inf.elast[9, 8], inf.elast[10, 9])
#   adsurv.elast.inf[i] <- sum(inf.elast[4, 3], inf.elast[5, 4], inf.elast[6, 5], inf.elast[7, 6], inf.elast[8, 7], inf.elast[9, 8], inf.elast[10, 9], inf.elast[11, 10], inf.elast[12, 11], inf.elast[13, 12], inf.elast[14, 13], inf.elast[15, 14], inf.elast[16, 15], inf.elast[17, 16], inf.elast[18, 17])
# }
# 
# # 95% intervals on lambda 
# healthy.cred<- quantile(healthy.eigenval1, c(0.25, 0.5, 0.75))
# spillover.cred<- quantile(spillover.eigenval1, c(0.25, 0.5, 0.75))
# infected.cred<- quantile(infected.eigenval1, c(0.25, 0.5, 0.75))
# 
# # 95% intervals for fecundity, juvenile survival, and adult survival in each environment (healthy and infected)
# fecund.he <- quantile(fecund.elast.he, c(0.25, 0.5, 0.75))
# fecund.inf <- quantile(fecund.elast.inf, c(.25, 0.5, 0.75))
# juvsurv.he <- quantile(juvsurv.elast.he, c(0.25, 0.5, 0.75))
# juvsurv.inf <- quantile(juvsurv.elast.inf, c(0.25, 0.5, 0.75))
# adsurv.he <- quantile(adsurv.elast.he, c(0.25, 0.5, 0.75))
# adsurv.inf <- quantile(adsurv.elast.inf, c(0.25, 0.5, 0.75))
# 
# par(mfrow = c(1, 1), mar = c(4, 8, 2, 2), las = 1, cex.lab = 1.0)
# plot(c(1, 1) ~ c(fecund.he[1], fecund.he[3]), lty = 1, xlim = c(0, 1), ylim = c(0, 7), type = "l", xlab = expression(paste("Elasticity of ", lambda, " to rate", sep = "")), ylab = "", yaxt = "n", lwd = 2)
# lines(c(2, 2) ~ c(fecund.inf[1], fecund.inf[3]), lty = 2, col = "red", lwd = 2)
# lines(c(3, 3) ~ c(juvsurv.he[1], juvsurv.he[3]), lty = 1, col = "black", lwd = 2)
# lines(c(4, 4) ~ c(juvsurv.inf[1], juvsurv.inf[3]), lty = 2, col = "red", lwd = 2)
# lines(c(5, 5) ~ c(adsurv.he[1], adsurv.he[3]), lty = 1, col = "black", lwd = 2)
# lines(c(6, 6) ~ c(adsurv.inf[1], adsurv.inf[3]), lty = 2, col = "red", lwd = 2)
# lines(c(0.75, 1.25) ~ c(fecund.he[2], fecund.he[2]), lty = 1, col = "black", lwd = 2)
# lines(c(1.75, 2.25) ~ c(fecund.inf[2], fecund.inf[2]), lty = 2, col = "red", lwd = 2)
# lines(c(2.75, 3.25) ~ c(juvsurv.he[2], juvsurv.he[2]), lty = 1, col = "black", lwd = 2)
# lines(c(3.75, 4.25) ~ c(juvsurv.inf[2], juvsurv.inf[2]), lty = 2, col = "red", lwd = 2)
# lines(c(4.75, 5.25) ~ c(adsurv.he[2], adsurv.he[2]), lty = 1, col = "black", lwd = 2)
# lines(c(5.75, 6.25) ~ c(adsurv.inf[2], adsurv.inf[2]), lty = 2, col = "red", lwd = 2)
# axis(side = 2, at = c(1:6), cex.axis = 1.0, labels = c("Fecundity (he)", "Fecundity (pers)", "Juvenile survival (he)", "Juvenile surv (pers)", "Adult survival (he)", "Adult survival (pers)"))
# 
# par(mfrow = c(1, 3), cex.axis = 1.2, cex.lab = 1.5)
# hist(healthy.eigenval1, xlim = c(min(min(healthy.eigenval1), min(spillover.eigenval1), min(infected.eigenval1)), max(max(healthy.eigenval1), max(spillover.eigenval1), max(infected.eigenval1))), breaks = 15, xlab = expression(lambda), main = "Healthy", col = "grey80")
# abline(v = 1, col = "red", lwd = 3)
# hist(spillover.eigenval1, xlim = c(min(min(healthy.eigenval1), min(spillover.eigenval1), min(infected.eigenval1)), max(max(healthy.eigenval1), max(spillover.eigenval1), max(infected.eigenval1))), breaks = 15, ylab = "", xlab = expression(lambda), main = "Introduction", col = "grey80")
# abline(v = 1, col = "red", lwd = 3)
# hist(infected.eigenval1, xlim = c(min(min(healthy.eigenval1), min(spillover.eigenval1), min(infected.eigenval1)), max(max(healthy.eigenval1), max(spillover.eigenval1), max(infected.eigenval1))), breaks = 15, ylab = "", xlab = expression(lambda), main = "Persistence", col = "grey80")
# abline(v = 1, col = "red", lwd = 3)
# 
# # environ-state-specific elasticities
# he.sens <- sensitivity(healthy.leslie.list[[1]])
# he.elast <- (1 / Re(eigen(healthy.leslie.list[[1]])$values[1])) * he.sens * healthy.leslie.list[[1]]
# round(he.sens, 2)
# round(he.elast, 2)
# fecund.elast.he <- sum(he.elast[1, ])
# juvsurv.elast.he <- sum(he.elast[2, 1], he.elast[3, 2])
# adsurv.elast.he <- sum(he.elast[4, 3], he.elast[5, 4], he.elast[6, 5], he.elast[7, 6], he.elast[8, 7], he.elast[9, 8], he.elast[10, 9], he.elast[11, 10], he.elast[12, 11], he.elast[13, 12], he.elast[14, 13], he.elast[15, 14], he.elast[16, 15], he.elast[17, 16], he.elast[18, 17])
# 
# inf.sens <- sensitivity(infected.leslie.list[[1]])
# inf.elast <- (1 / Re(eigen(infected.leslie.list[[1]])$values[1])) * inf.sens * infected.leslie.list[[1]]
# round(inf.sens, 2)
# round(inf.elast, 2)
# fecund.elast.inf <- sum(inf.elast[1, ])
# juvsurv.elast.inf <- sum(inf.elast[2, 1], inf.elast[3, 2])
# #sum(inf.elast[5, 4], inf.elast[6, 5], inf.elast[7, 6], inf.elast[8, 7], inf.elast[9, 8], inf.elast[10, 9])
# adsurv.elast.inf <- sum(inf.elast[4, 3], inf.elast[5, 4], inf.elast[6, 5], inf.elast[7, 6], inf.elast[8, 7], inf.elast[9, 8], inf.elast[10, 9], inf.elast[11, 10], inf.elast[12, 11], inf.elast[13, 12], inf.elast[14, 13], inf.elast[15, 14], inf.elast[16, 15], inf.elast[17, 16], inf.elast[18, 17])
# 
# # intervals on age-structure in healthy and infected environments 
# age.bounds.he <- age.bounds.inf <- matrix(NA, nrow = 18, ncol = 3)
# for(i in 1:18){
#   age.bounds.he[i, ] <- quantile(abs(age.struct.he[i, ]), c(0.025, 0.5, 0.975))
#   age.bounds.inf[i, ] <- quantile(abs(age.struct.inf[i, ]), c(0.025, 0.5, 0.975))
# }
# 
# par(mfrow = c(1, 1))
# plot(age.bounds.he[1, c(1, 3)] ~ c(1, 1), type = "l", xlim = c(0.5, 18.5), ylim = c(0, .5), lwd = 2, xlab = "age (years)", ylab = "Expected proportion of pop")
# segments(x0 = 1 + .3, x1 = 1 + .3, y0 = age.bounds.inf[1, 1], y1 = age.bounds.inf[1, 2], lwd = 2, col = "red")
# for(i in 2: dim(age.bounds.he)[1]){
#   segments(x0 = i, x1 = i, y0 = age.bounds.he[i, 1], y1 = age.bounds.he[i, 2], lwd = 2)
#   segments(x0 = i + .3, x1 = i + .3, y0 = age.bounds.inf[i, 1], y1 = age.bounds.inf[i, 2], lwd = 2, col = "red")
# }
# leg.text <- c("healthy", "persistently infected")
# legend("topright", leg.text, col = c("black", "red"), lwd = c(2, 2), bty = "n")
# 
# #--------------------------------------------------------#
# #-- Sensitivies and Elasticities using Vec-permutation --#
# #--------------------------------------------------------#
# 
# vec.permut.fun <- function(reps, alpha, gamma){
#   fade.out.elast <- intro.elast <- rep(NA, reps)
#   for(i in 1:reps){
# # P is the vec-permutation matrix.
# P <- matrix(NA, nrow = 18 * 3, ncol = 18 * 3)
# for(j in 1:18){
#   odd.vec <- c(rep(0, j - 1), 1, rep(0, 18 - j))
#     P[j * 3 - 2 , ] <- c(odd.vec, rep(0, 18), rep(0, 18))
#     P[j * 3 - 1, ] <- c(rep(0, 18), odd.vec, rep(0, 18))
#     P[j * 3, ] <- c(rep(0, 18), rep(0, 18), odd.vec)
# }
# 
# # B is block diagonal, with 3 17x17 blocks for the 3 environmental states.
# healthy.leslie <- update.leslie.fun(current.state = "healthy", he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
# spillover.leslie <- update.leslie.fun(current.state = "spillover", he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
# endemic.leslie <- update.leslie.fun(current.state = "infected", he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
# leslie.list <- list(healthy.leslie, spillover.leslie, endemic.leslie)
# B <- bdiag(leslie.list)
# 
# # M is block diagonal with 17 3x3 blocks for the 17 demographic states
# #alpha <- 1/10
# #gamma <- 1/20
# small.M <- rbind(c(1 - alpha, alpha, 0), c(gamma, 0, 1 - gamma), c(gamma, 0, 1 - gamma))
# M <- bdiag(list(small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M))
# 
# # A is population projection matrix with environmental stochasticity
# A <- t(P) %*% M %*% P %*% B
# A_eigens <- eigen(A) # complex????
# S_a <- sensitivity(A)
# 
# # Sensitivity (environmental transitions)
# # assume disease status is updated before demography.
# S_m <- P %*% t(B) %*% S_a %*% t(P)
# round(S_m, 2)[1:10, 1:10]
# #E_m <- (1 / A_eigens$value[1]) * M * S_m # regular * because it's a Hadamard production
# E_m <- (1 / .92) * M * S_m # regular * because it's a Hadamard production
# round(E_m, 2)[1:10, 1:10]
# # compare elasticities of fade-out to elasticity of reintroduction
# fade.out.elast[i] <- sum(E_m[3, 3], E_m[6, 6], E_m[9, 9], E_m[12, 12], E_m[15, 15], E_m[18, 18], E_m[21, 21], E_m[24, 24], E_m[27, 27], E_m[30, 30], E_m[33, 33], E_m[36, 36], E_m[39, 39], E_m[42, 42], E_m[45, 45], E_m[48, 48], E_m[51, 51])
# intro.elast[i] <- sum(E_m[1, 1], E_m[4, 4], E_m[7, 7], E_m[10, 10], E_m[13, 13], E_m[16, 16], E_m[19, 19], E_m[22, 22], E_m[25, 25], E_m[28, 28], E_m[31, 31], E_m[34, 34], E_m[37, 37], E_m[40, 40], E_m[43, 43], E_m[46, 46], E_m[49, 49])
# }
# return(list(fade.out.elast = fade.out.elast, intro.elast = intro.elast))
# }
# 
# vec.permut.test1 <- vec.permut.fun(reps = 1000, alpha = 0.2, gamma = 0.1)
# vec.permut.test2 <- vec.permut.fun(reps = 1000, alpha = 0.1, gamma = 0.1)
# vec.permut.test3 <- vec.permut.fun(reps = 1000, alpha = 0.05, gamma = 0.1)
# 
# fade.out.elasts <- quantile(vec.permut.test$fade.out.elast, c(0.025, 0.975))
# intro.elasts <- quantile(vec.permut.test$intro.elast, c(0.025, 0.975))
# 
# par(mfrow = c(2, 3), cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.4)
# hist(vec.permut.test3$fade.out.elast, ylab = "frequency", xlab = "fade-out elasticity", col = "grey80", main = "E(time to reintro) = 20yrs")
# hist(vec.permut.test2$fade.out.elast, ylab = "frequency", xlab = "fade-out elasticity", col = "grey80", main = "E(time to reintro) = 10yrs")
# hist(vec.permut.test1$fade.out.elast, ylab = "frequency", xlab = "fade-out elasticity", col = "grey80", main = "E(time to reintro) = 5yrs")
# hist(vec.permut.test3$intro.elast, main = "", ylab = "frequency", xlab = "introduction elasticity", col = "grey80")
# hist(vec.permut.test2$intro.elast, main = "", ylab = "frequency", xlab = "introduction elasticity", col = "grey80")
# hist(vec.permut.test1$intro.elast, main = "", ylab = "frequency", xlab = "introduction elasticity", col = "grey80")
# 
# # Sensitivity (demographic transitions)
# S_b <- S_a %*% t(P) %*% t(M) %*% P
# round(S_b, 2)[1:20, 1:20]
# E_b <- (1 / .92) * B * S_b # regular * because it's a Hadamard production
# round(E_b, 2)

#-------------------------------------------------------#
#-- Figure 2: Leslie Matrix posterior param estimates --#
#-------------------------------------------------------#

#-- Summer Lamb Survival --#
par(cex.main = .8, mfrow = c(1, 2), cex.main = 1.0, las = 1)
plot(xlim = c(1.5, 5.5), xaxt = "n", ylim = c(0, 1), x = -1, y = 1, main = "Age-specific survival", xlab = "age class", ylab = "P(survives)")
for(i in 2:5){ 
  segments(x0 = i, x1 = i, y0 = summary(coda.samples.agespecsurv.pn)[[2]][i, 1], y1 = summary(coda.samples.agespecsurv.pn)[[2]][i, 5], col = "red", lwd = 2)
  segments(x0 = i + 0.25, x1 = i + 0.25, y0 = summary(coda.samples.agespecsurv.he)[[2]][i, 1], y1 = summary(coda.samples.agespecsurv.he)[[2]][i, 5], col = "grey30", lwd = 2, lty = 2)
  #  segments(x0 = i + 0.15, x1 = i + 0.15, y0 = summary(coda.samples.sheep.agespecrepro.premovi)[[2]][i, 1], y1 = summary(coda.samples.agespecsurv.premovi)[[2]][i, 5], col = "grey60", lwd = 2, lty = 3)
}
axis(side = 1, at = 2:5, labels = c("<2.5", "2.5-7", "8-13", ">13"))
legend("bottomright", cex = .8, c("persistent years", "healthy years"), lty = c(1, 2), col = c("red", "grey30"), lwd = c(2, 2), bty = "n")

plot(xlim = c(1.5, 5.5), xaxt = "n", ylim = c(0, 1), x = -1, y = 1, main = "Age-specific reproduction", xlab = "age class", ylab = "P(weans a lamb)")
for(i in 2:5){ 
  segments(x0 = i, x1 = i, y0 = summary(coda.samples.sheep.agespecrepro.pn)[[2]][i, 1], y1 = summary(coda.samples.sheep.agespecrepro.pn)[[2]][i, 5], col = "red", lwd = 2)
  segments(x0 = i + 0.25, x1 = i + 0.25, y0 = summary(coda.samples.sheep.agespecrepro.he)[[2]][i, 1], y1 = summary(coda.samples.sheep.agespecrepro.he)[[2]][i, 5], col = "grey30", lwd = 2, lty = 2)
# segments(x0 = i + 0.15, x1 = i + 0.15, y0 = summary(coda.samples.sheep.agespecrepro.premovi)[[2]][i, 1], y1 = summary(coda.samples.sheep.agespecrepro.premovi)[[2]][i, 5], col = "grey60", lwd = 2, lty = 3)
}
axis(side = 1, at = 2:5, labels = c("<2.5", "2.5-7", "8-13", ">13"))
#legend("topright", c("lamb disease years", "lamb healthy years"), lty = c(1, 2), col = c("red", "grey30"), lwd = c(2, 2), bty = "n")

#-- Age-specific ewe survival --#


#-----------------------------------------------------#
#-- Figure 3:Simulations with various fadeout times --# 
#-----------------------------------------------------#


par(mfrow = c(1, 5), cex.lab = 1.4)
plot(popsize.30.gamma.05[1, ] ~ seq(1, timesteps), type = "l", xaxt = "n", ylim = c(0, 1500), xlab = "year", ylab = "population size", main = "20 years", bty = "n")
for(i in 2:reps){
  lines(popsize.30.gamma.05[i, ] ~ seq(1, timesteps), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
for(j in 1:4){
  lines(popsize.30.gamma.05[gamma.05.quants[j], ] ~ seq(1, timesteps), type = "l", col = "black", lwd = 2)  
}
lines(popsize.30.gamma.05[gamma.05.meds, ] ~ seq(1, timesteps), type = "l", col = "red", lwd = 2)  
axis(side = 1, at = (seq(1:7) * 10 - 10), labels = c("-10","0", "10","20", "30", "40", "50"))

plot(popsize.30.gamma.1[1, ] ~ seq(1, timesteps), type = "l", xaxt = "n", ylim = c(0, 1500), xlab = "year", ylab = "", main = "10 years", bty = "n")
for(i in 2:reps){
  lines(popsize.30.gamma.1[i, ] ~ seq(1, timesteps), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
for(j in 1:4){
  lines(popsize.30.gamma.1[gamma.1.quants[j], ] ~ seq(1, timesteps), type = "l", col = "black", lwd = 2)  
}
lines(popsize.30.gamma.1[gamma.1.meds, ] ~ seq(1, timesteps), type = "l", col = "red", lwd = 2)  
axis(side = 1, at = (seq(1:7) * 10 - 10), labels = c("-10","0", "10","20", "30", "40", "50"))

plot(popsize.30.gamma.2[1, ] ~ seq(1, timesteps), type = "l", xaxt = "n", ylim = c(0, 1500), xlab = "year", ylab = "", main = "5 years", bty = "n")
for(i in 2:reps){
  lines(popsize.30.gamma.2[i, ] ~ seq(1, timesteps), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
for(j in 1:4){
  lines(popsize.30.gamma.2[gamma.2.quants[j], ] ~ seq(1, timesteps), type = "l", col = "black", lwd = 2)  
}
lines(popsize.30.gamma.2[gamma.2.meds, ] ~ seq(1, timesteps), type = "l", col = "red", lwd = 2) 
axis(side = 1, at = (seq(1:7) * 10 - 10), labels = c("-10","0", "10","20", "30", "40", "50"))

plot(popsize.30.gamma1[1, ] ~ seq(1, timesteps), type = "l", xaxt = "n", ylim = c(0, 1500), xlab = "year", ylab = "", main = "1 year", bty = "n")
for(i in 2:reps){
  lines(popsize.30.gamma1[i, ] ~ seq(1, timesteps), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
for(j in 1:4){
  lines(popsize.30.gamma1[gamma1.quants[j], ] ~ seq(1, timesteps), type = "l", col = "black", lwd = 2)  
}
lines(popsize.30.gamma1[gamma1.meds, ] ~ seq(1, timesteps), type = "l", col = "red", lwd = 2) 
axis(side = 1, at = (seq(1:7) * 10 - 10), labels = c("-10","0", "10","20", "30", "40", "50"))

plot(popsize.he[1, ] ~ seq(1, timesteps), type = "l", xaxt = "n", ylim = c(0, 1000), xlab = "year", ylab = "", main = "No disease", bty = "n")
for(i in 2:reps){
  lines(popsize.he[i, ] ~ seq(1, timesteps), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
for(j in 1:3){
  lines(popsize.he[he.quants[j], ] ~ seq(1, timesteps), type = "l", col = "black", lwd = 2)  
}
lines(popsize.he[he.med, ] ~ seq(1, timesteps), type = "l", col = "red", lwd = 2)  
axis(side = 1, at = (seq(1:7) * 10 - 10), labels = c("-10","0", "10","20", "30", "40", "50"))

#--------------------------------------------------#
#-- Figure 4: pop size after 30 years -------------#
#--------------------------------------------------#
par(mfrow = c(1, 2))
plot(x = 0, y = 0, xlim = c(0, 16), ylim = c(0, 2000), xlab = "Expected years to fade-out", ylab = "Simulated Pop Size After 30 years")
abline(h = 1000, col = "grey80", lty = 2)
abline(h = 1500, col = "grey80", lty = 2)
abline(h = 500, col = "grey80", lty = 2)
for(i in 1:length(yearstofadeout.seq)){
  segments(x0 = yearstofadeout.seq[i], x1 = yearstofadeout.seq[i], y0 = popsize.30.quants[i, 1], y1 = popsize.30.quants[i, 5])
  segments(x0 = yearstofadeout.seq[i] - 0.1, x1 = yearstofadeout.seq[i] + 0.1, y0 = popsize.30.quants[i, 3], y1 = popsize.30.quants[i, 3])
  segments(x0 = yearstofadeout.seq[i] + .2, x1 = yearstofadeout.seq[i] + .2, y0 = popsize.30.alpha.33.quants[i, 1], y1 = popsize.30.alpha.33.quants[i, 5], col = "red", lwd = 2)
  segments(x0 = yearstofadeout.seq[i] + .1, x1 = yearstofadeout.seq[i] + .3, y0 = popsize.30.alpha.33.quants[i, 3], y1 = popsize.30.alpha.33.quants[i, 3], col = "red", lwd = 2)
}
leg.text1 <- c("Expect 10 years between introductions", "Expect 3 years between introductions")
legend("topright", bty = "n", leg.text1, lty = c(1, 1), col = c("black", "red"), cex = .6)

plot(x = 0, y = 0, xlim = c(0, 16), ylim = c(0, 2000), xlab = "Expected years to reintroduction", ylab = "Simulated Pop Size After 30 years")
abline(h = 1000, col = "grey80", lty = 2)
abline(h = 1500, col = "grey80", lty = 2)
abline(h = 500, col = "grey80", lty = 2)
for(i in 1:length(yearstoreintro)){
  segments(x0 = yearstoreintro[i], x1 = yearstoreintro[i], y0 = popsize.alpha.gamma.33.30.quants[i, 1], y1 = popsize.alpha.gamma.33.30.quants[i, 5])
  segments(x0 = yearstoreintro[i] - 0.1, x1 = yearstoreintro[i] + 0.1, y0 = popsize.alpha.gamma.33.30.quants[i, 3], y1 = popsize.alpha.gamma.33.30.quants[i, 3])
  segments(x0 = yearstoreintro[i] + .2, x1 = yearstoreintro[i] + .2, y0 = popsize.alpha.gamma.1.30.quants[i, 1], y1 = popsize.alpha.gamma.1.30.quants[i, 5], col = "red", lwd = 2)
  segments(x0 = yearstoreintro[i] + .1, x1 = yearstoreintro[i] + .3, y0 = popsize.alpha.gamma.1.30.quants[i, 3], y1 = popsize.alpha.gamma.1.30.quants[i, 3], col = "red", lwd = 2)
}
leg.text2 <- c("Expect 3 years to fade-out", "Expect 10 years to fade-out")
legend("topleft", bty = "n", leg.text2, lty = c(1, 1), col = c("black", "red"), cex = .6)

