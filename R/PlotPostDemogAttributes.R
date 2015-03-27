PlotPostDemogAttributes <- function(reps, joint.posterior.coda, samples.to.draw, 
                                    tot.chains, posterior.names, sex.ratio, 
                                    age.class.ind)
{

  healthy.leslie.list <- infected.leslie.list <- spillover.leslie.list.1  <- spillover.leslie.list.5 <- vector("list", length = reps)
  healthy.eigenval1 <- infected.eigenval1 <- spillover.eigenval1.1 <- spillover.eigenval1.5 <-  rep(NA, reps)
  juvsurv.elast.he <- juvsurv.elast.inf <- juvsurv.elast.sp.1 <- juvsurv.elast.sp.5 <- rep(NA, reps)
  adsurv.elast.he <- adsurv.elast.inf <- adsurv.elast.sp.1  <- adsurv.elast.sp.5 <- rep(NA, reps)
  fecund.elast.he <- fecund.elast.inf <- fecund.elast.sp.1 <- fecund.elast.sp.5 <- rep(NA, reps)
  oldage.elast.he <- oldage.elast.inf <- oldage.elast.sp.1 <- oldage.elast.sp.5 <- rep(NA, reps)
  age.struct.he <- age.struct.inf <- age.struct.sp.1 <- age.struct.sp.5 <- matrix(NA, nrow = 19, ncol = reps)
  # joint.posterior.coda <- ipm11.coda
  # age.class.ind <- c(1, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6) 

  for(i in 1:reps){
    healthy.leslie.list[[i]] <- UpdateLeslieFun(current.state = "healthy", samples.to.draw = samples.to.draw, tot.chains = tot.chains, 
                                                joint.posterior.coda = joint.posterior.coda, 
                                                posterior.name = posterior.names, sex.ratio = sex.ratio, intro.cost = 0)
    infected.leslie.list[[i]] <- UpdateLeslieFun(current.state = "infected", samples.to.draw = samples.to.draw, tot.chains = tot.chains, 
                                                 joint.posterior.coda = joint.posterior.coda, 
                                                 posterior.name = posterior.names, sex.ratio = sex.ratio, intro.cost = 0)
    spillover.leslie.list.1[[i]] <- UpdateLeslieFun(current.state = "spillover", samples.to.draw = samples.to.draw, tot.chains = tot.chains, 
                                                    joint.posterior.coda = joint.posterior.coda, 
                                                    posterior.name = posterior.names, sex.ratio = sex.ratio, intro.cost = .3)
    spillover.leslie.list.5[[i]] <- UpdateLeslieFun(current.state = "spillover", samples.to.draw = samples.to.draw, tot.chains = tot.chains, 
                                                    joint.posterior.coda = joint.posterior.coda, 
                                                    posterior.name = posterior.names, sex.ratio = sex.ratio, intro.cost = .5)
    healthy.eigenval1[i] <- Re(eigen(healthy.leslie.list[[i]])$values[1]) # strip off only real part of eigenvalue
    infected.eigenval1[i] <- Re(eigen(infected.leslie.list[[i]])$values[1])
    spillover.eigenval1.1[i] <- Re(eigen(spillover.leslie.list.1[[i]])$values[1])
    spillover.eigenval1.5[i] <- Re(eigen(spillover.leslie.list.5[[i]])$values[1])
    he.eigen.rescale <- sum(Re(eigen(healthy.leslie.list[[i]])$vectors[, 1]))
    inf.eigen.rescale <- sum(Re(eigen(infected.leslie.list[[i]])$vectors[, 1]))
    sp.eigen.rescale.1 <- sum(Re(eigen(spillover.leslie.list.1[[i]])$vectors[, 1]))
    sp.eigen.rescale.5 <- sum(Re(eigen(spillover.leslie.list.5[[i]])$vectors[, 1]))
    age.struct.he[, i] <- Re(eigen(healthy.leslie.list[[i]])$vectors[, 1]) / he.eigen.rescale
    age.struct.inf[, i] <- Re(eigen(infected.leslie.list[[i]])$vectors[, 1]) / inf.eigen.rescale
    age.struct.sp.1[, i] <- Re(eigen(spillover.leslie.list.1[[i]])$vectors[, 1]) / sp.eigen.rescale.1
    age.struct.sp.5[, i] <- Re(eigen(spillover.leslie.list.5[[i]])$vectors[, 1]) / sp.eigen.rescale.5
  
    he.sens <- sensitivity(healthy.leslie.list[[i]])
    he.elast <- (1 / Re(eigen(healthy.leslie.list[[i]])$values[1])) * he.sens * healthy.leslie.list[[i]]
    fecund.elast.he[i] <- sum(he.elast[1, ])
    juvsurv.elast.he[i] <- sum(he.elast[2, 1], he.elast[3, 2])
    adsurv.elast.he[i] <- sum(he.elast[4, 3], he.elast[5, 4], he.elast[6, 5], he.elast[7, 6], 
                              he.elast[8, 7], he.elast[9, 8], he.elast[10, 9], he.elast[11, 10], 
                              he.elast[12, 11], he.elast[13, 12], he.elast[14, 13], 
                              he.elast[15, 14], he.elast[16, 15], he.elast[17, 16], he.elast[18, 17])
    oldage.elast.he[i] <- sum(he.elast[13, 12], he.elast[14, 13], he.elast[15, 14], 
                              he.elast[16, 15], he.elast[17, 16], he.elast[18, 17])
  
  inf.sens <- sensitivity(infected.leslie.list[[i]])
  inf.elast <- (1 / Re(eigen(infected.leslie.list[[i]])$values[1])) * inf.sens * infected.leslie.list[[i]]
  fecund.elast.inf[i] <- sum(inf.elast[1, ])
  juvsurv.elast.inf[i] <- sum(inf.elast[2, 1], inf.elast[3, 2])
  adsurv.elast.inf[i] <- sum(inf.elast[4, 3], inf.elast[5, 4], inf.elast[6, 5], inf.elast[7, 6], 
                             inf.elast[8, 7], inf.elast[9, 8], inf.elast[10, 9], inf.elast[11, 10],
                             inf.elast[12, 11], inf.elast[13, 12], inf.elast[14, 13], 
                             inf.elast[15, 14], inf.elast[16, 15], inf.elast[17, 16], inf.elast[18, 17])
  oldage.elast.inf[i] <- sum(inf.elast[13, 12], inf.elast[14, 13], inf.elast[15, 14], 
                             inf.elast[16, 15], inf.elast[17, 16], inf.elast[18, 17])
  
  sp.sens.1 <- sensitivity(spillover.leslie.list.1[[i]])
  sp.elast.1 <- (1 / Re(eigen(spillover.leslie.list.1[[i]])$values[1])) * sp.sens.1 * spillover.leslie.list.1[[i]]
  fecund.elast.sp.1[i] <- sum(sp.elast.1[1, ])
  juvsurv.elast.sp.1[i] <- sum(sp.elast.1[2, 1], sp.elast.1[3, 2])
  adsurv.elast.sp.1[i] <- sum(sp.elast.1[4, 3], sp.elast.1[5, 4], sp.elast.1[6, 5], 
                              sp.elast.1[7, 6], sp.elast.1[8, 7], sp.elast.1[9, 8], 
                              sp.elast.1[10, 9], sp.elast.1[11, 10], sp.elast.1[12, 11], 
                              sp.elast.1[13, 12], sp.elast.1[14, 13], sp.elast.1[15, 14], 
                              sp.elast.1[16, 15], sp.elast.1[17, 16], sp.elast.1[18, 17])
  oldage.elast.sp.1[i] <- sum(sp.elast.1[13, 12], sp.elast.1[14, 13], sp.elast.1[15, 14], 
                              sp.elast.1[16, 15], sp.elast.1[17, 16], sp.elast.1[18, 17])
  
  sp.sens.5 <- sensitivity(spillover.leslie.list.5[[i]])
  sp.elast.5 <- (1 / Re(eigen(spillover.leslie.list.5[[i]])$values[1])) * sp.sens.5 * spillover.leslie.list.5[[i]]
  fecund.elast.sp.5[i] <- sum(sp.elast.5[1, ])
  juvsurv.elast.sp.5[i] <- sum(sp.elast.5[2, 1], sp.elast.5[3, 2])
  adsurv.elast.sp.5[i] <- sum(sp.elast.5[4, 3], sp.elast.5[5, 4], sp.elast.5[6, 5], 
                              sp.elast.5[7, 6], sp.elast.5[8, 7], sp.elast.5[9, 8], 
                              sp.elast.5[10, 9], sp.elast.5[11, 10], sp.elast.5[12, 11], 
                              sp.elast.5[13, 12], sp.elast.5[14, 13], sp.elast.5[15, 14], 
                              sp.elast.5[16, 15], sp.elast.5[17, 16], sp.elast.5[18, 17])
  oldage.elast.sp.5[i] <- sum(sp.elast.5[13, 12], sp.elast.5[14, 13], sp.elast.5[15, 14], 
                              sp.elast.5[16, 15], sp.elast.5[17, 16], sp.elast.5[18, 17])
  
  print(i)
}

# 95% intervals on lambda 
healthy.cred <- quantile(healthy.eigenval1, c(0.25, 0.5, 0.75))
infected.cred <- quantile(infected.eigenval1, c(0.25, 0.5, 0.75))
spillover.cred.1 <- quantile(spillover.eigenval1.1, c(0.25, 0.5, 0.75))

# 95% intervals for fecundity, juvenile survival, and adult survival in each environment (healthy and infected)
fecund.he <- quantile(fecund.elast.he, c(0.05, 0.5, 0.975))
fecund.inf <- quantile(fecund.elast.inf, c(.05, 0.5, 0.975))
fecund.sp <- quantile(fecund.elast.sp.1, c(.05, 0.5, 0.975))
fecund.sp.5 <- quantile(fecund.elast.sp.5, c(.05, 0.5, 0.975))
juvsurv.he <- quantile(juvsurv.elast.he, c(0.05, 0.5, 0.975))
juvsurv.inf <- quantile(juvsurv.elast.inf, c(0.05, 0.5, 0.975))
juvsurv.sp <- quantile(juvsurv.elast.sp.1, c(0.05, 0.5, 0.975))
juvsurv.sp.5 <- quantile(juvsurv.elast.sp.5, c(0.05, 0.5, 0.975))
adsurv.he <- quantile(adsurv.elast.he, c(0.05, 0.5, 0.975))
adsurv.inf <- quantile(adsurv.elast.inf, c(0.05, 0.5, 0.975))
adsurv.sp <- quantile(adsurv.elast.sp.1, c(0.05, 0.5, 0.975))
adsurv.sp.5 <- quantile(adsurv.elast.sp.5, c(0.05, 0.5, 0.975))

# 95% intervals for ratio of old-age elasticity to all adsurv elasticity 
old.ad.ratio.he <- quantile(oldage.elast.he / adsurv.elast.he, c(0.05, 0.25, 0.5, 0.75, 0.975))
old.ad.ratio.inf <- quantile(oldage.elast.inf / adsurv.elast.inf, c(0.05, 0.25, 0.5, 0.75, 0.975))
old.ad.ratio.sp1 <- quantile(oldage.elast.sp.1 / adsurv.elast.sp.1, c(0.05, 0.25, 0.5, 0.75, 0.975))
old.ad.ratio.sp5 <- quantile(oldage.elast.sp.5 / adsurv.elast.sp.5, c(0.05, 0.25, 0.5, 0.75, 0.975))

# par(mfrow = c(1, 1), mar = c(4, 12, 2, 2), las = 1, cex.lab = 1.0)
# plot(c(1, 1) ~ c(fecund.he[1], fecund.he[3]), lty = 1, xlim = c(0, 1), ylim = c(0, 13), type = "l", xlab = expression(paste("Elasticity of ", lambda, " to rate", sep = "")), ylab = "", yaxt = "n", lwd = 2)
# lines(c(2, 2) ~ c(fecund.sp[1], fecund.sp[3]), lty = 2, col = "red", lwd = 2)
# lines(c(3, 3) ~ c(fecund.sp.5[1], fecund.sp.5[3]), lty = 2, col = "purple", lwd = 2)
# lines(c(4, 4) ~ c(fecund.inf[1], fecund.inf[3]), lty = 2, col = "blue", lwd = 2)
# lines(c(5, 5) ~ c(juvsurv.he[1], juvsurv.he[3]), lty = 1, col = "black", lwd = 2)
# lines(c(6, 6) ~ c(juvsurv.sp[1], juvsurv.sp[3]), lty = 2, col = "red", lwd = 2)
# lines(c(7, 7) ~ c(juvsurv.sp.5[1], juvsurv.sp.5[3]), lty = 2, col = "purple", lwd = 2)
# lines(c(8, 8) ~ c(juvsurv.inf[1], juvsurv.inf[3]), lty = 2, col = "blue", lwd = 2)
# lines(c(9, 9) ~ c(adsurv.he[1], adsurv.he[3]), lty = 1, col = "black", lwd = 2)
# lines(c(10, 10) ~ c(adsurv.sp[1], adsurv.sp[3]), lty = 2, col = "red", lwd = 2)
# lines(c(11, 11) ~ c(adsurv.sp.5[1], adsurv.sp.5[3]), lty = 2, col = "purple", lwd = 2)
# lines(c(12, 12) ~ c(adsurv.inf[1], adsurv.inf[3]), lty = 2, col = "blue", lwd = 2)
# lines(c(0.75, 1.25) ~ c(fecund.he[2], fecund.he[2]), lty = 1, col = "black", lwd = 2)
# lines(c(1.75, 2.25) ~ c(fecund.sp[2], fecund.sp[2]), lty = 2, col = "red", lwd = 2)
# lines(c(2.75, 3.25) ~ c(fecund.sp.5[2], fecund.sp.5[2]), lty = 2, col = "purple", lwd = 2)
# lines(c(3.75, 4.25) ~ c(fecund.inf[2], fecund.inf[2]), lty = 2, col = "blue", lwd = 2)
# lines(c(4.75, 5.25) ~ c(juvsurv.he[2], juvsurv.he[2]), lty = 1, col = "black", lwd = 2)
# lines(c(5.75, 6.25) ~ c(juvsurv.sp[2], juvsurv.sp[2]), lty = 2, col = "red", lwd = 2)
# lines(c(6.75, 7.25) ~ c(juvsurv.sp.5[2], juvsurv.sp.5[2]), lty = 2, col = "purple", lwd = 2)
# lines(c(7.75, 8.25) ~ c(juvsurv.inf[2], juvsurv.inf[2]), lty = 2, col = "blue", lwd = 2)
# lines(c(8.75, 9.25) ~ c(adsurv.he[2], adsurv.he[2]), lty = 1, col = "black", lwd = 2)
# lines(c(9.75, 10.25) ~ c(adsurv.sp[2], adsurv.sp[2]), lty = 2, col = "red", lwd = 2)
# lines(c(10.75, 11.25) ~ c(adsurv.sp.5[2], adsurv.sp.5[2]), lty = 2, col = "purple", lwd = 2)
# lines(c(11.75, 12.25) ~ c(adsurv.inf[2], adsurv.inf[2]), lty = 2, col = "blue", lwd = 2)
# axis(side = 2, at = c(1:12), cex.axis = 1.0, labels = c("Fecundity (he)", "Fecundity (sp .5)", "Fecundity (sp .1)", "Fecundity (pers)", "Juvenile surv (he)", "Juvenile surv (sp .5)", "Juvenile surv (sp .1)", "Juvenile surv (pers)", "Adult surv (he)", "Adult surv (sp .5)", "Adult surv (sp .1)", "Adult surv (pers)"))
# 
# k <- hist(c(healthy.eigenval1, infected.eigenval1), breaks = 15, plot = F)
# # par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.5)
# hist(healthy.eigenval1, xlim = c(min(min(healthy.eigenval1), min(infected.eigenval1)), max(max(healthy.eigenval1), max(infected.eigenval1))), breaks = k$breaks, xlab = expression(lambda), main = "Healthy", col = "grey80")
# abline(v = 1, col = "red", lwd = 3)
# abline(v = 1, col = "red", lwd = 3)
# hist(infected.eigenval1, xlim = c(min(min(healthy.eigenval1), min(infected.eigenval1)), max(max(healthy.eigenval1), max(infected.eigenval1))), breaks = k$breaks, ylab = "", xlab = expression(lambda), main = "Persistence", col = "grey80")
# abline(v = 1, col = "red", lwd = 3)

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
age.bounds.he <- age.bounds.inf <- age.bounds.sp.1 <- age.bounds.sp.5 <- matrix(NA, nrow = 18, ncol = 3)
for(i in 1:18){
  age.bounds.he[i, ] <- quantile(abs(age.struct.he[i, ]), c(0.025, 0.5, 0.975))
  age.bounds.inf[i, ] <- quantile(abs(age.struct.inf[i, ]), c(0.025, 0.5, 0.975))
  age.bounds.sp.1[i, ] <- quantile(abs(age.struct.sp.1[i, ]), c(0.025, 0.5, 0.975))
  age.bounds.sp.5[i, ] <- quantile(abs(age.struct.sp.5[i, ]), c(0.025, 0.5, 0.975))
}


# gen times
# require(MASS)
# reps <- 100
gen.time.test.he <- rep(NA, reps)
gen.time.test.inf <- rep(NA, reps)
for(i in 1:reps){
  gen.time.test.he[i] <- GenTimeFun(current.state = "healthy", sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost = .3)$gen.time
  gen.time.test.inf[i] <- GenTimeFun(current.state = "infected", sex.ratio, samples.to.draw, tot.chains, joint.posterior.coda, posterior.names, intro.cost = .3)$gen.time
}

#par(mfrow = c(1, 4), cex.axis = 1, cex.lab = 1, oma = c(3, 7, 0, 2), mar = c(2, 4, 1, 1), las = 2)
par(mfrow = c(2, 2), cex.axis = 1, cex.lab = 1.2, oma = c(3, 7, 0, 2), mar = c(4, 4, 4, 1), las = 1)
plot(c(1, 1) ~ c(fecund.he[1], fecund.he[3]), lty = 1, xlim = c(0, 1), ylim = c(0, 7), type = "l", xlab = expression(paste("Elasticity of ", lambda, " to rate", sep = "")), ylab = "", yaxt = "n", lwd = 2)
lines(c(2, 2) ~ c(fecund.inf[1], fecund.inf[3]), lty = 1, col = "grey60", lwd = 2)
lines(c(3, 3) ~ c(juvsurv.he[1], juvsurv.he[3]), lty = 1, col = "black", lwd = 2)
lines(c(4, 4) ~ c(juvsurv.inf[1], juvsurv.inf[3]), lty = 1, col = "grey60", lwd = 2)
lines(c(5, 5) ~ c(adsurv.he[1], adsurv.he[3]), lty = 1, col = "black", lwd = 2)
lines(c(6, 6) ~ c(adsurv.inf[1], adsurv.inf[3]), lty = 1, col = "grey60", lwd = 2)
lines(c(0.75, 1.25) ~ c(fecund.he[2], fecund.he[2]), lty = 1, col = "black", lwd = 2)
lines(c(1.75, 2.25) ~ c(fecund.inf[2], fecund.inf[2]), lty = 1, col = "grey60", lwd = 2)
lines(c(2.75, 3.25) ~ c(juvsurv.he[2], juvsurv.he[2]), lty = 1, col = "black", lwd = 2)
lines(c(3.75, 4.25) ~ c(juvsurv.inf[2], juvsurv.inf[2]), lty = 1, col = "grey60", lwd = 2)
lines(c(4.75, 5.25) ~ c(adsurv.he[2], adsurv.he[2]), lty = 1, col = "black", lwd = 2)
lines(c(5.75, 6.25) ~ c(adsurv.inf[2], adsurv.inf[2]), lty = 1, col = "grey60", lwd = 2)
axis(side = 2, at = c(1:6), cex.axis = 1.2, labels = c("Fecundity (he)", "Fecundity (pers)", "Juvenile surv (he)", "Juvenile surv (pers)", "Adult surv (he)",  "Adult surv (pers)"))
leg.text <- c("healthy", "infected")
legend("bottomright", leg.text, col = c("black", "grey60"), lwd = c(2, 2), lty = c(1, 1), bty = "n", cex = 1.2)
#mtext(side = 1, line = 4, outer = F, cex = 1, expression(paste("Elasticity of ", lambda, " to vital rate", sep = "")), las = 0)

plot(density(healthy.eigenval1), main = "", ylim = c(0, 20), xlim = c(0.7, 1.3), lwd = 2, lty = 1, xlab = expression(paste(lambda)))
lines(density(infected.eigenval1), col = "grey60", lwd = 2, lty = 2)
#  lines(density(spillover.eigenval1.1), col = "purple", lwd = 2, lty = 1)
#  lines(density(spillover.eigenval1.5), col = "red", lwd = 2, lty = 1)
abline(v = 1, col = "grey40", lty = 3)
leg.text <- c("healthy", "infected")
legend(y = 20, x = .7, leg.text, col = c("black", "grey60"), lwd = c(2, 2), lty = c(1, 2), bty = "n", cex = 1.2)
#mtext(side = 1, line = 4, outer = F, cex = 1, expression(paste(lambda)), las = 1)

plot(age.bounds.he[1, c(1, 3)] ~ c(1, 1), lty = 1, type = "l", xlim = c(0.5, 18.5), ylim = c(0, .2), lwd = 2, xlab = "age (years)", ylab = "Proportion of population")
segments(x0 = 1 + .3, x1 = 1 + .3, y0 = age.bounds.inf[1, 1], lty = 1, y1 = age.bounds.inf[1, 2], lwd = 2, col = "grey60")
for(i in 2: dim(age.bounds.he)[1]){
  segments(x0 = i, x1 = i, y0 = age.bounds.he[i, 1], y1 = age.bounds.he[i, 2], lwd = 2, lty = 1)
  segments(x0 = i + .3, x1 = i + .3, y0 = age.bounds.inf[i, 1], y1 = age.bounds.inf[i, 2], lwd = 2, col = "grey60", lty = 1)
  #    segments(x0 = i + .5, x1 = i + .5, y0 = age.bounds.sp.1[i, 1], y1 = age.bounds.sp.1[i, 2], lwd = 2, col = "purple")
  #    segments(x0 = i + .7, x1 = i + .7, y0 = age.bounds.sp.5[i, 1], y1 = age.bounds.sp.5[i, 2], lwd = 2, col = "red")
}
leg.text <- c("healthy", "infected")
legend(y = 0.2, x = 12.5, leg.text, col = c("black", "grey60"), lwd = c(2, 2), lty = c(1, 1), bty = "n", cex = 1.2)
#mtext(side = 1, line = 4, outer = F, cex = 1, expression("Ewe age"), las = 1)

plot(x = 1, y = 1, pch = 0, cex = 0, xlim = c(0.5, 9), ylim = c(0, 1.2), xlab = "Generation time (yrs)", ylab = "Density")
lines(density(gen.time.test.he), lty = 1, lwd = 2)
lines(density(gen.time.test.inf), col = "grey60", lwd = 2, lty = 2)
leg.text <- c("healthy", "infected")
legend(y = 1.2, x = 7.5, leg.text, col = c("black", "grey60"), lwd = c(2, 2), lty = c(1, 2), bty = "n", cex = 1.2)
#mtext(side = 1, line = 4, outer = F, cex = 1, expression("Generation time (yrs)"), las = 1)

}
