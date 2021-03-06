PlotRedWhiteNoise <- function(RedWhiteNoiseSimsOut, 
                              timesteps, 
                              reps)
{
# extract 2.5th and 97.5th quantiles --#
gam.2.alph.05.quants <- which(popsize.gam.2.alph.05[, timesteps - 1] %in% as.numeric(quantile(popsize.gam.2.alph.05[, timesteps - 1], c(0.025, 0.25, 0.5, 0.75, 0.975), type = 3, na.rm = T)))
gam.2.alph.2.quants <- which(popsize.gam.2.alph.2[, timesteps - 1] %in% as.numeric(quantile(popsize.gam.2.alph.2[, timesteps - 1], c(0.025, 0.25, 0.5, 0.75, 0.975), type = 3, na.rm = T)))
gam.2.alph1.quants <- which(popsize.gam.2.alph1[, timesteps - 1] %in% as.numeric(quantile(popsize.gam.2.alph1[, timesteps - 1], c(0.025, 0.25, 0.5, 0.75, 0.975), type = 3, na.rm = T)))

gam.2.alph.05.meds <- which(popsize.gam.2.alph.05[, timesteps - 1] %in% as.numeric(quantile(popsize.gam.2.alph.05[, timesteps - 1], .5, type = 3, na.rm = T)))
gam.2.alph.2.meds <- which(popsize.gam.2.alph.2[, timesteps - 1] %in% as.numeric(quantile(popsize.gam.2.alph.2[, timesteps - 1], .5, type = 3, na.rm = T)))
gam.2.alph1.meds <- which(popsize.gam.2.alph1[, timesteps - 1] %in% as.numeric(quantile(popsize.gam.2.alph1[, timesteps - 1], .5, type = 3, na.rm = T)))

unstr.gam.2.alph.05.quants <- which(unstr.popsize.gam.2.alph.05[, timesteps - 1] %in% as.numeric(quantile(unstr.popsize.gam.2.alph.05[, timesteps - 1], c(0.025, 0.25, 0.5, 0.75, 0.975), type = 3, na.rm = T)))
unstr.gam.2.alph.2.quants <- which(unstr.popsize.gam.2.alph.2[, timesteps - 1] %in% as.numeric(quantile(unstr.popsize.gam.2.alph.2[, timesteps - 1], c(0.025, 0.25, 0.5, 0.75, 0.975), type = 3, na.rm = T)))
unstr.gam.2.alph1.quants <- which(unstr.popsize.gam.2.alph1[, timesteps - 1] %in% as.numeric(quantile(unstr.popsize.gam.2.alph1[, timesteps - 1], c(0.025, 0.25, 0.5, 0.75, 0.975), type = 3, na.rm = T)))

unstr.gam.2.alph.05.meds <- which(unstr.popsize.gam.2.alph.05[, timesteps - 1] %in% as.numeric(quantile(unstr.popsize.gam.2.alph.05[, timesteps - 1], .5, type = 3, na.rm = T)))
unstr.gam.2.alph.2.meds <- which(unstr.popsize.gam.2.alph.2[, timesteps - 1] %in% as.numeric(quantile(unstr.popsize.gam.2.alph.2[, timesteps - 1], .5, type = 3, na.rm = T)))
unstr.gam.2.alph1.meds <- which(unstr.popsize.gam.2.alph1[, timesteps - 1] %in% as.numeric(quantile(unstr.popsize.gam.2.alph1[, timesteps - 1], .5, type = 3, na.rm = T)))

unstr.intros.gam.2.alph.05.quants <- which(unstr.intros.popsize.gam.2.alph.05[, timesteps - 1] %in% as.numeric(quantile(unstr.intros.popsize.gam.2.alph.05[, timesteps - 1], c(0.025, 0.25, 0.5, 0.75, 0.975), type = 3, na.rm = T)))
unstr.intros.gam.2.alph.2.quants <- which(unstr.intros.popsize.gam.2.alph.2[, timesteps - 1] %in% as.numeric(quantile(unstr.intros.popsize.gam.2.alph.2[, timesteps - 1], c(0.025, 0.25, 0.5, 0.75, 0.975), type = 3, na.rm = T)))
unstr.intros.gam.2.alph1.quants <- which(unstr.intros.popsize.gam.2.alph1[, timesteps - 1] %in% as.numeric(quantile(unstr.intros.popsize.gam.2.alph1[, timesteps - 1], c(0.025, 0.25, 0.5, 0.75, 0.975), type = 3, na.rm = T)))

unstr.intros.gam.2.alph.05.meds <- which(unstr.intros.popsize.gam.2.alph.05[, timesteps - 1] %in% as.numeric(quantile(unstr.intros.popsize.gam.2.alph.05[, timesteps - 1], .5, type = 3, na.rm = T)))
unstr.intros.gam.2.alph.2.meds <- which(unstr.intros.popsize.gam.2.alph.2[, timesteps - 1] %in% as.numeric(quantile(unstr.intros.popsize.gam.2.alph.2[, timesteps - 1], .5, type = 3, na.rm = T)))
unstr.intros.gam.2.alph1.meds <- which(unstr.intros.popsize.gam.2.alph1[, timesteps - 1] %in% as.numeric(quantile(unstr.intros.popsize.gam.2.alph1[, timesteps - 1], .5, type = 3, na.rm = T)))

# build boxplot data
boxgrps <- rep(c(1, 2, 3), each = reps)
endsize.alph.05 <- c(popsize.gam.2.alph.05[, timesteps - 1], unstr.popsize.gam.2.alph.05[, timesteps - 1], unstr.intros.popsize.gam.2.alph.05[, timesteps - 1])
endsize.alph.2 <- c(popsize.gam.2.alph.2[, timesteps - 1], unstr.popsize.gam.2.alph.2[, timesteps - 1], unstr.intros.popsize.gam.2.alph.2[, timesteps - 1])
endsize.alph1 <- c(popsize.gam.2.alph1[, timesteps - 1], unstr.popsize.gam.2.alph1[, timesteps - 1], unstr.intros.popsize.gam.2.alph1[, timesteps - 1])



plot.reps <- reps

layout(matrix(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12), nrow = 3, byrow = T))
par(oma = c(2, 5, 2, 1), mar = c(2, 2, 2, 2))
  if(log.y == T)
  {
  plot((popsize.gam.2.alph.05[1, 1:59] + 1) ~ seq(1:59), type = "l", log = "y", ylim = c(1, 7000), xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((popsize.gam.2.alph.05[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in c(1, 2, 4, 5)){
    lines((popsize.gam.2.alph.05[gam.2.alph.05.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((popsize.gam.2.alph.05[gam.2.alph.05.meds, 1:59] + 1) ~ seq(1:59), type = "l", col = "red", lwd = 2)  
  mtext(side = 3, outer = F, line = 1, "Structured disease", las = 1)
  mtext(side = 2, outer = F, line = 5, "Intro each 20 years", las = 0)
  
  plot((unstr.popsize.gam.2.alph.05[1, 1:59] + 1) ~ seq(1:59), type = "l", ylim = c(1, 7000), log = "y", xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((unstr.popsize.gam.2.alph.05[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in c(1, 2, 4, 5)){
    lines((unstr.popsize.gam.2.alph.05[gam.2.alph.05.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((unstr.popsize.gam.2.alph.05[gam.2.alph.05.meds, 1:59] + 1) ~ seq(1:59), type = "l", col = "red", lwd = 2)  
  mtext(side = 3, outer = F, line = 1, "Random disease", las = 1)
  
  plot((unstr.intros.popsize.gam.2.alph.05[1, 1:59] + 1) ~ seq(1:59), log = "y", type = "l", ylim = c(1, 7000), xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((unstr.intros.popsize.gam.2.alph.05[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in c(1, 2, 4, 5)){
    lines((unstr.intros.popsize.gam.2.alph.05[gam.2.alph.05.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((unstr.intros.popsize.gam.2.alph.05[gam.2.alph.05.meds, 1:59] + 1) ~ seq(1:59), type = "l", col = "red", lwd = 2)  
  mtext(side = 3, outer = F, line = 1, "Spillovers only", las = 1)
  
  boxplot((endsize.alph.05 + 1)~ boxgrps, col = c("white", "black", "grey60"), log = "y", ylim = c(1, 7000), bty = "n")
  leg.text <- c("Structured", "Random", "Spillover")
  legend("topleft", leg.text, fill = c("white", "black", "grey60"), bty = "n", pt.cex = 2, cex = .6)
  mtext(side = 3, outer = F, line = 1, "N at year 60", las = 1)
  
  
  ### END FIRST ROW
  
  plot((popsize.gam.2.alph.2[1, 1:59] + 1) ~ seq(1:59), type = "l", log = "y", ylim = c(1, 7000), xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((popsize.gam.2.alph.2[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in c(1, 2, 4, 5)){
    lines((popsize.gam.2.alph.2[gam.2.alph.2.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((popsize.gam.2.alph.2[gam.2.alph.2.meds, 1:59] + 1)~ seq(1:59), type = "l", col = "red", lwd = 2)  
  mtext(side = 2, outer = F, line = 5, "Intro each 5 years", las = 0)
  
  plot((unstr.popsize.gam.2.alph.2[1, 1:59] + 1) ~ seq(1:59), log = "y", type = "l", ylim = c(1, 7000), xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((unstr.popsize.gam.2.alph.2[i, 1:59] + 1) ~ seq(1:59), log = "y", type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in c(1, 2, 4, 5)){
    lines((unstr.popsize.gam.2.alph.2[gam.2.alph.2.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((unstr.popsize.gam.2.alph.2[gam.2.alph.2.meds, 1:59] + 1) ~ seq(1:59), type = "l", col = "red", lwd = 2)  
  
  plot((unstr.intros.popsize.gam.2.alph.2[1, 1:59] + 1) ~ seq(1:59), log = "y", type = "l", ylim = c(1, 7000), xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((unstr.intros.popsize.gam.2.alph.2[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in c(1, 2, 4, 5)){
    lines((unstr.intros.popsize.gam.2.alph.2[gam.2.alph.2.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((unstr.intros.popsize.gam.2.alph.2[gam.2.alph.2.meds, 1:59] + 1) ~ seq(1:59), type = "l", col = "red", lwd = 2)  
  
  boxplot((endsize.alph.2 + 1)~ boxgrps, col = c("white", "black", "grey60"), log = "y", ylim = c(1, 7000))
  leg.text <- c("Structured", "Random", "Spillover")
  legend("top", leg.text, fill = c("white", "black", "grey60"), bty = "n", cex = .8)
  
  ### END SECOND ROW
  
  plot((popsize.gam.2.alph1[1, 1:59] + 1) ~ seq(1:59), type = "l", ylim = c(1, 7000), log = "y", xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((popsize.gam.2.alph1[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in c(1, 2, 4, 5)){
    lines((popsize.gam.2.alph1[gam.2.alph1.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((popsize.gam.2.alph1[gam.2.alph1.meds, 1:59] + 1) ~ seq(1:59), type = "l", col = "red", lwd = 2)  
  mtext(side = 2, outer = F, line = 5, "Intro each year", las = 0)
  
  
  plot((unstr.popsize.gam.2.alph1[1, 1:59] + 1) ~ seq(1:59), type = "l", log = "y", ylim = c(1, 7000), xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((unstr.popsize.gam.2.alph1[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in c(1, 2, 4, 5)){
    lines((unstr.popsize.gam.2.alph1[gam.2.alph1.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((unstr.popsize.gam.2.alph1[gam.2.alph1.meds, 1:59] + 1) ~ seq(1:59), type = "l", col = "red", lwd = 2)  
  
  plot((unstr.intros.popsize.gam.2.alph1[1, 1:59] + 1) ~ seq(1:59), type = "l", log = "y", ylim = c(1, 7000), xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((unstr.intros.popsize.gam.2.alph1[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25), log = "y")
  }
  for(j in c(1, 2, 4, 5)){
    lines((unstr.intros.popsize.gam.2.alph1[gam.2.alph1.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((unstr.intros.popsize.gam.2.alph1[gam.2.alph1.meds, 1:59] + 1) ~ seq(1:59), type = "l", col = "red", lwd = 2)  
  
  boxplot((endsize.alph1 + 1)~ boxgrps, col = c("white", "black", "grey60"), log = "y", ylim = c(1, 7000))
  leg.text <- c("Structured", "Random", "Spillover")
  legend("top", leg.text, fill = c("white", "black", "grey60"), bty = "n", cex = .8)
  
}

else 
  { # y not logged
    
  plot((popsize.gam.2.alph.05[1, 1:59] + 1) ~ seq(1:59), type = "l", ylim = c(1, 1600), xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((popsize.gam.2.alph.05[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in c(1, 2, 4, 5)){
    lines((popsize.gam.2.alph.05[gam.2.alph.05.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((popsize.gam.2.alph.05[gam.2.alph.05.meds, 1:59] + 1) ~ seq(1:59), type = "l", col = "red", lwd = 2)  
  mtext(side = 3, outer = F, line = 1, "Structured disease", las = 1)
  mtext(side = 2, outer = F, line = 5, "Intro each 20 years", las = 0)
  
  plot((unstr.popsize.gam.2.alph.05[1, 1:59] + 1) ~ seq(1:59), type = "l", ylim = c(1, 1600), xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((unstr.popsize.gam.2.alph.05[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in c(1, 2, 4, 5)){
    lines((unstr.popsize.gam.2.alph.05[gam.2.alph.05.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((unstr.popsize.gam.2.alph.05[gam.2.alph.05.meds, 1:59] + 1) ~ seq(1:59), type = "l", col = "red", lwd = 2)  
  mtext(side = 3, outer = F, line = 1, "Random disease", las = 1)
  
  plot((unstr.intros.popsize.gam.2.alph.05[1, 1:59] + 1) ~ seq(1:59), type = "l", ylim = c(1, 1600), xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((unstr.intros.popsize.gam.2.alph.05[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in c(1, 2, 4, 5)){
    lines((unstr.intros.popsize.gam.2.alph.05[gam.2.alph.05.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((unstr.intros.popsize.gam.2.alph.05[gam.2.alph.05.meds, 1:59] + 1) ~ seq(1:59), type = "l", col = "red", lwd = 2)  
  mtext(side = 3, outer = F, line = 1, "Spillovers only", las = 1)
  
  boxplot((endsize.alph.05 + 1)~ boxgrps, col = c("white", "black", "grey60"), ylim = c(1, 1600), bty = "n")
  leg.text <- c("Structured", "Random", "Spillover")
  legend("topleft", leg.text, fill = c("white", "black", "grey60"), bty = "n", pt.cex = 2, cex = .6)
  mtext(side = 3, outer = F, line = 1, "N at year 60", las = 1)
  
  
  ### END FIRST ROW
  
  plot((popsize.gam.2.alph.2[1, 1:59] + 1) ~ seq(1:59), type = "l", ylim = c(1, 1600), xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((popsize.gam.2.alph.2[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in c(1, 2, 4, 5)){
    lines((popsize.gam.2.alph.2[gam.2.alph.2.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((popsize.gam.2.alph.2[gam.2.alph.2.meds, 1:59] + 1)~ seq(1:59), type = "l", col = "red", lwd = 2)  
  mtext(side = 2, outer = F, line = 5, "Intro each 5 years", las = 0)
  
  plot((unstr.popsize.gam.2.alph.2[1, 1:59] + 1) ~ seq(1:59),  type = "l", ylim = c(1, 1600), xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((unstr.popsize.gam.2.alph.2[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in c(1, 2, 4, 5)){
    lines((unstr.popsize.gam.2.alph.2[gam.2.alph.2.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((unstr.popsize.gam.2.alph.2[gam.2.alph.2.meds, 1:59] + 1) ~ seq(1:59), type = "l", col = "red", lwd = 2)  
  
  plot((unstr.intros.popsize.gam.2.alph.2[1, 1:59] + 1) ~ seq(1:59), type = "l", ylim = c(1, 1600), xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((unstr.intros.popsize.gam.2.alph.2[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in c(1, 2, 4, 5)){
    lines((unstr.intros.popsize.gam.2.alph.2[gam.2.alph.2.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((unstr.intros.popsize.gam.2.alph.2[gam.2.alph.2.meds, 1:59] + 1) ~ seq(1:59), type = "l", col = "red", lwd = 2)  
  
  boxplot((endsize.alph.2 + 1)~ boxgrps, col = c("white", "black", "grey60"), ylim = c(1, 1600))
  leg.text <- c("Structured", "Random", "Spillover")
  legend("top", leg.text, fill = c("white", "black", "grey60"), bty = "n", cex = .8)
  
  ### END SECOND ROW
  
  plot((popsize.gam.2.alph1[1, 1:59] + 1) ~ seq(1:59), type = "l", ylim = c(1, 1600), xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((popsize.gam.2.alph1[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in c(1, 2, 4, 5)){
    lines((popsize.gam.2.alph1[gam.2.alph1.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((popsize.gam.2.alph1[gam.2.alph1.meds, 1:59] + 1) ~ seq(1:59), type = "l", col = "red", lwd = 2)  
  mtext(side = 2, outer = F, line = 5, "Intro each year", las = 0)
  
  
  plot((unstr.popsize.gam.2.alph1[1, 1:59] + 1) ~ seq(1:59), type = "l", ylim = c(1, 1600), xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((unstr.popsize.gam.2.alph1[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
  }
  for(j in c(1, 2, 4, 5)){
    lines((unstr.popsize.gam.2.alph1[gam.2.alph1.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((unstr.popsize.gam.2.alph1[gam.2.alph1.meds, 1:59] + 1) ~ seq(1:59), type = "l", col = "red", lwd = 2)  
  
  plot((unstr.intros.popsize.gam.2.alph1[1, 1:59] + 1) ~ seq(1:59), type = "l", ylim = c(1, 1600), xlab = "year", ylab = "Population size", main = "", bty = "n")
  for(i in 2:plot.reps){
    lines((unstr.intros.popsize.gam.2.alph1[i, 1:59] + 1) ~ seq(1:59), type = "l", col = rgb(.35, .35, .35, alpha = .25), log = "y")
  }
  for(j in c(1, 2, 4, 5)){
    lines((unstr.intros.popsize.gam.2.alph1[gam.2.alph1.quants[j], 1:59] + 1) ~ seq(1:59), type = "l", col = "black", lwd = 2)  
  }
  lines((unstr.intros.popsize.gam.2.alph1[gam.2.alph1.meds, 1:59] + 1) ~ seq(1:59), type = "l", col = "red", lwd = 2)  
  
  boxplot((endsize.alph1 + 1)~ boxgrps, col = c("white", "black", "grey60"), ylim = c(1, 1600))
  leg.text <- c("Structured", "Random", "Spillover")
  legend("top", leg.text, fill = c("white", "black", "grey60"), bty = "n", cex = .8)
  
  }

#   # Time to extinction
#   gam.2.alph.05.fadeout <- gam.2.alph.2.fadeout <- gam.2.alph1.fadeout <- rep(NA, reps)
#   unstr.gam.2.alph.05.fadeout <- unstr.gam.2.alph.2.fadeout <- unstr.gam.2.alph1.fadeout <- rep(NA, reps)
# unstr.intros.gam.2.alph.05.fadeout <- unstr.intros.gam.2.alph.2.fadeout <- unstr.intros.gam.2.alph1.fadeout <- rep(NA, reps)
# for(i in 1:reps){
#     gam.2.alph.05.fadeout[i] <- min(which(popsize.gam.2.alph.05[i, -timesteps] < 1))
#     gam.2.alph.2.fadeout[i] <- min(which(popsize.gam.2.alph.2[i, -timesteps] < 1))
#     gam.2.alph1.fadeout[i] <- min(which(popsize.gam.2.alph1[i, -timesteps] < 1))
#  
#     unstr.gam.2.alph.05.fadeout[i] <- min(which(unstr.popsize.gam.2.alph.05[i, -timesteps] < 1))
#     unstr.gam.2.alph.2.fadeout[i] <- min(which(unstr.popsize.gam.2.alph.2[i, -timesteps] < 1))
#     unstr.gam.2.alph1.fadeout[i] <- min(which(unstr.popsize.gam.2.alph1[i, -timesteps] < 1))
# 
#     unstr.intros.gam.2.alph.05.fadeout[i] <- min(which(unstr.intros.popsize.gam.2.alph.05[i, -timesteps] < 1))
#     unstr.intros.gam.2.alph.2.fadeout[i] <- min(which(unstr.intros.popsize.gam.2.alph.2[i, -timesteps] < 1))
#     unstr.intros.gam.2.alph1.fadeout[i] <- min(which(unstr.intros.popsize.gam.2.alph1[i, -timesteps] < 1))
# }
# #  hist(gam.2.alph.05.fadeout)
# #  hist(gam.2.alph.2.fadeout)
#   par(mfrow = c(3, 1))
#   hist(gam.2.alph1.fadeout, xlim = c(0, 60))
# 
# #  hist(unstr.gam.2.alph.05.fadeout)
# #  hist(unstr.gam.2.alph.2.fadeout)
#   hist(unstr.gam.2.alph1.fadeout, xlim = c(0, 60))
# 
# #  hist(unstr.intros.gam.2.alph.05.fadeout)
# #  hist(unstr.intros.gam.2.alph.2.fadeout)
#   hist(unstr.intros.gam.2.alph1.fadeout, xlim = c(0, 60))
# 
#   unstr.intros.gam.2.alph.05.quants <- which(unstr.intros.popsize.gam.2.alph.05[, timesteps - 1] %in% as.numeric(quantile(unstr.intros.popsize.gam.2.alph.05[, timesteps - 1], c(0.025, 0.25, 0.5, 0.75, 0.975), type = 3, na.rm = T)))
#   unstr.intros.gam.2.alph.2.quants <- which(unstr.intros.popsize.gam.2.alph.2[, timesteps - 1] %in% as.numeric(quantile(unstr.intros.popsize.gam.2.alph.2[, timesteps - 1], c(0.025, 0.25, 0.5, 0.75, 0.975), type = 3, na.rm = T)))
#   unstr.intros.gam.2.alph1.quants <- which(unstr.intros.popsize.gam.2.alph1[, timesteps - 1] %in% as.numeric(quantile(unstr.intros.popsize.gam.2.alph1[, timesteps - 1], c(0.025, 0.25, 0.5, 0.75, 0.975), type = 3, na.rm = T)))
# 
#   unstr.intros.gam.2.alph.05.meds <- which(unstr.intros.popsize.gam.2.alph.05[, timesteps - 1] %in% as.numeric(quantile(unstr.intros.popsize.gam.2.alph.05[, timesteps - 1], .5, type = 3, na.rm = T)))
#   unstr.intros.gam.2.alph.2.meds <- which(unstr.intros.popsize.gam.2.alph.2[, timesteps - 1] %in% as.numeric(quantile(unstr.intros.popsize.gam.2.alph.2[, timesteps - 1], .5, type = 3, na.rm = T)))
#   unstr.intros.gam.2.alph1.meds <- which(unstr.intros.popsize.gam.2.alph1[, timesteps - 1] %in% as.numeric(quantile(unstr.intros.popsize.gam.2.alph1[, timesteps - 1], .5, type = 3, na.rm = T)))

}