#---------------------------------------------------------#
#-- checks for changepoints in bighorn pop trajectories --#
#---------------------------------------------------------#
require(SiZer)
require(bfast)

compd.data <- read.csv("./Data/compiled_data_summary_130919.csv", header = T, sep = ",")
aso <- subset(compd.data, Pop == "Asotin", select = c("year", "PopEst"))
aso <- aso[complete.cases(aso), ]
bb <- subset(compd.data, Pop == "BlackButte", select = c("year", "PopEst"))
bb <- bb[complete.cases(bb), ]
bc <- subset(compd.data, Pop == "BigCanyon", select = c("year", "PopEst"))
bc <- bc[complete.cases(bc), ]
im <- subset(compd.data, Pop == "Imnaha", select = c("year", "PopEst"))
im <- im[complete.cases(im), ]
lo <- subset(compd.data, Pop == "Lostine", select = c("year", "PopEst"))
lo <- lo[complete.cases(lo), ]
lhc <- subset(compd.data, Pop == "LowerHells", select = c("year", "PopEst"))
lhc <- lhc[complete.cases(lhc), ]
mv <- subset(compd.data, Pop == "MountainView" & year >= 1994, select = c("year", "PopEst"))
mv <- mv[complete.cases(mv), ]
mu <- subset(compd.data, Pop == "MuirCreek", select = c("year", "PopEst"))
mu <- mu[complete.cases(mu), ]
my <- subset(compd.data, Pop == "MyersCreek", select = c("year", "PopEst"))
my <- my[complete.cases(my), ]
rb <- subset(compd.data, Pop == "Redbird" & year >= 1992, select = c("year", "PopEst"))
rb <- rb[complete.cases(rb), ]
uhc <- subset(compd.data, Pop == "UHCOR")
uhc <- subset(compd.data, Pop == "UHCOR", select = c("year", "PopEst", "r_est"))
uhc <- uhc[complete.cases(uhc), ]
we <- subset(compd.data, Pop == "Wenaha", select = c("year", "PopEst"))
we <- we[complete.cases(we), ]

aso.pw <- piecewise.linear(x = aso$year, y = aso$PopEst, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
bb.pw <- piecewise.linear(x = bb$year, y = bb$PopEst, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
bc.pw <- piecewise.linear(x = bc$year, y = bc$PopEst, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
im.pw <- piecewise.linear(x = im$year, y = im$PopEst, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
lo.pw <- piecewise.linear(x = lo$year, y = lo$PopEst, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
lhc.pw <- piecewise.linear(x = lhc$year, y = lhc$PopEst, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
mv.pw <- piecewise.linear(x = mv$year, y = mv$PopEst, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
mu.pw <- piecewise.linear(x = mu$year, y = mu$PopEst, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
my.pw <- piecewise.linear(x = my$year, y = my$PopEst, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
rb.pw <- piecewise.linear(x = rb$year, y = rb$PopEst, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
uhc.pw <- piecewise.linear(x = uhc$year, y = uhc$PopEst, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
we.pw <- piecewise.linear(x = we$year, y = we$PopEst, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)

#-- look at slope changes in pops that have at least three years of data prior to Movi intro --#
par(mfrow = c(2, 1), oma = c(0, 3, 0, 0), mex = .8)
plot(x = 0, y = 0, xlim = c(-20, 20), ylim = c(0, 12), yaxt = "n", xlab = "", ylab = "")
segments(y0 = 12, y1 = 12, x0 = aso.pw$intervals[1, 3], x1 = aso.pw$intervals[2, 3])
segments(y0 = 11, y1 = 11, x0 = bb.pw$intervals[1, 3], x1 = bb.pw$intervals[2, 3])
segments(y0 = 10, y1 = 10, x0 = bc.pw$intervals[1, 3], x1 = bc.pw$intervals[2, 3])
segments(y0 = 9, y1 = 9, x0 = im.pw$intervals[1, 3], x1 = im.pw$intervals[2, 3])
segments(y0 = 8, y1 = 8, x0 = lo.pw$intervals[1, 3], x1 = lo.pw$intervals[2, 3])
segments(y0 = 7, y1 = 7, x0 = lhc.pw$intervals[1, 3], x1 = lhc.pw$intervals[2, 3])
segments(y0 = 6, y1 = 6, x0 = mv.pw$intervals[1, 3], x1 = mv.pw$intervals[2, 3])
segments(y0 = 5, y1 = 5, x0 = mu.pw$intervals[1, 3], x1 = mu.pw$intervals[2, 3])
segments(y0 = 4, y1 = 4, x0 = my.pw$intervals[1, 3], x1 = my.pw$intervals[2, 3])
segments(y0 = 3, y1 = 3, x0 = rb.pw$intervals[1, 3], x1 = rb.pw$intervals[2, 3])
segments(y0 = 2, y1 = 2, x0 = uhc.pw$intervals[1, 3], x1 = uhc.pw$intervals[2, 3])
segments(y0 = 1, y1 = 1, x0 = we.pw$intervals[1, 3], x1 = we.pw$intervals[2, 3])
abline(v = 0, lwd = 2, lty = 2, col = "red")
axis(side = 2, at = 1:12, labels = c( "Wenaha", "UHC", "Redbird", "Myers", "Muir", "Mount V", "LHC", "Lostine", "Imnaha", "Big Canyon", "Black Butte", "Asotin"), las = 2)

plot(x = 0, y = 0, xlim = c(-20, 20), ylim = c(0, 12), yaxt = "n", xlab = "", ylab = "")
segments(y0 = 12, y1 = 12, x0 = aso.pw$intervals[1, 1], x1 = aso.pw$intervals[2, 1])
segments(y0 = 11, y1 = 11, x0 = bb.pw$intervals[1, 1] - 1995, x1 = bb.pw$intervals[2, 1] - 1995)
segments(y0 = 10, y1 = 10, x0 = bc.pw$intervals[1, 1] - 2000, x1 = bc.pw$intervals[2, 1] - 2000)
segments(y0 = 9, y1 = 9, x0 = im.pw$intervals[1, 1] - 2000, x1 = im.pw$intervals[2, 1] - 2000)
segments(y0 = 8, y1 = 8, x0 = lo.pw$intervals[1, 1] - 1986, x1 = lo.pw$intervals[2, 1] - 1986)
segments(y0 = 7, y1 = 7, x0 = lhc.pw$intervals[1, 1] - 1995, x1 = lhc.pw$intervals[2, 1] - 1995)
segments(y0 = 6, y1 = 6, x0 = mv.pw$intervals[1, 1] - 1995, x1 = mv.pw$intervals[2, 1] - 1995)
segments(y0 = 5, y1 = 5, x0 = mu.pw$intervals[1, 1] - 2000, x1 = mu.pw$intervals[2, 1] - 2000)
segments(y0 = 4, y1 = 4, x0 = my.pw$intervals[1, 1] - 2002, x1 = my.pw$intervals[2, 1] - 2002)
segments(y0 = 3, y1 = 3, x0 = rb.pw$intervals[1, 1] - 1995, x1 = rb.pw$intervals[2, 1] - 1995)
segments(y0 = 2, y1 = 2, x0 = uhc.pw$intervals[1, 1] - 1997, x1 = uhc.pw$intervals[2, 1] - 1997)
segments(y0 = 1, y1 = 1, x0 = we.pw$intervals[1, 1] - 1996, x1 = we.pw$intervals[2, 1] - 1996)
abline(v = 0, lwd = 2, lty = 2, col = "red")
axis(side = 2, at = 1:12, labels = c( "Wenaha", "UHC", "Redbird", "Myers", "Muir", "Mount V", "LHC", "Lostine", "Imnaha", "Big Canyon", "Black Butte", "Asotin"), las = 2)


aso.pw.r <- piecewise.linear(x = aso$year, y = aso$r_est, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
bb.pw.r <- piecewise.linear(x = bb$year, y = bb$r_est, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
bc.pw.r <- piecewise.linear(x = bc$year, y = bc$r_est, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
im.pw.r <- piecewise.linear(x = im$year, y = im$r_est, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
lo.pw.r <- piecewise.linear(x = lo$year, y = lo$r_est, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
lhc.pw.r <- piecewise.linear(x = lhc$year, y = lhc$r_est, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
mv.pw.r <- piecewise.linear(x = mv$year, y = mv$r_est, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
mu.pw.r <- piecewise.linear(x = mu$year, y = mu$r_est, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
my.pw.r <- piecewise.linear(x = my$year, y = my$r_est, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
rb.pw.r <- piecewise.linear(x = rb$year, y = rb$r_est, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
we.pw.r <- piecewise.linear(x = we$year, y = we$r_est, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)

# Plot estimated changepoints 
plot(x = 0, y = 0, xlim = c(1990, 2012), ylim = c(0, 12), yaxt = "n", xlab = "", ylab = "")
segments(y0 = 12, y1 = 12, x0 = aso.pw$intervals[1, 1], x1 = aso.pw$intervals[2, 1])
segments(y0 = 11, y1 = 11, x0 = bb.pw$intervals[1, 1], x1 = bb.pw$intervals[2, 1])
segments(y0 = 10, y1 = 10, x0 = bc.pw$intervals[1, 1], x1 = bc.pw$intervals[2, 1])
segments(y0 = 9, y1 = 9, x0 = im.pw$intervals[1, 1], x1 = im.pw$intervals[2, 1])
segments(y0 = 8, y1 = 8, x0 = lo.pw$intervals[1, 1], x1 = lo.pw$intervals[2, 1])
segments(y0 = 7, y1 = 7, x0 = rb.pw$intervals[1, 1], x1 = rb.pw$intervals[2, 1])
segments(y0 = 6, y1 = 6, x0 = we.pw$intervals[1, 1], x1 = we.pw$intervals[2, 1])
axis(side = 2, at = 1:12, labels = c("", "", "", "", "", "Wenaha", "Redbird", "Lostine", "Imnaha", "Big Canyon", "Black Butte", "Asotin"), las = 2)

# bootstrap stepwise growth intervals and refit changepoint model
aso$Diffs <- rep(NA, dim(aso)[1])

PopTrajectoryBoot <- function(popdata, nboot)
{
  boot.trajectory <- matrix(NA, nrow = nboot, ncol = dim(popdata)[1])
  changepoint.lb <- changepoint.ub <- rep(NA, nboot)
  slope1.est <- slope2.est <- rep(NA, nboot)
  slope1.lb <- slope1.ub <- rep(NA, nboot)
  slope2.lb <- slope2.ub <- rep(NA, nboot)
  popdata$Diffs <- rep(NA, dim(popdata)[1])
  for(i in 2:dim(popdata)[1]){
    popdata$Diffs[i] <- popdata$PopEst[i] - popdata$PopEst[i - 1]
  }
  for(j in 1:nboot){
  boot.diffs <- popdata$Diffs[sample(2:dim(popdata)[1], rep = T)]
  boot.trajectory[j, 1] <- popdata$PopEst[1]
    for(i in 2:dim(popdata)[1]){
      boot.trajectory[j, i] <- max(boot.trajectory[j, i - 1] + boot.diffs[i], 1)
    }
  pw.fit <- piecewise.linear(x = seq(min(popdata$year), max(popdata$year - 1)), y = boot.trajectory[j, -dim(popdata)[1]], middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
  changepoint.lb[j] <- pw.fit$intervals[1, 1]
  changepoint.ub[j] <- pw.fit$intervals[2, 1]
  slope1.lb[j] <- pw.fit$intervals[1, 2]
  slope1.ub[j] <- pw.fit$intervals[2, 2]
  slope1.est[j] <- pw.fit$model$coef[2]
  slope2.est[j] <- pw.fit$model$coef[3]
  slope2.lb[j] <- pw.fit$intervals[1, 4]
  slope2.ub[j] <- pw.fit$intervals[2, 4]
  print(j)
  }
  return(list(changepoint.lb, changepoint.ub, slope1.est, slope2.est, slope1.lb, slope1.ub, slope2.lb, slope2.ub))
}

aso.boot.test <- PopTrajectoryBoot(popdata = aso, nboot = 20)
bb.boot.test <- PopTrajectoryBoot(popdata = bb, nboot = 20)
bc.boot.test <- PopTrajectoryBoot(popdata = bc, nboot = 20)
im.boot.test <- PopTrajectoryBoot(popdata = im, nboot = 20)
lo.boot.test <- PopTrajectoryBoot(popdata = lo, nboot = 20)
mv.boot.test <- PopTrajectoryBoot(popdata = mv, nboot = 20)
mu.boot.test <- PopTrajectoryBoot(popdata = mu, nboot = 20)
my.boot.test <- PopTrajectoryBoot(popdata = my, nboot = 20)
rb.boot.test <- PopTrajectoryBoot(popdata = rb, nboot = 20)
we.boot.test <- PopTrajectoryBoot(popdata = we, nboot = 20)

par(mfrow = c(3, 4))
hist(aso.boot.test[[2]] - aso.boot.test[[1]], xlim = c(0, 22), col = "grey80", main = "Asotin")
abline(v = aso.pw$intervals[2, 1] - aso.pw$intervals[1, 1], col = "red", lwd = 2)
hist(bb.boot.test[[2]] - bb.boot.test[[1]], xlim = c(0, 22), col = "grey80", main = "Black Butte")
abline(v = bb.pw$intervals[2, 1] - bb.pw$intervals[1, 1], col = "red", lwd = 2)
hist(bc.boot.test[[2]] - bc.boot.test[[1]], xlim = c(0, 22), col = "grey80", main = "Big Canyon")
abline(v = bc.pw$intervals[2, 1] - bc.pw$intervals[1, 1], col = "red", lwd = 2)
hist(im.boot.test[[2]] - im.boot.test[[1]], xlim = c(0, 22), col = "grey80", main = "Imnaha")
abline(v = im.pw$intervals[2, 1] - im.pw$intervals[1, 1], col = "red", lwd = 2)
hist(lo.boot.test[[2]] - lo.boot.test[[1]], xlim = c(0, 22), col = "grey80", main = "Lostine")
abline(v = lo.pw$intervals[2, 1] - mv.pw$intervals[1, 1], col = "red", lwd = 2)
hist(mv.boot.test[[2]] - mv.boot.test[[1]], xlim = c(0, 22), col = "grey80", main = "Mountain View")
abline(v = mv.pw$intervals[2, 1] - mv.pw$intervals[1, 1], col = "red", lwd = 2)
hist(mu.boot.test[[2]] - mu.boot.test[[1]], xlim = c(0, 22), col = "grey80", main = "Muir Creek")
abline(v = mu.pw$intervals[2, 1] - mu.pw$intervals[1, 1], col = "red", lwd = 2)
hist(my.boot.test[[2]] - my.boot.test[[1]], xlim = c(0, 22), col = "grey80", main = "Myers Creek")
abline(v = my.pw$intervals[2, 1] - my.pw$intervals[1, 1], col = "red", lwd = 2)
hist(rb.boot.test[[2]] - rb.boot.test[[1]], xlim = c(0, 22), col = "grey80", main = "Redbird")
abline(v = rb.pw$intervals[2, 1] - rb.pw$intervals[1, 1], col = "red", lwd = 2)
hist(we.boot.test[[2]] - we.boot.test[[1]], xlim = c(0, 22), col = "grey80", main = "Wenaha")
abline(v = we.pw$intervals[2, 1] - we.pw$intervals[1, 1], col = "red", lwd = 2)

par(mfrow = c(3, 4))
hist(aso.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Asotin")
abline(v = aso.pw$model$coef[3], col = "red", lwd = 2)
hist(bb.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Black Butte")
abline(v = bb.pw$model$coef[3], col = "red", lwd = 2)
hist(bc.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Big Canyon")
abline(v = bc.pw$model$coef[3], col = "red", lwd = 2)
hist(im.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Imnaha")
abline(v = im.pw$model$coef[3], col = "red", lwd = 2)
hist(lo.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Lostine")
abline(v = lo.pw$model$coef[3], col = "red", lwd = 2)
hist(mv.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Mountain View")
abline(v = mv.pw$model$coef[3], col = "red", lwd = 2)
hist(mu.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Muir Creek")
abline(v = mu.pw$model$coef[3], col = "red", lwd = 2)
hist(my.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Myers Creek")
abline(v = my.pw$model$coef[3], col = "red", lwd = 2)
hist(rb.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Redbird")
abline(v = rb.pw$model$coef[3], col = "red", lwd = 2)
hist(we.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Wenaha")
abline(v = we.pw$model$coef[3], col = "red", lwd = 2)

par(mfrow = c(3, 4))
hist(aso.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Asotin")
abline(v = aso.pw$model$coef[3], col = "red", lwd = 2)
hist(bb.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Black Butte")
abline(v = bb.pw$model$coef[3], col = "red", lwd = 2)
hist(bc.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Big Canyon")
abline(v = bc.pw$model$coef[3], col = "red", lwd = 2)
hist(im.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Imnaha")
abline(v = im.pw$model$coef[3], col = "red", lwd = 2)
hist(lo.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Lostine")
abline(v = lo.pw$model$coef[3], col = "red", lwd = 2)
hist(mv.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Mountain View")
abline(v = mv.pw$model$coef[3], col = "red", lwd = 2)
hist(mu.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Muir Creek")
abline(v = mu.pw$model$coef[3], col = "red", lwd = 2)
hist(my.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Myers Creek")
abline(v = my.pw$model$coef[3], col = "red", lwd = 2)
hist(rb.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Redbird")
abline(v = rb.pw$model$coef[3], col = "red", lwd = 2)
hist(we.boot.test[[4]], xlim = c(-50, 50), col = "grey80", main = "Wenaha")
abline(v = we.pw$model$coef[3], col = "red", lwd = 2)

par(mfrow = c(3, 4))
hist(aso.boot.test[[4]] - aso.boot.test[[3]], xlim = c(-50, 50), col = "grey80", main = "Asotin")
abline(v = (aso.pw$model$coef[3] - aso.pw$model$coef[2]), col = "red", lwd = 2)
hist(bb.boot.test[[4]] - bb.boot.test[[3]], breaks = 40, xlim = c(-50, 50), col = "grey80", main = "Black Butte")
abline(v = (bb.pw$model$coef[3] - bb.pw$model$coef[2]), col = "red", lwd = 2)
hist(bc.boot.test[[4]] - bc.boot.test[[3]], breaks = 40, xlim = c(-50, 50), col = "grey80", main = "Big Canyon")
abline(v = (bc.pw$model$coef[3] - bc.pw$model$coef[2]), col = "red", lwd = 2)
hist(im.boot.test[[4]] - im.boot.test[[3]], breaks = 40, xlim = c(-50, 50), col = "grey80", main = "Imnaha")
abline(v = (im.pw$model$coef[3] - im.pw$model$coef[2]), col = "red", lwd = 2)
hist(lo.boot.test[[4]] - lo.boot.test[[3]], breaks = 40, xlim = c(-50, 50), col = "grey80", main = "Lostine")
abline(v = (lo.pw$model$coef[3] - lo.pw$model$coef[2]), col = "red", lwd = 2)
hist(mv.boot.test[[4]] - mv.boot.test[[3]], breaks = 40, xlim = c(-50, 50), col = "grey80", main = "Mountain View")
abline(v = (mv.pw$model$coef[3] - mv.pw$model$coef[2]), col = "red", lwd = 2)
hist(mu.boot.test[[4]] - mu.boot.test[[3]], breaks = 40, xlim = c(-50, 50), col = "grey80", main = "Muir Creek")
abline(v = (mu.pw$model$coef[3] - mu.pw$model$coef[2]), col = "red", lwd = 2)
hist(my.boot.test[[4]] - my.boot.test[[3]], breaks = 40, xlim = c(-50, 50), col = "grey80", main = "Myers Creek")
abline(v = (my.pw$model$coef[3] - my.pw$model$coef[2]), col = "red", lwd = 2)
hist(rb.boot.test[[4]] - rb.boot.test[[3]], breaks = 40, xlim = c(-50, 50), col = "grey80", main = "Redbird")
abline(v = (rb.pw$model$coef[3] - rb.pw$model$coef[2]), col = "red", lwd = 2)
hist(we.boot.test[[4]] - we.boot.test[[3]], breaks = 40, xlim = c(-50, 50), col = "grey80", main = "Wenaha")
abline(v = (we.pw$model$coef[3] - we.pw$model$coef[2]), col = "red", lwd = 2)

#-- classify years as pre, intro, post and compare via mixed-effects model with RE on year --#
require(lme4)
pop.level.data <- vector("list", length(levels(factor(compd.data$Pop))))
for(i in 2:length(levels(factor(compd.data$Pop)))){
  k <- subset(compd.data, Pop == levels(factor(compd.data$Pop))[i])
  include <- rep(NA, dim(k)[1])
  for(j in 2:length(include)){
    include[j] <- ifelse(k$NoSheepRel[j - 1] == 0 & k$NoSheepRem[j - 1] == 0, 1, 0)
  }
  pop.level.data[[i]] <- k[include == 1, ]
}
compd.data.notranslocation <- do.call("rbind", pop.level.data)
compd.data.lmer <- subset(compd.data.notranslocation, select = c("r_est", "PreIntroPost", "Pop"))

compd.data.lmer <- compd.data.lmer[complete.cases(compd.data.lmer), ]
compd.data.lmer$PrePost <- ifelse(compd.data.lmer$PreIntroPost == "Pre", "Pre", "Post")
r_est_model <- lmer(compd.data.lmer$r_est ~ factor(compd.data.lmer$PreIntroPost) + (1 | compd.data.lmer$Pop))
r_est_model2 <- lmer(compd.data.lmer$r_est ~ factor(compd.data.lmer$PrePost) + (1 | compd.data.lmer$Pop))
anova(r_est_model, r_est_model2)
compd.data.nointro <- subset(compd.data.lmer, PreIntroPost != "Intro")
r_est_model3 <- lmer(compd.data.nointro$r_est ~ factor(compd.data.nointro$PreIntroPost) + (1 | compd.data.nointro$Pop))
summary(r_est_model3)
compd.data.nointro$PreIntroPost <- factor(compd.data.nointro$PreIntroPost)
compd.data.nointro$PreIntroPost <- factor(compd.data.nointro$PreIntroPost, levels(compd.data.nointro$PreIntroPost)[c(2, 1)])

par(mfrow = c(1, 1))
boxplot(compd.data.nointro$r_est ~ compd.data.nointro$PreIntroPost, xaxt = "n", col = "grey80", ylab = "estimated r", names = c("Post all-age \n die-off", "Pre all-age \n die-off"))
abline(h = 0, lty = 2, lwd = 2, col = "grey40")
axis(side = 1, at = c(1, 2), labels = c("Before\n die-off", "After\n die-off"))

compd.data.nointro[which.max(compd.data.nointro$r_est), ]

# ewe count timeseries
plot(xlim = c(1970, 2012), ylim = c(0, 250), x = 0, y = 0, cex = 0)


# sims where we treat disease as one-year phenomenon (build time series with as many pres and as many posts as empirically observed, but not blocked)
pop.level.data <- new.ts <- vector("list", length(levels(factor(compd.data$Pop))))
for(i in 2:length(levels(factor(compd.data$Pop)))){
  k <- subset(compd.data, Pop == levels(factor(compd.data$Pop))[i])
  include <- rep(NA, dim(k)[1])
  for(j in 2:length(include)){
    include[j] <- ifelse(k$NoSheepRel[j - 1] == 0 & k$NoSheepRem[j - 1] == 0, 1, 0)
  }
  pop.level.data[[i]] <- k[include == 1, ]
  new.r.est <- sample(na.omit(pop.level.data[[i]]$r_est), dim(pop.level.data[[i]])[1], rep = T)
  new.ts[[i]] <- rep(NA, length(new.r.est))
  new.ts[[i]][1] <- 40
  for(j in 2:length(new.ts[[i]])){
    new.ts[[i]][j] <- new.ts[[i]][j - 1] + new.ts[[i]][j - 1] * new.r.est[j]
  }
}
compd.data.notranslocation <- do.call("rbind", pop.level.data)
compd.data.lmer <- subset(compd.data.notranslocation, select = c("r_est", "PreIntroPost", "Pop"))


layout(matrix(c(1, 1, 4, 2, 3, 4), nrow = 2, byrow = T))
par(mar = c(3, 6, 1, 1))
plot(xlim = c(0, 35), ylim = c(0, 450), x = 0, y = 0, cex = 0, ylab = "Randomly ordered r values", xlab = "Years since initial restoration")
for(i in 2:length(levels(factor(compd.data$Pop)))){
  lines(new.ts[[i]] ~ seq(1 : length(new.ts[[i]])), col = "grey30", lty = 1)
}
mtext(side = 2, line = 6, inner = T, "Simulated \npopulation estimates", cex = .8)
plot(xlim = c(0, 26), ylim = c(0, 250), x = 0, y = 0, cex = 0, ylab = "Pre-die-off", xlab = "Years since initial restoration")
for(i in 2:length(levels(factor(compd.data$Pop)))){
  k <- subset(compd.data, Pop == levels(factor(compd.data$Pop))[i] & PreIntroPost == "Pre")
  k$std.year <- k$year - min(na.omit(k$year))
  lines(k$PopEst ~ k$std.year, col = "grey30", lty = 1)
}
mtext(side = 2, line = 6, inner = T, "Observed \npopulation estimates", cex = .8)

plot(xlim = c(0, 26), ylim = c(0, 250), x = 0, y = 0, cex = 0, ylab = "Post-die-off", xlab = "Years since first all-age die-off")
for(i in 2:length(levels(factor(compd.data$Pop)))){
  l <- subset(compd.data, Pop == levels(factor(compd.data$Pop))[i] & PreIntroPost == "Post")
  l$std.year <- l$year - min(na.omit(l$year))
  lines(l$PopEst ~ l$std.year, col = "grey30", lty = 1)
}
boxplot(compd.data.nointro$r_est ~ compd.data.nointro$PreIntroPost, xaxt = "n", col = "grey80", ylab = "estimated r", names = c("Post all-age \n die-off", "Pre all-age \n die-off"))
abline(h = 0, lty = 2, lwd = 2, col = "grey40")
axis(side = 1, at = c(1, 2), labels = c("Before\n die-off", "After\n die-off"))
