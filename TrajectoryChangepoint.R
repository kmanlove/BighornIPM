#---------------------------------------------------------#
#-- checks for changepoints in bighorn pop trajectories --#
#---------------------------------------------------------#
require(SiZer)
require(bfast)

compd.data <- read.csv("./Data/compiled_data_summary_130919.csv", header = T, sep = "")
compd.data <- subset(compd.data, year %in% 1990:2012)
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
lhc <- subset(compd.data, Pop == "LHCOR", select = c("year", "PopEst"))
lhc <- lhc[complete.cases(lhc), ]
mv <- subset(compd.data, Pop == "MountainView", select = c("year", "PopEst"))
mv <- mv[complete.cases(mv), ]
mu <- subset(compd.data, Pop == "MuirCreek", select = c("year", "PopEst"))
mu <- mu[complete.cases(mu), ]
my <- subset(compd.data, Pop == "MyersCreek", select = c("year", "PopEst"))
my <- my[complete.cases(my), ]
rb <- subset(compd.data, Pop == "Redbird", select = c("year", "PopEst"))
rb <- rb[complete.cases(rb), ]
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
we.pw <- piecewise.linear(x = we$year, y = we$PopEst, middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)

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
  slope2.lb[j] <- pw.fit$intervals[1, 4]
  slope2.ub[j] <- pw.fit$intervals[2, 4]
  print(j)
  }
  return(list(changepoint.lb, changepoint.ub, slope1.lb, slope1.ub, slope2.lb, slope2.ub))
}

aso.boot.test <- PopTrajectoryBoot(popdata = aso, nboot = 20)
hist(aso.boot.test[[2]] - aso.boot.test[[1]], xlim = c(0, 22), col = "grey80")
abline(v = aso.pw$intervals[2, 1] - aso.pw$intervals[1, 1], col = "red", lwd = 2)

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
