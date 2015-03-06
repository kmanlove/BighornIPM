#-----------------------------------#
#-- Motivating simulation example --#
#-----------------------------------#

# assume stationary series
# use Gaussian noise filter
require(astsa)

x <- filter(rnorm(100), filter=rep(1,3), circular=TRUE)
plot(x ~ seq(1:100), type = "l")

AR <- .8

x1 = arima.sim(list(order = c(1, 0, 0), ar = .8), n = 100)
x2 = arima.sim(list(order = c(1, 0, 0), ar = -.8), n = 100)
par(mfrow=c(2,1))
plot(x1, main=(expression(AR(1)~~~phi==+.9)))  # ~ is a space and == is equal  
plot(x2, main=(expression(AR(1)~~~phi==-.9)))

# let xs be "r"s
typical.var <- .05

ts.sim <- function(reps, ar.param, typical.var, max.pop.size)
{
  N.1.100 <- N.2.100 <- N.3.100 <- rep(NA, reps)
  x1 <- x2 <- x3 <- vector("list", reps)
  for(j in 1:reps){
  x1[[j]] = arima.sim(list(order = c(1, 0, 0), ar = ar.param), n = 100)
  x2[[j]] = arima.sim(list(order = c(1, 0, 0), ar = (ar.param * -1)), n = 100)
  x3[[j]] = rnorm(0, 1, n = 100)

  N1 <- N2 <- N3 <- rep(NA, 100)
  N1[1] <- N2[1] <- N3[1] <- 50
  for(i in 1:99){
    N1[i + 1] <- max(min(max.pop.size, N1[i] + N1[i] * x1[[j]][i] * typical.var), 1)
    N2[i + 1] <- max(min(max.pop.size, N2[i] + N2[i] * x2[[j]][i] * typical.var), 1)
    N3[i + 1] <- max(min(max.pop.size, N3[i] + N3[i] * x3[[j]][i] * typical.var), 1)
  }
  N.1.100[j] <- N1[100]
  N.2.100[j] <- N2[100]
  N.3.100[j] <- N3[100]
  }
  
  return(list(N.1.100, N.2.100, N.3.100))
}

reps <- 1000
ts.sim.test <- ts.sim(reps = reps, ar.param = .9, typical.var = .05, max.pop.size = 500)
ts.sim.out <- cbind(c(ts.sim.test[[1]], ts.sim.test[[2]], ts.sim.test[[3]]), c(rep(1, reps), rep(2, reps), rep(3, reps)))
var(ts.sim.test[[1]])
var(ts.sim.test[[2]])
var(ts.sim.test[[3]])

par(mfrow = c(1, 1), oma = c(2, 2, 2, 2), mar = c(2, 4, 1, 1), las = 1)
boxplot(ts.sim.out[, 1] ~ ts.sim.out[ , 2], log = "y", ylab = "Population size after 100 years", xaxt = "n")
boxplot(ts.sim.out[, 1] ~ ts.sim.out[ , 2], ylab = "Population size after 100 years", xaxt = "n")
axis(side = 1, at = c(1, 2, 3), labels = c("Positive", "Negative", "None"), las = 1, cex.axis = .8)
mtext(side = 1, line = 3, outer = F, "Autocorrelation structure")
sim.lm <- lm(ts.sim.out[, 1] ~ factor(ts.sim.out[, 2]))
anova(sim.lm)
fligner.test.out <- fligner.test(ts.sim.out[, 1] ~ factor(ts.sim.out[, 2]))
