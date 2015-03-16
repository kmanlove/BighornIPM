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

#------------------------------#
#-- Simulations for figure 0 --#
#------------------------------#
require(pomp)
pompExample(gillespie.sir)
plot(gillespie.sir)

sir.step <- function (x, t, params, delta.t, ...) {
  ## unpack the parameters
  N <- params["N"] # population size
  gamma <- params["gamma"] # recovery rate
  mu <- params["mu"] # birth rate = death rate
  beta <- params["beta"] # contact rate
  foi <- beta * x["I"]/N # the force of infection
  trans <- c(
    ## births are assumed to be Poisson:
    rpois(n=1,lambda=mu * N *  delta.t),
    ## exits from S
    reulermultinom(n=1,size=x["S"],rate=c(foi,mu),dt=delta.t),
    ## exits from I:
    reulermultinom(n=1,size=x["I"],rate=c(gamma,mu),dt=delta.t),
    ## exits from R:
    reulermultinom(n=1,size=x["R"],rate=c(mu),dt=delta.t)
  )
  ## now connect the compartments
  x["S"] <- x["S"]+trans[1]-trans[2]-trans[3]
  x["I"] <- x["I"]+trans[2]-trans[4]-trans[5]
  x["R"] <- x["R"]+trans[4]-trans[6]
  x["cases"] <- x["cases"]+trans[4]
  x
}

sir.step <- '
double rate[6]; // transition rates
double trans[6]; // transition numbers
// compute the transition rates
rate[0] = mu
*
N; // birth into susceptible class
rate[1] = beta
*
I/N; // force of infection
rate[2] = mu; // death from susceptible class
rate[3] = gamma; // recovery
rate[4] = mu; // death from infectious class
rate[5] = mu; // death from recovered class
// compute the transition numbers
trans[0] = rpois(rate[0]
*
dt); // births are Poisson
reulermultinom(2,S,&rate[1],dt,&trans[1]);
reulermultinom(2,I,&rate[3],dt,&trans[3]);
reulermultinom(1,R,&rate[5],dt,&trans[5]);
// balance the equations
S += trans[0]-trans[1]-trans[2];
I += trans[1]-trans[3]-trans[4];
R += trans[3]-trans[5];
incid += trans[3]; // incidence is cumulative recoveries;
'

rmeas <- 'cases = rnbinom_mu(theta,rho*incid);'
dmeas <- 'lik = dnbinom_mu(cases,theta,rho * incid,give_log);'

pomp(
  data = data.frame(
    cases = NA,
    time = seq(0, 10, by = 1 / 52)
  ),
  times = "time",
  t0 = -1 / 52,
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),
  rprocess = euler.sim(
    step.fun = Csnippet(sir.step),
    delta.t=1/52/20
  ),
  statenames=c("S","I","R","incid"),
  paramnames=c(
    "gamma","mu","theta","beta",
    "N","rho",
    "S.0","I.0","R.0"
  ),
  zeronames=c("incid"),
  initializer=function(params, t0, ...) {
    x0 <- c(S=0,I=0,R=0,incid=0)
    fracs <- params[c("S.0","I.0","R.0")]
    x0[1:3] <- round(params['N'] * fracs/sum(fracs))
    x0
  },
  params=c(
    N=500000,beta=400,
    gamma=26,mu=1/50,rho=0.1,theta=100,
    S.0=26/400,I.0=0.002,R.0=1
  )
) -> sir

simulate(sir,seed=1914679908L) -> sir

plot(sir)

si.step <- '
double rate[4]; // transition rates
double trans[4]; // transition numbers
// compute the transition rates
rate[0] = mu
*
N; // birth into susceptible class
rate[1] = beta
*
I/N; // force of infection
rate[2] = mu; // death from susceptible class
rate[3] = mu; // death from infectious class
// compute the transition numbers
trans[0] = rpois(rate[0]
*
dt); // births are Poisson
reulermultinom(2,S,&rate[1],dt,&trans[1]);
reulermultinom(1,I,&rate[3],dt,&trans[3]);
// balance the equations
S += trans[0]-trans[1]-trans[2];
I += trans[1]-trans[3];
incid += trans[3]; // incidence is cumulative recoveries;
'

si <- pomp(
  data = data.frame(
    cases = NA,
    time = seq(0, 10, by = 1 / 52)
  ),
  times = "time",
  t0 = -1 / 52,
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),
  rprocess = euler.sim(
    step.fun = Csnippet(si.step),
    delta.t=1/52/20
  ),
  statenames=c("S","I","incid"),
  paramnames=c(
    "mu","theta","beta",
    "N","rho",
    "S.0","I.0"
  ),
  zeronames=c("incid"),
  initializer=function(params, t0, ...) {
    x0 <- c(S=0,I=0,incid=0)
    fracs <- params[c("S.0","I.0")]
    x0[1:2] <- round(params['N'] * fracs/sum(fracs))
    x0
  },
  params=c(
    N=500000,beta=1.0001,mu=1/50,rho=0.1,theta=100,
    S.0=200/400,I.0=0.002
  )
)

simulate(si, seed=1914679908L) -> si

plot(si)

si.data <- as.data.frame(si)
sir.data <-as.data.frame(sir)

par(mfrow = c(2, 2))
plot(si.data$I ~ si.data$time, type = "l", ylab = "Prevalence", xlab = "time")
acf(si.data$I, lag.max = 520, main = "")
plot(sir.data$I ~ sir.data$time, type = "l", ylab = "Prevalence", xlab = "time")
acf(sir.data$I / (sir.data$S + sir.data$I + sir.data$R), lag.max = 520, main = "")

# generate prevalence function
y <- rep(NA, 1000)
y[1] <- 1
k <- -.1
L <- 1
x <- seq(1:1000)
for(i in 2:length(y)){
#  y[i] <- r * y[i - 1] * (1 - (y[i - 1] / K))
  y[i] <- L / (1 + exp(-1 * k * (i - 50)))
}

x.in <- x / 1000
plot(y ~ x, type = "l", ylab = "Host vital rate", xlab = "Pathogen prevalence")
  
# generalized logistic function
y[1] <- 1.1
A <- .9
K <- 1.1
B <- -.02
M <- .2
v <- Q <- .01

for(i in 2:length(y)){
  #  y[i] <- r * y[i - 1] * (1 - (y[i - 1] / K))
  y[i] <- A + ((K - A) / ((1 + (Q) * exp(-1 * B *(i - M * 1000))) ^ (1 / v)))
}

x.in <- x / 1000
plot(y ~ x, type = "l", ylab = "Host vital rate", xlab = "Pathogen prevalence")

prev.si <- round(si.data$I / 500000, 3) * 1000 + 1
prev.sir <- round(sir.data$I / 500000, 3) * 1000 + 1

rate.through.time.si <- y[prev.si]
rate.through.time.sir <- y[prev.sir]

plot(rate.through.time.si ~ seq(1:521), type = "l")
plot(rate.through.time.sir ~ seq(1:521), type = "l")


# figure 0
layout(matrix(c(1, 2, 3, 4, 1, 5, 6, 7), byrow = T, nrow = 2))
par(oma = c(0, 1, 0, 0), mex = .8)
plot(y ~ x, type = "l", ylab = "Host vital rate", xlab = "Pathogen prevalence", xaxt = "n")
axis(side = 1, at = c(0, 250, 500, 750, 1000), labels = c(0, .25, .5, .75, 1.0))
mtext("(a)", side = 3, adj = 0.15, line = -1.5, cex = .8)
plot(sir.data$I ~ sir.data$time, type = "l", ylim = c(0, 1300), ylab = "Prevalence", xlab = "Time (years)")
mtext("(b)", side = 3, adj = 0.15, line = -1.5, cex = .8)
acf(sir.data$I / (sir.data$S + sir.data$I + sir.data$R), lag.max = 520, main = "", xaxt = "n")
axis(side = 1, at = c(0, 200, 400, 600, 800, 1000), labels = c(0, 2, 4, 6, 8, 10))
mtext("(c)", side = 3, adj = 0.15, line = -1.5, cex = .8)
plot(rate.through.time.sir ~ seq(1:521), type = "l", xaxt = "n", ylim = c(1.096, 1.101), xlab = "Time (years)", ylab = "Vital rate through time")
axis(side = 1, at = c(0, 200, 400, 600, 800, 1000), labels = c(0, 2, 4, 6, 8, 10))
mtext("(d)", side = 3, adj = 0.15, line = -1.5, cex = .8)
plot(si.data$I ~ si.data$time, type = "l", ylab = "Prevalence", xlab = "Time (years)")
mtext("(e)", side = 3, adj = 0.15, line = -1.5, cex = .8)
acf(si.data$I, lag.max = 520, main = "", xaxt = "n")
axis(side = 1, at = c(0, 200, 400, 600, 800, 1000), labels = c(0, 2, 4, 6, 8, 10))
mtext("(f)", side = 3, adj = 0.15, line = -1.5, cex = .8)
plot(rate.through.time.si ~ seq(1:521), type = "l", xaxt = "n", xlab = "Time (years)", ylim = c(0.9, 1.13), ylab = "Vital rate through time")
axis(side = 1, at = c(0, 200, 400, 600, 800, 1000), labels = c(0, 2, 4, 6, 8, 10))
mtext("(g)", side = 3, adj = 0.15, line = -1.5, cex = .8)
