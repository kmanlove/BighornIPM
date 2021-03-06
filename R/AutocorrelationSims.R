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

sir <- pomp(
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
    delta.t = 1 / 52 / 20
  ),
  statenames=c("S","I","R","incid"),
  paramnames=c(
    "gamma",
    "mu",
    "theta",
    "beta",
    "N",
    "rho",
    "S.0",
    "I.0",
    "R.0"
  ),
  zeronames = c("incid"),
  initializer= function(params, t0, ...) {
    x0 <- c(S = 0, I = 0, R = 0, incid = 0)
    fracs <- params[c("S.0", "I.0", "R.0")]
    x0[1:3] <- round(params['N'] * fracs / sum(fracs))
    x0
  },
  params=c(
    N = 500000, 
#    N = 5000, 
    beta = 400,
#    beta = 350,
#    gamma = 36,
    gamma = 26,
    mu = 1/50,
#    mu = 1/5,
    rho = 0.1, 
    theta = 100,
#    S.0 = 26/400,
 S.0 = 30/400,
#S.0 = 1/4,
    I.0 = 0.002,
    R.0 = 1
  )
)

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
    N=5000,beta=1.0001,mu=1/50,rho=0.1,theta=100,
    S.0=200/400,I.0=0.002
#     #    N = 500000, 
#     N = 5000, 
#     #    beta = 400,
#     beta = 1.00010,
# #    gamma = 26,
#     #    gamma = 26,
#     #    mu = 1/50,
#     mu = 1/5,
#     rho = 0.1, 
#     theta = 100,
#     #    S.0 = 26/400,
#     S.0 = 26/40,
#     I.0 = 0.002,
#     R.0 = 1
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
prev.si <- round(si.data$I / (si.data$I + si.data$S), 3) * 1000 + 1
prev.sir <- round(sir.data$I / (sir.data$S + sir.data$I + sir.data$R), 3) * 1000 + 1
prev.si.y <- A + ((K - A) / ((1 + (Q) * exp(-1 * B *(prev.si - M * 1000))) ^ (1 / v)))
prev.sir.y <- A + ((K - A) / ((1 + (Q) * exp(-1 * B *(prev.sir - M * 1000))) ^ (1 / v)))
prev.si.y.t.plus.1 <- prev.si.y[-1]
prev.sir.y.t.plus.1 <- prev.sir.y[-1]
  
rate.through.time.si <- y[prev.si]
rate.through.time.sir <- y[prev.sir]

plot(rate.through.time.si ~ seq(1:521), type = "l")
plot(rate.through.time.sir ~ seq(1:521), type = "l")

rate.through.time.si.t.plus.1 <- rate.through.time.si[-1]
rate.through.time.sir.t.plus.1 <- rate.through.time.sir[-1]


# figure 0
layout(matrix(c(1, 2, 3, 4, 5, 1, 6, 7, 8, 9), byrow = T, nrow = 2))
par(oma = c(0, 1, 0, 0), mex = .8)
plot(y ~ x, type = "l", ylab = "Host vital rate", xlab = "Pathogen prevalence", xaxt = "n")
axis(side = 1, at = c(0, 250, 500, 750, 1000), labels = c(0, .25, .5, .75, 1.0))
mtext("(a)", side = 3, adj = 0.15, line = -1.5, cex = .8)
plot(sir.data$I ~ sir.data$time, type = "l", ylim = c(0, 1300), ylab = "Prevalence", xlab = "Time (years)")
mtext("(b)", side = 3, adj = 0.15, line = -1.5, cex = .8)
acf(sir.data$I / (sir.data$S + sir.data$I + sir.data$R), lag.max = 520, main = "", xaxt = "n")
axis(side = 1, at = c(0, 200, 400, 600, 800, 1000), labels = c(0, 2, 4, 6, 8, 10))
mtext("(c)", side = 3, adj = 0.15, line = -1.5, cex = .8)
plot(prev.sir.y.t.plus.1 ~ prev.sir.y[-length(prev.sir.y)], xlab = "Vital rate (t)", ylab = "Vital rate (t + 1)", type = "l", xlim = c(1.096, 1.0965), ylim = c(1.096, 1.0965))
mtext("(d)", side = 3, adj = 0.15, line = -1.5, cex = .8)
plot(prev.sir.y ~ seq(1:521), type = "l", xaxt = "n", ylim = c(1.096, 1.0965), xlab = "Time (years)", ylab = "Vital rate")
axis(side = 1, at = c(0, 200, 400, 600, 800, 1000), labels = c(0, 2, 4, 6, 8, 10))
mtext("(e)", side = 3, adj = 0.15, line = -1.5, cex = .8)
plot(si.data$I ~ si.data$time, type = "l", ylab = "Prevalence", xlab = "Time (years)")
mtext("(f)", side = 3, adj = 0.15, line = -1.5, cex = .8)
acf(si.data$I, lag.max = 520, main = "", xaxt = "n")
axis(side = 1, at = c(0, 200, 400, 600, 800, 1000), labels = c(0, 2, 4, 6, 8, 10))
mtext("(g)", side = 3, adj = 0.15, line = -1.5, cex = .8)
plot(prev.si.y.t.plus.1 ~ prev.si.y[-length(prev.si.y)], type = "l", xlab = "Vital rate (t)", ylab = " Vital rate (t + 1)")
mtext("(h)", side = 3, adj = 0.15, line = -1.5, cex = .8)
plot(prev.si.y ~ seq(1:521), type = "l", xaxt = "n", xlab = "Time (years)", ylim = c(0.9, 1.13), ylab = "Vital rate")
axis(side = 1, at = c(0, 200, 400, 600, 800, 1000), labels = c(0, 2, 4, 6, 8, 10))
mtext("(i)", side = 3, adj = 0.15, line = -1.5, cex = .8)


#------------------------------------#
#-- SIR sims for Figure 2C ----------#
#------------------------------------#
require(pomp)

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
# test

N.levels <- exp(seq(log(10), log(100000), length.out = 15))
Ro.levels <- exp(seq(log(1), log(18), length.out = 10))
gamma.levels <- c(1/10, 1/100, 1/1000)
nu <- .7
mu <- 1 / (20 * 365)
collar.prop <- .25
N <- 1000
delta <- 1 / (.15 * collar.prop * N)
p <- .1 # probability of dying from infection

beta.levels  <- (Ro.levels * (mu + gamma.levels) * mu) / ((1 - p) * nu)
param.mat <- as.data.frame(expand.grid(list(N.levels, beta.levels, gamma.levels)))
names(param.mat) <- c("N", "beta", "gamma")

rmeas <- 'cases = rnbinom_mu(theta, rho*incid);'
dmeas <- 'lik = dnbinom_mu(cases, theta, rho * incid, give_log);'


sir.sim <- function(param.mat, reps){
  sir.test <- sir.test.sim <- sir.test.data <- fade.out <- vector("list", dim(param.mat)[1])
  for(i in 1:dim(param.mat)[1]){
    sir.test[[i]] <- pomp(
      data = data.frame(
        cases = NA,
        time = seq(0, 100, by = 1 / 52)
      ),
      times = "time",
      t0 = -1 / 52,
      dmeasure = Csnippet(dmeas),
      rmeasure = Csnippet(rmeas),
      rprocess = euler.sim(
        step.fun = Csnippet(sir.step),
        delta.t = 1 / 52 / 20
      ),
      statenames=c("S","I","R","incid"),
      paramnames=c(
        "gamma",
        "mu",
        "theta",
        "beta",
        "N",
        "rho",
        "S.0",
        "I.0",
        "R.0"
      ),
      zeronames = c("incid"),
      initializer= function(params, t0, ...) {
        x0 <- c(S = 0, I = 0, R = 0, incid = 0)
        fracs <- params[c("S.0", "I.0", "R.0")]
        x0[1:3] <- round(params['N'] * fracs / sum(fracs))
        x0
      },
      params=c(
        N = param.mat$N[i], 
        beta = param.mat$beta[i],
        gamma = param.mat$gamma[i],
        mu = 1/20,
        rho = 1, 
        theta = N,
        S.0 = (param.mat$N[i] - 1) / param.mat$N[i],
        I.0 = (1 / param.mat$N[i]),
        R.0 = 0
      )
    )
    sir.test.sim[[i]] <- vector("list", reps)
    sir.test.data[[i]] <- vector("list", reps)
    fade.out[[i]] <- rep(NA, reps)
    for(j in 1:reps){
#      sir.test.sim[[i]][[j]] <- simulate(sir.test[[i]], seed = paste(i, j, "1914679908", sep = ""))
      seed.in <- as.numeric(paste(i, j, sep = ""))
      sir.test.sim[[i]][[j]] <- simulate(sir.test[[i]], seed = seed.in)
      sir.test.data[[i]][[j]] <- as.data.frame(sir.test.sim[[i]][[j]])
      fade.out[[i]][j] <- min(which(sir.test.data[[i]][[j]]$I == 0))
    } # reps
    print(i)
  } # param mat
  
  return(list(fade.out = fade.out, sir.test.sim = sir.test.sim, sir.test.data = sir.test.data))
  
}

reps.in <- 100
sir.sim.batch <- sir.sim(param.mat[1:50, ], reps = reps.in)

param.mat.with.sims <- matrix(NA, nrow = 50, ncol = dim(param.mat)[2] + reps.in)
for(i in 1:dim(param.mat.with.sims)[1]){
  param.mat.with.sims[i, ] <- unlist(c(param.mat[i, ], as.vector(unlist(sir.sim.batch$fade.out[[i]]))))
}

param.mat.with.sims

write.csv(param.mat.with.sims, "./Data/SimmedSIRs/ParamMatSims3_1_50.csv")

plot(sir.sim.batch$sir.test.sim[[1]][[1]])
plot(sir.sim.batch$sir.test.sim[[2]][[1]])
plot(sir.sim.batch$sir.test.sim[[2]][[2]])
plot(sir.sim.batch$sir.test.sim[[3]][[1]])

#----------------------#
#-- plot simmed data --#
#----------------------#
# read in simmed data
sims_1_50 <- read.csv("./Data/SimmedSIRs/ParamMatSims2_1_50.csv")
sims_51_100 <- read.csv("./Data/SimmedSIRs/ParamMatSims2_51_100.csv")
sims_101_150 <- read.csv("./Data/SimmedSIRs/ParamMatSims2_101_150.csv")
sims_151_200 <- read.csv("./Data/SimmedSIRs/ParamMatSims2_151_200.csv")
sims_201_250 <- read.csv("./Data/SimmedSIRs/ParamMatSims2_201_250.csv")
sims_251_300 <- read.csv("./Data/SimmedSIRs/ParamMatSims2_251_300.csv")
sims_301_350 <- read.csv("./Data/SimmedSIRs/ParamMatSims2_301_350.csv")
sims_351_400 <- read.csv("./Data/SimmedSIRs/ParamMatSims2_351_400.csv")
sims_401_450 <- read.csv("./Data/SimmedSIRs/ParamMatSims2_401_450.csv")

full.simmed.data <- rbind(sims_1_50[1:50, ], 
                          sims_51_100[1:50, ], 
                          sims_101_150[1:50, ], 
                          sims_151_200[1:50, ], 
                          sims_201_250[1:50, ],
                          sims_251_300[1:50, ], 
                          sims_301_350[1:50, ],
                          sims_351_400[1:50, ], 
                          sims_401_450[1:50, ]
                          )
full.simmed.data$Q025 <- full.simmed.data$Q25 <- full.simmed.data$Q50 <- full.simmed.data$Q75 <- full.simmed.data$Q975 <- rep(NA, dim(full.simmed.data)[1])

full.simmed.data[is.infinite(as.matrix(full.simmed.data))] <- 100000

for(i in 1:dim(full.simmed.data)[1]){
  # extract quantiles
  full.simmed.data$Q025[i] <- quantile(full.simmed.data[i, 4:103], 0.025)
  full.simmed.data$Q25[i] <- quantile(full.simmed.data[i, 4:103], 0.25)
  full.simmed.data$Q50[i] <- quantile(full.simmed.data[i, 4:103], 0.5)
  full.simmed.data$Q75[i] <- quantile(full.simmed.data[i, 4:103], 0.75)
  full.simmed.data$Q975[i] <- quantile(full.simmed.data[i, 4:103], 0.975)
}

names(full.simmed.data) <- c("ParamSet", "N", "Ro", "gamma", seq(1:100), "Q025", "Q25", "Q50", "Q75", "Q975")

gamma.1 <- subset(full.simmed.data, full.simmed.data$gamma == 0.1)
gamma.01 <- subset(full.simmed.data, full.simmed.data$gamma == 0.01)
gamma.001 <- subset(full.simmed.data, full.simmed.data$gamma == 0.001)

get.ccs <- function(data, gamma.in){
  # split data on levels of Ro
  CCS <- min.n.pos <- rep(NA, length(levels(factor(data$Ro))))
  Nvals <- levels(factor(data$N))
  for(j in 1:length(levels(factor(data$Ro)))){
    k <- subset(data, Ro == levels(factor(data$Ro))[j])
    min.n.pos[j] <- which.min(abs((unlist(k$Q50)) - 365))
    CCS[j] <- k$N[min.n.pos[j]]
  }
  plot.data <- as.data.frame(cbind(CCS, levels(factor(data$Ro)), rep(gamma.in, length(CCS)), min.n.pos))
  names(plot.data) <- c("CCS", "Ro", "gamma", "min.n.pos")
  return(plot.data)
}

gamma.1.ccs <- get.ccs(data = gamma.1, gamma.in = .1)
gamma.01.ccs <- get.ccs(data = gamma.01, gamma.in = .01)
gamma.001.ccs <- get.ccs(data = gamma.001, gamma.in = .001)

plot(as.numeric(as.character(gamma.1.ccs$CCS)) ~ as.numeric(as.character(gamma.1.ccs$Ro)), ylim = c(10, 10000), log = "xy", pch = 16, col = "black", type = "b")
points(as.numeric(as.character(gamma.001.ccs$CCS)) ~ as.numeric(as.character(gamma.001.ccs$Ro)), pch = 16, col = "red", type = "b")
points(as.numeric(as.character(gamma.01.ccs$CCS)) ~ as.numeric(as.character(gamma.01.ccs$Ro)), pch = 16, col = "blue", type = "b")

#-- persistence time, Ro by gamma (no CCS) --#
gamma.1.ro1 <- subset(gamma.1, Ro == levels(factor(gamma.1$Ro))[1])
gamma.1.ro1.perstime <- as.vector(unlist(gamma.1.ro1[, 5:104]))
gamma.1.ro1.N <- rep(gamma.1.ro1$N, times = 100)

gamma.1.ro2 <- subset(gamma.1, Ro == levels(factor(gamma.1$Ro))[2])
gamma.1.ro2.perstime <- as.vector(unlist(gamma.1.ro2[, 5:104]))
gamma.1.ro2.N <- rep(gamma.1.ro2$N, times = 100)

gamma.1.ro3 <- subset(gamma.1, Ro == levels(factor(gamma.1$Ro))[3])
gamma.1.ro3.perstime <- as.vector(unlist(gamma.1.ro3[, 5:104]))
gamma.1.ro3.N <- rep(gamma.1.ro3$N, times = 100)

gamma.1.ro4 <- subset(gamma.1, Ro == levels(factor(gamma.1$Ro))[4])
gamma.1.ro4.perstime <- as.vector(unlist(gamma.1.ro4[, 5:104]))
gamma.1.ro4.N <- rep(gamma.1.ro4$N, times = 100)

gamma.1.ro5 <- subset(gamma.1, Ro == levels(factor(gamma.1$Ro))[5])
gamma.1.ro5.perstime <- as.vector(unlist(gamma.1.ro5[, 5:104]))
gamma.1.ro5.N <- rep(gamma.1.ro5$N, times = 100)

gamma.1.ro6 <- subset(gamma.1, Ro == levels(factor(gamma.1$Ro))[6])
gamma.1.ro6.perstime <- as.vector(unlist(gamma.1.ro6[, 5:104]))
gamma.1.ro6.N <- rep(gamma.1.ro6$N, times = 100)

gamma.1.ro10 <- subset(gamma.1, Ro == levels(factor(gamma.1$Ro))[10])
gamma.1.ro10.perstime <- as.vector(unlist(gamma.1.ro10[, 5:104]))
gamma.1.ro10.N <- rep(gamma.1.ro10$N, times = 100)


# gamma = 0.01
gamma.01.ro1 <- subset(gamma.01, Ro == levels(factor(gamma.01$Ro))[1])
gamma.01.ro1.perstime <- as.vector(unlist(gamma.01.ro1[, 5:104]))
gamma.01.ro1.N <- rep(gamma.01.ro1$N, times = 100)

gamma.01.ro2 <- subset(gamma.01, Ro == levels(factor(gamma.01$Ro))[2])
gamma.01.ro2.perstime <- as.vector(unlist(gamma.01.ro2[, 5:104]))
gamma.01.ro2.N <- rep(gamma.01.ro2$N, times = 100)

gamma.01.ro3 <- subset(gamma.01, Ro == levels(factor(gamma.01$Ro))[3])
gamma.01.ro3.perstime <- as.vector(unlist(gamma.01.ro3[, 5:104]))
gamma.01.ro3.N <- rep(gamma.01.ro3$N, times = 100)

gamma.01.ro4 <- subset(gamma.01, Ro == levels(factor(gamma.01$Ro))[4])
gamma.01.ro4.perstime <- as.vector(unlist(gamma.01.ro4[, 5:104]))
gamma.01.ro4.N <- rep(gamma.01.ro4$N, times = 100)

gamma.01.ro5 <- subset(gamma.01, Ro == levels(factor(gamma.01$Ro))[5])
gamma.01.ro5.perstime <- as.vector(unlist(gamma.01.ro5[, 5:104]))
gamma.01.ro5.N <- rep(gamma.01.ro5$N, times = 100)

gamma.01.ro6 <- subset(gamma.01, Ro == levels(factor(gamma.01$Ro))[6])
gamma.01.ro6.perstime <- as.vector(unlist(gamma.01.ro6[, 5:104]))
gamma.01.ro6.N <- rep(gamma.01.ro6$N, times = 100)

gamma.01.ro10 <- subset(gamma.01, Ro == levels(factor(gamma.01$Ro))[10])
gamma.01.ro10.perstime <- as.vector(unlist(gamma.01.ro10[, 5:104]))
gamma.01.ro10.N <- rep(gamma.01.ro10$N, times = 100)


# gamma = 0.001
gamma.001.ro1 <- subset(gamma.001, Ro == levels(factor(gamma.001$Ro))[1])
gamma.001.ro1.perstime <- as.vector(unlist(gamma.001.ro1[, 5:104]))
gamma.001.ro1.N <- rep(gamma.001.ro1$N, times = 100)
cbind(gamma.001.ro1.perstime, gamma.001.ro1.N)

gamma.001.ro2 <- subset(gamma.001, Ro == levels(factor(gamma.001$Ro))[2])
gamma.001.ro2.perstime <- as.vector(unlist(gamma.001.ro2[, 5:104]))
gamma.001.ro2.N <- rep(gamma.001.ro2$N, times = 100)

gamma.001.ro3 <- subset(gamma.001, Ro == levels(factor(gamma.001$Ro))[3])
gamma.001.ro3.perstime <- as.vector(unlist(gamma.001.ro3[, 5:104]))
gamma.001.ro3.N <- rep(gamma.001.ro3$N, times = 100)

gamma.001.ro4 <- subset(gamma.001, Ro == levels(factor(gamma.001$Ro))[4])
gamma.001.ro4.perstime <- as.vector(unlist(gamma.001.ro4[, 5:104]))
gamma.001.ro4.N <- rep(gamma.001.ro4$N, times = 100)

gamma.001.ro5 <- subset(gamma.001, Ro == levels(factor(gamma.001$Ro))[5])
gamma.001.ro5.perstime <- as.vector(unlist(gamma.001.ro5[, 5:104]))
gamma.001.ro5.N <- rep(gamma.001.ro5$N, times = 100)

gamma.001.ro6 <- subset(gamma.001, Ro == levels(factor(gamma.001$Ro))[6])
gamma.001.ro6.perstime <- as.vector(unlist(gamma.01.ro6[, 5:104]))
gamma.001.ro6.N <- rep(gamma.001.ro6$N, times = 100)

gamma.001.ro10 <- subset(gamma.001, Ro == levels(factor(gamma.001$Ro))[10])
gamma.001.ro10.perstime <- as.vector(unlist(gamma.01.ro10[, 5:104]))
gamma.001.ro10.N <- rep(gamma.001.ro10$N, times = 100)

# Figure 
par(mfrow = c(1, 3))
plot(gamma.1.ro1.perstime ~ gamma.1.ro1.N, main = "gamma = .1", log = "xy", ylab = "Persistence time", xlab = "Population size", ylim = c(2, 5200), xlim = c(100, 100000))
points(gamma.1.ro2.perstime ~ gamma.1.ro2.N, col = "blue")
points(gamma.1.ro3.perstime ~ gamma.1.ro3.N, col = "purple")
points(gamma.1.ro4.perstime ~ gamma.1.ro4.N, col = "red")
points(gamma.1.ro5.perstime ~ gamma.1.ro5.N, col = "orange")
points(gamma.1.ro6.perstime ~ gamma.1.ro6.N, col = "yellow")
points(gamma.1.ro10.perstime ~ gamma.1.ro10.N, col = "green")
abline(h = 365)
abline(v = 500, col = "grey40")
abline(v = 1000, col = "grey40")
abline(v = 10000, col = "grey40")

plot(gamma.01.ro1.perstime ~ gamma.01.ro1.N, main = "gamma = .01", log = "xy", ylab = "Persistence time", xlab = "Population size", ylim = c(2, 5200), xlim = c(100, 100000))
points(gamma.01.ro2.perstime ~ gamma.01.ro2.N, col = "blue")
points(gamma.01.ro3.perstime ~ gamma.01.ro3.N, col = "purple")
points(gamma.01.ro4.perstime ~ gamma.01.ro4.N, col = "red")
points(gamma.01.ro5.perstime ~ gamma.01.ro5.N, col = "orange")
points(gamma.01.ro6.perstime ~ gamma.01.ro6.N, col = "yellow")
points(gamma.01.ro10.perstime ~ gamma.01.ro10.N, col = "green")
abline(h = 365)
abline(v = 500, col = "grey40")
abline(v = 1000, col = "grey40")
abline(v = 10000, col = "grey40")

plot(gamma.001.ro1.perstime ~ gamma.001.ro1.N, main = "gamma = .001", log = "xy", ylab = "Persistence time", xlab = "Population size", ylim = c(2, 5200), xlim = c(100, 100000))
points(gamma.001.ro2.perstime ~ gamma.001.ro2.N, col = "blue")
points(gamma.001.ro3.perstime ~ gamma.001.ro3.N, col = "purple")
points(gamma.001.ro4.perstime ~ gamma.001.ro4.N, col = "red")
points(gamma.001.ro5.perstime ~ gamma.001.ro5.N, col = "orange")
points(gamma.001.ro6.perstime ~ gamma.001.ro6.N, col = "yellow")
points(gamma.001.ro10.perstime ~ gamma.001.ro10.N, col = "green")
abline(h = 365)
abline(v = 500, col = "grey40")
abline(v = 1000, col = "grey40")
abline(v = 10000, col = "grey40")

