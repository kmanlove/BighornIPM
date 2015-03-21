#-----------------------------------------------------------------#
#-- empirical autocorrelation in disease status at Hells Canyon --#
#-----------------------------------------------------------------#

compd.data <- read.csv("./Data/compiled_data_summary_130919.csv", header = T, sep = ",")
compd.data <- subset(compd.data, !(Pop == "Imnaha" & year <= 1999))
compd.data$PNIndLambs <- with(compd.data, 
                              ifelse(CLASS == c("HEALTHY"), 1, 
                                     ifelse(CLASS %in% c("ADULTS", "ALL_AGE", "ALL_AGE_SUSP", "LAMBS", "LAMBS_SUSP"), 2, NA))
)
compd.data$PNIndEwes <- with(compd.data, 
                             PNIndLambs <- ifelse(CLASS == c("HEALTHY"), 1, 
                                                  ifelse(CLASS %in% c("ALL_AGE", "ALL_AGE_SUSP", "LAMBS", "LAMBS_SUSP", "ADULTS"), 2, NA))
)

compd.data <- compd.data[-476, ]
compd.data$Pop <- factor(compd.data$Pop)
compd.data <- subset(compd.data, year >= 1995 & year <= 2012)

suspected.classes <- c("", "ALL_AGE_SUSP", "LAMBS_SUSP", "LAMBS_SUS{", "LAMB_SUSP")
known.status <- subset(compd.data, !(CLASS_SUSP %in% suspected.classes))

known.status$Binary.class <- ifelse(is.na(known.status$CLASS_SUSP), "NA", 
                                    ifelse(known.status$CLASS_SUSP == "HEALTHY", "Healthy", "Pneumonia"))

par(oma = c(1, 3, 0, 0), mar = c(2, 3, 1, 1), cex.axis = .7)
plot(as.numeric(factor(known.status$Pop)) ~ known.status$year, ylab = "", col = ifelse(known.status$Binary.class == "NA", "white", ifelse(known.status$Binary.class == "Healthy", "black", "red")), pch = 15, yaxt = "n", cex = 1)
axis(side = 2, at = seq(1:length(levels(factor(known.status$Pop)))), labels = levels(factor(known.status$Pop)), las = 2, cex.axis = .7)

#-----------------------------------------------#
#-- get CIs on state-transition probabilities --#
#-----------------------------------------------#
# state-transitions with no autocorr structure would have random draws of vals
# point-est with no structure is 
he.prob <- 54 / (54 + 72) # 0.4286
pn.prob <- 72 / (54 + 72) # 0.5714

# bootstrap pops and recalculate
nboot <- 100
n.pops <- length(levels(known.status$Pop))

empautocorr.boot <- function(nboot){
  he.prob.nostructure <- rep(NA, nboot)
  hehe.nostruct <- hepn.nostruct <- pnhe.nostruct <- pnpn.nostruct <- rep(NA, nboot)
  hehe.struct <- hepn.struct <- pnhe.struct <- pnpn.struct <- rep(NA, nboot)
  resample.pops <- resample.list <- resample.dat <- no.nas <- matches <- structure.match.tab <- structure.match.vec <- vector("list", nboot)
  for(i in 1:nboot){
    resample.pops[[i]] <- sample(1:n.pops, n.pops, replace = T)
    resample.list[[i]] <- vector("list", n.pops)
    for(j in 1:n.pops)
    {
      resample.list[[i]][[j]] <- subset(known.status, Pop == levels(known.status$Pop)[resample.pops[[i]][j]])
    }
    resample.dat[[i]] <- do.call("rbind", resample.list[[i]])
    
    # no structure
    dis.status.nostructure <- table(resample.dat[[i]]$Binary.class)
    he.prob.nostructure[i] <- dis.status.nostructure["Healthy"] / (dis.status.nostructure["Healthy"] + dis.status.nostructure["Pneumonia"])
    hehe.nostruct[i] <- he.prob.nostructure[i] * he.prob.nostructure[i]
    hepn.nostruct[i] <- he.prob.nostructure[i] * (1 - he.prob.nostructure[i])
    pnhe.nostruct[i] <- he.prob.nostructure[i] * (1 - he.prob.nostructure[i])
    pnpn.nostruct[i] <- (1 - he.prob.nostructure[i]) * (1 - he.prob.nostructure[i])
    
    
    # structure
    matches[[i]] <- vector("list", n.pops)
    for(j in 1:length(resample.list[[i]])){
      no.nas[[i]] <- subset(resample.list[[i]][[j]], !Binary.class == "NA")
      matches[[i]][[j]] <- rep(NA, dim(no.nas[[i]])[1])
      if(dim(no.nas[[i]])[1] <= 1){
        matches[[i]][[j]] <- NA
      } else {
        matches[[i]][[j]] <- rep(NA, dim(no.nas[[i]])[1])
      for(k in 2:dim(no.nas[[i]])[1]){
        matches[[i]][[j]][k] <- ifelse(no.nas[[i]]$year[k] != (no.nas[[i]]$year[k-1] + 1), "NA",
                           ifelse(no.nas[[i]]$Binary.class[k] == "Healthy" & no.nas[[i]]$Binary.class[k - 1] == "Healthy", "HeHe", 
                           ifelse(no.nas[[i]]$Binary.class[k] == "Healthy" & no.nas[[i]]$Binary.class[k - 1] == "Pneumonia", "HePn",
                           ifelse(no.nas[[i]]$Binary.class[k] == "Pneumonia" & no.nas[[i]]$Binary.class[k - 1] == "Pneumonia", "PnPn",
                                  "PnHe"
                                  ))))
        
      }
      }
    }
    structure.match.vec[[i]] <- as.vector(unlist(matches[[i]]))
    structure.match.vec[[i]] <- structure.match.vec[[i]][!(structure.match.vec[[i]] == "NA")]
    structure.match.tab[[i]] <- table(structure.match.vec[[i]])
    hehe.struct[i] <- structure.match.tab[[i]]["HeHe"] / sum(structure.match.tab[[i]])
    pnpn.struct[i] <- structure.match.tab[[i]]["PnPn"] / sum(structure.match.tab[[i]])
    hepn.struct[i] <- structure.match.tab[[i]]["HePn"] / sum(structure.match.tab[[i]])
    pnhe.struct[i] <- structure.match.tab[[i]]["PnHe"] / sum(structure.match.tab[[i]])
    
    print(i)
  }
  
  return(list(hehe.nostruct = hehe.nostruct, 
              hepn.nostruct = hepn.nostruct, 
              pnpn.nostruct = pnpn.nostruct, 
              pnhe.nostruct = pnhe.nostruct, 
              hehe.struct = hehe.struct, 
              hepn.struct = hepn.struct, 
              pnhe.struct = pnhe.struct, 
              pnpn.struct = pnpn.struct))
  
}

empautocorr.test <- empautocorr.boot(1000)

hehe.quants.nostr <- quantile(empautocorr.test$hehe.nostruct, c(0.05, 0.95))
hepn.quants.nostr <- quantile(empautocorr.test$hepn.nostruct, c(0.05, 0.95))
pnpn.quants.nostr <- quantile(empautocorr.test$pnpn.nostruct, c(0.05, 0.95))
pnhe.quants.nostr <- quantile(empautocorr.test$pnhe.nostruct, c(0.05, 0.95))
hehe.quants.str <- quantile(empautocorr.test$hehe.struct, c(0.05, 0.95))
hepn.quants.str <- quantile(na.omit(empautocorr.test$hepn.struct), c(0.05, 0.95))
pnpn.quants.str <- quantile(empautocorr.test$pnpn.struct, c(0.05, 0.95))
pnhe.quants.str <- quantile(empautocorr.test$pnhe.struct, c(0.05, 0.95))
same.state.quants.str <- quantile(c(empautocorr.test$pnpn.struct, empautocorr.test$hehe.struct), c(0.05, 0.95))
cross.state.quants.str <- quantile(c(na.omit(empautocorr.test$hepn.struct), empautocorr.test$pnhe.struct), c(0.05, 0.95))
same.state.quants.nostr <- quantile(c(empautocorr.test$pnpn.nostruct, empautocorr.test$hehe.nostruct), c(0.05, 0.95))
cross.state.quants.nostr <- quantile(c(na.omit(empautocorr.test$hepn.nostruct), empautocorr.test$pnhe.nostruct), c(0.05, 0.95))

par(oma = c(2, 5, 0, 0))
plot(x = 0, y = 0, cex = 0, ylim = c(0, 8), yaxt = "n", xlim = c(0, .7))
#segments(x0 = same.state.quants.nostr[1], x1 = same.state.quants.nostr[2], y0 = 12, y1 =12, lty = 2)
#segments(x0 = same.state.quants.str[1], x1 = same.state.quants.str[2], y0 = 11, y1 = 11)
#segments(x0 = cross.state.quants.nostr[1], x1 = cross.state.quants.nostr[2], y0 = 10, y1 = 10, lty = 2)
#segments(x0 = cross.state.quants.str[1], x1 = cross.state.quants.str[2], y0 = 9, y1 = 9)
segments(x0 = pnpn.quants.nostr[1], x1 = pnpn.quants.nostr[2], y0 = 8, y1 = 8, lty = 2)
segments(x0 = hehe.quants.nostr[1], x1 = hehe.quants.nostr[2], y0 = 7, y1 = 7, lty = 2)
segments(x0 = hepn.quants.nostr[1], x1 = hepn.quants.nostr[2], y0 = 6, y1 = 6, lty = 2)
segments(x0 = pnhe.quants.nostr[1], x1 = pnhe.quants.nostr[2], y0 = 5, y1 = 5, lty = 2)
segments(x0 = pnpn.quants.str[1], x1 = pnpn.quants.str[2], y0 = 4, y1 = 4)
segments(x0 = hehe.quants.str[1], x1 = hehe.quants.str[2], y0 = 3, y1 = 3)
segments(x0 = hepn.quants.str[1], x1 = hepn.quants.str[2], y0 = 2, y1 = 2)
segments(x0 = pnhe.quants.str[1], x1 = pnhe.quants.str[2], y0 = 1, y1 = 1)
axis(side = 2, at = seq(1:8), labels = c("PnHe Struct", "HePn Struct", 
                                         "HeHe Struct", "PnPn Struct",
                                         "PnHe Not Struct", "HePn Not Struct",
                                         "HeHe Not Struct", "PnPn Not Struct"), las = 2)

par(mfrow = c(1, 2), oma = c(0, 0, 0, 0), mar = c(4, 5, 1, 1), cex.axis = .7, cex.lab = .8)
plot(as.numeric(factor(known.status$Pop)) ~ known.status$year, ylim = c(-2.5, 16), xlab = "year", ylab = "", col = ifelse(known.status$Binary.class == "NA", "white", ifelse(known.status$Binary.class == "Healthy", "black", "red")), pch = 15, yaxt = "n", cex = 1)
leg.text <- c("Disease", "No Disease")
legend("bottomright", leg.text, cex = .7, pch = c(15, 15), pt.cex = c(1, 1), col = c("red", "black"), bty = "n")
axis(side = 2, at = seq(1:length(levels(factor(known.status$Pop)))), labels = levels(factor(known.status$Pop)), las = 2, cex.axis = .7)
plot(x = 0, y = 0, cex = 0, ylim = c(-1, 8), yaxt = "n", xlim = c(0, .8), ylab = "", xlab = "P(state-transition)")
#segments(x0 = same.state.quants.nostr[1], x1 = same.state.quants.nostr[2], y0 = 12, y1 =12, lty = 2)
#segments(x0 = same.state.quants.str[1], x1 = same.state.quants.str[2], y0 = 11, y1 = 11)
#segments(x0 = cross.state.quants.nostr[1], x1 = cross.state.quants.nostr[2], y0 = 10, y1 = 10, lty = 2)
#segments(x0 = cross.state.quants.str[1], x1 = cross.state.quants.str[2], y0 = 9, y1 = 9)
segments(x0 = pnpn.quants.nostr[1], x1 = pnpn.quants.nostr[2], y0 = 8, y1 = 8, lty = 2)
segments(x0 = hehe.quants.nostr[1], x1 = hehe.quants.nostr[2], y0 = 7, y1 = 7, lty = 2)
segments(x0 = hepn.quants.nostr[1], x1 = hepn.quants.nostr[2], y0 = 6, y1 = 6, lty = 2)
segments(x0 = pnhe.quants.nostr[1], x1 = pnhe.quants.nostr[2], y0 = 5, y1 = 5, lty = 2)
segments(x0 = pnpn.quants.str[1], x1 = pnpn.quants.str[2], y0 = 4, y1 = 4)
segments(x0 = hehe.quants.str[1], x1 = hehe.quants.str[2], y0 = 3, y1 = 3)
segments(x0 = hepn.quants.str[1], x1 = hepn.quants.str[2], y0 = 2, y1 = 2)
segments(x0 = pnhe.quants.str[1], x1 = pnhe.quants.str[2], y0 = 1, y1 = 1)
leg.text2 <- c("Unstructured", "Structured")
legend("bottomright", leg.text2, cex = .7, lty = c(2, 1), bty = "n")
axis(side = 2, at = seq(1:8), labels = c("Persist - Healthy", "Healthy - Persist", 
                                         "Healthy - Healthy", "Persist - Persist",
                                         "Persist - Healthy", "Healthy - Persist",
                                         "Healthy - Healthy", "Persist - Persist"), las = 2)

#-----------------------------------------------------#
#-- CCS/Ro/gamma figure: Keeling and Rohani pg. 208 --#
#-----------------------------------------------------#
# Proportion of animals collared 
collar.prop <- .1

# Host birth rate 
  # (needs to be pulsed... Gaussian over 1 month)
nu <- .7
  
# Host life expectancy (in days)
mu <- 1 / (20 * 365)

# Movement rate: 
    # assume one collared animal moves between two HC pops every year
    # we see collar.prop proportion of moves
    # so we see 1 move, there are actually 1/collar.prop moves occurring
    # each pop has chance 1/15 of getting move
    # so delta = 1 / (15 * collar.prop) * propI * sqrt(N)
collar.prop <- .25
N <- 1000
delta <- 1 / (.15 * collar.prop * N)

# P(disease-induced mort | infection) # see Keeling and Rohani p. 34
p <- .1 # probability of dying from infection

# Ro (Keeling and Rohani pg. 35)
Ro <- (beta * (1 - p) * nu) / ((mu + gamma) * mu) 

#----------------#
#-- Sim Design --#
#----------------#

# Goal: 
  # find pop size that experiences one fade-out per year
  # return persistence period

# Structure:
  # daily timestep for 100 years
  # put intro 6 months prior to birthpulse (so that intro is in rut)
Tmax <- 365 * 10

  # simulate over panel of 
    # 1. N's; run 100 sims/N
    # 2. Ro's: run systematically over logged range of Ro from 1 to 20
    # 3. 1/gamma: run at 2-wk; 3 mon; 6 mon; 1 year

N.levels <- c(10, 100, 1000, 10000)
Ro.levels <- c(1.25, 1.5, 3.0, 6, 9, 18)
gamma.levels <- c(1/10, 1/100, 1/1000)
beta.levels  <- (Ro.levels * (mu + gamma.levels) * mu) / ((1 - p) * nu)
param.mat <- as.data.frame(expand.grid(list(N.levels, Ro.levels, gamma.levels)))
names(param.mat) <- c("N", "Ro", "gamma")


sir.sim <- function(param.mat, nu, mu, p, delta, Tmax){
  X <- Y <- Z <- matrix(NA, nrow = dim(param.mat)[1], ncol = Tmax)
  pers.time <- matrix(0, nrow = dim(param.mat)[1], ncol = Tmax)
  pers.time.vec <- rep(NA, dim(param.mat)[1])
  
  # loop over parameter space
  for(j in 1:dim(param.mat)[1]){
    # initialize X, Y, Z (= S, I, R respectively)
    X[j, 1] <- param.mat$N[j] - 1
    Y[j, 1] <- 1
    Z[j, 1] <- 0
      
    # determine beta
    beta <- (param.mat$Ro[j] * (mu + param.mat$gamma[j]) * mu) / ((1 - p) * nu)
    
    # loop over timesteps
    for(i in 2:Tmax){
      # stop loop if Y[i - 1] == 0
      if(Y[j, i - 1] == 0){
        X[j, i] <- X[j, i - 1]
        Y[j, i] <- Y[j, i - 1]
        Z[j, i] <- Z[j, i - 1]
        pers.time[j, i] <- i
      } else { # if Y[i - 1] != 0
      # 1. Is timestep in birth pulse?
        # bp is truncated Gaussian for one month centered at 182 (from 167 to 197)
        if((i %% 365 < 197) & (i %% 365 > 167 )){ # if in birthpulse, do this
          births <- (nu * (X[j, i - 1] + Y[j, i - 1] + Z[j, i - 1])) / (30 * 2) # "*2" because only females reproduce
          X[j, i] <- max(0, floor(X[j, i - 1] + births - mu * X[j, i - 1] - beta * X[j, i - 1] * Y[j, i - 1]))
        } else { # if not in birthpulse, do this
          X[j, i] <- max(0, floor(X[j, i - 1] - mu * X[j, i - 1] - beta * X[j, i - 1] * Y[j, i - 1]))         
        }
        Y[j, i] <- max(0, floor(Y[j, i - 1] + beta * X[j, i - 1] * Y[j, i - 1]  + delta * (X[j, i - 1] + Y[j, i - 1] + Z[j, i - 1]) - ((param.mat$gamma[j] + mu) / (1 - p)) * Y[j, i - 1]))
        Z[j, i] <- max(0, floor(Z[j, i - 1] + param.mat$gamma[j] * Y[j, i - 1] - mu * Z[j, i - 1] - delta * (X[j, i - 1] + Y[j, i - 1] + Z[j, i - 1])))
      }
      print(paste(j, "_", i, sep = ""))
    }
    pers.time.vec[j] <- which(pers.time[j, ] != 0)[1]
  }
  return(list(X = X, Y = Y, Z = Z, pers.time = pers.time, pers.time.vec = pers.time.vec, param.mat = param.mat))
}

# test
Tmax <- 1000
N.levels <- c(10, 100, 1000)
Ro.levels <- c(1.25, 3.0,  9)
gamma.levels <- c(1/100, 1/1000)
param.mat <- as.data.frame(expand.grid(list(N.levels, Ro.levels, gamma.levels)))
names(param.mat) <- c("N", "Ro", "gamma")


sir.test <- sir.sim(param.mat = param.mat, nu = nu, mu = mu, p = p, delta = delta, Tmax = Tmax)
sir.test$pers.time.vec
param.mat.post <- cbind(sir.test$param.mat, sir.test$pers.time.vec)

sir.stoch.sim <- function(param.mat, nu, mu, p, delta, Tmax){
  X <- Y <- Z <- matrix(NA, nrow = dim(param.mat)[1], ncol = Tmax)
  pers.time <- matrix(0, nrow = dim(param.mat)[1], ncol = Tmax)
  pers.time.vec <- rep(NA, dim(param.mat)[1])
  
  # loop over parameter space
  for(j in 1:dim(param.mat)[1]){
    # initialize X, Y, Z (= S, I, R respectively)
    X[j, 1] <- param.mat$N[j] - 1
    Y[j, 1] <- 1
    Z[j, 1] <- 0
    
    # determine beta
    beta <- (param.mat$Ro[j] * (mu + param.mat$gamma[j]) * mu) / ((1 - p) * nu)
    
    # loop over timesteps
    for(i in 2:Tmax){
      # stop loop if Y[i - 1] == 0
      if(Y[j, i - 1] == 0){
        X[j, i] <- X[j, i - 1]
        Y[j, i] <- Y[j, i - 1]
        Z[j, i] <- Z[j, i - 1]
        pers.time[j, i] <- i
      } else { # if Y[i - 1] != 0
        # 1. Is timestep in birth pulse?
        # bp is truncated Gaussian for one month centered at 182 (from 167 to 197)
        if((i %% 365 < 197) & (i %% 365 > 167 )){ # if in birthpulse, do this
          births <- (nu * (X[j, i - 1] + Y[j, i - 1] + Z[j, i - 1])) / (30 * 2) # "*2" because only females reproduce
          X[j, i] <- max(0, floor(X[j, i - 1] + births - mu * X[j, i - 1] - beta * X[j, i - 1] * Y[j, i - 1]))
        } else { # if not in birthpulse, do this
          X[j, i] <- max(0, floor(X[j, i - 1] - mu * X[j, i - 1] - beta * X[j, i - 1] * Y[j, i - 1]))         
        }
        Y[j, i] <- max(0, floor(Y[j, i - 1] + beta * X[j, i - 1] * Y[j, i - 1]  + delta * (X[j, i - 1] + Y[j, i - 1] + Z[j, i - 1]) - ((param.mat$gamma[j] + mu) / (1 - p)) * Y[j, i - 1]))
        Z[j, i] <- max(0, floor(Z[j, i - 1] + param.mat$gamma[j] * Y[j, i - 1] - mu * Z[j, i - 1] - delta * (X[j, i - 1] + Y[j, i - 1] + Z[j, i - 1])))
      }
      print(paste(j, "_", i, sep = ""))
    }
    pers.time.vec[j] <- which(pers.time[j, ] != 0)[1]
  }
  return(list(X = X, Y = Y, Z = Z, pers.time = pers.time, pers.time.vec = pers.time.vec, param.mat = param.mat))
}

