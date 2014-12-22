# BHS IPM Round 11 -- no marrays; ewe-age- and disease-status-specific reproduction; ewe-age- and disease-status-specific sls
# November 11, 2014

# test branch. 

# 0. Load required packages
require(rjags)
require(runjags)

# 1. Load data.

studysheep <- read.csv("~/work/Kezia/Research/EcologyPapers/ClustersAssocations_V2/ClustersAssociations/Data/RevisedData_11Sept2013/Study_sheep_toothage_original_012612.csv", header = T, sep = "\t")
lambs <- read.csv("~/work/Kezia/Research/EcologyPapers/ClustersAssocations_V2/ClustersAssociations/Data/MergedLambData_26Mar2013.csv", header = T)
compd.data <- read.csv("~/work/Kezia/Research/EcologyPapers/ClustersAssocations_V2/ClustersAssociations/Data/compiled_data_summary_130919.csv", header = T, sep = "")
#compd.data$PNInd <- ifelse(compd.data$CLASS == c("HEALTHY", "ADULTS"), 1, ifelse(compd.data$CLASS %in% c("ALL_AGE", "ALL_AGE_SUSP", "LAMBS", "LAMBS_SUSP"), 2, NA))
compd.data$PNInd <- ifelse(compd.data$CLASS == c("HEALTHY"), 1, ifelse(compd.data$CLASS %in% c("ALL_AGE", "ALL_AGE_SUSP", "LAMBS", "LAMBS_SUSP", "ADULTS"), 2, NA))
compd.data <- compd.data[-476, ]
compd.data$Pop <- factor(compd.data$Pop)
compd.data <- subset(compd.data, year >= 1997 & year <= 2012)

# extract ewes with tooth ages
ewes.with.teeth <- subset(studysheep, SEX == "F" & is.na(Tooth_Age) == F)

# extract lambs born to ewes with tooth ages
lambs.with.dam.age <- subset(lambs, EWEID %in% levels(factor(ewes.with.teeth$ID)))

# rename populations
ewes.with.teeth$END_Population <- ifelse(
  ewes.with.teeth$END_Population == "AS", "Asotin", ifelse(
    ewes.with.teeth$END_Population == "BC", "BigCanyon", ifelse(
      ewes.with.teeth$END_Population == "BB", "BlackButte", ifelse(
        ewes.with.teeth$END_Population == "IM", "Imnaha", ifelse(
          ewes.with.teeth$END_Population == "LO", "Lostine", ifelse(
            ewes.with.teeth$END_Population == "UHCOR", "UHCOR", ifelse(
              ewes.with.teeth$END_Population == "RB", "Redbird", ifelse(
                ewes.with.teeth$END_Population == "WE", "Wenaha", ifelse(
                  ewes.with.teeth$END_Population == "SM", "SheepMountain", ifelse(
                    ewes.with.teeth$END_Population == "SC", "UpperSaddle", ifelse(
                      ewes.with.teeth$END_Population == "BRC", "BearCreek", ifelse(
                        ewes.with.teeth$END_Population == "MU", "MuirCreek", ifelse(
                          ewes.with.teeth$END_Population == "MY", "MyersCreek", NA)))))))))))))
ewes.with.teeth <- subset(ewes.with.teeth, is.na(END_Population) == F)

#-----------------------#
#-- INDEX DEFINITIONS --#
#-----------------------#
# i = (1, ..., nindividuals) (sampled ewes)
# j = (1, ..., npops) (pops)
# t = (1, ..., 2012 - 1997) (years)
# a = (1, ..., 18) (ewe ages in years)
# c = (1, ..., 5) (ewe ages in Festa-Bianchet age-classes)

# state-space data
factor.list <- list(compd.data$Pop, compd.data$year)
Oad <- tapply(compd.data$Ewes, factor.list, sum)
Ojuv <- tapply(compd.data$Lambs, factor.list, sum)
compd.data$NoFemRem.nonas <- ifelse(is.na(compd.data$NoFemRem) == T, 0, compd.data$NoFemRem)
Osls <- round(tapply(compd.data$RadEwesWLambs * compd.data$SumLambSurv, factor.list, sum))
RadEwes <- round(tapply(compd.data$RadEwes, factor.list, sum))
#RadEwesWLambsTable <- tapply(compd.data$RadEwesWLambs, factor.list, sum)
#SLSTable <- tapply(compd.data$SumLambSurv, factor.list, sum)
RemovedEwes <- tapply(compd.data$NoFemRem.nonas, factor.list, sum)
for(j in 1: 16){
  RemovedEwes[j, ] <- ifelse(is.na(RemovedEwes[j, ]) == T, 0, RemovedEwes[j, ])
  Osls[j, ] <- ifelse(is.na(Osls[j, ]) == T, 0, Osls[j, ])
  RadEwes[j, ] <- ifelse(is.na(RadEwes[j, ]) == T, 0, RadEwes[j, ])
}

n.years <- 2012 - 1997

# pop counts 
y <- Ojuv + Oad

# build j x t matrix of pop-year disease statuses
popyr.dis.status <- matrix(NA, nrow = length(levels(factor(ewe.pop.ind))), ncol = 2012 - 1997)
for(i in 1:length(levels(factor(ewe.pop.ind)))){
  for(j in 1 : (2012 - 1997)){
    k <- subset(compd.data, as.character(Pop) == as.character(levels(factor(ewe.pop.ind)))[i] & year == (1997 + j))
    popyr.dis.status[i, j] <- ifelse(dim(k)[1] == 0, NA, k$PNInd)
  }
}

popyr.dis.status <- ifelse(is.na(popyr.dis.status) == T, 3, popyr.dis.status)

# CJS data 
# recode all entries prior to 1997 with 1997 as entry bioyr (to match with population counts strings in IPM)
ewes.with.teeth$ENTRY_BIOYR2 <- ifelse(ewes.with.teeth$ENTRY_BIOYR <= 1996, 1997, ewes.with.teeth$ENTRY_BIOYR)
# recode a few specific animals to have entry bioyrs no earlier than the first year with compd.data for their resident population
sheep.ind <- which(ewes.with.teeth$ID == "02MY33")
sheep.ind2 <- which(ewes.with.teeth$ID == "02MY39") # she needs to come out -- died on 04/20/02.
sheep.ind3 <- which(ewes.with.teeth$ID == "02MY41")
sheep.ind4 <- which(ewes.with.teeth$ID == "01QZ70")
ewes.with.teeth[c(sheep.ind, sheep.ind3, sheep.ind4), ]$ENTRY_BIOYR2 <- 2002
ewes.with.teeth <- subset(ewes.with.teeth, ID != "02MY39")

# limit to 1997 - 2012 
# build i x t matrix of observation status for each ewe in each year
# build i x t matrix of true ewe-ages during each observation year
# NOTE: "ch" stands for "capture history" for each ewe. 
ch.full <- ch.full2 <- ewe.age <- matrix(NA, nrow = dim(ewes.with.teeth)[1], ncol = 2012 - 1997)
ewe.pop.ind <- rep(NA, dim(ewes.with.teeth)[1])
for(i in 1:dim(ch.full)[1]){
  k <- subset(ewes.with.teeth, as.character(ID) == as.character(ewes.with.teeth$ID)[i])
  ewe.pop.ind[i] <- as.character(k$END_Population)
  ch.full[i, ] <- c(rep(0, (k$ENTRY_BIOYR2[1] - 1997)), rep(1, floor(k$END_BIOYR) - floor(k$ENTRY_BIOYR2)), rep(0, 2012-(floor(k$END_BIOYR[1]))))       
  if((floor(k$END_BIOYR) - floor(k$ENTRY_BIOYR2) - 1) >= 0){
    ewe.age[i, ] <- c(rep(0, (k$ENTRY_BIOYR2[1] - 1997)),  seq(floor(k$AENTRY), (floor(k$AENTRY) + (floor(k$END_BIOYR) - floor(k$ENTRY_BIOYR2) - 1))), rep(0, 2012-(floor(k$END_BIOYR[1]))))    
  } else ewe.age[i, ] <- rep(0, dim(ch.full)[2])
}

# ewe.age <- ewe.age + 1
ewe.age <- ifelse(ewe.age == 0, 20, ewe.age)
ewe.pop.ind.num <- as.numeric(as.factor(ewe.pop.ind))

# create vector with occasion of marking:
get.first <- function(x) min(which(x != 0))
f.init <- apply(ch.full, 1, get.first)

# for now, pull out ewes who didn't survive their collaring year (whose rows in ch are entirely 0's)
ch <- ch.full[-which(f.init == "Inf"), ]
f <- apply(ch, 1, get.first)

# build a x 1 vector of age-class specifications (maps 1:18 age in years to 1:5 age in age-class)
#age.class.ind <- c(1, 2, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6)
#age.class.ind <- c(1, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5)
#age.class.ind <- c(1, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6)
age.class.ind <- c(1, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6)



#---------------------------------------#
#-- CJS reproduction --------------#
#---------------------------------------#
# This is PRODUCTION, not WEANING --#
#-- fill in ewes who didn't reproduce --#
ewe.repro.list <- vector("list", dim(ewes.with.teeth)[1])
for(i in 1:length(ewe.repro.list)){
  k <- subset(ewes.with.teeth, ID == levels(factor(ewes.with.teeth$ID))[i])
  years <- seq(k$ENTRY_BIOYR2, k$END_BIOYR)
  ewe.age.prod <- seq(floor(as.numeric(as.character(k$AENTRY))), floor(as.numeric(as.character(k$AENTRY))) + length(years) - 1)
  lamb.status <- rep(NA, length(years))
  pop.name <- rep(as.character(k$END_Population[1]), length(years))
  for(j in 1:length(years)){
    lamb.year <- subset(lambs, EWEID == levels(factor(ewes.with.teeth$ID))[i] & YEAR == years[j])
    lamb.status[j] <- ifelse(dim(lamb.year)[1] == 0, 0, 1)
  }
  ewe.repro.list[[i]] <- data.frame(cbind(years, pop.name, ewe.age.prod, lamb.status))
}

age.spec.ewe.prod <- do.call(rbind, ewe.repro.list)
age.spec.ewe.prod$years <- factor(age.spec.ewe.prod$years)
age.spec.ewe.prod$pop.name <- factor(age.spec.ewe.prod$pop.name)
ewe.prod.pop <- as.numeric(age.spec.ewe.prod$pop.name)
ewe.prod.age <- as.numeric(as.character(age.spec.ewe.prod$ewe.age.prod)) + 1
ewe.prod.year <- as.numeric(factor(age.spec.ewe.prod$years))
ewe.prod.success <- as.numeric(as.character(age.spec.ewe.prod$lamb.status))


#----------------------------------#
#-- 2.4. CJS weaning --------------#
#----------------------------------#
# This is WEANING, not PRODUCTION --#
# CENSOR2 == 0 implies lamb survived; CENSOR2 == 1 implies lamb died
# check with table for Aso vs. BB: table(lambs$CENSOR2, lambs$PN, lambs$HERD)
ewe.wean.list <- vector("list", dim(ewes.with.teeth)[1])
for(i in 1:length(ewe.wean.list)){
  k <- subset(ewes.with.teeth, ID == levels(factor(ewes.with.teeth$ID))[i])
  years <- seq(k$ENTRY_BIOYR2, k$END_BIOYR)
  ewe.age.wean <- seq(floor(as.numeric(as.character(k$AENTRY))), floor(as.numeric(as.character(k$AENTRY))) + length(years) - 1)
  wean.status <- rep(NA, length(years))
  pop.name <- rep(as.character(k$END_Population[1]), length(years))
  for(j in 1:length(years)){
    wean.year <- subset(lambs, EWEID == levels(factor(ewes.with.teeth$ID))[i] & YEAR == years[j])
    #    wean.given.lambed.status[j] <- ifelse(dim(wean.given.lambed.year)[1] == 0, NA, ifelse(wean.given.lambed.year$CENSOR2 == 0, 1, 0))
    wean.status[j] <- ifelse(dim(wean.year)[1] == 0, NA, ifelse(wean.year$CENSOR2 == 0, 1, 0))
#    wean.status[j] <- ifelse(dim(wean.year)[1] == 0, 0, ifelse(wean.year$CENSOR2 == 0, 1, 0))
  }
  ewe.wean.list[[i]] <- data.frame(cbind(years, pop.name, ewe.age.wean, wean.status))
}

#-- eliminate ewes who didn't reproduce --#
age.spec.ewe.wean <- do.call(rbind, ewe.wean.list)
age.spec.ewe.wean <- subset(age.spec.ewe.wean, is.na(wean.status) == F)
age.spec.ewe.wean$years <- factor(age.spec.ewe.wean$years)
age.spec.ewe.wean$pop.name <- factor(age.spec.ewe.wean$pop.name)
ewe.wean.pop <- as.numeric(age.spec.ewe.wean$pop.name)
ewe.wean.age <- as.numeric(as.character(age.spec.ewe.wean$ewe.age.wean)) + 1
ewe.wean.year <- as.numeric(factor(age.spec.ewe.wean$years))
ewe.wean.success <- as.numeric(as.character(age.spec.ewe.wean$wean.status))


#------------------------------------#
#-- IPM -----------------------------#
#------------------------------------#

# bug model
sink("ipm11.bug")
cat("
    model{
    
    #--------------------------#
    #-- Define priors ---------#
    #--------------------------#
    #-- state-space --#
    # observation error
    tauy <- pow(sigma.y, -2)
    sigma.y ~ dunif(0, 50)
    sigma2.y <- pow(sigma.y, 2)
    
    # initial population sizes. switch to Poisson to get integer-valued group sizes
    for(j in 1:n.pops){ #N[j = pop, t = year, a = age.class]. leave Njuv, Nad, Ntot as 2-D.
    for(a in 1:18){ # 
    N[j, 1, a] ~ dpois(10)  # all ages.
    } #a
    Nad[j, 1] <- sum(N[j, 1, 2:18])
    Njuv[j, 1] <- sum(N[j, 1, 1])
    Ntot[j, 1] <- Nad[j, 1] + Njuv[j, 1]
    } #j
    
    #-- survival and recapture probabilities for CJS --#
    #-- NOTE: for CJS, need phi calculated specifically for each individual during each year. 
    #-- So, CJS phis need to loop over nind and n.years. --#
    # get phi estimates for each individual (i) 
    for(i in 1:nind){
    for(t in f[i] : (n.years)){
#    logit(phi.individ.adsurv[i, t]) <- beta.adsurv[popyr.dis.status[ewe.pop.ind.num[i], t], age.class.ind[ewe.age[i, t]]] + time.re.adsurv[t]
    logit(phi.individ.adsurv[i, t]) <- beta.adsurv[popyr.dis.status[ewe.pop.ind.num[i], t], age.class.ind[ewe.age[i, t]]]
    #-- this logit pulls out the effect for ewe i's survival prob given her age class and her pop's current disease status
    } #t
    } #i
    
    for(r in 1:n.repros){
#    logit(phi.individ.repro[r]) <- beta.repro[popyr.dis.status[ewe.prod.pop[r], ewe.prod.year[r]], age.class.ind[ewe.prod.age[r]]] + time.re.repro[ewe.prod.year[r]]
    logit(phi.individ.repro[r]) <- beta.repro[popyr.dis.status[ewe.prod.pop[r], ewe.prod.year[r]], age.class.ind[ewe.prod.age[r]]] 
    # reproduction needs to be in its own loop, over number of repros (not number of ewes)
    } #r
    
    for(w in 1:n.weans){
    logit(phi.individ.wean[w]) <- beta.wean[popyr.dis.status[ewe.wean.pop[w], ewe.wean.year[w]], age.class.ind[ewe.wean.age[w]]] + time.re.wean[ewe.wean.year[w]]
#    logit(phi.individ.wean[w]) <- beta.wean[popyr.dis.status[ewe.wean.pop[w], ewe.wean.year[w]], age.class.ind[ewe.wean.age[w]]]
    } #w
    
    # get phi estimates for each popyr (j) 
    for(t in 1:n.years){
    for(j in 1:n.pops){
#    logit(phi.popyr.overwinter[j, t]) <- beta.overwinter[popyr.dis.status[j, t]] + time.re.overwinter[t]
#    logit(phi.popyr.overwinter[j, t]) <- beta.overwinter[popyr.dis.status[j, t]] 
    for(a in 1:n.ages){
    logit(phi.popyr.adsurv[j, t, a]) <- beta.adsurv[popyr.dis.status[j, t], age.class.ind[a]]
    logit(phi.popyr.repro[j, t, a]) <- beta.repro[popyr.dis.status[j, t], age.class.ind[a]] 
#    logit(phi.popyr.wean[j, t, a]) <- beta.wean[popyr.dis.status[j, t], age.class.ind[a]] 
#     logit(phi.popyr.adsurv[j, t, a]) <- beta.adsurv[popyr.dis.status[j, t], age.class.ind[a]] + time.re.adsurv[t]
#     logit(phi.popyr.repro[j, t, a]) <- beta.repro[popyr.dis.status[j, t], age.class.ind[a]] + time.re.repro[t]
     logit(phi.popyr.wean[j, t, a]) <- beta.wean[popyr.dis.status[j, t], age.class.ind[a]] + time.re.wean[t]
#     #-- this logit pulls out the effect for each age-class in pop-year i, using pop-year i's disease status and the time re
    } #a
    } #j
    } #t
    
#     # Specificy priors on the 2-d matrix of betas (called in the CJS logit survival function) and the time re.
     for(t in 1:(n.years)){
#     time.re.adsurv[t] ~ dnorm(0, tau.time.adsurv) #-- random system-wide year effect
#     time.re.repro[t] ~ dnorm(0, tau.time.repro) #-- random system-wide year effect
     time.re.wean[t] ~ dnorm(0, tau.time.wean) #-- random system-wide year effect
#     time.re.overwinter[t] ~ dnorm(0, tau.time.overwinter) #-- random system-wide year effect
     }
    
    for(d in 1:n.dis.states){
    beta.overwinter[d] ~ dnorm(0, 0.01)T(-10, 10) # overwinter survival isn't mapped to ewe age. 
    for(a in 1:n.age.classes){
    beta.adsurv[d , a] ~ dnorm(0, 0.01)T(-10, 10)
    #}
    #    for(a in 2:n.age.classes){
    beta.repro[d , a] ~ dnorm(0, 0.01)T(-10, 10)
    beta.wean[d , a] ~ dnorm(0, 0.01)T(-10, 10)
    }
    }
    
#     # hyperpriors for time.re's
#     sigma.time.adsurv ~ dunif(0, 10)
#     tau.time.adsurv <- pow(sigma.time.adsurv, -2)
#     sigma.time2.adsurv <- pow(sigma.time.adsurv, 2)
#     
#     sigma.time.repro ~ dunif(0, 10)
#     tau.time.repro <- pow(sigma.time.repro, -2)
#     sigma.time2.repro <- pow(sigma.time.repro, 2)   
#     
    sigma.time.wean ~ dunif(0, 10)
    tau.time.wean <- pow(sigma.time.wean, -2)
    sigma.time2.wean <- pow(sigma.time.wean, 2)   
#     
#     sigma.time.overwinter ~ dunif(0, 10)
#     tau.time.overwinter <- pow(sigma.time.overwinter, -2)
#     sigma.time2.overwinter <- pow(sigma.time.overwinter, 2)   
    
    #----------------------------------------#
    #-- Likelihoods of the single datasets --#
    #----------------------------------------#
    
    #-------------------------------------------------------------------#
    #-- 3.1. Likelihood for population count data (state-space model) --#
    #-------------------------------------------------------------------#
    # 3.1.1. System process
    for(j in 1:n.pops){
    for(t in 2 : (n.years - 1)){
    # repro probs need to be checked....
    # loop to get number of offspring produced by each age-class last year
    for(a in 2:18){
    Nrepro[j, t, a] ~ dbin(phi.popyr.repro[j, t - 1, a - 1], N[j, t - 1, a - 1])
#    Nrepro[j, t, a] ~ dbin(phi.popyr.repro[j, t, a - 1], N[j, t - 1, a - 1])
    Nwean[j, t, a] ~ dbin(phi.popyr.wean[j, t, a], Nrepro[j, t, a]) 
#    Nwean[j, t, a] ~ dbin(phi.popyr.wean[j, t, a], Nrepro[j, t - 1, a - 1]) 
    # Note: Weaning updates are from last year in this version of the model
    }
    N[j, t, 1] <- sum(Nwean[j, t, 2:18])
    for(a in 2:18){
#    N[j, t, a] ~ dbin(phi.popyr.adsurv[j, t - 1, a - 1], N[j, t - 1, a - 1])
    N[j, t, a] ~ dbin(phi.popyr.adsurv[j, t, a - 1], N[j, t - 1, a - 1])
    Njuv.agespec[j, t, a] ~ dbin(phi.popyr.wean[j, t, a], N[j, t, a])
} #a
    Nad[j, t] <- max(sum(N[j, t, 2:18]) - RemovedEwes[j, t], 1) 
    # subtract (known number of) removed ewes from pop count before doing observation (Oad)
    #    Nad[j, t] <- sum(N[j, t, 2:18])
    Nfall[j, t] <- sum(Nwean[j, t, 2:18])
#    Njuv[j, t] ~ dbin(phi.popyr.overwinter[j, t], Nfall[j, t])
    Njuv[j, t] <- sum(Njuv.agespec[j, t, 2:18])
    Ntot[j, t] <- Nad[j, t] + Njuv[j, t]
    } #t
    
    # 3.1.2. Observation process
    for (t in 1 : (n.years - 1)){
    # add in my observation data
    #    Osls[j, t] ~ dbin(phi.popyr.wean[j, t, a], RadEwes[j, t])
    # need to modify model for beta.wean. Right now, all age-classes estimated independently. Need
    # single fixed effect to apply here (since I don't have ages on all RadEwes; need adjustment based solely on pop and year.)
    Ojuv[j, t] ~ dpois(Njuv[j, t] + 1)
    Oad[j, t] ~ dpois(Nad[j, t] + 1)
    y[j, t] <- Ojuv[j, t] + Oad[j, t]
    } #t
    } #j
    
    #-----------------------------------------------------------#
    #-- 3.2. Likelihood for adult survival data: CJS model --#
    #-----------------------------------------------------------#
    for(i in 1 : nind){
    # define state at first cpature
    z[i, f[i]] ~ dbern(1)
    for(t in (f[i] + 1) : n.years){
    # state process
    z[i, t] ~ dbern(mu.adsurv[i, t])
    mu.adsurv[i, t] <- phi.individ.adsurv[i, t - 1] * z[i, t - 1] 
    } #t
    } #i
    
    #-----------------------------------------------------------#
    #-- 3.3. Likelihood for reproduction data: CJS model --#
    #-----------------------------------------------------------#
    for(r in 1 : n.repros){
    # state process
    z.repro[r] ~ dbern(mu.repro[r])
    mu.repro[r] <- phi.individ.repro[r]
    # might be smart to add observation error in here...
    } #r
    
    #-----------------------------------------------------------#
    #-- 3.4. Likelihood for weaning data: CJS model --#
    #-----------------------------------------------------------#
    for(w in 1 : n.weans){
    # state process
    z.wean[w] ~ dbern(mu.wean[w])
    mu.wean[w] <- phi.individ.wean[w]
    # might be smart to add observation error in here...
    } #w
    
    }", fill = T)
sink()



#---------------------------#
#-- Bundle IPM data --------#
#---------------------------#
ipm11.data <- list(z = ch,
                   f = f, 
                   nind = dim(ch)[1], 
                   n.years = n.years,
                   n.pops = dim(popyr.dis.status)[1],
                   #                   n.age.classes = (length(levels(factor(age.class.ind))) + 1),
                   n.age.classes = (length(levels(factor(age.class.ind)))),
                   n.ages = 18,
                   n.dis.states = length(levels(factor(popyr.dis.status))),
                   ewe.age = ewe.age,
                   age.class.ind = age.class.ind,
                   popyr.dis.status = popyr.dis.status,
                   ewe.pop.ind.num = ewe.pop.ind.num,
                   ewe.prod.pop = ewe.prod.pop,
                   ewe.prod.age = ewe.prod.age,
                   ewe.prod.year = ewe.prod.year,
                   z.repro = ewe.prod.success,
                   n.repros = dim(age.spec.ewe.prod)[1],
                   ewe.wean.pop = ewe.wean.pop,
                   ewe.wean.age = ewe.wean.age,
                   ewe.wean.year = ewe.wean.year,
                   z.wean = ewe.wean.success,
                   n.weans = dim(age.spec.ewe.wean)[1],
                   Ojuv = Ojuv,
                   Oad = Oad,
#                   Osls = Osls,
#                   RadEwes = RadEwes,
                   RemovedEwes = RemovedEwes
)



# initial values
n.occasions <- 15

# function to create ch.inits
ch.init <- function(ch, f){
  for(i in 1:dim(ch)[1]){
    ch[i, 1 : f[i]] <- NA
  }
  return(ch)
}

ipm11.inits <- function(){
  list(
#     sigma.time.adsurv = runif(1, 0, 10),
#        sigma.time.repro = runif(1, 0, 10),
        sigma.time.wean = runif(1, 0, 10), 
#        sigma.time.overwinter = runif(1, 0, 10),
       #       mean.p.repro = runif(1, 0, 1),
       sigma.y = runif(1, 0, 10)
  )
}


# parameters to monitor
#ipm11.params <- c("beta.adsurv", "beta.repro", "beta.wean", "beta.overwinter")
#ipm11.params <- c("beta.adsurv", "beta.repro", "beta.wean", "beta.overwinter", "sigma.time.overwinter")
ipm11.params <- c("beta.adsurv", "beta.repro", "beta.wean", "sigma.time.wean")
#ipm11.params <- c("beta.adsurv", "beta.repro", "beta.wean", "beta.overwinter", "sigma.time.adsurv", "sigma.time.repro", "sigma.time.wean", "sigma.time.overwinter")

# mcmc settings
ni <- 1000
nt <- 3
nb <- 500
nc <- 3

# MAY NEED TO REINITIALIZE A NUMBER OF TIMES TO GET APPROPRIATE INITIAL VALUES. 

# call JAGS from R
ipm11.call <- jags.model("ipm11.bug",
                         data = ipm11.data,
                         inits = ipm11.inits,
                         n.chains = nc,
                         n.adapt = nb
)


update(ipm11.call, ni)

ipm11.coda <- coda.samples(ipm11.call,
                           ipm11.params,
                           ni)

summary(ipm11.coda)
gelman.diag(ipm11.coda)

coda.summary.obj.11 <- summary(ipm11.coda)
row.names(coda.summary.obj.11[[2]])
# beta.posts.adsurv <- coda.summary.obj.11[[2]][1:21, ]
# beta.posts.overwinter <- coda.summary.obj.11[[2]][22:24, ]
# beta.posts.repro <- coda.summary.obj.11[[2]][25:45, ]
# beta.posts.wean <- coda.summary.obj.11[[2]][46:66, ]
# 
# beta.posts.adsurv <- coda.summary.obj.11[[2]][1:18, ]
# beta.posts.overwinter <- coda.summary.obj.11[[2]][19:21, ]
# beta.posts.repro <- coda.summary.obj.11[[2]][22:39, ]
# beta.posts.wean <- coda.summary.obj.11[[2]][40:57, ]

beta.posts.adsurv <- coda.summary.obj.11[[2]][1:18, ]
beta.posts.repro <- coda.summary.obj.11[[2]][19:36, ]
beta.posts.wean <- coda.summary.obj.11[[2]][37:54, ]
# beta.posts.adsurv <- coda.summary.obj.11[[2]][1:15, ]
# beta.posts.overwinter <- coda.summary.obj.11[[2]][16:18, ]
# beta.posts.repro <- coda.summary.obj.11[[2]][19:33, ]
# beta.posts.wean <- coda.summary.obj.11[[2]][34:48, ]

# plot betas out
plot.cols <- c("white", "black", "red")
par(mfrow = c(2, 2))
plot(-1, -1, ylim = c(0, 1), xlim = c(3, 15), ylab = "Probability of survival", xlab = "Class")
for(i in 3:15){
  segments(x0 = i, x1 = i, y0 = (exp(beta.posts.adsurv[i, 1])) / (1 + exp(beta.posts.adsurv[i, 1])), y1 = exp(beta.posts.adsurv[i, 5]) / (1 + exp(beta.posts.adsurv[i, 5])), col = plot.cols[(i %% 3 + 1)], lty = (i %% 3 + 1), lwd = 2)
  segments(x0 = i - 0.25, x1 = i + 0.25, y0 = (exp(beta.posts.adsurv[i, 3])) / (1 + exp(beta.posts.adsurv[i, 3])), y1 = exp(beta.posts.adsurv[i, 3]) / (1 + exp(beta.posts.adsurv[i, 3])), col = plot.cols[(i %% 3 + 1)], lty = (i %% 3 + 1), lwd = 2)
}

plot(-1, -1, ylim = c(0, 1), xlim = c(3, 15), ylab = "Probability of reproducing", xlab = "Class")
for(i in 3:15){
  segments(x0 = i, x1 = i, y0 = (exp(beta.posts.repro[i, 1])) / (1 + exp(beta.posts.repro[i, 1])), y1 = (exp(beta.posts.repro[i, 5])) / (1 + exp(beta.posts.repro[i, 5])), col = plot.cols[(i %% 3 + 1)], lty = (i %% 3 + 1), lwd = 2)
  segments(x0 = i - 0.25, x1 = i + 0.25, y0 = (exp(beta.posts.repro[i, 3])) / (1 + exp(beta.posts.repro[i, 3])), y1 = exp(beta.posts.repro[i, 3]) / (1 + exp(beta.posts.repro[i, 3])), col = plot.cols[(i %% 3 + 1)], lty = (i %% 3 + 1), lwd = 2)
}

plot(-1, -1, ylim = c(0, 1), xlim = c(3, 15), ylab = "Probability of weaning", xlab = "Class")
for(i in 3:15){
  segments(x0 = i, x1 = i, y0 = (exp(beta.posts.wean[i, 1])) / (1 + exp(beta.posts.wean[i, 1])), y1 = (exp(beta.posts.wean[i, 5])) / (1 + exp(beta.posts.wean[i, 5])), col = plot.cols[(i %% 3 + 1)], lty = (i %% 3 + 1), lwd = 2)
  segments(x0 = i - 0.25, x1 = i + 0.25, y0 = (exp(beta.posts.wean[i, 3])) / (1 + exp(beta.posts.wean[i, 3])), y1 = exp(beta.posts.wean[i, 3]) / (1 + exp(beta.posts.wean[i, 3])), col = plot.cols[(i %% 3 + 1)], lty = (i %% 3 + 1), lwd = 2)
}

# plot(-1, -1, ylim = c(0, 1), xlim = c(0.5, 3.5), ylab = "Probability of surviving overwinter", xlab = "Class")
# for(i in 1:3){
#   segments(x0 = i, x1 = i, y0 = ((exp(beta.posts.overwinter[i, 1])) / (1 + exp(beta.posts.overwinter[i, 1]))), y1 = (exp(beta.posts.overwinter[i, 5])) / (1 + exp(beta.posts.overwinter[i, 5])), col = plot.cols[(i %% 3 + 1)], lty = (i %% 3 + 1), lwd = 2)
#   segments(x0 = i - 0.25, x1 = i + 0.25, y0 = (exp(beta.posts.overwinter[i, 3])) / (1 + exp(beta.posts.overwinter[i, 3])), y1 = exp(beta.posts.overwinter[i, 3]) / (1 + exp(beta.posts.overwinter[i, 3])), col = plot.cols[(i %% 3 + 1)], lty = (i %% 3 + 1), lwd = 2)
# }

