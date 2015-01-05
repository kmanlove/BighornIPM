# BHS IPM Round 11 -- no marrays; ewe-age- and disease-status-specific reproduction; ewe-age- and disease-status-specific sls
# November 11, 2014
# 0. Load required packages.
require(rjags)
require(runjags)

# 1. Load data.

studysheep <- read.csv("~/work/Kezia/Research/EcologyPapers/ClustersAssocations_V2/ClustersAssociations/Data/RevisedData_11Sept2013/Study_sheep_toothage_original_012612.csv", header = T, sep = "\t")
lambs <- read.csv("~/work/Kezia/Research/EcologyPapers/ClustersAssocations_V2/ClustersAssociations/Data/MergedLambData_26Mar2013.csv", header = T)
compd.data <- read.csv("~/work/Kezia/Research/EcologyPapers/ClustersAssocations_V2/ClustersAssociations/Data/compiled_data_summary_130919.csv", header = T, sep = "")
compd.data <- subset(compd.data, !(Pop == "Imnaha" & year <= 1999))
compd.data$PNIndLambs <- ifelse((compd.data$Pop == "Asotin" & compd.data$year <= 2010) | (compd.data$Pop == "BigCanyon" & compd.data$year <= 2000) | (compd.data$Pop == "MuirCreek" & compd.data$year <= 2000), 1, ifelse(is.na(compd.data$CLASS) == T, NA, 2))
compd.data$PNIndEwes <- ifelse((compd.data$Pop == "Asotin" & compd.data$year <= 2010) | (compd.data$Pop == "BigCanyon" & compd.data$year <= 1999) | (compd.data$Pop == "MuirCreek" & compd.data$year <= 1999), 1, ifelse(is.na(compd.data$CLASS) == T, NA, 2))

compd.data <- compd.data[-476, ]
compd.data$Pop <- factor(compd.data$Pop)
compd.data <- subset(compd.data, year >= 1997 & year <= 2012)

# extract ewes with tooth ages
ewes.with.teeth <- subset(studysheep, SEX == "F" & is.na(Tooth_Age) == F)

# extract lambs born to ewes with tooth ages
lambs.with.dam.age <- subset(lambs, EWEID %in% levels(factor(ewes.with.teeth$ID)))

# rename populations in ewes.with.teeth
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
# cut compd.data down to just the pops in ewes.with.teeth, and just years 1997 and forward
compd.data <- subset(compd.data, Pop %in% levels(factor(ewes.with.teeth$END_Population)) & year >= 1997)
compd.data$Pop <- factor(compd.data$Pop)
factor.list <- list(compd.data$Pop, compd.data$year)
Oad <- tapply(compd.data$Ewes, factor.list, sum) # table observed adults in each pop-year
Ojuv <- tapply(compd.data$Lambs, factor.list, sum) # table observed lambs in each pop-year
compd.data$NoFemRem.nonas <- ifelse(is.na(compd.data$NoFemRem) == T, 0, compd.data$NoFemRem)
compd.data$NoFemRel.nonas <- ifelse(is.na(compd.data$NoFemRel) == T, 0, compd.data$NoFemRel)
#Osls <- round(tapply(compd.data$RadEwesWLambs * compd.data$SumLambSurv, factor.list, sum))
#RadEwes <- round(tapply(compd.data$RadEwes, factor.list, sum))
RemovedEwes <- tapply(compd.data$NoFemRem.nonas, factor.list, sum)
AddedEwes <- tapply(compd.data$NoFemRel.nonas, factor.list, sum)
for(j in 1:12){
  RemovedEwes[j, ] <- ifelse(is.na(RemovedEwes[j, ]) == T, 0, RemovedEwes[j, ])
  AddedEwes[j, ] <- ifelse(is.na(AddedEwes[j, ]) == T, 0, AddedEwes[j, ])
  #  Osls[j, ] <- ifelse(is.na(Osls[j, ]) == T, 0, Osls[j, ])
  #  RadEwes[j, ] <- ifelse(is.na(RadEwes[j, ]) == T, 0, RadEwes[j, ])
}

n.years <- 2012 - 1996

# pop counts 
y <- Ojuv + Oad # total ewes and lambs observed in each pop-year

# build j x t matrix of pop-year disease statuses (separate for lambs and adults, to separate out winter outbreak events)
popyr.dis.status.lambs <- popyr.dis.status.adults <- matrix(NA, nrow = length(levels(factor(ewe.pop.ind))), ncol = 2012 - 1996)
for(i in 1:length(levels(factor(ewe.pop.ind)))){
  for(j in 1 : (2012 - 1996)){
    k <- subset(compd.data, as.character(Pop) == as.character(levels(factor(ewe.pop.ind)))[i] & year == (1996 + j))
    popyr.dis.status.lambs[i, j] <- ifelse(dim(k)[1] == 0, NA, k$PNIndLambs)
    popyr.dis.status.adults[i, j] <- ifelse(dim(k)[1] == 0, NA, k$PNIndEwes)
  }
}

#popyr.dis.status <- ifelse(is.na(popyr.dis.status) == T, 3, popyr.dis.status)
popyr.dis.status.lambs <- ifelse(is.na(popyr.dis.status.lambs) == T, 3, popyr.dis.status.lambs)
popyr.dis.status.adults <- ifelse(is.na(popyr.dis.status.adults) == T, 3, popyr.dis.status.adults)

# CJS data 
# recode all entries prior to 1997 with 1997 as entry bioyr (to match with population counts strings in IPM)
ewes.with.teeth$ENTRY_BIOYR2 <- ifelse(ewes.with.teeth$ENTRY_BIOYR <= 1996, 1996, ewes.with.teeth$ENTRY_BIOYR)
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
ewe.surv.list <- vector("list", dim(ewes.with.teeth)[1])
ch.full <- ch.full2 <- ewe.age <- matrix(NA, nrow = dim(ewes.with.teeth)[1], ncol = 2012 - 1996)
ewe.pop.ind <- rep(NA, dim(ewes.with.teeth)[1])
for(i in 1:dim(ch.full)[1]){
  k <- subset(ewes.with.teeth, as.character(ID) == as.character(ewes.with.teeth$ID)[i])
  years <- seq(k$ENTRY_BIOYR2, k$END_BIOYR)
  ewe.pop.ind[i] <- as.character(k$END_Population)
  ewe.surv.status <- rep(NA, length(years))
  ch.full[i, ] <- c(rep(0, (k$ENTRY_BIOYR2[1] - 1996)), rep(1, floor(k$END_BIOYR) - floor(k$ENTRY_BIOYR2)), rep(0, 2012-(floor(k$END_BIOYR[1]))))         
  ewe.age[i, ] <- c(rep(0, (k$ENTRY_BIOYR2[1] - 1996)),  seq(floor(k$AENTRY), (floor(k$AENTRY) + (floor(k$END_BIOYR) - floor(k$ENTRY_BIOYR2)))), rep(0, 2012-(floor(k$END_BIOYR[1]) + 1)))    
  for(j in 1:length(years)){
    ewe.surv.status[j] <- ifelse(as.character(years[j]) == as.character(k$END_BIOYR[1]), "died", "survived")
  }
  ewe.surv.list[[i]] <- data.frame(cbind(years, rep(ewe.pop.ind[i], length(years)), ewe.surv.status, seq(floor(k$AENTRY), (floor(k$AENTRY) + (floor(k$END_BIOYR) - floor(k$ENTRY_BIOYR2))))))
}

ewe.age <- ifelse(ewe.age == 0, 20, ewe.age)
ewe.pop.ind.num <- as.numeric(as.factor(ewe.pop.ind))

age.spec.ewe.surv <- do.call(rbind, ewe.surv.list)
names(age.spec.ewe.surv) <- c("years", "pop", "ewe.surv.status", "ewe.age")
age.spec.ewe.surv$years <- as.numeric(as.character(age.spec.ewe.surv$years))
age.spec.ewe.surv$ewe.age <- as.numeric(as.character(age.spec.ewe.surv$ewe.age))
age.spec.ewe.surv <- subset(age.spec.ewe.surv, years >= 1997 & ewe.age >= 1)
ewe.surv.age <- as.numeric(age.spec.ewe.surv$ewe.age)
ewe.surv.status <- ifelse(as.numeric(age.spec.ewe.surv$ewe.surv.status) == 2, 1, 0)
ewe.surv.year <- age.spec.ewe.surv$years - 1996
ewe.surv.pop <- age.spec.ewe.surv$pop

he.ewe.yrs <- subset(age.spec.ewe.surv, (pop == "Asotin" & years <= 2010) | (pop == "MuirCreek" & years <= 1999) | (pop == "BigCanyon" & years <= 1999)  | (pop == "Imnaha" & years <= 1999))
he.ewe.yrs$age.class <- ifelse(he.ewe.yrs$ewe.age >= 1 & he.ewe.yrs$ewe.age < 3, 2, 
                               ifelse(he.ewe.yrs$ewe.age >= 3 & he.ewe.yrs$ewe.age < 7, 3,
                                      ifelse(he.ewe.yrs$ewe.age >= 7 & he.ewe.yrs$ewe.age < 16, 4,
                                             ifelse(he.ewe.yrs$ewe.age >= 16, 5, 1))))
table(he.ewe.yrs$age.class, he.ewe.yrs$ewe.surv.status)



# create vector with occasion of marking:
get.first <- function(x) min(which(x != 0))
f.init <- apply(ch.full, 1, get.first)

# for now, pull out ewes who didn't survive their collaring year (whose rows in ch are entirely 0's)
ch <- ch.full[-which(f.init == "Inf"), ]
f <- apply(ch, 1, get.first)

# build a x 1 vector of age-class specifications (maps 1:18 age in years to 1:5 age in age-class)
# age class 6 is NOT CURRENTLY MARKED. 
age.class.ind <- c(1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6) 

#---------------------------------------------#
#-- ewes followed by age, pop, health class --#
#---------------------------------------------#
ewes.he.aso <- subset(ewes.with.teeth, (Population == "Asotin" & ENTRY_BIOYR <= 2010 & DEAD == 1))
ewes.he.bcan <- subset(ewes.with.teeth, (Population == "BigCanyon" & ENTRY_BIOYR <= 1999 & DEAD == 1))
ewes.he.muir <- subset(ewes.with.teeth, (Population == "MuirCreek" & ENTRY_BIOYR <= 1999 & DEAD == 1))
ewes.he.imn <- subset(ewes.with.teeth,  (Population == "Imnaha" & ENTRY_BIOYR <= 1999 & DEAD == 1))



#----------------------------------#
#-- 2.4. CJS weaning --------------#
#----------------------------------#
# This is weaning without production: get a 1 for weaning, 0 for censoring, NA for no record in lamb data  --#
# CENSOR2 == 0 implies lamb survived; CENSOR2 == 1 implies lamb died
ewe.wean.list <- vector("list", dim(ewes.with.teeth)[1])
for(i in 1:length(ewe.wean.list)){
  k <- subset(ewes.with.teeth, ID == levels(factor(ewes.with.teeth$ID))[i])
  years <- seq(k$ENTRY_BIOYR2, k$END_BIOYR)
  ewe.age.wean <- seq(floor(as.numeric(as.character(k$AENTRY))), floor(as.numeric(as.character(k$AENTRY))) + length(years) - 1)
  wean.status <- rep(NA, length(years))
  pop.name <- rep(as.character(k$END_Population[1]), length(years))
  for(j in 1:length(years)){
    wean.year <- subset(lambs, EWEID == levels(factor(ewes.with.teeth$ID))[i] & YEAR == years[j])
    wean.status[j] <- ifelse(dim(wean.year)[1] == 0, NA, ifelse(wean.year$CENSOR2 == 0, 1, 0))
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
ewe.wean.year <- as.numeric(as.factor(as.numeric(as.character(age.spec.ewe.wean$years))))
ewe.wean.success <- as.numeric(as.character(age.spec.ewe.wean$wean.status))
table(age.spec.ewe.wean$pop.name, age.spec.ewe.wean$wean.status)


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
    for(s in 1:n.adsurvs){
    logit(phi.individ.adsurv[s]) <- beta.adsurv[popyr.dis.status.adults[ewe.surv.pop[s], ewe.surv.year[s]], age.class.ind[ewe.surv.age[s]]] + time.re.adsurv[ewe.surv.year[s]]
    } #s
    
    for(w in 1:n.weans){
    logit(phi.individ.wean[w]) <- beta.wean[popyr.dis.status.lambs[ewe.wean.pop[w], ewe.wean.year[w]], age.class.ind[ewe.wean.age[w]]] + time.re.wean[ewe.wean.year[w]]
    } #w
    
    # get phi estimates for each popyr (j) 
    for(t in 1:n.years){
    for(j in 1:n.pops){
    for(a in 1:n.ages){
    logit(phi.popyr.adsurv[j, t, a]) <- beta.adsurv[popyr.dis.status.adults[j, t], age.class.ind[a]] + time.re.adsurv[t]
    logit(phi.popyr.wean[j, t, a]) <- beta.wean[popyr.dis.status.lambs[j, t], age.class.ind[a]] + time.re.wean[t]
    #-- this logit pulls out the effect for each age-class in pop-year i, using pop-year i's disease status and the time re
    } #a
    } #j
    } #t
    
    # Specificy priors on the 2-d matrix of betas (called in the CJS logit survival function) and the time re.
    for(t in 1:(n.years)){
    time.re.adsurv[t] ~ dnorm(0, tau.time.adsurv) #-- random system-wide year effect
    time.re.wean[t] ~ dnorm(0, tau.time.wean) #-- random system-wide year effect
    }
    
    for(d in 1:n.dis.states){
    for(a in 1:n.age.classes){
    beta.adsurv[d , a] ~ dnorm(0, 0.01)T(-10, 10)
    beta.wean[d , a] ~ dnorm(0, 0.01)T(-10, 10)
    } #a
    } #d
    
    # hyperpriors for time.re's
    sigma.time.adsurv ~ dunif(0, 10)
    tau.time.adsurv <- pow(sigma.time.adsurv, -2)
    sigma.time2.adsurv <- pow(sigma.time.adsurv, 2)
    
    sigma.time.wean ~ dunif(0, 10)
    tau.time.wean <- pow(sigma.time.wean, -2)
    sigma.time2.wean <- pow(sigma.time.wean, 2)   
    
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
    Nwean[j, t, a] ~ dbin(phi.popyr.wean[j, t, a], N[j, t, a])  # ewes age first, reproduce second. 
    }
    N[j, t, 1] <- sum(Nwean[j, t, 2:18])
    for(a in 2:18){ # loop for age-specific survival and reproduction numbers
    N[j, t, a] ~ dbin(phi.popyr.adsurv[j, t, a - 1], N[j, t - 1, a - 1])
    Njuv.agespec[j, t, a] ~ dbin(phi.popyr.wean[j, t, a], N[j, t, a])
    } #a
    Nad[j, t] <- max(sum(N[j, t, 2:18]) - RemovedEwes[j, t] + AddedEwes[j, t], 1) 
    # subtract (known number of) removed ewes from pop count before doing observation (Oad)
    Njuv[j, t] <- sum(Njuv.agespec[j, t, 2:18])
    Ntot[j, t] <- Nad[j, t] + Njuv[j, t]
    } #t
    
    # 3.1.2. Observation process
    for (t in 1 : (n.years - 1)){
    # add in observation data
    Ojuv[j, t] ~ dpois(Njuv[j, t] + 1)
    Oad[j, t] ~ dpois(Nad[j, t] + 1)
    y[j, t] <- Ojuv[j, t] + Oad[j, t]
    } #t
    } #j
    
    #---------------------------------------------------------------------------------#
    #-- 3.2. Likelihood for adult survival data: CJS model: see Kery Schaub pg. 180 --#
    #---------------------------------------------------------------------------------#
    
    for(s in 1 : n.adsurvs){
    z.adsurv[s] ~ dbern(mu.adsurv[s])
    mu.adsurv[s] <- phi.individ.adsurv[s]
    } #s
    
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
ipm11.data <- list(
  ewe.surv.pop = ewe.surv.pop,
  ewe.surv.age = ewe.surv.age,
  ewe.surv.year = ewe.surv.year,
  z.adsurv = ewe.surv.status,
  n.adsurvs = length(ewe.surv.status),
  n.years = n.years,
  n.pops = dim(popyr.dis.status.adults)[1],
  n.age.classes = (length(levels(factor(age.class.ind)))),
  n.ages = 18,
  n.dis.states = length(levels(factor(popyr.dis.status.adults))),
  age.class.ind = age.class.ind,
  popyr.dis.status.adults = popyr.dis.status.adults,
  popyr.dis.status.lambs = popyr.dis.status.lambs,
  ewe.wean.pop = ewe.wean.pop,
  ewe.wean.age = ewe.wean.age,
  ewe.wean.year = ewe.wean.year,
  z.wean = ewe.wean.success,
  n.weans = dim(age.spec.ewe.wean)[1],
  Ojuv = Ojuv,
  Oad = Oad,
  AddedEwes = AddedEwes,
  RemovedEwes = RemovedEwes
)



# initial values
n.occasions <- 16

# function to create ch.inits
ch.init <- function(ch, f){
  for(i in 1:dim(ch)[1]){
    ch[i, 1 : f[i]] <- NA
  }
  return(ch)
}

ipm11.inits <- function(){
  list(
    sigma.time.adsurv = runif(1, 0, 10),
    sigma.time.wean = runif(1, 0, 10), 
    sigma.y = runif(1, 0, 10)
  )
}


# parameters to monitor
ipm11.params <- c("beta.adsurv", "beta.wean", "sigma.time.wean", "sigma.time.adsurv")

# mcmc settings
ni <- 2000
nt <- 3
nb <- 1000
nc <- 3

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
convg.diags <- gelman.diag(ipm11.coda)
#write.csv(convg.diags[[1]], "~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Posteriors/IPM/MoviDef/GelmanRubinDiags_26Dec2014.csv")

coda.summary.obj.11 <- summary(ipm11.coda)
#write.csv(coda.summary.obj.11[[2]], "~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Posteriors/IPM/MoviDef/PosteriorQuantiles_26Dec2014.csv")
row.names(coda.summary.obj.11[[2]])

beta.posts.adsurv <- coda.summary.obj.11[[2]][1:18, ]
beta.posts.wean <- coda.summary.obj.11[[2]][19:36, ]

# plot betas out
plot.cols <- c("white", "black", "red")
par(mfrow = c(1, 2))
plot(-1, -1, ylim = c(0, 1), xlim = c(3, 15), ylab = "Probability of survival", xlab = "Class")
for(i in 3:15){
  segments(x0 = i, x1 = i, y0 = (exp(beta.posts.adsurv[i, 1])) / (1 + exp(beta.posts.adsurv[i, 1])), y1 = exp(beta.posts.adsurv[i, 5]) / (1 + exp(beta.posts.adsurv[i, 5])), col = plot.cols[(i %% 3 + 1)], lty = (i %% 3 + 1), lwd = 2)
  segments(x0 = i - 0.25, x1 = i + 0.25, y0 = (exp(beta.posts.adsurv[i, 3])) / (1 + exp(beta.posts.adsurv[i, 3])), y1 = exp(beta.posts.adsurv[i, 3]) / (1 + exp(beta.posts.adsurv[i, 3])), col = plot.cols[(i %% 3 + 1)], lty = (i %% 3 + 1), lwd = 2)
}

plot(-1, -1, ylim = c(0, 1), xlim = c(3, 15), ylab = "Probability of weaning", xlab = "Class")
for(i in 3:15){
  segments(x0 = i, x1 = i, y0 = (exp(beta.posts.wean[i, 1])) / (1 + exp(beta.posts.wean[i, 1])), y1 = (exp(beta.posts.wean[i, 5])) / (1 + exp(beta.posts.wean[i, 5])), col = plot.cols[(i %% 3 + 1)], lty = (i %% 3 + 1), lwd = 2)
  segments(x0 = i - 0.25, x1 = i + 0.25, y0 = (exp(beta.posts.wean[i, 3])) / (1 + exp(beta.posts.wean[i, 3])), y1 = exp(beta.posts.wean[i, 3]) / (1 + exp(beta.posts.wean[i, 3])), col = plot.cols[(i %% 3 + 1)], lty = (i %% 3 + 1), lwd = 2)
}
