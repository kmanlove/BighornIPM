#---------------------------------------------#
#-- Code for recruitment vs. adult survival --#
#---------------------------------------------#
require(rjags)
require(runjags)
require(popbio)
require(Matrix)


studysheep <- read.csv("~/work/Kezia/Research/EcologyPapers/ClustersAssocations_V2/ClustersAssociations/Data/RevisedData_11Sept2013/Study_sheep_toothage_original_012612.csv", header = T, sep = "\t")
# check for new studysheep
lambs <- read.csv("~/work/Kezia/Research/EcologyPapers/ClustersAssocations_V2/ClustersAssociations/Data/MergedLambData_26Mar2013.csv", header = T)
compd.data <- read.csv("~/work/Kezia/Research/EcologyPapers/ClustersAssocations_V2/ClustersAssociations/Data/compiled_data_summary_130919.csv", header = T, sep = "")
# call compd.data from Dropbox -> BHS -> Data -> Compiled Data -> Current. 15 Sept 2014 vsn!!

#---------------------------------------------#
#-- 1) Age-specific reproduction -------------#
#---------------------------------------------#
#-- 1) build data. Use Festa-Bianchet 2006 stages: lamb, yearling, 2-7, 8-13, >13 --#
#-- keep each year cohort distinct from the others --#
#----- m-array runs from 1995 - 2010 --#
#-- data --#

ewes.with.teeth <- subset(studysheep, SEX == "F" & is.na(Tooth_Age) == F)
lambs.with.dam.age <- subset(lambs, EWEID %in% levels(factor(ewes.with.teeth$ID)))
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


# #-- fill in ewes who didn't reproduce --#
# ewe.list <- vector("list", dim(ewes.with.teeth)[1])
# for(i in 1:length(ewe.list)){
#   k <- subset(ewes.with.teeth, ID == levels(factor(ewes.with.teeth$ID))[i])
#   years <- seq(k$ENTRY_BIOYR, k$END_BIOYR)
#   ewe.age <- seq(floor(as.numeric(as.character(k$AENTRY))), k$Tooth_Age)
#   lamb.status <- pn.status <- rep(NA, length(years))
#   pop.name <- rep(as.character(k$Population[1]), length(years))
#   for(j in 1:length(years)){
#     #    pop <- subset(compd.data, as.character(Pop) == as.character(k$Population)[1] & year == years[j])
#     pop <- subset(compd.data, as.character(Pop) == as.character(k$END_Population)[1] & year == years[j])
#     lamb.year <- subset(lambs, EWEID == levels(factor(ewes.with.teeth$ID))[i] & YEAR == years[j])
#     lamb.status[j] <- ifelse(dim(lamb.year)[1] == 0, 2, lamb.year$CENSOR2 + 1)
#     #  CENSOR2 == 1 when lamb DIES/ CENSOR2 == 0 when lamb survives
#     #  Code no lamb as dead lamb. 
#     pn.status[j] <- ifelse(dim(pop)[1] == 0, NA, as.character(pop$CLASS[1]))
#   }
#   ewe.list[[i]] <- data.frame(cbind(years, ewe.age, lamb.status, pn.status, pop.name))
# }


#-- fill in ewes who didn't reproduce --#
ewe.list <- vector("list", dim(ewes.with.teeth)[1])
for(i in 1:length(ewe.list)){
  k <- subset(ewes.with.teeth, ID == levels(factor(ewes.with.teeth$ID))[i])
  years <- seq(k$ENTRY_BIOYR, k$END_BIOYR)
  ewe.age <- seq(floor(as.numeric(as.character(k$AENTRY))), k$Tooth_Age)
  lamb.status <- pn.status <- rep(NA, length(years))
  pop.name <- rep(as.character(k$Population[1]), length(years))
  for(j in 1:length(years)){
    #    pop <- subset(compd.data, as.character(Pop) == as.character(k$Population)[1] & year == years[j])
    pop <- subset(compd.data, as.character(Pop) == as.character(k$END_Population)[1] & year == years[j])
    lamb.year <- subset(lambs, EWEID == levels(factor(ewes.with.teeth$ID))[i] & YEAR == years[j])
    lamb.status[j] <- ifelse(dim(lamb.year)[1] == 0, NA, lamb.year$CENSOR2 + 1)
    #  CENSOR2 == 1 when lamb DIES/ CENSOR2 == 0 when lamb survives
    #  Code no lamb as dead lamb. 
    pn.status[j] <- ifelse(dim(pop)[1] == 0, NA, as.character(pop$CLASS[1]))
  }
  ewe.list[[i]] <- data.frame(cbind(years, ewe.age, lamb.status, pn.status, pop.name))
}

#compd.data.subset <- subset(compd.data, select = c("Pop", "year", "CLASS", "CLASS_SUSP", "SumLambSurv", "Recr", "RadEwes", "RadEwesWLambs"), year >= 1995)
#compd.data.subset[1:10, ]

age.spec.ewe.prod.data.orig <- do.call(rbind, ewe.list)
age.spec.ewe.prod.dat.small <- subset(age.spec.ewe.prod.data.orig, is.na(pn.status) == F & is.na(lamb.status) == F)
#table(age.spec.ewe.prod.dat.small$pop.name, age.spec.ewe.prod.dat.small$pn.status, age.spec.ewe.prod.dat.small$years)
repro.array <- table(as.numeric(as.character(age.spec.ewe.prod.dat.small$ewe.age)), age.spec.ewe.prod.dat.small$lamb.status)
age.spec.prod.pn <- subset(age.spec.ewe.prod.data.orig, is.na(pn.status) == F & (pn.status == "LAMBS" | pn.status == "ALL_AGE" | pn.status == "ADULTS"))
repro.array.pn <- rbind(table(as.numeric(as.character(age.spec.prod.pn$ewe.age)), age.spec.prod.pn$lamb.status)[-1, ], c(0, 0), c(0, 0))

#age.spec.prod.he <- subset(age.spec.ewe.prod.data.orig, is.na(pn.status) == F & (pn.status == "ADULTS" | pn.status == "HEALTHY"))
age.spec.prod.he <- subset(age.spec.ewe.prod.data.orig, is.na(pn.status) == F & (pn.status == "HEALTHY"))
repro.array.he <- rbind(table(as.numeric(as.character(age.spec.prod.he$ewe.age)), age.spec.prod.he$lamb.status), c(0, 0))

#-- how many ewe-years in each age-class in each PN state? --#
young.prime.pn <- table(as.numeric(as.character(age.spec.prod.pn$ewe.age)) <= 7, age.spec.prod.pn$lamb.status)
prime.pn <- table(as.numeric(as.character(age.spec.prod.pn$ewe.age)) >= 8 & as.numeric(as.character(age.spec.prod.pn$ewe.age)) <= 12, age.spec.prod.pn$lamb.status)
old.pn <- table(as.numeric(as.character(age.spec.prod.pn$ewe.age)) >= 13, age.spec.prod.pn$lamb.status)

young.prime.he <- table(as.numeric(as.character(age.spec.prod.he$ewe.age)) <= 7, age.spec.prod.he$lamb.status)
prime.he <- table(as.numeric(as.character(age.spec.prod.he$ewe.age)) >= 8 & as.numeric(as.character(age.spec.prod.he$ewe.age)) <= 12, age.spec.prod.he$lamb.status)
old.he <- table(as.numeric(as.character(age.spec.prod.he$ewe.age)) >= 13, age.spec.prod.he$lamb.status)

aso.repro.dat <- subset(age.spec.ewe.prod.dat.small, pop.name == "Asotin")
imn.repro.dat <- subset(age.spec.ewe.prod.dat.small, pop.name == "Imnaha" & years %in% c("1998", "1999", "2000", "2001","2002", "1997"))
bcan.repro.dat <- subset(age.spec.ewe.prod.dat.small, pop.name == "BigCanyon" & years %in% c("1998", "1999"))
muir.repro.dat <- subset(age.spec.ewe.prod.dat.small, pop.name == "MuirCreek" & years %in% c("1998", "1999"))

movi.neg.repro <- as.data.frame(rbind(aso.repro.dat, imn.repro.dat, bcan.repro.dat, muir.repro.dat))
repro.array.premovi <- rbind(c(0, 0), table(movi.neg.repro$ewe.age, movi.neg.repro$lamb.status))
table(movi.neg.repro$ewe.age, movi.neg.repro$pn.status)

movi.young.prime.he <- table(as.numeric(as.character(movi.neg.repro$ewe.age)) <= 7, movi.neg.repro$lamb.status)
movi.prime.he <- table(as.numeric(as.character(movi.neg.repro$ewe.age)) >= 8 & as.numeric(as.character(movi.neg.repro$ewe.age)) <= 12, movi.neg.repro$lamb.status)
movi.old.he <- table(as.numeric(as.character(movi.neg.repro$ewe.age)) >= 13, movi.neg.repro$lamb.status)

#-- calculate number of lambs produced (in aggregate) by each age class --#
n.ages <- dim(repro.array)[1]
r <- r.pn <- r.he <- r.premovi <- rep(NA, n.ages)
for(a in 1:n.ages){
  r[a] <- sum(repro.array[a, ])
  r.pn[a] <- sum(repro.array.pn[a, ])
  r.he[a] <- sum(repro.array.he[a, ])
  r.premovi[a] <- sum(repro.array.premovi[a, ])
}

#-- MODEL --#


sink("sheep.agespecrepro.binom.bug")
cat("
    model{
    
    #----------------#
    #-- Parameters --#
    #----------------#
    #-- wean.<i> = ewes in ith age class who weaned a lamb
    
    #----------------#
    #-- Priors ------#
    #----------------#    
    
    wean.1 ~ dunif(0, 1)
    wean.2 ~ dunif(0, 1)
    wean.3 ~ dunif(0, 1)
    wean.4 ~ dunif(0, 1)
    wean.5 ~ dunif(0, 1)
    
    #------------------#
    #-- Likelihood ----#  
    #------------------#
    for(a in 1:n.ages){
      repro.vec[a] ~ dbinom(lamb.pr[a], r[a])
    }
    
    #-- define cell probabilities --#
    # yearlings
    lamb.pr[1] <-  wean.1     # probability lamb is weaned
    # lamb.pr[1, 3] <-  nolamb.1   # probability no lamb is observed
    # 2-year-olds
    lamb.pr[2] <-  wean.2     # probability lamb is weaned
    # 3-7 year-olds
    for(a in 3:7){
      lamb.pr[a] <-  wean.3     # probability lamb is weaned
    } #a
    # 8 - 13 year-olds
    for(a in 8:13){
      lamb.pr[a] <-  wean.4     # probability lamb is weaned
    } #a  
    # > 13 year-olds
    for(a in 14:19){
      lamb.pr[a] <-  wean.5     # probability lamb is weaned
    } #a  
    }
    ", fill = T)
sink()

#-- pneumonia-years only 
sheep.agespecrepro.data.pn <- list(n.ages = n.ages,
                                   r = r.pn, 
                                   repro.vec = repro.array.pn[ , 1]
)

sheep.agespecrepro.data.he <- list(n.ages = n.ages,
                                   r = r.he, 
                                   repro.vec = repro.array.he[ , 1]
)

# initial values
sheep.agespecrepro.inits <- function(){
  list(
    wean.1 = runif(0, 1),
    wean.2 = runif(0, 1),
    wean.3 = runif(0, 1),
    wean.4 = runif(0, 1),
    wean.5 = runif(0, 1)
  )
}

# parameters monitored
sheep.agespecrepro.parameters <- c(
  "wean.1", "wean.2", "wean.3", "wean.4", "wean.5"
)

# MCMC settings
ni <- 20000
nt <- 3
nb <- 10000
nc <- 3

#-- Pneumonia-year model
# call JAGS from R
sheep.agespecrepro.pn <- jags.model("sheep.agespecrepro.binom.bug",
                                    data = sheep.agespecrepro.data.pn,
                                    inits = sheep.agespecrepro.inits,
                                    n.chains = nc,
                                    n.adapt = nb
)

update(sheep.agespecrepro.pn, ni)

coda.samples.sheep.agespecrepro.pn <- coda.samples(sheep.agespecrepro.pn,
                                                   sheep.agespecrepro.parameters,
                                                   ni)

#-- Healthy-year model
# call JAGS from R
sheep.agespecrepro.he <- jags.model("sheep.agespecrepro.binom.bug",
                                    data = sheep.agespecrepro.data.he,
                                    inits = sheep.agespecrepro.inits,
                                    n.chains = nc,
                                    n.adapt = nb
)

update(sheep.agespecrepro.he, ni)

coda.samples.sheep.agespecrepro.he <- coda.samples(sheep.agespecrepro.he,
                                                   sheep.agespecrepro.parameters,
                                                   ni)

#-- Pre-Movi-year model
# call JAGS from R
sheep.agespecrepro.data.premovi <- list(n.ages = n.ages,
                                        r = r.premovi, 
                                        repro.vec = repro.array.premovi[ , 1]
)

sheep.agespecrepro.premovi <- jags.model("sheep.agespecrepro.binom.bug",
                                         data = sheep.agespecrepro.data.premovi,
                                         inits = sheep.agespecrepro.inits,
                                         n.chains = nc,
                                         n.adapt = nb
)

update(sheep.agespecrepro.premovi, ni)

coda.samples.sheep.agespecrepro.premovi <- coda.samples(sheep.agespecrepro.premovi,
                                                        sheep.agespecrepro.parameters,
                                                        ni)

premovi.repro <- rbind(coda.samples.sheep.agespecrepro.premovi[[1]][1001:2000, ], coda.samples.sheep.agespecrepro.premovi[[2]][1001:2000, ], coda.samples.sheep.agespecrepro.premovi[[3]][1001:2000, ])
he.repro <- rbind(coda.samples.sheep.agespecrepro.he[[1]][1001:2000, ], coda.samples.sheep.agespecrepro.he[[2]][1001:2000, ], coda.samples.sheep.agespecrepro.he[[3]][1001:2000, ])
pn.repro <- rbind(coda.samples.sheep.agespecrepro.pn[[1]][1001:2000, ], coda.samples.sheep.agespecrepro.pn[[2]][1001:2000, ], coda.samples.sheep.agespecrepro.pn[[3]][1001:2000, ])

#write.csv(he.repro, "~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Posteriors/Reproduction/HealthyReproPost_30Sept2014.csv")
#write.csv(pn.repro, "~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Posteriors/Reproduction/PNReproPost_30Sept2014.csv")

#-- summaries of MCMC outputs --#

summary(coda.samples.sheep.agespecrepro.pn)[[2]]
summary(coda.samples.sheep.agespecrepro.he)[[2]]
summary(coda.samples.sheep.agespecrepro.premovi)[[2]]

gelman.diag(coda.samples.sheep.agespecrepro.pn)
gelman.diag(coda.samples.sheep.agespecrepro.he)
gelman.diag(coda.samples.sheep.agespecrepro.premovi)

#-- plot intervals over time --#
par(cex.main = .8, mfrow = c(1, 1))
plot(xlim = c(1.5, 5.5), xaxt = "n", ylim = c(0, 1), x = -1, y = 1, main = "95% credible intervals for age-specific reproduction", xlab = "age class", ylab = "P(weans a lamb)")
for(i in 2:5){ 
  segments(x0 = i, x1 = i, y0 = summary(coda.samples.sheep.agespecrepro.pn)[[2]][i, 1], y1 = summary(coda.samples.sheep.agespecrepro.pn)[[2]][i, 5], col = "red", lwd = 2)
  segments(x0 = i + 0.25, x1 = i + 0.25, y0 = summary(coda.samples.sheep.agespecrepro.he)[[2]][i, 1], y1 = summary(coda.samples.sheep.agespecrepro.he)[[2]][i, 5], col = "grey30", lwd = 2, lty = 2)
  segments(x0 = i + 0.15, x1 = i + 0.15, y0 = summary(coda.samples.sheep.agespecrepro.premovi)[[2]][i, 1], y1 = summary(coda.samples.sheep.agespecrepro.premovi)[[2]][i, 5], col = "grey60", lwd = 2, lty = 3)
}
axis(side = 1, at = 2:5, labels = c("<2.5", "2.5-7", "8-13", ">13"))
legend("topright", c("lamb disease years", "lamb healthy years", "pre-Movi years"), lty = c(1, 2, 3), col = c("red", "grey30", "grey60"), lwd = c(2, 2, 2), bty = "n")



# #-----------------------------------------------#
# #-- 2) Age-specific survival PN vs. Healthy ----#
# #-----------------------------------------------#
# #-- 1) build data. Use Festa-Bianchet 2006 stages: lamb, yearling, 2-7, 8-13, >13 --#
# #-- keep each year cohort distinct from the others --#
# #----- m-array runs from 1995 - 2010 --#
# ewes <- subset(studysheep, SEX == "F" & is.na(Tooth_Age) == F)
# ewe.m.array <- matrix(NA, nrow = length(levels(factor(ewes$ID))), ncol = 2011 - 1995)
# stage.entry <- rep(NA, length(levels(factor(ewes$ID))))
# for(i in 1:dim(ewe.m.array)[1]){
#   k <- subset(ewes, ID == levels(factor(ewes$ID))[i])
#   ewe.m.array[i, 1 : (k$ENTRY_BIOYR - 1994)] <- 0
#   ewe.m.array[i, (k$ENTRY_BIOYR - 1994)] <- k$AENTRY
#   ewe.m.array[i, (min(k$ENTRY_BIOYR - 1994 + 1, dim(ewe.m.array)[2]) : min(k$END_BIOYR - 1994 + 1, dim(ewe.m.array)[2]))] <- k$AENTRY + (1 : length((min(k$ENTRY_BIOYR - 1994 + 1, dim(ewe.m.array)[2]) : min(k$END_BIOYR - 1994 + 1, dim(ewe.m.array)[2]))))
#   ewe.m.array[i, (min(k$END_BIOYR - 1994 + 1, dim(ewe.m.array)[2]) : min(k$END_BIOYR - 1994 + 1, dim(ewe.m.array)[2]))] <- 0
# }
# 
# ewe.m.array.stage.v2 <- ifelse(ewe.m.array <= 1, 1, 
#                                floor(ewe.m.array)
# )
# 
# ewe.m.array.stage.v2[is.na(ewe.m.array.stage.v2) == T] <- 19
# for(i in 1:dim(ewe.m.array.stage.v2)[1]){
#   for(j in 2:dim(ewe.m.array.stage.v2)[2]){
#     ewe.m.array.stage.v2[i, j] <- ifelse(ewe.m.array.stage.v2[i, j - 1] >= 2 & ewe.m.array.stage.v2[i, j] == 1, 19, ewe.m.array.stage.v2[i, j])
#   }
# }
# 
# ewe.m.array.stage.v2[1:10, ]
# table(ewe.m.array.stage.v2)
# 
# # f is a vector containing year of first capture for each sheep
# get.first <- function(x){
#   min(which(x >= 2))
# }
# f.v2 <- apply(ewe.m.array.stage.v2, 1, get.first)
# 
# # states are 1, 2, 3, 4, 5, 6 #
# 
# sink("sheep.agespecsurv.v2.bug")
# cat("
#     model{
#     
#     #----------------#
#     #-- Parameters --#
#     #----------------#
#     #-- alpha.<i> = aging (stage-transition)
#     #-- s.<i> = survival (within stage; no transition)
#     
#     #----------------#
#     #-- Priors ------#
#     #----------------#
#     for(t in 1:(n.years - 1)){
#     s.2[t] <- mean.s2
#     s.3[t] <- mean.s3
#     s.4[t] <- mean.s4
#     s.5[t] <- mean.s5
#     collar.prop[t] <- mean.collar.prop
#     }
#     mean.s2 ~ dunif(0, 1)
#     mean.s3 ~ dunif(0, 1)
#     mean.s4 ~ dunif(0, 1)
#     mean.s5 ~ dunif(0, 1)
#     mean.collar.prop ~ dunif(0, 1)
#     
#     # Define state-transitions and observation matrices
#     for(i in 1:nind){
#     # State Process: Define probabilities of state S(t + 1) | S(t)
#     for(t in f[i] :(n.years - 1)){
#     # read matrix as (left-hand state number into right-hand state number)
#     # Leslie matrix row 1 (uncollared)
#     ps[1, i, t, 1] <- 1 - (collar.prop[t])
#     ps[1, i, t, 2] <- collar.prop[t]
#     ps[1, i, t, 3] <- collar.prop[t]
#     ps[1, i, t, 4] <- collar.prop[t]
#     ps[1, i, t, 5] <- collar.prop[t]
#     ps[1, i, t, 6] <- collar.prop[t]
#     ps[1, i, t, 7] <- collar.prop[t]
#     ps[1, i, t, 8] <- collar.prop[t]
#     ps[1, i, t, 9] <- collar.prop[t]
#     ps[1, i, t, 10] <- collar.prop[t]
#     ps[1, i, t, 11] <- collar.prop[t]
#     ps[1, i, t, 12] <- collar.prop[t]
#     ps[1, i, t, 13] <- collar.prop[t]
#     ps[1, i, t, 14] <- collar.prop[t]
#     ps[1, i, t, 15] <- collar.prop[t]
#     ps[1, i, t, 16] <- collar.prop[t]
#     ps[1, i, t, 17] <- collar.prop[t]
#     ps[1, i, t, 18] <- collar.prop[t]
#     ps[1, i, t, 19] <- 0
#     # Leslie matrix row 2 (yearlings)
#     ps[2, i, t, 1] <- 0
#     ps[2, i, t, 2] <- 0
#     ps[2, i, t, 3] <- s.2[t]
#     ps[2, i, t, 4] <- 0
#     ps[2, i, t, 5] <- 0
#     ps[2, i, t, 6] <- 0
#     ps[2, i, t, 7] <- 0
#     ps[2, i, t, 8] <- 0
#     ps[2, i, t, 9] <- 0
#     ps[2, i, t, 10] <- 0
#     ps[2, i, t, 11] <- 0
#     ps[2, i, t, 12] <- 0
#     ps[2, i, t, 13] <- 0
#     ps[2, i, t, 14] <- 0
#     ps[2, i, t, 15] <- 0
#     ps[2, i, t, 16] <- 0
#     ps[2, i, t, 17] <- 0
#     ps[2, i, t, 18] <- 0
#     ps[2, i, t, 19] <- (1 - s.2[t])
#     # Leslie matrix row 3 (2-year-olds)
#     ps[3, i, t, 1] <- 0
#     ps[3, i, t, 2] <- 0
#     ps[3, i, t, 3] <- 0
#     ps[3, i, t, 4] <- s.3[t]
#     ps[3, i, t, 5] <- 0
#     ps[3, i, t, 6] <- 0
#     ps[3, i, t, 7] <- 0
#     ps[3, i, t, 8] <- 0
#     ps[3, i, t, 9] <- 0
#     ps[3, i, t, 10] <- 0
#     ps[3, i, t, 11] <- 0
#     ps[3, i, t, 12] <- 0
#     ps[3, i, t, 13] <- 0
#     ps[3, i, t, 14] <- 0
#     ps[3, i, t, 15] <- 0
#     ps[3, i, t, 16] <- 0
#     ps[3, i, t, 17] <- 0
#     ps[3, i, t, 18] <- 0
#     ps[3, i, t, 19] <- (1 - s.3[t])
#     # Leslie matrix row 4 (3-year-olds)
#     ps[4, i, t, 1] <- 0
#     ps[4, i, t, 2] <- 0
#     ps[4, i, t, 3] <- 0
#     ps[4, i, t, 4] <- 0
#     ps[4, i, t, 5] <- s.3[t]
#     ps[4, i, t, 6] <- 0
#     ps[4, i, t, 7] <- 0
#     ps[4, i, t, 8] <- 0
#     ps[4, i, t, 9] <- 0
#     ps[4, i, t, 10] <- 0
#     ps[4, i, t, 11] <- 0
#     ps[4, i, t, 12] <- 0
#     ps[4, i, t, 13] <- 0
#     ps[4, i, t, 14] <- 0
#     ps[4, i, t, 15] <- 0
#     ps[4, i, t, 16] <- 0
#     ps[4, i, t, 17] <- 0
#     ps[4, i, t, 18] <- 0
#     ps[4, i, t, 19] <- (1 - s.3[t])
#     # Leslie matrix row 5 (4-year-olds)
#     ps[5, i, t, 1] <- 0
#     ps[5, i, t, 2] <- 0
#     ps[5, i, t, 3] <- 0
#     ps[5, i, t, 4] <- 0
#     ps[5, i, t, 5] <- 0
#     ps[5, i, t, 6] <- s.3[t]
#     ps[5, i, t, 7] <- 0
#     ps[5, i, t, 8] <- 0
#     ps[5, i, t, 9] <- 0
#     ps[5, i, t, 10] <- 0
#     ps[5, i, t, 11] <- 0
#     ps[5, i, t, 12] <- 0
#     ps[5, i, t, 13] <- 0
#     ps[5, i, t, 14] <- 0
#     ps[5, i, t, 15] <- 0
#     ps[5, i, t, 16] <- 0
#     ps[5, i, t, 17] <- 0
#     ps[5, i, t, 18] <- 0
#     ps[5, i, t, 19] <- (1 - s.3[t])
#     # Leslie matrix row 6 (5-year-olds)
#     ps[6, i, t, 1] <- 0
#     ps[6, i, t, 2] <- 0
#     ps[6, i, t, 3] <- 0
#     ps[6, i, t, 4] <- 0
#     ps[6, i, t, 5] <- 0
#     ps[6, i, t, 6] <- 0
#     ps[6, i, t, 7] <- s.3[t]
#     ps[6, i, t, 8] <- 0
#     ps[6, i, t, 9] <- 0
#     ps[6, i, t, 10] <- 0
#     ps[6, i, t, 11] <- 0
#     ps[6, i, t, 12] <- 0
#     ps[6, i, t, 13] <- 0
#     ps[6, i, t, 14] <- 0
#     ps[6, i, t, 15] <- 0
#     ps[6, i, t, 16] <- 0
#     ps[6, i, t, 17] <- 0
#     ps[6, i, t, 18] <- 0
#     ps[6, i, t, 19] <- (1 - s.3[t])
#     # Leslie matrix row 7 (6-year-olds)
#     ps[7, i, t, 1] <- 0
#     ps[7, i, t, 2] <- 0
#     ps[7, i, t, 3] <- 0
#     ps[7, i, t, 4] <- 0
#     ps[7, i, t, 5] <- 0
#     ps[7, i, t, 6] <- 0
#     ps[7, i, t, 7] <- 0
#     ps[7, i, t, 8] <- s.3[t]
#     ps[7, i, t, 9] <- 0
#     ps[7, i, t, 10] <- 0
#     ps[7, i, t, 11] <- 0
#     ps[7, i, t, 12] <- 0
#     ps[7, i, t, 13] <- 0
#     ps[7, i, t, 14] <- 0
#     ps[7, i, t, 15] <- 0
#     ps[7, i, t, 16] <- 0
#     ps[7, i, t, 17] <- 0
#     ps[7, i, t, 18] <- 0
#     ps[7, i, t, 19] <- (1 - s.3[t])
#     # Leslie matrix row 8 (7-year-olds)
#     ps[8, i, t, 1] <- 0
#     ps[8, i, t, 2] <- 0
#     ps[8, i, t, 3] <- 0
#     ps[8, i, t, 4] <- 0
#     ps[8, i, t, 5] <- 0
#     ps[8, i, t, 6] <- 0
#     ps[8, i, t, 7] <- 0
#     ps[8, i, t, 8] <- 0
#     ps[8, i, t, 9] <- s.3[t]
#     ps[8, i, t, 10] <- 0
#     ps[8, i, t, 11] <- 0
#     ps[8, i, t, 12] <- 0
#     ps[8, i, t, 13] <- 0
#     ps[8, i, t, 14] <- 0
#     ps[8, i, t, 15] <- 0
#     ps[8, i, t, 16] <- 0
#     ps[8, i, t, 17] <- 0
#     ps[8, i, t, 18] <- 0
#     ps[8, i, t, 19] <- (1 - s.3[t])
#     # Leslie matrix row 9 (8-year-olds)
#     ps[9, i, t, 1] <- 0
#     ps[9, i, t, 2] <- 0
#     ps[9, i, t, 3] <- 0
#     ps[9, i, t, 4] <- 0
#     ps[9, i, t, 5] <- 0
#     ps[9, i, t, 6] <- 0
#     ps[9, i, t, 7] <- 0
#     ps[9, i, t, 8] <- 0
#     ps[9, i, t, 9] <- 0
#     ps[9, i, t, 10] <- s.4[t]
#     ps[9, i, t, 11] <- 0
#     ps[9, i, t, 12] <- 0
#     ps[9, i, t, 13] <- 0
#     ps[9, i, t, 14] <- 0
#     ps[9, i, t, 15] <- 0
#     ps[9, i, t, 16] <- 0
#     ps[9, i, t, 17] <- 0
#     ps[9, i, t, 18] <- 0
#     ps[9, i, t, 19] <- (1 - s.4[t])
#     # Leslie matrix row 10 (9-year-olds)
#     ps[10, i, t, 1] <- 0
#     ps[10, i, t, 2] <- 0
#     ps[10, i, t, 3] <- 0
#     ps[10, i, t, 4] <- 0
#     ps[10, i, t, 5] <- 0
#     ps[10, i, t, 6] <- 0
#     ps[10, i, t, 7] <- 0
#     ps[10, i, t, 8] <- 0
#     ps[10, i, t, 9] <- 0
#     ps[10, i, t, 10] <- 0
#     ps[10, i, t, 11] <- s.4[t]
#     ps[10, i, t, 12] <- 0
#     ps[10, i, t, 13] <- 0
#     ps[10, i, t, 14] <- 0
#     ps[10, i, t, 15] <- 0
#     ps[10, i, t, 16] <- 0
#     ps[10, i, t, 17] <- 0
#     ps[10, i, t, 18] <- 0
#     ps[10, i, t, 19] <- (1 - s.4[t])
#     # Leslie matrix row 11 (10-year-olds)
#     ps[11, i, t, 1] <- 0
#     ps[11, i, t, 2] <- 0
#     ps[11, i, t, 3] <- 0
#     ps[11, i, t, 4] <- 0
#     ps[11, i, t, 5] <- 0
#     ps[11, i, t, 6] <- 0
#     ps[11, i, t, 7] <- 0
#     ps[11, i, t, 8] <- 0
#     ps[11, i, t, 9] <- 0
#     ps[11, i, t, 10] <- 0
#     ps[11, i, t, 11] <- 0
#     ps[11, i, t, 12] <- s.4[t]
#     ps[11, i, t, 13] <- 0
#     ps[11, i, t, 14] <- 0
#     ps[11, i, t, 15] <- 0
#     ps[11, i, t, 16] <- 0
#     ps[11, i, t, 17] <- 0
#     ps[11, i, t, 18] <- 0
#     ps[11, i, t, 19] <- (1 - s.4[t])
#     # Leslie matrix row 12 (11-year-olds)
#     ps[12, i, t, 1] <- 0
#     ps[12, i, t, 2] <- 0
#     ps[12, i, t, 3] <- 0
#     ps[12, i, t, 4] <- 0
#     ps[12, i, t, 5] <- 0
#     ps[12, i, t, 6] <- 0
#     ps[12, i, t, 7] <- 0
#     ps[12, i, t, 8] <- 0
#     ps[12, i, t, 9] <- 0
#     ps[12, i, t, 10] <- 0
#     ps[12, i, t, 11] <- 0
#     ps[12, i, t, 12] <- 0
#     ps[12, i, t, 13] <- s.4[t]
#     ps[12, i, t, 14] <- 0
#     ps[12, i, t, 15] <- 0
#     ps[12, i, t, 16] <- 0
#     ps[12, i, t, 17] <- 0
#     ps[12, i, t, 18] <- 0
#     ps[12, i, t, 19] <- (1 - s.4[t])
#     # Leslie matrix row 13 (12-year-olds)
#     ps[13, i, t, 1] <- 0
#     ps[13, i, t, 2] <- 0
#     ps[13, i, t, 3] <- 0
#     ps[13, i, t, 4] <- 0
#     ps[13, i, t, 5] <- 0
#     ps[13, i, t, 6] <- 0
#     ps[13, i, t, 7] <- 0
#     ps[13, i, t, 8] <- 0
#     ps[13, i, t, 9] <- 0
#     ps[13, i, t, 10] <- 0
#     ps[13, i, t, 11] <- 0
#     ps[13, i, t, 12] <- 0
#     ps[13, i, t, 13] <- 0
#     ps[13, i, t, 14] <- s.4[t]
#     ps[13, i, t, 15] <- 0
#     ps[13, i, t, 16] <- 0
#     ps[13, i, t, 17] <- 0
#     ps[13, i, t, 18] <- 0
#     ps[13, i, t, 19] <- (1 - s.4[t])
#     # Leslie matrix row 14 (13-year-olds)
#     ps[14, i, t, 1] <- 0
#     ps[14, i, t, 2] <- 0
#     ps[14, i, t, 3] <- 0
#     ps[14, i, t, 4] <- 0
#     ps[14, i, t, 5] <- 0
#     ps[14, i, t, 6] <- 0
#     ps[14, i, t, 7] <- 0
#     ps[14, i, t, 8] <- 0
#     ps[14, i, t, 9] <- 0
#     ps[14, i, t, 10] <- 0
#     ps[14, i, t, 11] <- 0
#     ps[14, i, t, 12] <- 0
#     ps[14, i, t, 13] <- 0
#     ps[14, i, t, 14] <- 0
#     ps[14, i, t, 15] <- s.4[t]
#     ps[14, i, t, 16] <- 0
#     ps[14, i, t, 17] <- 0
#     ps[14, i, t, 18] <- 0 
#     ps[14, i, t, 19] <- (1 - s.4[t])
#     # Leslie matrix row 15 (14-year-olds)
#     ps[15, i, t, 1] <- 0
#     ps[15, i, t, 2] <- 0
#     ps[15, i, t, 3] <- 0
#     ps[15, i, t, 4] <- 0
#     ps[15, i, t, 5] <- 0
#     ps[15, i, t, 6] <- 0
#     ps[15, i, t, 7] <- 0
#     ps[15, i, t, 8] <- 0
#     ps[15, i, t, 9] <- 0
#     ps[15, i, t, 10] <- 0
#     ps[15, i, t, 11] <- 0
#     ps[15, i, t, 12] <- 0
#     ps[15, i, t, 13] <- 0
#     ps[15, i, t, 14] <- 0
#     ps[15, i, t, 15] <- 0
#     ps[15, i, t, 16] <- s.5[t]
#     ps[15, i, t, 17] <- 0
#     ps[15, i, t, 18] <- 0
#     ps[15, i, t, 19] <- (1 - s.5[t])
#     # Leslie matrix row 16 (15-year-olds)
#     ps[16, i, t, 1] <- 0
#     ps[16, i, t, 2] <- 0
#     ps[16, i, t, 3] <- 0
#     ps[16, i, t, 4] <- 0
#     ps[16, i, t, 5] <- 0
#     ps[16, i, t, 6] <- 0
#     ps[16, i, t, 7] <- 0
#     ps[16, i, t, 8] <- 0
#     ps[16, i, t, 9] <- 0
#     ps[16, i, t, 10] <- 0
#     ps[16, i, t, 11] <- 0
#     ps[16, i, t, 12] <- 0
#     ps[16, i, t, 13] <- 0
#     ps[16, i, t, 14] <- 0
#     ps[16, i, t, 15] <- 0
#     ps[16, i, t, 16] <- 0
#     ps[16, i, t, 17] <- s.5[t]
#     ps[16, i, t, 18] <- 0
#     ps[16, i, t, 19] <- (1 - s.5[t])
#     # Leslie matrix row 17 (16-year-olds)
#     ps[17, i, t, 1] <- 0
#     ps[17, i, t, 2] <- 0
#     ps[17, i, t, 3] <- 0
#     ps[17, i, t, 4] <- 0
#     ps[17, i, t, 5] <- 0
#     ps[17, i, t, 6] <- 0
#     ps[17, i, t, 7] <- 0
#     ps[17, i, t, 8] <- 0
#     ps[17, i, t, 9] <- 0
#     ps[17, i, t, 10] <- 0
#     ps[17, i, t, 11] <- 0
#     ps[17, i, t, 12] <- 0
#     ps[17, i, t, 13] <- 0
#     ps[17, i, t, 14] <- 0
#     ps[17, i, t, 15] <- 0
#     ps[17, i, t, 16] <- 0
#     ps[17, i, t, 17] <- 0
#     ps[17, i, t, 18] <- s.5[t]
#     ps[17, i, t, 19] <- (1 - s.5[t])
#     # Leslie matrix row 18 (dead)
#     ps[18, i, t, 1] <- 0
#     ps[18, i, t, 2] <- 0
#     ps[18, i, t, 3] <- 0
#     ps[18, i, t, 4] <- 0
#     ps[18, i, t, 5] <- 0
#     ps[18, i, t, 6] <- 0
#     ps[18, i, t, 7] <- 0
#     ps[18, i, t, 8] <- 0
#     ps[18, i, t, 9] <- 0
#     ps[18, i, t, 10] <- 0
#     ps[18, i, t, 11] <- 0
#     ps[18, i, t, 12] <- 0
#     ps[18, i, t, 13] <- 0
#     ps[18, i, t, 14] <- 0
#     ps[18, i, t, 15] <- 0
#     ps[18, i, t, 16] <- 0
#     ps[18, i, t, 17] <- 0
#     ps[18, i, t, 18] <- 0
#     ps[18, i, t, 19] <- 1
#     # Leslie matrix row 18 (dead)
#     ps[19, i, t, 1] <- 0
#     ps[19, i, t, 2] <- 0
#     ps[19, i, t, 3] <- 0
#     ps[19, i, t, 4] <- 0
#     ps[19, i, t, 5] <- 0
#     ps[19, i, t, 6] <- 0
#     ps[19, i, t, 7] <- 0
#     ps[19, i, t, 8] <- 0
#     ps[19, i, t, 9] <- 0
#     ps[19, i, t, 10] <- 0
#     ps[19, i, t, 11] <- 0
#     ps[19, i, t, 12] <- 0
#     ps[19, i, t, 13] <- 0
#     ps[19, i, t, 14] <- 0
#     ps[19, i, t, 15] <- 0
#     ps[19, i, t, 16] <- 0
#     ps[19, i, t, 17] <- 0
#     ps[19, i, t, 18] <- 0
#     ps[19, i, t, 19] <- 1
#     } #t
#     } #i
#     #------------------#
#     #-- Likelihood ----#
#     #------------------#
#       for(i in 1:nind){
#        for(t in (f[i] + 1): n.years){
#          z[i, t] ~ dcat(ps[z[i, t - 1], i, t - 1, ])
#         } #t
#       } #i
#     }
#     ", fill = T)
# sink()
# 
# # bundle data
# sheep.agespecsurv.data.v2 <- list(z = ewe.m.array.stage.v2,
#                                   #                              y = ewe.m.array.stage,
#                                   f = f.v2,
#                                   n.years = dim(ewe.m.array.stage.v2)[2],
#                                   nind = dim(ewe.m.array.stage.v2)[1]
# )
# # #-- pneumonia-years only 
# # sheep.agespecrepro.data.pn <- list(n.ages = n.ages,
# #                                    r = r.pn, 
# #                                    repro.vec = repro.array.pn[ , 1]
# # )
# 
# 
# # initial values
# sheep.agespecsurv.inits.v2 <- function(){
#   list(
#     mean.s2 = runif(1, 0, 1),
#     mean.s3 = runif(1, 0, 1),
#     mean.s4 = runif(1, 0, 1),
#     mean.s5 = runif(1, 0, 1),
#     mean.collar.prop = runif(1, 0, 1)
#   )
# }
# 
# # parameters monitored
# sheep.agespecsurv.parameters.v2 <- c(
#   "mean.s2",
#   "mean.s3",
#   "mean.s4",
#   "mean.s5",
#   "mean.collar.prop")
# 
# # MCMC settings
# ni <- 2000
# nt <- 3
# nb <- 1000
# nc <- 3
# 
# # call JAGS from R
# sheep.agespecsurv.v2 <- jags.model("sheep.agespecsurv.v2.bug",
#                                    data = sheep.agespecsurv.data.v2,
#                                    inits = sheep.agespecsurv.inits.v2,
#                                    n.chains = nc,
#                                    n.adapt = nb
# )
# 
# update(sheep.agespecsurv.v2, ni)
# 
# coda.samples.sheep.agespecsurv.v2 <- coda.samples(sheep.agespecsurv.v2,
#                                                   sheep.agespecsurv.parameters.v2,
#                                                   ni)
# 
# summary(coda.samples.sheep.agespecsurv.v2)
# gelman.diag(coda.samples.sheep.agespecsurv.v2)
# 
# he.surv.post <- rbind(coda.samples.sheep.agespecsurv.v2[[1]][1001:2000, ], coda.samples.sheep.agespecsurv.v2[[2]][1001:2000, ], coda.samples.sheep.agespecsurv.v2[[3]][1001:2000, ])
  
  
#----------------------------------------------------#
#-- Adult survival PN/Healthy Revised ---------------#
#----------------------------------------------------#
# data format: one row per ewe-year. 
ewes.with.teeth <- subset(studysheep, SEX == "F" & is.na(Tooth_Age) == F)
ewe.years <- rep(NA, dim(ewes.with.teeth)[1])
ewe.dat <- vector("list", dim(ewes.with.teeth)[1])
for(i in 1:dim(ewes.with.teeth)[1]){
  ewe.years[i] <- ewes.with.teeth$END_BIOYR[i] - ewes.with.teeth$ENTRY_BIOYR[i] + 1
  ewe.dat[[i]] <- matrix(NA, nrow = ewe.years[i], ncol = 5)
  loop.years <- ewe.years[i]
  for(j in 1:loop.years){
    year.status <- subset(compd.data, as.character(Pop) == as.character(ewes.with.teeth$Population)[i] & year == (ewes.with.teeth$ENTRY_BIOYR[i] + (j - 1)))
    ewe.dat[[i]][j, 1] <- as.character(ewes.with.teeth$ID[i]) 
    ewe.dat[[i]][j, 2] <- seq(ewes.with.teeth$ENTRY_BIOYR[i], ewes.with.teeth$END_BIOYR[i])[j]
    ewe.dat[[i]][j, 3] <- ewes.with.teeth$AENTRY[i] + (j - 1)
    ewe.dat[[i]][j, 4] <- as.character(year.status$CLASS)[1]
    ewe.dat[[i]][j, 5] <- ifelse(j != loop.years, 0, ifelse(ewes.with.teeth$DEAD == 1, 1, NA))
  }
}  

ewe.rowperyear.dat <- as.data.frame(do.call("rbind", ewe.dat))
names(ewe.rowperyear.dat) <- c("ID", "Year", "EstAge", "PNClass", "Dead")
ewe.rowperyear.dat$EstAge <- as.numeric(as.character(ewe.rowperyear.dat$EstAge))
ewe.rowperyear.dat$AgeClass <- ifelse(ewe.rowperyear.dat$EstAge <= 2, 2,  ifelse(ewe.rowperyear.dat$EstAge > 2 & ewe.rowperyear.dat$EstAge <= 7, 3, ifelse(ewe.rowperyear.dat$EstAge > 7 & ewe.rowperyear.dat$EstAge <= 13, 4, 5)))

#-- calculate number of ewe-years survived (in aggregate) by each age class in each disease state --#
pn.eweyrs <- subset(ewe.rowperyear.dat, PNClass %in% c("ADULTS", "ALL_AGE", "LAMBS"))
he.eweyrs <- subset(ewe.rowperyear.dat, !(PNClass %in% c("ADULTS", "ALL_AGE", "LAMBS")))

surv.array <- table(round(ewe.rowperyear.dat$EstAge), ewe.rowperyear.dat$Dead)
pn.surv.array <- rbind(c(0, 0), c(0, 0), table(round(pn.eweyrs$EstAge), pn.eweyrs$Dead), c(0, 0))
he.surv.array <- table(round(he.eweyrs$EstAge), he.eweyrs$Dead)

n.eweages <- 19
s <- s.pn <- s.he <- s.premovi <- rep(NA, n.eweages)
for(a in 1:n.eweages){
  s[a] <- sum(surv.array[a, ])
  s.pn[a] <- sum(pn.surv.array[a, ])
  s.he[a] <- sum(he.surv.array[a, ])
#  r.premovi[a] <- sum(repro.array.premovi[a, ])
}

#-- MODEL --#
sink("agespecsurv.binom.bug")
cat("
    model{
    
    #----------------#
    #-- Parameters --#
    #----------------#
    #-- wean.<i> = ewes in ith age class who weaned a lamb
    
    #----------------#
    #-- Priors ------#
    #----------------#    
    
    surv.1 ~ dunif(0, 1)
    surv.2 ~ dunif(0, 1)
    surv.3 ~ dunif(0, 1)
    surv.4 ~ dunif(0, 1)
#    wean.5 ~ dunif(0, 1)
    surv.5 ~ dunif(0, 1)
    
    #------------------#
    #-- Likelihood ----#  
    #------------------#
    for(a in 1:n.eweages){
#    repro.vec[a] ~ dbinom(lamb.pr[a], r[a])
      surv.vec[a] ~ dbinom(ewe.pr[a], s[a])
    }
    
    #-- define cell probabilities --#
    # yearlings
#    lamb.pr[1] <-  wean.1     # probability lamb is weaned
    ewe.pr[1] <-  surv.1   # probability no lamb is observed
    # 2-year-olds
#    lamb.pr[2] <-  wean.2     # probability lamb is weaned
    ewe.pr[2] <-  surv.2     # probability lamb is weaned
    # 3-7 year-olds
    for(a in 3:7){
#    lamb.pr[a] <-  wean.3     # probability lamb is weaned
    ewe.pr[a] <-  surv.3     # probability lamb is weaned
    } #a
    # 8 - 13 year-olds
    for(a in 8:13){
#    lamb.pr[a] <-  wean.4     # probability lamb is weaned
    ewe.pr[a] <-  surv.4     # probability lamb is weaned
    } #a  
    # > 13 year-olds
    for(a in 14:19){
#    lamb.pr[a] <-  wean.5     # probability lamb is weaned
    ewe.pr[a] <-  surv.5     # probability lamb is weaned
    } #a  
    }
    ", fill = T)
sink()

#-- pneumonia-years only 
agespecsurv.data.pn <- list(n.eweages = n.eweages,
                                   s = s.pn, 
                                   surv.vec = pn.surv.array[ , 1]
)

agespecsurv.data.he <- list(n.eweages = n.eweages,
                                   s = s.he, 
                                   surv.vec = he.surv.array[ , 1]
)

# initial values
agespecsurv.inits <- function(){
  list(
    surv.1 = runif(0, 1),
    surv.2 = runif(0, 1),
    surv.3 = runif(0, 1),
    surv.4 = runif(0, 1),
    surv.5 = runif(0, 1)
  )
}

# parameters monitored
agespecsurv.parameters <- c(
  "surv.1", "surv.2", "surv.3", "surv.4", "surv.5"
)

# MCMC settings
ni <- 20000
nt <- 3
nb <- 10000
nc <- 3

#-- Pneumonia-year model
# call JAGS from R
agespecsurv.pn <- jags.model("agespecsurv.binom.bug",
                                    data = agespecsurv.data.pn,
                                    inits = agespecsurv.inits,
                                    n.chains = nc,
                                    n.adapt = nb
)

update(agespecsurv.pn, ni)

coda.samples.agespecsurv.pn <- coda.samples(agespecsurv.pn,
                                                   agespecsurv.parameters,
                                                   ni)

#-- Healthy-year model
# call JAGS from R
agespecsurv.he <- jags.model("agespecsurv.binom.bug",
                                    data = agespecsurv.data.he,
                                    inits = agespecsurv.inits,
                                    n.chains = nc,
                                    n.adapt = nb
)

update(agespecsurv.he, ni)

coda.samples.agespecsurv.he <- coda.samples(agespecsurv.he,
                                                   agespecsurv.parameters,
                                                   ni)

#coda.samples.agespecsurv.he <- read.csv("~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Posteriors/Survival/HealthyReproPost_30Sept2014.csv")
#coda.samples.agespecsurv.pn <- read.csv("~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Posteriors/Survival/PNReproPost_30Sept2014.csv")

#-- plot intervals over time --#
par(cex.main = .8, mfrow = c(1, 1))
plot(xlim = c(1.5, 5.5), xaxt = "n", ylim = c(0, 1), x = -1, y = 1, main = "95% credible intervals for age-specific survival", xlab = "age class", ylab = "P(survives)")
for(i in 2:5){ 
  segments(x0 = i, x1 = i, y0 = summary(coda.samples.agespecsurv.pn)[[2]][i, 1], y1 = summary(coda.samples.agespecsurv.pn)[[2]][i, 5], col = "red", lwd = 2)
  segments(x0 = i + 0.25, x1 = i + 0.25, y0 = summary(coda.samples.agespecsurv.he)[[2]][i, 1], y1 = summary(coda.samples.agespecsurv.he)[[2]][i, 5], col = "grey30", lwd = 2, lty = 2)
#  segments(x0 = i + 0.15, x1 = i + 0.15, y0 = summary(coda.samples.sheep.agespecrepro.premovi)[[2]][i, 1], y1 = summary(coda.samples.agespecsurv.premovi)[[2]][i, 5], col = "grey60", lwd = 2, lty = 3)
}
axis(side = 1, at = 2:5, labels = c("<2.5", "2.5-7", "8-13", ">13"))
legend("topright", c("disease years", "healthy years"), lty = c(1, 2), col = c("red", "grey30"), lwd = c(2, 2), bty = "n")

#premovi.surv <- rbind(coda.samples.sheep.agespecrepro.premovi[[1]][1001:2000, ], coda.samples.sheep.agespecrepro.premovi[[2]][1001:2000, ], coda.samples.sheep.agespecrepro.premovi[[3]][1001:2000, ])
he.surv <- rbind(coda.samples.agespecsurv.he[[1]][1:10000, ], coda.samples.agespecsurv.he[[2]][1:10000, ], coda.samples.agespecsurv.he[[3]][1:10000, ])
pn.surv <- rbind(coda.samples.agespecsurv.pn[[1]][1:10000, ], coda.samples.agespecsurv.pn[[2]][1:10000, ], coda.samples.agespecsurv.pn[[3]][1:10000, ])
#write.csv(he.surv, "~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Posteriors/Survival/HealthyReproPost_30Sept2014.csv")
#write.csv(pn.surv, "~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Posteriors/Survival/PNReproPost_30Sept2014.csv")

#-------------------------------------------------------------#
#-- Recruitment in pneumonia and healthy years ---------------#
#-------------------------------------------------------------#
compd.data.recr <- subset(compd.data, select = c("RadEwes", "RadEwesWLambs", "SumLambSurv", "Recr"), year >= 1995)
pn.yrs <- subset(compd.data, !(CLASS %in% c("HEALTHY")) & year >= 1995)
# SumLambSurv references survival only for those ewes observed at least once with a lamb
# 1) calculate radiocollared ewes who reproduced.
pnyear.prop.repro <- pn.yrs$RadEwesWLambs / pn.yrs$RadEwes
# 2) multiply the proportion who reproduced by SLS to get estimated end-of-summer ewe:lamb ratio
pnyear.endofsummer.ewelambrat <- pn.yrs$SumLambSurv * pnyear.prop.repro
# 3) calculate multiplicative change between end of summer ewe-lamb ratio and recr
pnyear.sls.recr.ratio <- 1 - (pn.yrs$SumLambSurv - pn.yrs$Recr)

he.yrs <- subset(compd.data, CLASS %in% c("HEALTHY") & year >= 1995)
# 1) calculate radiocollared ewes who reproduced.
heyear.prop.repro <- he.yrs$RadEwesWLambs / he.yrs$RadEwes
# 2) multiply the proportion who reproduced by SLS to get estimated end-of-summer ewe:lamb ratio
heyear.endofsummer.ewelambrat <- he.yrs$SumLambSurv * heyear.prop.repro
# 3) calculate multiplicative change between end of summer ewe-lamb ratio and recr
heyear.sls.recr.ratio <- 1 - (he.yrs$SumLambSurv - he.yrs$Recr)
  
par(mfrow = c(2, 1))
hist(pnyear.sls.recr.ratio, col = "grey80", xlim = c(0, 4), main = "Pneumonia years for lambs", xlab = "Estimate of total overwinter lamb mortality")
hist(heyear.sls.recr.ratio, col = "grey80", xlim = c(0, 4), main = "Healthy years (healthy + adult-only)", xlab = "Estimate of total overwinter lamb mortality")




#-------------------------------------------------------------------------#
#-- 3) Population trajectory simulations for diseased and healthy years --#
#-------------------------------------------------------------------------#
#-- basic sim structure: 
#----------------------- 1) build state vector with current pop age structure
#----------------------- 2) draw environmental state based on alpha and gamma
#----------------------- 3) draw Leslie matrix paramters from posteriors associated with environmental state
#----------------------- 4) project population forward and record new age structure


#-- generate inital age-structure based on limiting distribution for healthy-year states
#-- (for now, make it up and bias it low) --#
# projection function

he.surv <- read.csv("~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Posteriors/Survival/HealthyReproPost_30Sept2014.csv")
pn.surv <- read.csv("~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Posteriors/Survival/PNReproPost_30Sept2014.csv")
he.repro <- read.csv("~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Posteriors/Reproduction/HealthyReproPost_30Sept2014.csv")
pn.repro <- read.csv("~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Posteriors/Reproduction/PNReproPost_30Sept2014.csv")

sp.repro.post <- pn.repro
sp.surv.post <- he.surv * runif(1, .85, 1)
inf.surv.post <- pn.surv
he.surv.post <- he.surv
he.repro.post <- he.repro
pn.repro.post <- pn.repro

pn.recr <- na.omit(pnyear.sls.recr.ratio) * .5
he.recr <- na.omit(pnyear.sls.recr.ratio) * .5 # multiply by .5 to eliminate ram-lambs at recruitment.

repros.med <- c(0, rep(median(he.repro.post[1:3000, 2]), 2), rep(median(he.repro.post[1:3000, 3]), 4), rep(median(he.repro.post[1:3000, 4]), 6), rep(median(he.repro.post[1:3000, 5]), 5))    
#    survs <- c(he.repro.post[sample(1:3000, 1), 1], rep(he.surv.post[sample(1:3000, 1), 2], 2), rep(he.surv.post[sample(1:3000, 1), 3], 4), rep(he.surv.post[sample(1:3000, 1), 4], 6), rep(he.surv.post[sample(1:3000, 1), 5], 4))
survs.med <- c(median(he.recr[(1:length(he.recr))]), rep(median(he.surv.post[1:3000, 2]), 2), rep(median(he.surv.post[1:3000, 3]), 4), rep(median(he.surv.post[1:3000, 4]), 6), rep(median(he.surv.post[1:3000, 5]), 4))
leslie.init <- rbind(repros.med, cbind(diag(c(survs.med)), rep(0, length(survs.med))))
eigen(leslie.init)$vectors[, 1]

stable.age.structure <- abs(eigen(leslie.init)$vectors[, 1]) / sum(abs(eigen(leslie.init)$vectors[, 1]))

N.init <- 100
ages.init <- floor(N.init * stable.age.structure)
#ages.init.orig <- c(100, 50, 40, 40, 40, 30, 30, 30, 30, 30, 30, 20, 20, 20, 20, 10, 10, 10)  
#ages.init <- round(c(sum(ages.init.orig[2:18]) / 4 * .4, ages.init.orig[2:18] / 4) )
# build vector to keep track of states
current.state <- rep(NA, 30)
current.state[1] <- "healthy"

# function to update environmental state
update.status.fun <- function(alpha, gamma, current.state){
  current.state.new <- rep(NA, 1)
  if(current.state == "healthy"){
    gets.infected <- rbinom(1, 1, alpha)
    current.state.new[1] <- ifelse(gets.infected == 0, "healthy", "spillover")
  }
  else if(current.state == "spillover"){
    current.state.new[1] <- "infected"
    } else if(current.state == "infected"){
      fade.out <- rbinom(1, 1, gamma)
      current.state.new[1] <- ifelse(fade.out == 1, "healthy", "infected")
    }
  return(list(current.state.new = current.state.new))
}

# function to update Leslie matrix parameters
  #  need to expand Leslie to be 18x18... also, need individuals to age....
update.leslie.fun <- function(current.state, he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr){
  leslie.out <- rep(NA, 6, 6)
  if(current.state == "healthy"){
    repros <- c(0, rep(he.repro.post[sample(1:3000, 1), 2], 2), rep(he.repro.post[sample(1:3000, 1), 3], 4), rep(he.repro.post[sample(1:3000, 1), 4], 6), rep(he.repro.post[sample(1:3000, 1), 5], 5))    
#    survs <- c(he.repro.post[sample(1:3000, 1), 1], rep(he.surv.post[sample(1:3000, 1), 2], 2), rep(he.surv.post[sample(1:3000, 1), 3], 4), rep(he.surv.post[sample(1:3000, 1), 4], 6), rep(he.surv.post[sample(1:3000, 1), 5], 4))
    survs <- c(he.recr[sample(1:length(he.recr), 1)], rep(he.surv.post[sample(1:3000, 1), 2], 2), rep(he.surv.post[sample(1:3000, 1), 3], 4), rep(he.surv.post[sample(1:3000, 1), 4], 6), rep(he.surv.post[sample(1:3000, 1), 5], 4))
    leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs))))
  }
  else {
    repros <- c(0, rep(inf.repro.post[sample(1:3000, 1), 2], 2), rep(inf.repro.post[sample(1:3000, 1), 3], 4), rep(inf.repro.post[sample(1:3000, 1), 4], 6), rep(inf.repro.post[sample(1:3000, 1), 5], 5))    
#    survs <- c(inf.repro.post[sample(1:3000, 1), 1], rep(inf.surv.post[sample(1:3000, 1), 2], 2), rep(inf.surv.post[sample(1:3000, 1), 3], 4), rep(inf.surv.post[sample(1:3000, 1), 4], 6), rep(inf.surv.post[sample(1:3000, 1), 5], 4))
    survs <- c(pn.recr[sample(1:length(pn.recr), 1)], rep(inf.surv.post[sample(1:3000, 1), 2], 2), rep(inf.surv.post[sample(1:3000, 1), 3], 4), rep(inf.surv.post[sample(1:3000, 1), 4], 6), rep(inf.surv.post[sample(1:3000, 1), 5], 4))
    leslie <- rbind(repros, cbind(diag(c(survs)), rep(0, length(survs))))
  }
  return(leslie)
}


project.fun <- function(timesteps, ages.init, alpha, gamma, he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr){
  N <- matrix(NA, nrow = length(ages.init), ncol = timesteps)
  N[, 1] <- ages.init
  tot.pop.size <- log.lambda.s <- rep(NA, length = timesteps)
  disease.status <- rep(NA, timesteps)
#  disease.status[1] <- "healthy"
  disease.status[1:11] <- c(rep("healthy", 10), "spillover")
  for(i in 1:10){
    new.leslie <- update.leslie.fun(current.state = disease.status[i], he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
    N[, i + 1] <- round(t(N[ , i]) %*% new.leslie) 
    tot.pop.size[i] <- sum(N[ , i])
    if(i == 1){
      log.lambda.s[i] <- NA
    } else {
      log.lambda.s[i] <- ifelse(tot.pop.size[i] == 0, NA, log(tot.pop.size[i] / tot.pop.size[i - 1]))
    }
  }
  for(i in 11:(timesteps - 1)){
    new.leslie <- update.leslie.fun(current.state = disease.status[i], he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
    N[, i + 1] <- round(t(N[ , i]) %*% new.leslie)
    disease.status[i + 1] <- update.status.fun(alpha, gamma, current.state = disease.status[i])$current.state.new[1]
    tot.pop.size[i] <- sum(N[ , i])
    if(i == 1){
      log.lambda.s[i] <- NA
    } else {
      log.lambda.s[i] <- ifelse(tot.pop.size[i] == 0, NA, log(tot.pop.size[i] / tot.pop.size[i - 1]))
    }
  }
#  tot.pop.size <- apply(N, 2, sum)
#  out.list <- list(N = N, disease.status = disease.status, tot.pop.size = tot.pop.size)
  out.list <- list(N = N, disease.status = disease.status, tot.pop.size = tot.pop.size, log.lambda.s = log.lambda.s)
  return(out.list)
}

project.fun.out <- project.fun(timesteps = 20, ages.init = ages.init, alpha = .01, gamma = 1, he.repro = he.repro.post, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)

#--------------------------------------------------------#
#-- Model check: population trajectory with no disease --#
#--------------------------------------------------------#
he.project.fun <- function(timesteps, ages.init, alpha, gamma, he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr){
  N <- matrix(NA, nrow = length(ages.init), ncol = timesteps)
  N[, 1] <- ages.init
  tot.pop.size <- log.lambda.s <- rep(NA, length = timesteps)
  disease.status <- rep("healthy", timesteps)
  #  disease.status[1] <- "healthy"
#  disease.status[1:11] <- c(rep("healthy", 10), "spillover")
  for(i in 1:(timesteps - 1)){
    new.leslie <- update.leslie.fun(current.state = disease.status[i], he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
    N[, i + 1] <- t(N[ , i]) %*% new.leslie  
    tot.pop.size[i] <- sum(N[ , i])
    if(i == 1){
      log.lambda.s[i] <- NA
    } else {
      log.lambda.s[i] <- log(tot.pop.size[i] / tot.pop.size[i - 1])
    }
  }
#   for(i in 11:(timesteps - 1)){
#     new.leslie <- update.leslie.fun(current.state = disease.status[i], he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
#     N[, i + 1] <- t(N[ , i]) %*% new.leslie
#     disease.status[i + 1] <- update.status.fun(alpha, gamma, current.state = disease.status[i])$current.state.new[1]
#   }
#  tot.pop.size <- apply(N, 2, sum)
  out.list <- list(N = N, disease.status = disease.status, tot.pop.size = tot.pop.size, log.lambda.s = log.lambda.s)
  return(out.list)
}

sp.repro.post <- pn.repro
sp.surv.post <- he.surv * runif(1, .3, 1)
inf.surv.post <- pn.surv
he.surv.post <- he.surv

he.project.fun.out <- he.project.fun(timesteps = 20, ages.init = ages.init, alpha = .01, gamma = 1, he.repro = he.repro.post, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)

timesteps <- 60
reps <- 100
popsize.he <- log.lambda.s.he <- matrix(NA, ncol = timesteps, nrow = reps)

for(i in 1:reps){
  he.project <- he.project.fun(timesteps = timesteps, ages.init = ages.init, alpha = .01, gamma = .1, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
  popsize.he[i, ] <- he.project$tot.pop.size 
  log.lambda.s.he[i, ] <- he.project$log.lambda.s
}  

#write.csv(popsize.he, "~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Simulations/Healthy/HealthyPopsize_30Sept2014.csv", row.names = F)
#popsize.he <- as.matrix(read.csv("~/work/Kezia/Research/EcologyPapers/RecruitmentVsAdultSurv/Data/Simulations/Healthy/HealthyPopsize_30Sept2014.csv"))

he.quants <- which(popsize.he[, 40] %in% as.numeric(quantile(popsize.he[, 40], c(0.025, 0.975), type = 3)))
he.med <- which(popsize.he[, 40] %in% as.numeric(quantile(popsize.he[, 40], .5, type = 3)))

#par(mfrow = c(2, 1))
layout(matrix(c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3), nrow = 2, byrow = T))
plot(popsize.he[1, -c(1)] ~ seq(2, timesteps), type = "l", ylim = c(0, 3000), xlab = "year", ylab = "population size")
for(i in 2:reps){
  lines(popsize.he[i, -c(1)] ~ seq(2, timesteps), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
for(j in 1:3){
  lines(popsize.he[he.quants[j], -c(1, 2)] ~ seq(3, timesteps), type = "l", col = "black", lwd = 2)  
}
lines(popsize.he[he.med, -c(1, 2)] ~ seq(3, timesteps), type = "l", col = "red", lwd = 2)

plot(log.lambda.s.he[1, -c(1, 2)] ~ seq(3, timesteps), type = "l", ylim = c(-.5, .5), xlab = "year", ylab = expression(paste("log(", lambda, "s)", sep = "")))
for(i in 2:reps){
  lines(log.lambda.s.he[i, -c(1, 2)] ~ seq(3, timesteps), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
abline( h = 0, lty = 2, col = "red", lwd = 2)
boxplot(as.vector(log.lambda.s.he[, -c(1, 2)]), col = "grey80", ylim = c(-.5, .5), ylab = expression(paste("log(", lambda, "s)", sep = "")))
abline(h = 0, lty = 2, col = "red", lwd = 2)

#-----------------------------------------------------#
#-- Environmental model 1: Movi presence definition --#
#-----------------------------------------------------#

#------------------------------------------------------------------#
#-- project population sizes in loop over gammas of .1, .5, .9 ----#
#------------------------------------------------------------------#
popgrowth.sim.fun <- function(timesteps, reps, alpha.range, gamma.range, alpha.steps, gamma.steps){
  alphas <- 1 / seq(min(alpha.range), max(alpha.range), length.out = alpha.steps)  
  gammas <- 1 / seq(min(gamma.range), max(gamma.range), length.out = gamma.steps)
  alpha.gamma.frame <- expand.grid(alphas, gammas)
  popsize.ij <- loglambda.ij <- vector("list", dim(alpha.gamma.frame)[1])
  mean.lnlambda <- rep(NA, dim(alpha.gamma.frame)[1])
  for(i in 1:length(popsize.ij)){
    popsize.ij[[i]] <- loglambda.ij[[i]] <- matrix(NA, nrow = timesteps, ncol = reps)
    for(j in 1:reps){
      popgrowthsim.ij <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = alpha.gamma.frame[i, 1], gamma = alpha.gamma.frame[i, 2], he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
      popsize.ij[[i]][, j] <- popgrowthsim.ij$tot.pop.size
      loglambda.ij[[i]][, j] <- popgrowthsim.ij$log.lambda.s
    }
    mean.lnlambda[i] <- mean(na.omit(unlist(as.vector(loglambda.ij[[i]][-c(11, 12), ]))))
  }
  outlist <- list(alpha.gamma.frame = alpha.gamma.frame, popsize.ij = popsize.ij, loglambda.ij = loglambda.ij, mean.lnlambda = mean.lnlambda)
  return(outlist)
}

timesteps <- 30
reps <- 100
alpha.range <- c(1, 100)
gamma.range <- c(1, 100)
alpha.steps <- 10
gamma.steps <- 10
popgrowth.test <- popgrowth.sim.fun(timesteps, reps, alpha.range, gamma.range, alpha.steps, gamma.steps)
popgrowth.test <- outlist
x.vals <- 1 / popgrowth.test$alpha.gamma.frame[, 1]
y.vals <- 1 / popgrowth.test$alpha.gamma.frame[, 2]

require(grDevices)
color.ramp <- colorRampPalette(c("red", "grey20", "blue"))
cols <- color.ramp(100)
color.val <- cols[round(100 * (popgrowth.test$mean.lnlambda - min(popgrowth.test$mean.lnlambda)) / (abs(min(popgrowth.test$mean.lnlambda )) - (min(popgrowth.test$mean.lnlambda)))) ]
#color.val2 <- color.val / max(abs(color.val))
pt.size <- 10 * (popgrowth.test$mean.lnlambda - min(popgrowth.test$mean.lnlambda)) / (abs(min(popgrowth.test$mean.lnlambda )) - (min(popgrowth.test$mean.lnlambda)))
#color.in <- paste("grey", round(100 * color.val), sep = "")

par(mfrow = c(1, 1))
plot(y.vals  ~ x.vals, cex = pt.size, xlab = expression(paste("expected waiting time until reintroduction (", alpha, ")", sep = "")), ylab = expression(paste("waiting time until expected fade-out (", gamma, ")", sep = "")), col = color.val, pch = 16)
axis(side = 1, at = as.numeric(levels(factor(popgrowth.test$alpha.gamma.frame[, 1]))), labels = round(1 / as.numeric(levels(factor(popgrowth.test$alpha.gamma.frame[, 1]))), 2))
axis(side = 2, at = as.numeric(levels(factor(popgrowth.test$alpha.gamma.frame[, 2]))), labels = round(1 / as.numeric(levels(factor(popgrowth.test$alpha.gamma.frame[, 2]))), 2))

timesteps <- 60
reps <- 1000
popsize.30.gamma.05 <- matrix(NA, ncol = timesteps, nrow = reps)
popsize.30.gamma.1 <- matrix(NA, ncol = timesteps, nrow = reps)
popsize.30.gamma.2 <- matrix(NA, ncol = timesteps, nrow = reps)
popsize.30.gamma.5 <- matrix(NA, ncol = timesteps, nrow = reps)
popsize.30.gamma.9 <- matrix(NA, ncol = timesteps, nrow = reps)
popsize.30.gamma1 <- matrix(NA, ncol = timesteps, nrow = reps)

loglambda.30.gamma.05 <- matrix(NA, ncol = timesteps, nrow = reps)
loglambda.30.gamma.1 <- matrix(NA, ncol = timesteps, nrow = reps)
loglambda.30.gamma.2 <- matrix(NA, ncol = timesteps, nrow = reps)
loglambda.30.gamma.5 <- matrix(NA, ncol = timesteps, nrow = reps)
loglambda.30.gamma.9 <- matrix(NA, ncol = timesteps, nrow = reps)
loglambda.30.gamma1 <- matrix(NA, ncol = timesteps, nrow = reps)

for(i in 1:reps){
  project.30.gamma.05 <-  project.fun(timesteps = timesteps, ages.init = ages.init, alpha = 0, gamma = .05, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
  project.30.gamma.1 <-  project.fun(timesteps = timesteps, ages.init = ages.init, alpha = 0, gamma = .1, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
  project.30.gamma.2 <-  project.fun(timesteps = timesteps, ages.init = ages.init, alpha = 0, gamma = .2, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
  project.30.gamma.5 <-  project.fun(timesteps = timesteps, ages.init = ages.init, alpha = 0, gamma = .5, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
  project.30.gamma.9 <-  project.fun(timesteps = timesteps, ages.init = ages.init, alpha = 0, gamma = .9, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
  project.30.gamma1 <-  project.fun(timesteps = timesteps, ages.init = ages.init, alpha = 0, gamma = 1, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
  
  popsize.30.gamma.05[i, ] <- project.30.gamma.05$tot.pop.size
  popsize.30.gamma.1[i, ] <-  project.30.gamma.1$tot.pop.size
  popsize.30.gamma.2[i, ] <-  project.30.gamma.2$tot.pop.size
  popsize.30.gamma.5[i, ] <-  project.30.gamma.5$tot.pop.size
  popsize.30.gamma.9[i, ] <-  project.30.gamma.9$tot.pop.size
  popsize.30.gamma1[i, ]  <-  project.30.gamma1$tot.pop.size
  
  loglambda.30.gamma.05[i, ] <- project.30.gamma.05$log.lambda.s
  loglambda.30.gamma.1[i, ] <-  project.30.gamma.1$log.lambda.s
  loglambda.30.gamma.2[i, ] <-  project.30.gamma.2$log.lambda.s
  loglambda.30.gamma.5[i, ] <-  project.30.gamma.5$log.lambda.s
  loglambda.30.gamma.9[i, ] <-  project.30.gamma.9$log.lambda.s
  loglambda.30.gamma1[i, ]  <-  project.30.gamma1$log.lambda.s
}  

#subset(popsize.30.gamma.1, (popsize.30.gamma.1[, 30]) %in% c(floor(quantile(popsize.30.gamma.1[, 30], 0.025)), ceiling(quantile(popsize.30.gamma.1[, 30], 0.025))))
#                                           , 0.975))))

addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

# extract 2.5th and 97.5th quantiles --#
gamma.05.quants <- which(popsize.30.gamma.05[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.05[, timesteps], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))
gamma.1.quants <- which(popsize.30.gamma.1[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.1[, timesteps], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))
gamma.2.quants <- which(popsize.30.gamma.2[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.2[, timesteps], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))
gamma.5.quants <- which(popsize.30.gamma.5[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.5[, timesteps], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))
gamma.9.quants <- which(popsize.30.gamma.9[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.9[, timesteps], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))
gamma1.quants <- which(popsize.30.gamma1[, timesteps] %in% as.numeric(quantile(popsize.30.gamma1[, timesteps], c(0.025, 0.25, 0.75, 0.975), type = 3, na.rm = T)))

gamma.05.meds <- which(popsize.30.gamma.05[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.05[, timesteps], .5, type = 3, na.rm = T)))
gamma.1.meds <- which(popsize.30.gamma.1[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.1[, timesteps], .5, type = 3, na.rm = T)))
gamma.2.meds <- which(popsize.30.gamma.2[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.2[, timesteps], .5, type = 3, na.rm = T)))
gamma.5.meds <- which(popsize.30.gamma.5[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.5[, timesteps], .5, type = 3, na.rm = T)))
gamma.9.meds <- which(popsize.30.gamma.9[, timesteps] %in% as.numeric(quantile(popsize.30.gamma.9[, timesteps], .5, type = 3, na.rm = T)))
gamma1.meds <- which(popsize.30.gamma1[, timesteps] %in% as.numeric(quantile(popsize.30.gamma1[, timesteps], .5, type = 3, na.rm = T)))

plot.reps <- 100
par(mfrow = c(2, 5))
plot(popsize.30.gamma.05[1, 11:59] ~ seq(11:59), type = "l", ylim = c(0, 1500), xlab = "year", ylab = "population size", main = "20 years", bty = "n")
for(i in 2:plot.reps){
  lines(popsize.30.gamma.05[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
# for(j in 1:4){
#   lines(popsize.30.gamma.05[gamma.05.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(popsize.30.gamma.05[gamma.05.meds, ] ~ seq(11:59), type = "l", col = "red", lwd = 2)  

plot(popsize.30.gamma.1[1, 11:59] ~ seq(11:59), type = "l", ylim = c(0, 1500), xlab = "year", ylab = "", main = "10 years", bty = "n")
for(i in 2:plot.reps){
  lines(popsize.30.gamma.1[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
# for(j in 1:4){
#   lines(popsize.30.gamma.1[gamma.1.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(popsize.30.gamma.1[gamma.1.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2)  

plot(popsize.30.gamma.2[1, 11:59] ~ seq(11:59), type = "l", ylim = c(0, 1500), xlab = "year", ylab = "", main = "5 years", bty = "n")
for(i in 2:plot.reps){
  lines(popsize.30.gamma.2[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
# for(j in 1:4){
#   lines(popsize.30.gamma.2[gamma.2.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(popsize.30.gamma.2[gamma.2.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2) 

plot(popsize.30.gamma1[1, 11:59] ~ seq(11:59), type = "l", ylim = c(0, 1500), xlab = "year", ylab = "", main = "1 year", bty = "n")
for(i in 2:plot.reps){
  lines(popsize.30.gamma1[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
# for(j in 1:4){
#   lines(popsize.30.gamma1[gamma1.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(popsize.30.gamma1[gamma1.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2) 

plot(popsize.he[1, 11:59] ~ seq(11:59), type = "l", ylim = c(0, 1000), xlab = "year", ylab = "", main = "No disease", bty = "n")
for(i in 2:plot.reps){
  lines(popsize.he[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
# for(j in 1:3){
#   lines(popsize.he[he.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(popsize.he[he.med, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2)  

# log-lambda row
plot(loglambda.30.gamma.05[1, 11:59] ~ seq(11:59), type = "l", ylim = c(-1, 1), xlab = "year", ylab = "population size", main = "20 years", bty = "n")
for(i in 2:plot.reps){
  lines(loglambda.30.gamma.05[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
abline(h = 0, col = "red", lwd = 2)
# for(j in 1:4){
#   lines(loglambda.30.gamma.05[gamma.05.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(loglambda.30.gamma.05[gamma.05.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2)  

plot(loglambda.30.gamma.1[1, 11:59] ~ seq(11:59), type = "l", ylim = c(-1, 1), xlab = "year", ylab = "", main = "10 years", bty = "n")
for(i in 2:plot.reps){
  lines(loglambda.30.gamma.1[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
abline(h = 0, col = "red", lwd = 2)
# for(j in 1:4){
#   lines(loglambda.30.gamma.1[gamma.1.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(loglambda.30.gamma.1[gamma.1.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2)  

plot(loglambda.30.gamma.2[1, 11:59] ~ seq(11:59), type = "l", ylim = c(-1, 1), xlab = "year", ylab = "", main = "5 years", bty = "n")
for(i in 2:plot.reps){
  lines(loglambda.30.gamma.2[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
abline(h = 0, col = "red", lwd = 2)
# for(j in 1:4){
#   lines(loglambda.30.gamma.2[gamma.2.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(loglambda.30.gamma.2[gamma.2.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2) 

plot(loglambda.30.gamma1[1, 11:59] ~ seq(11:59), type = "l", ylim = c(-1, 1), xlab = "year", ylab = "", main = "1 year", bty = "n")
for(i in 2:plot.reps){
  lines(loglambda.30.gamma1[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
abline(h = 0, col = "red", lwd = 2)
# for(j in 1:4){
#   lines(loglambda.30.gamma1[gamma1.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(loglambda.30.gamma1[gamma1.meds, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2) 

plot( log.lambda.s.he[1, 11:59] ~ seq(11:59), type = "l", ylim = c(-1, 1), xlab = "year", ylab = "", main = "No disease", bty = "n")
for(i in 2:plot.reps){
  lines( log.lambda.s.he[i, 11:59] ~ seq(11:59), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
abline(h = 0, col = "red", lwd = 2)
# for(j in 1:3){
#   lines(loglambda.he[he.quants[j], 11:59] ~ seq(11:59), type = "l", col = "black", lwd = 2)  
# }
#lines(loglambda.he[he.med, 11:59] ~ seq(11:59), type = "l", col = "red", lwd = 2)  

par(mfrow = c(1, 6))
boxplot(na.omit(as.vector(loglambda.30.gamma.05)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)
boxplot(na.omit(as.vector(loglambda.30.gamma.1)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)
boxplot(na.omit(as.vector(loglambda.30.gamma.2)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)
boxplot(na.omit(as.vector(loglambda.30.gamma.5)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)
boxplot(na.omit(as.vector(loglambda.30.gamma.9)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)
boxplot(na.omit(as.vector(loglambda.30.gamma1)), ylim = c(-1.5, 1.5))
abline(h = 0, col = "red", lwd = 2)

#-------------------------------------------------------#
#-- costs of altering alpha/gamma on 30-year pop size --#
#-------------------------------------------------------#
#-- better to think about gammas as expected time to fade-out.... --#
#-- F(x) = 1/beta * (exp (-x / beta))
#-- where beta = expected survival of system
#-- so, gammas of interest are 1, 1/2, 1/3, 1/4

timesteps <- 40
reps <- 1000
yearstofadeout.seq <- seq(1, 15)
gamma.seq <- 1 / yearstofadeout.seq
popsize.30.gamma.mat <- matrix(NA, reps, length(gamma.seq))
popsize.30.gamma.alpha.33.mat <- matrix(NA, reps, length(gamma.seq))

for(i in 1:reps){
  for(j in 1:length(gamma.seq)){
    popsize.30.gamma.mat[i, j] <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = .1, gamma = gamma.seq[j], he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)$tot.pop.size[timesteps]
    popsize.30.gamma.alpha.33.mat[i, j] <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = .33, gamma = gamma.seq[j], he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)$tot.pop.size[timesteps]
  }
}

popsize.30.quants <- matrix(NA, length(gamma.seq), 5)
popsize.30.alpha.33.quants <- matrix(NA, length(gamma.seq), 5)
for(j in 1:length(gamma.seq)){
  popsize.30.quants[j, ] <- quantile(popsize.30.gamma.mat[, j], c(0.025, 0.25, 0.5, 0.75, 0.975))
  popsize.30.alpha.33.quants[j, ] <- quantile(popsize.30.gamma.alpha.33.mat[, j], c(0.025, 0.25, 0.5, 0.75, 0.975))
}

par(mfrow = c(1, 1))
plot(x = 0, y = 0, xlim = c(0, 16), ylim = c(0, 200), xlab = "Expected years to fade-out (alpha = 0.1)", ylab = "Simulated Pop Size After 30 years")
for(i in 1:length(yearstofadeout.seq)){
  segments(x0 = yearstofadeout.seq[i], x1 = yearstofadeout.seq[i], y0 = popsize.30.quants[i, 1], y1 = popsize.30.quants[i, 5])
  segments(x0 = yearstofadeout.seq[i] - 0.1, x1 = yearstofadeout.seq[i] + 0.1, y0 = popsize.30.quants[i, 3], y1 = popsize.30.quants[i, 3])
  segments(x0 = yearstofadeout.seq[i] + .2, x1 = yearstofadeout.seq[i] + .2, y0 = popsize.30.alpha.33.quants[i, 1], y1 = popsize.30.alpha.33.quants[i, 5], col = "red")
  segments(x0 = yearstofadeout.seq[i] + .1, x1 = yearstofadeout.seq[i] + .3, y0 = popsize.30.alpha.33.quants[i, 3], y1 = popsize.30.alpha.33.quants[i, 3], col = "red")
}

#-- Same thing, but over levels of alpha, with gamma fixed at e(fade-out) = 3 and e(fade-out) = 10
timesteps <- 40
reps <- 1000
yearstoreintro<- seq(1, 15)
alpha.seq <- 1 / yearstoreintro
popsize.30.alpha.gamma.33.mat <- matrix(NA, reps, length(alpha.seq))
popsize.30.alpha.gamma.1.mat <- matrix(NA, reps, length(alpha.seq))

for(i in 1:reps){
  for(j in 1:length(alpha.seq)){
    popsize.30.alpha.gamma.33.mat[i, j] <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = alpha.seq[j], gamma = .33, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)$tot.pop.size[timesteps]
    popsize.30.alpha.gamma.1.mat[i, j] <- project.fun(timesteps = timesteps, ages.init = ages.init, alpha = alpha.seq[j], gamma = .1, he.repro.post = he.repro, sp.repro.post = sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post = sp.surv.post, inf.surv.post = inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)$tot.pop.size[timesteps]
  }
}

popsize.alpha.gamma.33.30.quants <- matrix(NA, length(alpha.seq), 5)
popsize.alpha.gamma.1.30.quants <- matrix(NA, length(alpha.seq), 5)
for(j in 1:length(gamma.seq)){
  popsize.alpha.gamma.33.30.quants[j, ] <- quantile(popsize.30.alpha.gamma.33.mat[, j], c(0.025, 0.25, 0.5, 0.75, 0.975))
  popsize.alpha.gamma.1.30.quants[j, ] <- quantile(popsize.30.alpha.gamma.1.mat[, j], c(0.025, 0.25, 0.5, 0.75, 0.975))
}

par(mfrow = c(1, 1))
plot(x = 0, y = 0, xlim = c(0, 16), ylim = c(0, 100), xlab = "Expected years to reintroduction", ylab = "Simulated Pop Size After 30 years")
for(i in 1:length(yearstoreintro)){
  segments(x0 = yearstoreintro[i], x1 = yearstoreintro[i], y0 = popsize.alpha.gamma.33.30.quants[i, 1], y1 = popsize.alpha.gamma.33.30.quants[i, 5])
  segments(x0 = yearstoreintro[i] - 0.1, x1 = yearstoreintro[i] + 0.1, y0 = popsize.alpha.gamma.33.30.quants[i, 3], y1 = popsize.alpha.gamma.33.30.quants[i, 3])
  segments(x0 = yearstoreintro[i] + .2, x1 = yearstoreintro[i] + .2, y0 = popsize.alpha.gamma.1.30.quants[i, 1], y1 = popsize.alpha.gamma.1.30.quants[i, 5], col = "red")
  segments(x0 = yearstoreintro[i] + .1, x1 = yearstoreintro[i] + .3, y0 = popsize.alpha.gamma.1.30.quants[i, 3], y1 = popsize.alpha.gamma.1.30.quants[i, 3], col = "red")
}

#--------------------------------------------------------------------------#
#-- Exploration of lambda and age structure in each environmental state ---#
#--------------------------------------------------------------------------#
healthy.leslie.list <- spillover.leslie.list <- infected.leslie.list <- vector("list", length = 1000)
healthy.eigenval1 <- spillover.eigenval1 <- infected.eigenval1 <- rep(NA, 1000)
juvsurv.elast.he <- juvsurv.elast.inf <- adsurv.elast.he <- adsurv.elast.inf <- fecund.elast.he <- fecund.elast.inf <- rep(NA, 1000)
age.struct.he <- age.struct.inf <- matrix(NA, nrow = 18, ncol = 1000)

for(i in 1:1000){
  healthy.leslie.list[[i]] <- update.leslie.fun(current.state = "healthy", he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
  spillover.leslie.list[[i]] <- update.leslie.fun(current.state = "spillover", he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
  infected.leslie.list[[i]] <- update.leslie.fun(current.state = "infected", he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
  healthy.eigenval1[i] <- Re(eigen(healthy.leslie.list[[i]])$values[1]) # strip off only real part of eigenvalue
  spillover.eigenval1[i] <- Re(eigen(spillover.leslie.list[[i]])$values[1])
  infected.eigenval1[i] <- Re(eigen(infected.leslie.list[[i]])$values[1])
  he.eigen.rescale <- sum(Re(eigen(healthy.leslie.list[[i]])$vectors[, 1]))
  inf.eigen.rescale <- sum(Re(eigen(infected.leslie.list[[i]])$vectors[, 1]))
  age.struct.he[, i] <- Re(eigen(healthy.leslie.list[[i]])$vectors[, 1]) / he.eigen.rescale
  age.struct.inf[, i] <- Re(eigen(infected.leslie.list[[i]])$vectors[, 1]) / inf.eigen.rescale
  
  he.sens <- sensitivity(healthy.leslie.list[[i]])
  he.elast <- (1 / Re(eigen(healthy.leslie.list[[i]])$values[1])) * he.sens * healthy.leslie.list[[i]]
#  round(he.sens, 2)
#  round(he.elast, 2)
  fecund.elast.he[i] <- sum(he.elast[1, ])
  juvsurv.elast.he[i] <- sum(he.elast[2, 1], he.elast[3, 2])
  adsurv.elast.he[i] <- sum(he.elast[4, 3], he.elast[5, 4], he.elast[6, 5], he.elast[7, 6], he.elast[8, 7], he.elast[9, 8], he.elast[10, 9], he.elast[11, 10], he.elast[12, 11], he.elast[13, 12], he.elast[14, 13], he.elast[15, 14], he.elast[16, 15], he.elast[17, 16], he.elast[18, 17])
  
  inf.sens <- sensitivity(infected.leslie.list[[i]])
  inf.elast <- (1 / Re(eigen(infected.leslie.list[[i]])$values[1])) * inf.sens * infected.leslie.list[[i]]
#  round(inf.sens, 2)
#  round(inf.elast, 2)
  fecund.elast.inf[i] <- sum(inf.elast[1, ])
  juvsurv.elast.inf[i] <- sum(inf.elast[2, 1], inf.elast[3, 2])
  #sum(inf.elast[5, 4], inf.elast[6, 5], inf.elast[7, 6], inf.elast[8, 7], inf.elast[9, 8], inf.elast[10, 9])
  adsurv.elast.inf[i] <- sum(inf.elast[4, 3], inf.elast[5, 4], inf.elast[6, 5], inf.elast[7, 6], inf.elast[8, 7], inf.elast[9, 8], inf.elast[10, 9], inf.elast[11, 10], inf.elast[12, 11], inf.elast[13, 12], inf.elast[14, 13], inf.elast[15, 14], inf.elast[16, 15], inf.elast[17, 16], inf.elast[18, 17])
}

# 95% intervals on lambda 
healthy.cred<- quantile(healthy.eigenval1, c(0.25, 0.5, 0.75))
spillover.cred<- quantile(spillover.eigenval1, c(0.25, 0.5, 0.75))
infected.cred<- quantile(infected.eigenval1, c(0.25, 0.5, 0.75))

# 95% intervals for fecundity, juvenile survival, and adult survival in each environment (healthy and infected)
fecund.he <- quantile(fecund.elast.he, c(0.25, 0.5, 0.75))
fecund.inf <- quantile(fecund.elast.inf, c(.25, 0.5, 0.75))
juvsurv.he <- quantile(juvsurv.elast.he, c(0.25, 0.5, 0.75))
juvsurv.inf <- quantile(juvsurv.elast.inf, c(0.25, 0.5, 0.75))
adsurv.he <- quantile(adsurv.elast.he, c(0.25, 0.5, 0.75))
adsurv.inf <- quantile(adsurv.elast.inf, c(0.25, 0.5, 0.75))

par(mfrow = c(1, 1), mar = c(4, 8, 2, 2), las = 1, cex.lab = 1.0)
plot(c(1, 1) ~ c(fecund.he[1], fecund.he[3]), lty = 1, xlim = c(0, 1), ylim = c(0, 7), type = "l", xlab = expression(paste("Elasticity of ", lambda, " to rate", sep = "")), ylab = "", yaxt = "n", lwd = 2)
lines(c(2, 2) ~ c(fecund.inf[1], fecund.inf[3]), lty = 2, col = "red", lwd = 2)
lines(c(3, 3) ~ c(juvsurv.he[1], juvsurv.he[3]), lty = 1, col = "black", lwd = 2)
lines(c(4, 4) ~ c(juvsurv.inf[1], juvsurv.inf[3]), lty = 2, col = "red", lwd = 2)
lines(c(5, 5) ~ c(adsurv.he[1], adsurv.he[3]), lty = 1, col = "black", lwd = 2)
lines(c(6, 6) ~ c(adsurv.inf[1], adsurv.inf[3]), lty = 2, col = "red", lwd = 2)
lines(c(0.75, 1.25) ~ c(fecund.he[2], fecund.he[2]), lty = 1, col = "black", lwd = 2)
lines(c(1.75, 2.25) ~ c(fecund.inf[2], fecund.inf[2]), lty = 2, col = "red", lwd = 2)
lines(c(2.75, 3.25) ~ c(juvsurv.he[2], juvsurv.he[2]), lty = 1, col = "black", lwd = 2)
lines(c(3.75, 4.25) ~ c(juvsurv.inf[2], juvsurv.inf[2]), lty = 2, col = "red", lwd = 2)
lines(c(4.75, 5.25) ~ c(adsurv.he[2], adsurv.he[2]), lty = 1, col = "black", lwd = 2)
lines(c(5.75, 6.25) ~ c(adsurv.inf[2], adsurv.inf[2]), lty = 2, col = "red", lwd = 2)
axis(side = 2, at = c(1:6), cex.axis = 1.0, labels = c("Fecundity (he)", "Fecundity (pers)", "Juvenile survival (he)", "Juvenile surv (pers)", "Adult survival (he)", "Adult survival (pers)"))

par(mfrow = c(1, 3), cex.axis = 1.2, cex.lab = 1.5)
hist(healthy.eigenval1, xlim = c(min(min(healthy.eigenval1), min(spillover.eigenval1), min(infected.eigenval1)), max(max(healthy.eigenval1), max(spillover.eigenval1), max(infected.eigenval1))), breaks = 15, xlab = expression(lambda), main = "Healthy", col = "grey80")
abline(v = 1, col = "red", lwd = 3)
hist(spillover.eigenval1, xlim = c(min(min(healthy.eigenval1), min(spillover.eigenval1), min(infected.eigenval1)), max(max(healthy.eigenval1), max(spillover.eigenval1), max(infected.eigenval1))), breaks = 15, ylab = "", xlab = expression(lambda), main = "Introduction", col = "grey80")
abline(v = 1, col = "red", lwd = 3)
hist(infected.eigenval1, xlim = c(min(min(healthy.eigenval1), min(spillover.eigenval1), min(infected.eigenval1)), max(max(healthy.eigenval1), max(spillover.eigenval1), max(infected.eigenval1))), breaks = 15, ylab = "", xlab = expression(lambda), main = "Persistence", col = "grey80")
abline(v = 1, col = "red", lwd = 3)

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
age.bounds.he <- age.bounds.inf <- matrix(NA, nrow = 18, ncol = 3)
for(i in 1:18){
  age.bounds.he[i, ] <- quantile(abs(age.struct.he[i, ]), c(0.025, 0.5, 0.975))
  age.bounds.inf[i, ] <- quantile(abs(age.struct.inf[i, ]), c(0.025, 0.5, 0.975))
}

par(mfrow = c(1, 1))
plot(age.bounds.he[1, c(1, 3)] ~ c(1, 1), type = "l", xlim = c(0.5, 18.5), ylim = c(0, .5), lwd = 2, xlab = "age (years)", ylab = "Expected proportion of pop")
segments(x0 = 1 + .3, x1 = 1 + .3, y0 = age.bounds.inf[1, 1], y1 = age.bounds.inf[1, 2], lwd = 2, col = "red")
for(i in 2: dim(age.bounds.he)[1]){
  segments(x0 = i, x1 = i, y0 = age.bounds.he[i, 1], y1 = age.bounds.he[i, 2], lwd = 2)
  segments(x0 = i + .3, x1 = i + .3, y0 = age.bounds.inf[i, 1], y1 = age.bounds.inf[i, 2], lwd = 2, col = "red")
}
leg.text <- c("healthy", "persistently infected")
legend("topright", leg.text, col = c("black", "red"), lwd = c(2, 2), bty = "n")

#--------------------------------------------------------#
#-- Sensitivies and Elasticities using Vec-permutation --#
#--------------------------------------------------------#

vec.permut.fun <- function(reps, alpha, gamma){
  fade.out.elast <- intro.elast <- rep(NA, reps)
  for(i in 1:reps){
# P is the vec-permutation matrix.
P <- matrix(NA, nrow = 18 * 3, ncol = 18 * 3)
for(j in 1:18){
  odd.vec <- c(rep(0, j - 1), 1, rep(0, 18 - j))
    P[j * 3 - 2 , ] <- c(odd.vec, rep(0, 18), rep(0, 18))
    P[j * 3 - 1, ] <- c(rep(0, 18), odd.vec, rep(0, 18))
    P[j * 3, ] <- c(rep(0, 18), rep(0, 18), odd.vec)
}

# B is block diagonal, with 3 17x17 blocks for the 3 environmental states.
healthy.leslie <- update.leslie.fun(current.state = "healthy", he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
spillover.leslie <- update.leslie.fun(current.state = "spillover", he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
endemic.leslie <- update.leslie.fun(current.state = "infected", he.repro.post = he.repro, sp.repro.post, inf.repro.post = pn.repro, he.surv.post = he.surv.post, sp.surv.post, inf.surv.post, he.recr = he.recr, pn.recr = pn.recr)
leslie.list <- list(healthy.leslie, spillover.leslie, endemic.leslie)
B <- bdiag(leslie.list)

# M is block diagonal with 17 3x3 blocks for the 17 demographic states
#alpha <- 1/10
#gamma <- 1/20
small.M <- rbind(c(1 - alpha, alpha, 0), c(gamma, 0, 1 - gamma), c(gamma, 0, 1 - gamma))
M <- bdiag(list(small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M, small.M))

# A is population projection matrix with environmental stochasticity
A <- t(P) %*% M %*% P %*% B
A_eigens <- eigen(A) # complex????
S_a <- sensitivity(A)

# Sensitivity (environmental transitions)
# assume disease status is updated before demography.
S_m <- P %*% t(B) %*% S_a %*% t(P)
round(S_m, 2)[1:10, 1:10]
#E_m <- (1 / A_eigens$value[1]) * M * S_m # regular * because it's a Hadamard production
E_m <- (1 / .92) * M * S_m # regular * because it's a Hadamard production
round(E_m, 2)[1:10, 1:10]
# compare elasticities of fade-out to elasticity of reintroduction
fade.out.elast[i] <- sum(E_m[3, 3], E_m[6, 6], E_m[9, 9], E_m[12, 12], E_m[15, 15], E_m[18, 18], E_m[21, 21], E_m[24, 24], E_m[27, 27], E_m[30, 30], E_m[33, 33], E_m[36, 36], E_m[39, 39], E_m[42, 42], E_m[45, 45], E_m[48, 48], E_m[51, 51])
intro.elast[i] <- sum(E_m[1, 1], E_m[4, 4], E_m[7, 7], E_m[10, 10], E_m[13, 13], E_m[16, 16], E_m[19, 19], E_m[22, 22], E_m[25, 25], E_m[28, 28], E_m[31, 31], E_m[34, 34], E_m[37, 37], E_m[40, 40], E_m[43, 43], E_m[46, 46], E_m[49, 49])
}
return(list(fade.out.elast = fade.out.elast, intro.elast = intro.elast))
}

vec.permut.test1 <- vec.permut.fun(reps = 1000, alpha = 0.2, gamma = 0.1)
vec.permut.test2 <- vec.permut.fun(reps = 1000, alpha = 0.1, gamma = 0.1)
vec.permut.test3 <- vec.permut.fun(reps = 1000, alpha = 0.05, gamma = 0.1)

fade.out.elasts <- quantile(vec.permut.test$fade.out.elast, c(0.025, 0.975))
intro.elasts <- quantile(vec.permut.test$intro.elast, c(0.025, 0.975))

par(mfrow = c(2, 3), cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.4)
hist(vec.permut.test3$fade.out.elast, ylab = "frequency", xlab = "fade-out elasticity", col = "grey80", main = "E(time to reintro) = 20yrs")
hist(vec.permut.test2$fade.out.elast, ylab = "frequency", xlab = "fade-out elasticity", col = "grey80", main = "E(time to reintro) = 10yrs")
hist(vec.permut.test1$fade.out.elast, ylab = "frequency", xlab = "fade-out elasticity", col = "grey80", main = "E(time to reintro) = 5yrs")
hist(vec.permut.test3$intro.elast, main = "", ylab = "frequency", xlab = "introduction elasticity", col = "grey80")
hist(vec.permut.test2$intro.elast, main = "", ylab = "frequency", xlab = "introduction elasticity", col = "grey80")
hist(vec.permut.test1$intro.elast, main = "", ylab = "frequency", xlab = "introduction elasticity", col = "grey80")

# Sensitivity (demographic transitions)
S_b <- S_a %*% t(P) %*% t(M) %*% P
round(S_b, 2)[1:20, 1:20]
E_b <- (1 / .92) * B * S_b # regular * because it's a Hadamard production
round(E_b, 2)

#-------------------------------------------------------#
#-- Figure 2: Leslie Matrix posterior param estimates --#
#-------------------------------------------------------#

#-- Summer Lamb Survival --#
par(cex.main = .8, mfrow = c(1, 2), cex.main = 1.0, las = 1)
plot(xlim = c(1.5, 5.5), xaxt = "n", ylim = c(0, 1), x = -1, y = 1, main = "Age-specific survival", xlab = "age class", ylab = "P(survives)")
for(i in 2:5){ 
  segments(x0 = i, x1 = i, y0 = summary(coda.samples.agespecsurv.pn)[[2]][i, 1], y1 = summary(coda.samples.agespecsurv.pn)[[2]][i, 5], col = "red", lwd = 2)
  segments(x0 = i + 0.25, x1 = i + 0.25, y0 = summary(coda.samples.agespecsurv.he)[[2]][i, 1], y1 = summary(coda.samples.agespecsurv.he)[[2]][i, 5], col = "grey30", lwd = 2, lty = 2)
  #  segments(x0 = i + 0.15, x1 = i + 0.15, y0 = summary(coda.samples.sheep.agespecrepro.premovi)[[2]][i, 1], y1 = summary(coda.samples.agespecsurv.premovi)[[2]][i, 5], col = "grey60", lwd = 2, lty = 3)
}
axis(side = 1, at = 2:5, labels = c("<2.5", "2.5-7", "8-13", ">13"))
legend("bottomright", cex = .8, c("persistent years", "healthy years"), lty = c(1, 2), col = c("red", "grey30"), lwd = c(2, 2), bty = "n")

plot(xlim = c(1.5, 5.5), xaxt = "n", ylim = c(0, 1), x = -1, y = 1, main = "Age-specific reproduction", xlab = "age class", ylab = "P(weans a lamb)")
for(i in 2:5){ 
  segments(x0 = i, x1 = i, y0 = summary(coda.samples.sheep.agespecrepro.pn)[[2]][i, 1], y1 = summary(coda.samples.sheep.agespecrepro.pn)[[2]][i, 5], col = "red", lwd = 2)
  segments(x0 = i + 0.25, x1 = i + 0.25, y0 = summary(coda.samples.sheep.agespecrepro.he)[[2]][i, 1], y1 = summary(coda.samples.sheep.agespecrepro.he)[[2]][i, 5], col = "grey30", lwd = 2, lty = 2)
# segments(x0 = i + 0.15, x1 = i + 0.15, y0 = summary(coda.samples.sheep.agespecrepro.premovi)[[2]][i, 1], y1 = summary(coda.samples.sheep.agespecrepro.premovi)[[2]][i, 5], col = "grey60", lwd = 2, lty = 3)
}
axis(side = 1, at = 2:5, labels = c("<2.5", "2.5-7", "8-13", ">13"))
#legend("topright", c("lamb disease years", "lamb healthy years"), lty = c(1, 2), col = c("red", "grey30"), lwd = c(2, 2), bty = "n")

#-- Age-specific ewe survival --#


#-----------------------------------------------------#
#-- Figure 3:Simulations with various fadeout times --# 
#-----------------------------------------------------#


par(mfrow = c(1, 5), cex.lab = 1.4)
plot(popsize.30.gamma.05[1, ] ~ seq(1, timesteps), type = "l", xaxt = "n", ylim = c(0, 1500), xlab = "year", ylab = "population size", main = "20 years", bty = "n")
for(i in 2:reps){
  lines(popsize.30.gamma.05[i, ] ~ seq(1, timesteps), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
for(j in 1:4){
  lines(popsize.30.gamma.05[gamma.05.quants[j], ] ~ seq(1, timesteps), type = "l", col = "black", lwd = 2)  
}
lines(popsize.30.gamma.05[gamma.05.meds, ] ~ seq(1, timesteps), type = "l", col = "red", lwd = 2)  
axis(side = 1, at = (seq(1:7) * 10 - 10), labels = c("-10","0", "10","20", "30", "40", "50"))

plot(popsize.30.gamma.1[1, ] ~ seq(1, timesteps), type = "l", xaxt = "n", ylim = c(0, 1500), xlab = "year", ylab = "", main = "10 years", bty = "n")
for(i in 2:reps){
  lines(popsize.30.gamma.1[i, ] ~ seq(1, timesteps), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
for(j in 1:4){
  lines(popsize.30.gamma.1[gamma.1.quants[j], ] ~ seq(1, timesteps), type = "l", col = "black", lwd = 2)  
}
lines(popsize.30.gamma.1[gamma.1.meds, ] ~ seq(1, timesteps), type = "l", col = "red", lwd = 2)  
axis(side = 1, at = (seq(1:7) * 10 - 10), labels = c("-10","0", "10","20", "30", "40", "50"))

plot(popsize.30.gamma.2[1, ] ~ seq(1, timesteps), type = "l", xaxt = "n", ylim = c(0, 1500), xlab = "year", ylab = "", main = "5 years", bty = "n")
for(i in 2:reps){
  lines(popsize.30.gamma.2[i, ] ~ seq(1, timesteps), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
for(j in 1:4){
  lines(popsize.30.gamma.2[gamma.2.quants[j], ] ~ seq(1, timesteps), type = "l", col = "black", lwd = 2)  
}
lines(popsize.30.gamma.2[gamma.2.meds, ] ~ seq(1, timesteps), type = "l", col = "red", lwd = 2) 
axis(side = 1, at = (seq(1:7) * 10 - 10), labels = c("-10","0", "10","20", "30", "40", "50"))

plot(popsize.30.gamma1[1, ] ~ seq(1, timesteps), type = "l", xaxt = "n", ylim = c(0, 1500), xlab = "year", ylab = "", main = "1 year", bty = "n")
for(i in 2:reps){
  lines(popsize.30.gamma1[i, ] ~ seq(1, timesteps), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
for(j in 1:4){
  lines(popsize.30.gamma1[gamma1.quants[j], ] ~ seq(1, timesteps), type = "l", col = "black", lwd = 2)  
}
lines(popsize.30.gamma1[gamma1.meds, ] ~ seq(1, timesteps), type = "l", col = "red", lwd = 2) 
axis(side = 1, at = (seq(1:7) * 10 - 10), labels = c("-10","0", "10","20", "30", "40", "50"))

plot(popsize.he[1, ] ~ seq(1, timesteps), type = "l", xaxt = "n", ylim = c(0, 1000), xlab = "year", ylab = "", main = "No disease", bty = "n")
for(i in 2:reps){
  lines(popsize.he[i, ] ~ seq(1, timesteps), type = "l", col = rgb(.35, .35, .35, alpha = .25))
}
for(j in 1:3){
  lines(popsize.he[he.quants[j], ] ~ seq(1, timesteps), type = "l", col = "black", lwd = 2)  
}
lines(popsize.he[he.med, ] ~ seq(1, timesteps), type = "l", col = "red", lwd = 2)  
axis(side = 1, at = (seq(1:7) * 10 - 10), labels = c("-10","0", "10","20", "30", "40", "50"))

#--------------------------------------------------#
#-- Figure 4: pop size after 30 years -------------#
#--------------------------------------------------#
par(mfrow = c(1, 2))
plot(x = 0, y = 0, xlim = c(0, 16), ylim = c(0, 2000), xlab = "Expected years to fade-out", ylab = "Simulated Pop Size After 30 years")
abline(h = 1000, col = "grey80", lty = 2)
abline(h = 1500, col = "grey80", lty = 2)
abline(h = 500, col = "grey80", lty = 2)
for(i in 1:length(yearstofadeout.seq)){
  segments(x0 = yearstofadeout.seq[i], x1 = yearstofadeout.seq[i], y0 = popsize.30.quants[i, 1], y1 = popsize.30.quants[i, 5])
  segments(x0 = yearstofadeout.seq[i] - 0.1, x1 = yearstofadeout.seq[i] + 0.1, y0 = popsize.30.quants[i, 3], y1 = popsize.30.quants[i, 3])
  segments(x0 = yearstofadeout.seq[i] + .2, x1 = yearstofadeout.seq[i] + .2, y0 = popsize.30.alpha.33.quants[i, 1], y1 = popsize.30.alpha.33.quants[i, 5], col = "red", lwd = 2)
  segments(x0 = yearstofadeout.seq[i] + .1, x1 = yearstofadeout.seq[i] + .3, y0 = popsize.30.alpha.33.quants[i, 3], y1 = popsize.30.alpha.33.quants[i, 3], col = "red", lwd = 2)
}
leg.text1 <- c("Expect 10 years between introductions", "Expect 3 years between introductions")
legend("topright", bty = "n", leg.text1, lty = c(1, 1), col = c("black", "red"), cex = .6)

plot(x = 0, y = 0, xlim = c(0, 16), ylim = c(0, 2000), xlab = "Expected years to reintroduction", ylab = "Simulated Pop Size After 30 years")
abline(h = 1000, col = "grey80", lty = 2)
abline(h = 1500, col = "grey80", lty = 2)
abline(h = 500, col = "grey80", lty = 2)
for(i in 1:length(yearstoreintro)){
  segments(x0 = yearstoreintro[i], x1 = yearstoreintro[i], y0 = popsize.alpha.gamma.33.30.quants[i, 1], y1 = popsize.alpha.gamma.33.30.quants[i, 5])
  segments(x0 = yearstoreintro[i] - 0.1, x1 = yearstoreintro[i] + 0.1, y0 = popsize.alpha.gamma.33.30.quants[i, 3], y1 = popsize.alpha.gamma.33.30.quants[i, 3])
  segments(x0 = yearstoreintro[i] + .2, x1 = yearstoreintro[i] + .2, y0 = popsize.alpha.gamma.1.30.quants[i, 1], y1 = popsize.alpha.gamma.1.30.quants[i, 5], col = "red", lwd = 2)
  segments(x0 = yearstoreintro[i] + .1, x1 = yearstoreintro[i] + .3, y0 = popsize.alpha.gamma.1.30.quants[i, 3], y1 = popsize.alpha.gamma.1.30.quants[i, 3], col = "red", lwd = 2)
}
leg.text2 <- c("Expect 3 years to fade-out", "Expect 10 years to fade-out")
legend("topleft", bty = "n", leg.text2, lty = c(1, 1), col = c("black", "red"), cex = .6)

#---------------------------------------------
#-- Figure5: alpha, gamma, 30 ye pop si

