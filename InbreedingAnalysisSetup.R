## load libraries
library(nimble)
library(tidyverse)
library(mvtnorm)
library(boot)
library(lamW)
library(GGally)
library(coda)
library(mclust)
library(parallel)
library(survminer)
library(survival)
library(coda)
library(mcmcplots)
library(MCMCvis)
library(scales)
library(data.table)

rm(list=ls())

## source necessary R functions
#source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_Gompertz.R")
#source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_GompertzNim.R")
#source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_GompertzMakeham.R")
#source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_GompertzMakehamNim.R")
source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_Siler.R")
source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_SilerNim.R")
#source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_Expo.R")
source("../SimulationStudy/FirstPaperFiles/ModelComparison_FUNCTIONS.R")

## load data
load("Data/inbreed_data.RData")
CH <- readRDS("Data/BadgersNewREADY300621.rds")

## combine data
CH <- inner_join(CH, mlh)
summary(CH)

############## ADJUST THIS SECTION TO SUIT REQUIRED ANALYSIS ###########################

## filter to required badgers - for Infected as Cubs vs Uninfected life ##

CH <- CH %>%
  filter(infected_lifetime == 0 | infected_as_cub == 1)

CH$infected_lifetime <- as.factor(CH$infected_lifetime)

## set variables
sex <- as.numeric(CH$sex) - 1 # (males = 1, females = 0)
infection <- CH$infected_as_cub # (infected as cub = 1, uninfected lifetime = 0)
inbr <- CH$hom
inbrCAT <- CH$hom
inbrCAT <- if_else(inbrCAT >= median(inbrCAT), 1, 0) # (inbred = 1, outbred = 0)

############################
############################

## filter to required badgers - for Infected at any point in life vs Uninfected life ##
## set variables
sex <- as.numeric(CH$sex) - 1 # (males = 1, females = 0)
infection <- rep(0, nrow(CH))
#infection[CH$infected_lifetime == 0] # (infected lifetime = 1, uninfected lifetime = 0)
infection[CH$infected_as_cub == 1] <- 1

inbr <- CH$hom
inbrCAT <- CH$hom
inbrCAT <- if_else(inbrCAT >= median(inbrCAT), 1, 0) # (inbred = 1, outbred = 0)

############################
############################

## filter to required badgers - for 3 level categorical variable for infection
## set variables
infection3 <- CH$infected_lifetime
infection3[CH$infected_as_cub == 1] <- 2
summary(as.factor(infection3))
infection3 <- infection3 + 1

infection <- matrix(0, nrow = nrow(CH), ncol = 3)

for(i in 1:length(infection3)){
    infection[i, (infection3[i])] <- 1
  }
colnames(infection) <- c("Uninfected", "InfectedAsAdult", "InfectedAsCub")
sex <- as.numeric(CH$sex) - 1 # (males = 1, females = 0)
inbr <- CH$hom
inbrCAT <- CH$hom
inbrCAT <- if_else(inbrCAT >= median(inbrCAT), 1, 0) # (inbred = 1, outbred = 0)

########################################################################################


## set and sample some random seeds for different models
set.seed(42)
seeds <- round(runif(16, 0, 100000000))

## create output folder
dir.create("outputs")

## read in data
tKD <- CH$death
tB <- CH$birth
tL <- CH$last_seen
y <- CH$captures

## extract max possible capture time
tM <- CH$max_captures

## set up censoring vector (1 = interval, 2 = right)
censored <- ifelse(is.na(CH$death), 2, 1)
CH$censored <- censored

## summaries
stopifnot(all(tL[!is.na(tKD)] <= tKD[!is.na(tKD)]))
stopifnot(all(tB <= tL)) #some individuals died on the same day they were born

## normalise to survival times
tKD <- tKD - tB
tL <- tL - tB

## define censoring matrices
cint <- cbind(tL, tKD)
colnames(cint) <- NULL
cint[censored == 2, 2] <- cint[censored == 2, 1] 
cint[censored == 2, 1] <- 0

## check censoring times
summary(apply(cbind(tL, tKD), 1, diff)[censored == 1])
summary(apply(cbind(tL, tKD), 1, diff)[censored == 2])

## set up latent death times
tD <- rep(NA, length(tKD))

## set up nind
nind <- nrow(cint)
dind <- rep(1, length(tKD))

## set up sex data (female = 1, male = 2, unknown = 3)
#sex <- as.numeric(as.factor(CH$sex)) - 1
#sex[sex == 2] <- NA
#g <- length(levels(as.factor(CH$sex))) - 1 #remove unknown group

## set up inbreeding vector
#inb <- CH$f_inbreed

## set up zL/zU input matrix
zL <- tL
zU <- tKD

## save data
#rm(tL, CH, tB)

save.image("Data/badgerSexInb_ICubvUninf.RData")
#save.image("Data/badgerSexInb_FullLifeInfvUninf.RData")
save.image("Data/badgerSexInb_AdultInfvCubInfvUninf.RData")
