##########################################################
##                                                      ###
##       Inbreeding, Sex, Infection and SG              ###
##                                                      ###
###########################################################

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

## load data
load("../Data/badgerSexInb_ICubvUninf_SG.RData")
CH <- readRDS("../Data/CH_AllInbreedSG.rds")
#load("BadgersNimbleSG.RData")

## load distributions
source("../../SimulationStudy/FirstPaperFiles/Distributions/Dist_Siler.R")
source("../../SimulationStudy/FirstPaperFiles/Distributions/Dist_SilerNim.R")
#source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_GompertzMakeham.R")
#source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_GompertzMakehamNim.R")
source("../../SimulationStudy/FirstPaperFiles/ModelComparison_FUNCTIONS.R")

## make infected male the 'intercept' category
#sex <- ifelse(sex == 1, 0, 1)
#infection <- ifelse(infection == 1, 0, 1)

levels(as.factor(CH$socg))
socg <- as.numeric(CH$socg)
socg <- ifelse(socg > 6, socg - 1, socg)
socg <- ifelse(socg > 26, socg - 1, socg)
socg <- ifelse(socg > 30, socg - 1, socg)
socg <- ifelse(socg > 34, socg - 1, socg)


## set seed
set.seed(seeds[11])

nSG <- 40
## set up plot output file
#pdf("outputs/Sex3Infection_AllParameters.pdf")

code <- nimbleCode({
  ## survival components for dead badgers
  for (i in 1:nind) {
    ## likelihood for interval-truncated siler
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1mult[i], a2mult[i], b1mult[i], b2mult[i], c1mult[i])
    
    log(a1mult[i]) <- log(a1) + betaSEX[1] * sex[i] + 
      betaINF[1] * infection[i] +
      betaINBR[1] * inbr[i] * zINBR[1] +
      betaINBRINF[1] * inbr[i] * infection[i] * zINBRINF[1] +
      betaINBRSEX[1] * inbr[i] * sex[i] * zINBRSEX[1] +
      betaINBRINFSEX[1] * inbr[i] * infection[i] * sex[i] * zINBRINFSEX[1] +
      gammaSG[SG[i]]
    
    log(a2mult[i]) <- log(a2) + betaSEX[2] * sex[i] + 
      betaINF[2] * infection[i] +
      betaINBR[2] * inbr[i] * zINBR[2] +
      betaINBRINF[2] * inbr[i] * infection[i] * zINBRINF[2] +
      betaINBRSEX[2] * inbr[i] * sex[i] * zINBRSEX[2] +
      betaINBRINFSEX[2] * inbr[i] * infection[i] * sex[i] * zINBRINFSEX[2] +
      gammaSG[SG[i]]
    
    log(b1mult[i]) <- log(b1) + betaSEX[3] * sex[i] + 
      betaINF[3] * infection[i] +
      betaINBR[3] * inbr[i] * zINBR[3] +
      betaINBRINF[3] * inbr[i] * infection[i] * zINBRINF[3] +
      betaINBRSEX[3] * inbr[i] * sex[i] * zINBRSEX[3] +
      betaINBRINFSEX[3] * inbr[i] * infection[i] * sex[i] * zINBRINFSEX[3]
    
    log(b2mult[i]) <- log(b2) + betaSEX[4] * sex[i] + 
      betaINF[4] * infection[i] +
      betaINBR[4] * inbr[i] * zINBR[4] +
      betaINBRINF[4] * inbr[i] * infection[i] * zINBRINF[4] +
      betaINBRSEX[4] * inbr[i] * sex[i] * zINBRSEX[4] +
      betaINBRINFSEX[4] * inbr[i] * infection[i] * sex[i] * zINBRINFSEX[4]
    
    log(c1mult[i]) <- log(c1) + betaSEX[5] * sex[i] + 
      betaINF[5] * infection[i] +
      betaINBR[5] * inbr[i] * zINBR[5] +
      betaINBRINF[5] * inbr[i] * infection[i] * zINBRINF[5] +
      betaINBRSEX[5] * inbr[i] * sex[i] * zINBRSEX[5] +
      betaINBRINFSEX[5] * inbr[i] * infection[i] * sex[i] * zINBRINFSEX[5] +
      gammaSG[SG[i]]
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  for (k in 1:5) {
    betaSEX[k] ~ dnorm(0, sd = 2)
    betaINF[k] ~ dnorm(0, sd = 2)
    betaINBR[k] ~ dnorm(0, sd = 2)
    betaINBRINF[k] ~ dnorm(0, sd = 2)
    betaINBRSEX[k] ~ dnorm(0, sd = 2)
    betaINBRINFSEX[k] ~ dnorm(0, sd = 2)
    zINBRINF[k] ~ dbern(0.5)
    zINBRSEX[k] ~ dbern(0.5)
    zINBRINFSEX[k] ~ dbern(0.5)
    zINBR[k] ~ dbern(0.5)
 
    ## constraints for interaction terms
    constraint_dataINBRINF[k] ~ dconstraint(zINBRINF[k] <= zINBR[k])
    constraint_dataINBRSEX[k] ~ dconstraint(zINBRSEX[k] <= zINBR[k])
    constraint_dataINBRINFSEX[k] ~ dconstraint(zINBRINFSEX[k] <= zINBR[k] * zINBRINF[k] * zINBRSEX[k])
  }  
  
  for (r in 1:nSG){
    gammaSG[r] ~ dnorm(0, tau)
  }
  
  tau <- pow(sigma, -2)
  a1 ~ dexp(0.1)
  a2 ~ dexp(0.1)
  b1 ~ dexp(0.1)
  b2 ~ dexp(0.1)
  c1 ~ dexp(0.1)
  mean.p ~ dunif(0, 1)
  sigma ~ dunif(0, 10)
})

## set up data
consts <- list(nind = nind, tM = tM, sex = sex, infection = infection, inbr = inbr, SG = socg)

data <- list(
  y = y, cint = cint, censored = censored, tD = tD, dind = dind, constraint_dataINBRINF = rep(1, times = 5),
  constraint_dataINBRSEX = rep(1, times = 5), constraint_dataINBRINFSEX = rep(1, times = 5))

## find overdispersed initial values
tinitFn <- function(cint, censored) {
  apply(cbind(cint, censored), 1, function(x) {
    if(x[3] == 2) {
      y <- x[2] + rexp(1, 1)
    } else {
      y <- runif(1, x[1], x[2])
    }
    y
  })
}

initFn <- function(cint, censored, model) {
  ## get ML estimates as initial values
  optFn <- function(pars, t) {
    # browser()
    if(any(pars[1:5] < 0)) {
      return(NA)
    }
    ll <- sum(dSiler(t, a1 = pars[1], a2 = pars[2], b1 = pars[3], b2 = pars[4], c1 = pars[5], log = TRUE))  
  }
  valid <- 0
  while(valid == 0) {
    pars <- list(convergence = 1)
    k <- 0
    while(pars$convergence != 0 & k < 20) {
      ## sample missing values
      tD <- tinitFn(cint, censored)
      pars <- optim(rexp(5, 10), optFn, t = tD, control = list(fnscale = -1))
      k <- k + 1
    }
    if(k == 20) {
      stop("Can't sample initial values")
    }
    pars <- pars$par
    
    ## output initial values
    inits <- list(
      tD = tD,
      a1 = pars[1],
      a2 = pars[2],
      b1 = pars[3],
      b2 = pars[4],
      c1 = pars[5],
      mean.p = runif(1, 0, 1),
      betaSEX = rnorm(5, 0, 1),
      betaINF = rnorm(5, 0, 1),
      betaINBR  = rnorm(5, 0, 1),
      gammaSG = rnorm(40, 0, 1),
      zINBR = rep(0, times = 5),
      sigma = runif(1, 0, 5),
      betaINBRINF = rnorm(5, 0, 1),
      betaINBRSEX = rnorm(5, 0, 1),
      betaINBRINFSEX = rnorm(5, 0, 1),
      zINBRINF = rep(0, times = 5),
      zINBRSEX = rep(0, times = 5),
      zINBRINFSEX = rep(0, times = 5)
    )
    model$setInits(inits)
    valid <- ifelse(!is.finite(model$calculate()), 0, 1)
  }
  return(inits)
}

## build the model without initial values
## (will throw an initialisation warning)
model <- nimbleModel(code, constants = consts, data = data)

## compile the model
cIndicatorModel <- compileNimble(model)

## find list of valid initial values (needs compiled model
## for some reason)
inits <- list()
for(k in 1:2) {
  inits[[k]] <- initFn(cint, censored, cIndicatorModel)
}

## configure MCMC
config <- configureMCMC(model)
config$removeSamplers(c("a1", "a2", "b1", "b2", "c1"))
#config$addSampler(target = c("a1", "a2", "b1", "b2", "c1"), type = 'AF_slice', control = list(sliceWidth = 2, adaptInterval = 30))

config$addSampler(target = c("a1"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))
config$addSampler(target = c("a2"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))
config$addSampler(target = c("b1"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))
config$addSampler(target = c("b2"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))
config$addSampler(target = c("c1"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))

## load in custom RJ-MCMC samplers
source("../ModelFitting/MCMC_RJ_multi.R")

## Add reversible jump
configureRJ_multi(conf = config,   ## model configuration
                  targetNodes = c("betaINBR", "betaINBRINF", "betaINBRSEX", "betaINBRINFSEX"),
                  indicatorNodes = c("zINBR", "zINBRINF", "zINBRSEX", "zINBRINFSEX"),
                  control = list(mean = 0, scale = 1))

config$addMonitors("betaSEX", "betaINF", "betaINBR", "betaINBRINF", "betaINBRSEX", "betaINBRINFSEX")
config

rIndicatorMCMC <- buildMCMC(config)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

system.time(run <- runMCMC(cIndicatorMCMC, 
                           niter = 50000, 
                           nburnin = 2400, 
                           nchains = 2, 
                           inits = inits[[1]],
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

## save mcmc ouput
saveRDS(run, "Outputs/SilerFullModel_InbrCONTNewGeneData_runsamples.rds")

## save mcmc ouput
saveRDS(run, "Outputs/GMFullModel_InbrCATtD_runsamples.rds")
run <- readRDS("Outputs/SilerFullModel_InbrCONT_runsamples.rds")

samples <- as.matrix(run$samples)
#saveRDS(samples, "outputs/Sex3Infection_AllParameters_samples.rds")

MCMCsummary(run$samples)
plot(run$samples)
MCMCtrace(run$samples)
samples <- as.matrix(postdf)
samples <- samples[sample.int(nrow(samples), ceiling(nrow(samples) * 0.1)), ]
samples %>%
  as.data.frame() %>%
  select(betaINBRa1, betaINBRa2, betaINBRb1, betaINBRb2, betaINBRc) %>%
  ggpairs()

library(data.table)
zNames <- model$expandNodeNames(c('zINBR', 'zINBRINF', 'zINBRINFSEX', 'zINBRSEX'))
zCols <- which(colnames(samples) %in% zNames)
binary <- as.data.table((samples[, zCols] != 0) + 0)
res <- binary[ , .N, by=names(binary)]
res <- res[order(N, decreasing = T)]
res <- res[, prob := N/dim(samples)[1]]

res

