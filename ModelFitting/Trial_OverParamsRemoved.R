##########################################################
##                                                      ###
##       Linear regression of Inb against b2            ###
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
#load("Data/badgerSexInb_FullLifeInfvUninf.RData")
load("Data/badgerSexInb_AdultInfvCubInfvUninf.RData")

## load distributions
source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_Siler.R")
source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_SilerNim.R")
source("../SimulationStudy/FirstPaperFiles/ModelComparison_FUNCTIONS.R")

## set seed
set.seed(seeds[15])

## set up plot output file
#pdf("outputs/Sex3Infection_AllParameters.pdf")

code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood for interval-truncated siler
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1 * a1mult[i], a2 * a2mult[i], b1 * b1mult[i], b2 * b2mult[i], c1 * c1mult[i])
    
    log(a1mult[i]) <- betaS[1] * sex[i] * zS[1] + 
      (betaINFCUB[1] * infection[i, 3] * zIC[1] + betaINFADULT[1] * infection[i, 2]) * zIA[1] + 
      sex[i] * (betaSexINFCUB[1] * infection[i, 3] * zSIC[1] + betaSexINFADULT[1] * infection[i, 2]) * zSIA[1]
    
    log(a2mult[i]) <- betaS[2] * sex[i] * zS[2] + 
      (betaINFCUB[2] * infection[i, 3] * zIC[2] + betaINFADULT[2] * infection[i, 2]) * zIA[2] + 
      sex[i] * (betaSexINFCUB[2] * infection[i, 3] * zSIC[2] + betaSexINFADULT[2] * infection[i, 2]) * zSIA[2]
    
    log(b1mult[i]) <- betaS[3] * sex[i] * zS[3] + 
      (betaINFCUB[3] * infection[i, 3] * zIC[3] + betaINFADULT[3] * infection[i, 2]) * zIA[3] + 
      sex[i] * (betaSexINFCUB[3] * infection[i, 3] * zSIC[1] + betaSexINFADULT[3] * infection[i, 2]) * zSIA[3]
    
    log(b2mult[i]) <- betaS[4] * sex[i] * zS[4] + 
      (betaINFCUB[4] * infection[i, 3] * zIC[4]+ betaINFADULT[4] * infection[i, 2]) * zIA[4] + 
      sex[i] * (betaSexINFCUB[4] * infection[i, 3] * zSIC[4] + betaSexINFADULT[4] * infection[i, 2]) * zSIA[4]
    
    log(c1mult[i]) <- betaS[5] * sex[i] * zS[5] + 
      (betaINFCUB[5] * infection[i, 3] * zIC[5] + betaINFADULT[5] * infection[i, 2]) * zIA[5] + 
      sex[i] * (betaSexINFCUB[5] * infection[i, 3] * zSIC[5] + betaSexINFADULT[5] * infection[i, 2]) * zSIA[5]
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  for (k in 1:5) {
    betaS[k] ~ dnorm(0, sd = 1)
    betaINFCUB[k] ~ dnorm(0, sd = 1)
    betaINFADULT[k] ~ dnorm(0, sd = 1)
    betaSexINFCUB[k] ~ dnorm(0, sd = 1)
    betaSexINFADULT[k] ~ dnorm(0, sd = 1)
    zS[k] ~ dbern(0.5)
    zIC[k] ~ dbern(0.5)
    zSIC[k] ~ dbern(0.5)
    zIA[k] ~ dbern(0.5)
    zSIA[k] ~ dbern(0.5)
  }
  
  a1 ~ dexp(10)
  a2 ~ dexp(10)
  b1 ~ dexp(10)
  b2 ~ dexp(10)
  c1 ~ dexp(10)
  mean.p ~ dunif(0, 1)
  constraint_dataSexInf1 ~ dconstraint(zSIC[1] <= zS[1] * zIC[1])
  constraint_dataSexInf2 ~ dconstraint(zSIC[2] <= zS[2] * zIC[2])
  constraint_dataSexInf3 ~ dconstraint(zSIC[3] <= zS[3] * zIC[3])
  constraint_dataSexInf4 ~ dconstraint(zSIC[4] <= zS[4] * zIC[4])
  constraint_dataSexInf5 ~ dconstraint(zSIC[5] <= zS[5] * zIC[5])
  constraint_Inf1 ~ dconstraint(zIC[1] + 1 / zIA[1] + 1 <= 1)
  constraint_Inf2 ~ dconstraint(zIC[2] + zIA[2] >= 2) 
  constraint_Inf3 ~ dconstraint(zIC[3] + zIA[3] >= 2) 
  constraint_Inf4 ~ dconstraint(zIC[4] + zIA[4] >= 2) 
  constraint_Inf5 ~ dconstraint(zIC[5] + zIA[5] >= 2) 
  constraint_Inf6 ~ dconstraint(zSIC[1] + zSIA[1] >= 2) 
  constraint_Inf7 ~ dconstraint(zSIC[2] + zSIA[2] >= 2) 
  constraint_Inf8 ~ dconstraint(zSIC[3] + zSIA[3] >= 2) 
  constraint_Inf9 ~ dconstraint(zSIC[4] + zSIA[4] >= 2) 
  constraint_Inf10 ~ dconstraint(zSIC[5] + zSIA[5] >= 2)
  
  constraint_Inf11 ~ dconstraint(zIC[1] + zIA[1] <= 0)
  constraint_Inf12 ~ dconstraint(zIC[2] + zIA[2] <= 0) 
  constraint_Inf13 ~ dconstraint(zIC[3] + zIA[3] <= 0) 
  constraint_Inf14 ~ dconstraint(zIC[4] + zIA[4] <= 0) 
  constraint_Inf15 ~ dconstraint(zIC[5] + zIA[5] <= 0) 
  constraint_Inf16 ~ dconstraint(zSIC[1] + zSIA[1] <= 0) 
  constraint_Inf17 ~ dconstraint(zSIC[2] + zSIA[2] <= 0) 
  constraint_Inf18 ~ dconstraint(zSIC[3] + zSIA[3] <= 0) 
  constraint_Inf19 ~ dconstraint(zSIC[4] + zSIA[4] <= 0) 
  constraint_Inf20 ~ dconstraint(zSIC[5] + zSIA[5] <= 0)
})

0/0



## set up data
consts <- list(nind = nind, tM = tM, sex = sex, infection = infection, inbr = inbr)

data <- list(y = y, cint = cint, 
             censored = censored, tD = tD, dind = dind, constraint_dataSexInf1 = 1,
             constraint_dataSexInf2 = 1,
             constraint_dataSexInf3 = 1,
             constraint_dataSexInf4 = 1,
             constraint_dataSexInf5 = 1,
             constraint_Inf1 = 1,
             constraint_Inf2 = 1,
             constraint_Inf3 = 1,
             constraint_Inf4 = 1,
             constraint_Inf5 = 1,
             constraint_Inf6 = 1,
             constraint_Inf7 = 1,
             constraint_Inf8 = 1,
             constraint_Inf9 = 1,
             constraint_Inf10 = 1,
             constraint_Inf11 = 1,
             constraint_Inf12 = 1,
             constraint_Inf13 = 1,
             constraint_Inf14 = 1,
             constraint_Inf15 = 1,
             constraint_Inf16 = 1,
             constraint_Inf17 = 1,
             constraint_Inf18 = 1,
             constraint_Inf19 = 1,
             constraint_Inf20 = 1
             )

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
initFn <- function(cint, censored, sex, infection3) {
  ## get ML estimates as initial values
  optFn <- function(pars, t, sex, infection3) {
    # browser()
    if(any(pars[1:5] < 0)) {
      return(NA)
    }
    
    ll <- sum(dSiler(t, a1 = pars[1], a2 = pars[2], b1 = pars[3], b2 = pars[4], c1 = pars[5], log = TRUE))
    
  }
  pars <- list(convergence = 1)
  k <- 0
  while(pars$convergence != 0 & k < 20) {
    ## sample missing values
    tD <- tinitFn(cint, censored)
    pars <- optim(rexp(5, 10), optFn, t = tD, sex = sex, infection3 = infection3, control = list(fnscale = -1))
    k <- k + 1
  }
  if(k == 20) {
    stop("Can't sample initial values")
  }
  pars <- pars$par
  
  ## output initial values
  list(
    tD = tD,
    a1 = pars[1],
    a2 = pars[2],
    b1 = pars[3],
    b2 = pars[4],
    c1 = pars[5],
    mean.p = runif(1, 0, 1),
    betaS = rnorm(5, 0, 1),
    betaINFCUB = rnorm(5, 0, 1),
    betaINFADULT = rnorm(5, 0, 1),
    betaSexINFCUB = rnorm(5, 0, 1),
    betaSexINFADULT = rnorm(5, 0, 1),
    zS = rep(0, times = 5),
    zIC = rep(0, times = 5),
    zSIC = rep(0, times = 5),
    zIA = rep(0, times = 5),
    zSIA = rep(0, times = 5)
  )
}


## build the model
model <- nimbleModel(code, constants = consts,
                     data = data, 
                     inits = initFn(cint, censored, sex, infection3))

## configure model
config <- configureMCMC(model)

config$removeSamplers(c("a1", "a2", "b1", "b2", "c1"))
config$addSampler(target = c("a1"), type = 'slice', control = list(sliceWidth = 0.5, adaptInterval = 50))
config$addSampler(target = c("a2"), type = 'slice', control = list(sliceWidth = 0.5, adaptInterval = 50))
config$addSampler(target = c("b1"), type = 'slice', control = list(sliceWidth = 0.5, adaptInterval = 50))
config$addSampler(target = c("b2"), type = 'slice', control = list(sliceWidth = 0.8, adaptInterval = 20))
config$addSampler(target = c("c1"), type = 'slice', control = list(sliceWidth = 0.5, adaptInterval = 50))
#config$addSampler(target = c("betaINFCUB[1]", "betaINFADULT[1]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))
#config$addSampler(target = c("betaSexINFCUB[1]", "betaSexINFADULT[1]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))
#config$addSampler(target = c("betaINFCUB[2]", "betaINFADULT[2]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))
#config$addSampler(target = c("betaSexINFCUB[2]", "betaSexINFADULT[2]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))
#config$addSampler(target = c("betaINFCUB[3]", "betaINFADULT[3]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))
#config$addSampler(target = c("betaSexINFCUB[3]", "betaSexINFADULT[3]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))
#config$addSampler(target = c("betaINFCUB[4]", "betaINFADULT[4]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))
#config$addSampler(target = c("betaSexINFCUB[4]", "betaSexINFADULT[4]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))
#config$addSampler(target = c("betaINFCUB[5]", "betaINFADULT[5]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))
#config$addSampler(target = c("betaSexINFCUB[5]", "betaSexINFADULT[5]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))

## Add reversible jump
configureRJ(conf = config,   ## model configuration
            targetNodes = c("betaS", "betaINFCUB", "betaINFADULT", "betaSexINFCUB", "betaSexINFADULT"),    ## coefficients for selection
            indicatorNodes = c("zS", "zIC", "zIA", "zSIC", "zSIA"),   ## indicators paired with coefficients
            control = list(mean = 0, scale = 1))

config$addMonitors("betaS", "betaINFCUB", "betaINFADULT", "betaSexINFCUB", "betaSexINFADULT", "zS", "zIC", "zSIC", "zIA", "zSIA" )

config$printSamplers(c("betaS", "betaINFCUB", "betaINFADULT", "betaSexINFCUB", "betaSexINFADULT", "zS", "zIC", "zSIC", "zIA", "zSIA", "a1", "a2", "b1", "b2", "c1", "mean.p"))

rIndicatorMCMC <- buildMCMC(config)
cIndicatorModel <- compileNimble(model)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

system.time(run <- runMCMC(cIndicatorMCMC, 
                           niter = 20000, 
                           nburnin = 5000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

## save mcmc ouput
saveRDS(run, "outputs/Sex3Infection_OverParamsRemoved_runsamples.rds")

samples <- as.matrix(run$samples)
#saveRDS(samples, "outputs/Sex3Infection_AllParameters_samples.rds")

MCMCsummary(run$samples)
#plot(run$samples)
#MCMCtrace(run$samples)

#samples <- samples[sample.int(nrow(samples), ceiling(nrow(samples) * 0.1)), ]
#samples %>%
#  as.data.frame() %>%
#  ggpairs()

#dev.off()
#MCMCplot(run$samples)
#MCMCtrace(run$samples, pdf = F)

## save MCMC
#samples <- readRDS("outputs/Sex3Infection_AllParameters_samples.rds")

## Marginal probabilities of inclusion for each variable
zNames <- model$expandNodeNames('z')
zCols <- which(colnames(samples) %in% zNames)
binary <- as.data.table((samples[, zCols] != 0) + 0)
res <- binary[ , .N, by=names(binary)]
res <- res[order(N, decreasing = T)]
res <- res[, prob := N/dim(samples)[1]]
res
res
saveRDS(res, "outputs/Sex3Infection_OverParamsRemoved_PosteriorModelProbs.rds")
#res <- readRDS("outputs/Sex3Infection_AllParameters_PosteriorModelProbs.rds")

#samples <- as.data.frame(samples)

#z_indicators <- samples %>%
#  select(c(42:56)) %>%
#  colSums()

#z_indicators <- data.frame(z_indicators/sum(res$N))
#z_indicators$parameter <- c("a1_Sex", "a1_Infection", "a1_SexInfection", 
#                            "a2_Sex", "a2_Infection", "a2SexInfection",
#                            "b1_Sex", "b1_Infection", "b1SexInfection",
#                            "b2_Sex", "b2_Infection", "b2SexInfection",
#                            "c_Sex", "c_Infection", "cSexInfection")
#colnames(z_indicators) <- c("Inclusion_Prob", "parameter")#

#z_indicators

#ggplot(z_indicators, aes(x = parameter, y = Inclusion_Prob)) +
#  geom_point()
