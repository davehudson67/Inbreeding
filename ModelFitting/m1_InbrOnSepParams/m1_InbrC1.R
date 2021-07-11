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
set.seed(seeds[12])

## set up plot output file
#pdf("outputs/Sex3InfectionInbr_AllParameters.pdf")

code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood for interval-truncated siler
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1 * a1mult[i], a2 * a2mult[i], b1 * b1mult[i], b2, c1 * c1mult[i])
    
    log(a1mult[i]) <- #betaINBR * inbr[i] * z +
      beta[1] * sex[i]  +
      beta[2] + (betaINFCUB[1] * infection[i, 3] + betaINFADULT[1] * infection[i, 2])
    #betaINBR[2] * inbr[i] * sex[i] * z[2] +
    #betaINBR[3] * inbr[i] * (betaInbrInfCUB[1] * infection[i, 3] + betaInbrInfADULT[1] * infection[i, 2]) * z[3]
    
    log(a2mult[i]) <- #betaINBR * inbr[i] * z +
      beta[3] * sex[i]  +
      beta[4] + (betaINFCUB[2] * infection[i, 3] + betaINFADULT[2] * infection[i, 2])
    #betaINBR[5] * inbr[i] * sex[i] * z[5] +
    #betaINBR[6] * inbr[i] * (betaInbrInfCUB[2] * infection[i, 3] + betaInbrInfADULT[2] * infection[i, 2]) * z[6]
    
    log(b1mult[i]) <- #betaINBR * inbr[i] * z +
      #beta[5] * sex[i]  +
      beta[5] + (betaINFCUB[3] * infection[i, 3] + betaINFADULT[3] * infection[i, 2])
    #betaINBR[8] * inbr[i] * sex[i] * z[8] +
    #betaINBR[9] * inbr[i] * (betaInbrInfCUB[3] * infection[i, 3] + betaInbrInfADULT[3] * infection[i, 2]) * z[9]
    
    #log(b2mult[i]) <- betaINBR[4] * inbr[i] * z[4]
    #beta[7] * sex[i]  +
    #beta[8] + (betaINFCUB[4] * infection[i, 3] + betaINFADULT[4] * infection[i, 2]) +
    #betaINBR[11] * inbr[i] * sex[i] * z[11] +
    #betaINBR[12] * inbr[i] * (betaInbrInfCUB[4] * infection[i, 3] + betaInbrInfADULT[4] * infection[i, 2]) * z[12]
    
    log(c1mult[i]) <- betaINBR * inbr[i] * z +
      beta[6] * sex[i]  +
      beta[7] + (betaINFCUB[4] * infection[i, 3] + betaINFADULT[4] * infection[i, 2])
    #betaINBR[14] * inbr[i] * sex[i] * z[14] +
    #betaINBR[15] * inbr[i] * (betaInbrInfCUB[5] * infection[i, 3] + betaInbrInfADULT[5] * infection[i, 2]) * z[15]
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  for (k in 1:7){
    beta[k] ~ dnorm(0, sd = 1)
  }
  
  for (m in 1:4) {
    betaINFCUB[m] ~ dnorm(0, sd = 1)
    betaINFADULT[m] ~ dnorm(0, sd = 1)
  }
  
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  b2 ~ dexp(1)
  c1 ~ dexp(1)
  mean.p ~ dunif(0, 1)
  betaINBR ~ dnorm(0, sd = 1)
  z ~ dbern(0.5)
})


## set up data
consts <- list(nind = nind, tM = tM, sex = sex, infection = infection, inbr = inbr)

data <- list(y = y, cint = cint, 
             censored = censored, tD = tD, dind = dind)

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
initFn <- function(cint, censored) {
  ## get ML estimates as initial values
  optFn <- function(pars, t) {
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
    pars <- optim(rexp(5, 10), optFn, t = tD, control = list(fnscale = -1))
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
    beta = rnorm(7, 0, 1),
    betaINFCUB = rnorm(4, 0, 1),
    betaINFADULT = rnorm(4, 0, 1),
    betaINBR = rnorm(1, 0, 1),
    z = 0
  )
}


## build the model
model <- nimbleModel(code, constants = consts,
                     data = data, 
                     inits = initFn(cint, censored))

## configure model
config <- configureMCMC(model)
#help(samplers)
config$removeSamplers(c("a1", "a2", "b1", "b2", "c1", "betaINFCUB", "betaINFADULT"))
config$addSampler(target = c("a1"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))
config$addSampler(target = c("a2"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))
config$addSampler(target = c("b1"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))
config$addSampler(target = c("b2"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 20))
config$addSampler(target = c("c1"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))
config$addSampler(target = c("betaINFCUB[1]", "betaINFADULT[1]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 100))
#config$addSampler(target = c("betaInbrInfCUB[1]", "betaInbrInfADULT[1]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 100))
config$addSampler(target = c("betaINFCUB[2]", "betaINFADULT[2]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 100))
#config$addSampler(target = c("betaInbrInfCUB[2]", "betaInbrInfADULT[2]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 100))
config$addSampler(target = c("betaINFCUB[3]", "betaINFADULT[3]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 100))
#config$addSampler(target = c("betaInbrInfCUB[3]", "betaInbrInfADULT[3]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 100))
config$addSampler(target = c("betaINFCUB[4]", "betaINFADULT[4]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 100))
#config$addSampler(target = c("betaInbrInfCUB[4]", "betaInbrInfADULT[4]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 100))

## Add reversible jump
configureRJ(conf = config,   ## model configuration
            targetNodes = c("betaINBR"),    ## coefficients for selection
            indicatorNodes = c("z"),   ## indicators paired with coefficients
            control = list(mean = 0, scale = 2))

config$addMonitors("beta", "z")

config$printSamplers(c("beta", "z", "a1", "a2", "b1", "b2", "c1", "mean.p", "betaINFCUB", "betaINFADULT"))

rIndicatorMCMC <- buildMCMC(config)
cIndicatorModel <- compileNimble(model)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

system.time(run <- runMCMC(cIndicatorMCMC, 
                           niter = 500000, 
                           nburnin = 20000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

## save mcmc ouput
saveRDS(run, "outputs/m1_InbrOnSepParams/m1_InbrC1_runsamples.rds")
#run <- readRDS("outputs/Sex3InfectionInbr_AllParameters_InbrB2runsamples.rds")

samples <- as.matrix(run$samples)
#saveRDS(samples, "outputs/Sex3Infection_AllParameters_samples.rds")

MCMCsummary(run$samples)
#mcmcplot(run$samples)
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
saveRDS(res, "outputs/m1_InbrOnSepParams/m1_InbrC1_PosteriorModelProbs.rds")

#res <- readRDS("outputs/Sex3Infection_AllParameters_PosteriorModelProbs.rds")
#samples <- as.matrix(samples)
#samples <- as.data.frame(samples)

#z_indicators <- samples %>%
#  select(c(27:31)) %>%
#  colSums()

#z_indicators <- data.frame(z_indicators/sum(res$N))
#z_indicators$parameter <- c("Inbr_a1",
#                            "Inbr_a2",
#                            "Inbr_b1",
#                            "Inbr_b2",
#                            "Inbr_c1")
#colnames(z_indicators) <- c("Inclusion_Prob", "parameter")

#z_indicators

#ggplot(z_indicators, aes(x = parameter, y = Inclusion_Prob)) +
#  geom_point()
