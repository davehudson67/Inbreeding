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
    
    log(a1mult[i]) <- beta[1] * sex[i] * z[1] + 
      beta[2] + (betaINFCUB[1] * infection[i, 3] + betaINFADULT[1] * infection[i, 2]) * z[2] + 
      beta[3] * sex[i] * (betaSexINFCUB[1] * infection[i, 3] + betaSexINFADULT[1] * infection[i, 2]) * z[3]
    
    log(a2mult[i]) <- beta[4] * sex[i] * z[4] + 
      beta[5] + (betaINFCUB[2] * infection[i, 3] + betaINFADULT[2] * infection[i, 2]) * z[5] + 
      beta[6] * sex[i] * (betaSexINFCUB[2] * infection[i, 3] + betaSexINFADULT[2] * infection[i, 2]) * z[6]
    
    log(b1mult[i]) <- beta[7] * sex[i] * z[7] + 
      beta[8] + (betaINFCUB[3] * infection[i, 3] + betaINFADULT[3] * infection[i, 2]) * z[8] + 
      beta[9] * sex[i] * (betaSexINFCUB[3] * infection[i, 3] + betaSexINFADULT[3] * infection[i, 2]) * z[9]
    
    log(b2mult[i]) <- beta[10] * sex[i] * z[10] + 
      beta[11] + (betaINFCUB[4] * infection[i, 3] + betaINFADULT[4] * infection[i, 2]) * z[11] + 
      beta[12] * sex[i] * (betaSexINFCUB[4] * infection[i, 3] + betaSexINFADULT[4] * infection[i, 2]) * z[12]
    
    log(c1mult[i]) <- beta[13] * sex[i] * z[13] + 
      beta[14] + (betaINFCUB[5] * infection[i, 3] + betaINFADULT[5] * infection[i, 2]) * z[14] + 
      beta[15] * sex[i] * (betaSexINFCUB[5] * infection[i, 3] + betaSexINFADULT[5] * infection[i, 2]) * z[15]
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  for (j in 1:15){
    beta[j] ~ dnorm(0, sd = 1)
    z[j] ~ dbern(0.5)
  }
  
  for (k in 1:5) {
    betaINFCUB[k] ~ dnorm(0, sd = 1)
    betaINFADULT[k] ~ dnorm(0, sd = 1)
    betaSexINFCUB[k] ~ dnorm(0, sd = 1)
    betaSexINFADULT[k] ~ dnorm(0, sd = 1)
  }
  
  a1 ~ dexp(10)
  a2 ~ dexp(10)
  b1 ~ dexp(10)
  b2 ~ dexp(10)
  c1 ~ dexp(10)
  mean.p ~ dunif(0, 1)
  constraint_dataSexInf1 ~ dconstraint(z[3] <= z[1] * z[2])
  constraint_dataSexInf2 ~ dconstraint(z[6] <= z[4] * z[5])
  constraint_dataSexInf3 ~ dconstraint(z[9] <= z[7] * z[8])
  constraint_dataSexInf4 ~ dconstraint(z[12] <= z[10] * z[11])
  constraint_dataSexInf5 ~ dconstraint(z[14] <= z[13] * z[14])
  
})


## set up data
consts <- list(nind = nind, tM = tM, sex = sex, infection = infection, inbr = inbr)

data <- list(y = y, cint = cint, 
             censored = censored, tD = tD, dind = dind, constraint_dataSexInf1 = 1,
             constraint_dataSexInf2 = 1,
             constraint_dataSexInf3 = 1,
             constraint_dataSexInf4 = 1,
             constraint_dataSexInf5 = 1)

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
    beta = rnorm(15, 0, 1),
    betaINFCUB = rnorm(5, 0, 1),
    betaINFADULT = rnorm(5, 0, 1),
    betaSexINFCUB = rnorm(5, 0, 1),
    betaSexINFADULT = rnorm(5, 0, 1),
    z = rep(0, times = 15)
  )
}


## build the model
model <- nimbleModel(code, constants = consts,
                     data = data, 
                     inits = initFn(cint, censored, sex, infection3))

## configure model
config <- configureMCMC(model)

config$removeSamplers(c("a1", "a2", "b1", "b2", "c1", "betaINFCUB", "betaSexINFCUB"))
#config$addSampler(target = c("a1", "a2", "b1", "b2", "c1"), type = 'AF_slice', control = 20)
#config$addSampler(target = c("a1", "b1", "c1"), type = 'AF_slice', control = 20)
config$addSampler(target = c("a1"), type = 'slice', control = list(sliceWidth = 0.5, adaptInterval = 50))
config$addSampler(target = c("a2"), type = 'slice', control = list(sliceWidth = 0.5, adaptInterval = 50))
config$addSampler(target = c("b1"), type = 'slice', control = list(sliceWidth = 0.5, adaptInterval = 50))
config$addSampler(target = c("b2"), type = 'slice', control = list(sliceWidth = 0.8, adaptInterval = 20))
config$addSampler(target = c("c1"), type = 'slice', control = list(sliceWidth = 0.5, adaptInterval = 50))
#config$addSampler(target = c("betaINFADULT[4]"), type = 'slice')
config$addSampler(target = c("a1", "betaINFCUB[1]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))
config$addSampler(target = c("a1", "betaSexINFCUB[1]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))
config$addSampler(target = c("a2", "betaINFCUB[2]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))
config$addSampler(target = c("a2", "betaSexINFCUB[2]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))
config$addSampler(target = c("b1", "betaINFCUB[3]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))
config$addSampler(target = c("b1", "betaSexINFCUB[3]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))
config$addSampler(target = c("b2", "betaINFCUB[4]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 100))
config$addSampler(target = c("b2", "betaSexINFCUB[4]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 100))
config$addSampler(target = c("c1", "betaINFCUB[5]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))
config$addSampler(target = c("c1", "betaSexINFCUB[5]"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 200))

## Add reversible jump
configureRJ(conf = config,   ## model configuration
            targetNodes = c("beta"),    ## coefficients for selection
            indicatorNodes = c("z"),   ## indicators paired with coefficients
            control = list(mean = 0, scale = 1))

config$addMonitors("beta", "z")

config$printSamplers(c("beta", "z", "a1", "a2", "b1", "b2", "c1", "mean.p"))

rIndicatorMCMC <- buildMCMC(config)
cIndicatorModel <- compileNimble(model)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

system.time(run <- runMCMC(cIndicatorMCMC, 
                           niter = 100000, 
                           nburnin = 15000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

## save mcmc ouput
saveRDS(run, "outputs/Sex3Infection_AllParameters_runsamples.rds")

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
saveRDS(res, "outputs/Sex3Infection_AllParameters_PosteriorModelProbs.rds")
#res <- readRDS("outputs/Sex3Infection_AllParameters_PosteriorModelProbs.rds")

samples <- as.data.frame(samples)

z_indicators <- samples %>%
  select(c(42:56)) %>%
  colSums()

z_indicators <- data.frame(z_indicators/sum(res$N))
z_indicators$parameter <- c("a1_Sex", "a1_Infection", "a1_SexInfection", 
                            "a2_Sex", "a2_Infection", "a2SexInfection",
                            "b1_Sex", "b1_Infection", "b1SexInfection",
                            "b2_Sex", "b2_Infection", "b2SexInfection",
                            "c_Sex", "c_Infection", "cSexInfection")
colnames(z_indicators) <- c("Inclusion_Prob", "parameter")

z_indicators

ggplot(z_indicators, aes(x = parameter, y = Inclusion_Prob)) +
  geom_point()
