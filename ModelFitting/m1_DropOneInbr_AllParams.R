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
inbr_md <- readRDS("inbr_markerdropped.rds")
inbr_md$ID <- rownames(inbr_md)

## reduce dataset to records
inb <- inner_join(CH, inbr_md)
inbr <- select(inb, hom1, hom2, hom3, hom4, hom5, hom6, hom7, hom8, hom9, hom10, 
               hom11, hom12, hom13, hom14, hom15, hom16, hom17, hom18, hom19, hom20,
               hom21, hom22)
inbr <- inbr[-1]
inbr$homALL <- CH$hom

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
    tD[i] ~ dsilerNim(a1 * a1mult[i], a2 * a2mult[i], b1 * b1mult[i], b2 * b2mult[i], c1 * c1mult[i])
    
    log(a1mult[i]) <- 
      betaa1INBR[1] * inbr[i, 1] * z[1] +
      betaa1INBR[2] * inbr[i, 2] * z[2] +
      betaa1INBR[3] * inbr[i, 3] * z[3] +
      betaa1INBR[4] * inbr[i, 4] * z[4] +
      betaa1INBR[5] * inbr[i, 5] * z[5] +
      betaa1INBR[6] * inbr[i, 6] * z[6] +
      betaa1INBR[7] * inbr[i, 7] * z[7] +
      betaa1INBR[8] * inbr[i, 8] * z[8] +
      betaa1INBR[9] * inbr[i, 9] * z[9] +
      betaa1INBR[10] * inbr[i, 10] * z[10] +
      betaa1INBR[11] * inbr[i, 11] * z[11] +
      betaa1INBR[12] * inbr[i, 12] * z[12] +
      betaa1INBR[13] * inbr[i, 13] * z[13] +
      betaa1INBR[14] * inbr[i, 14] * z[14] +
      betaa1INBR[15] * inbr[i, 15] * z[15] +
      betaa1INBR[16] * inbr[i, 16] * z[16] +
      betaa1INBR[17] * inbr[i, 17] * z[17] +
      betaa1INBR[18] * inbr[i, 18] * z[18] +
      betaa1INBR[19] * inbr[i, 19] * z[19] +
      betaa1INBR[20] * inbr[i, 20] * z[20] +
      betaa1INBR[21] * inbr[i, 21] * z[21] +
      betaa1INBR[22] * inbr[i, 22] * z[22] +
      betaa1INBR[23] * inbr[i, 23] * z[23] +
      beta[1] * sex[i]  +
      beta[2] + (betaINFCUB[1] * infection[i, 3] + betaINFADULT[1] * infection[i, 2])
    
    log(a2mult[i]) <-
      betaa2INBR[1] * inbr[i, 1] * z[1] +
      betaa2INBR[2] * inbr[i, 2] * z[2] +
      betaa2INBR[3] * inbr[i, 3] * z[3] +
      betaa2INBR[4] * inbr[i, 4] * z[4] +
      betaa2INBR[5] * inbr[i, 5] * z[5] +
      betaa2INBR[6] * inbr[i, 6] * z[6] +
      betaa2INBR[7] * inbr[i, 7] * z[7] +
      betaa2INBR[8] * inbr[i, 8] * z[8] +
      betaa2INBR[9] * inbr[i, 9] * z[9] +
      betaa2INBR[10] * inbr[i, 10] * z[10] +
      betaa2INBR[11] * inbr[i, 11] * z[11] +
      betaa2INBR[12] * inbr[i, 12] * z[12] +
      betaa2INBR[13] * inbr[i, 13] * z[13] +
      betaa2INBR[14] * inbr[i, 14] * z[14] +
      betaa2INBR[15] * inbr[i, 15] * z[15] +
      betaa2INBR[16] * inbr[i, 16] * z[16] +
      betaa2INBR[17] * inbr[i, 17] * z[17] +
      betaa2INBR[18] * inbr[i, 18] * z[18] +
      betaa2INBR[19] * inbr[i, 19] * z[19] +
      betaa2INBR[20] * inbr[i, 20] * z[20] +
      betaa2INBR[21] * inbr[i, 21] * z[21] +
      betaa2INBR[22] * inbr[i, 22] * z[22] +
      betaa2INBR[23] * inbr[i, 23] * z[23] +
      beta[3] * sex[i]  +
      beta[4] + (betaINFCUB[2] * infection[i, 3] + betaINFADULT[2] * infection[i, 2])
    
    
    log(b1mult[i]) <-
      betab1INBR[1] * inbr[i, 1] * z[1] +
      betab1INBR[2] * inbr[i, 2] * z[2] +
      betab1INBR[3] * inbr[i, 3] * z[3] +
      betab1INBR[4] * inbr[i, 4] * z[4] +
      betab1INBR[5] * inbr[i, 5] * z[5] +
      betab1INBR[6] * inbr[i, 6] * z[6] +
      betab1INBR[7] * inbr[i, 7] * z[7] +
      betab1INBR[8] * inbr[i, 8] * z[8] +
      betab1INBR[9] * inbr[i, 9] * z[9] +
      betab1INBR[10] * inbr[i, 10] * z[10] +
      betab1INBR[11] * inbr[i, 11] * z[11] +
      betab1INBR[12] * inbr[i, 12] * z[12] +
      betab1INBR[13] * inbr[i, 13] * z[13] +
      betab1INBR[14] * inbr[i, 14] * z[14] +
      betab1INBR[15] * inbr[i, 15] * z[15] +
      betab1INBR[16] * inbr[i, 16] * z[16] +
      betab1INBR[17] * inbr[i, 17] * z[17] +
      betab1INBR[18] * inbr[i, 18] * z[18] +
      betab1INBR[19] * inbr[i, 19] * z[19] +
      betab1INBR[20] * inbr[i, 20] * z[20] +
      betab1INBR[21] * inbr[i, 21] * z[21] +
      betab1INBR[22] * inbr[i, 22] * z[22] +
      betab1INBR[23] * inbr[i, 23] * z[23] +
      beta[5] * sex[i]  +
      beta[6] + (betaINFCUB[3] * infection[i, 3] + betaINFADULT[3] * infection[i, 2])
    
    log(b2mult[i]) <-
      betab2INBR[1] * inbr[i, 1] * z[1] +
      betab2INBR[2] * inbr[i, 2] * z[2] +
      betab2INBR[3] * inbr[i, 3] * z[3] +
      betab2INBR[4] * inbr[i, 4] * z[4] +
      betab2INBR[5] * inbr[i, 5] * z[5] +
      betab2INBR[6] * inbr[i, 6] * z[6] +
      betab2INBR[7] * inbr[i, 7] * z[7] +
      betab2INBR[8] * inbr[i, 8] * z[8] +
      betab2INBR[9] * inbr[i, 9] * z[9] +
      betab2INBR[10] * inbr[i, 10] * z[10] +
      betab2INBR[11] * inbr[i, 11] * z[11] +
      betab2INBR[12] * inbr[i, 12] * z[12] +
      betab2INBR[13] * inbr[i, 13] * z[13] +
      betab2INBR[14] * inbr[i, 14] * z[14] +
      betab2INBR[15] * inbr[i, 15] * z[15] +
      betab2INBR[16] * inbr[i, 16] * z[16] +
      betab2INBR[17] * inbr[i, 17] * z[17] +
      betab2INBR[18] * inbr[i, 18] * z[18] +
      betab2INBR[19] * inbr[i, 19] * z[19] +
      betab2INBR[20] * inbr[i, 20] * z[20] +
      betab2INBR[21] * inbr[i, 21] * z[21] +
      betab2INBR[22] * inbr[i, 22] * z[22] +
      betab2INBR[23] * inbr[i, 23] * z[23]
    
    
    log(c1mult[i]) <-
      betacINBR[1] * inbr[i, 1] * z[1] +
      betacINBR[2] * inbr[i, 2] * z[2] +
      betacINBR[3] * inbr[i, 3] * z[3] +
      betacINBR[4] * inbr[i, 4] * z[4] +
      betacINBR[5] * inbr[i, 5] * z[5] +
      betacINBR[6] * inbr[i, 6] * z[6] +
      betacINBR[7] * inbr[i, 7] * z[7] +
      betacINBR[8] * inbr[i, 8] * z[8] +
      betacINBR[9] * inbr[i, 9] * z[9] +
      betacINBR[10] * inbr[i, 10] * z[10] +
      betacINBR[11] * inbr[i, 11] * z[11] +
      betacINBR[12] * inbr[i, 12] * z[12] +
      betacINBR[13] * inbr[i, 13] * z[13] +
      betacINBR[14] * inbr[i, 14] * z[14] +
      betacINBR[15] * inbr[i, 15] * z[15] +
      betacINBR[16] * inbr[i, 16] * z[16] +
      betacINBR[17] * inbr[i, 17] * z[17] +
      betacINBR[18] * inbr[i, 18] * z[18] +
      betacINBR[19] * inbr[i, 19] * z[19] +
      betacINBR[20] * inbr[i, 20] * z[20] +
      betacINBR[21] * inbr[i, 21] * z[21] +
      betacINBR[22] * inbr[i, 22] * z[22] +
      betacINBR[23] * inbr[i, 23] * z[23] +
      beta[7] * sex[i]  +
      beta[8] + (betaINFCUB[4] * infection[i, 3] + betaINFADULT[4] * infection[i, 2])
    
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  for (j in 1:23){
    betaa1INBR[j] ~ dnorm(0, sd = 1)
    betaa2INBR[j] ~ dnorm(0, sd = 1)
    betab1INBR[j] ~ dnorm(0, sd = 1)
    betab2INBR[j] ~ dnorm(0, sd = 1)
    betacINBR[j] ~ dnorm(0, sd = 1)
    z[j] ~ dbern(0.5)
  }
  
  for (k in 1:8){
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
  dconstraint ~ dconstraint(sum(z[1:22]) <= 1)
  
})


## set up data
consts <- list(nind = nind, tM = tM, sex = sex, infection = infection, inbr = inbr)

data <- list(y = y, cint = cint, 
             censored = censored, tD = tD, dind = dind, dconstraint = 1)

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
    beta = rnorm(8, 0, 1),
    betaINFCUB = rnorm(4, 0, 1),
    betaINFADULT = rnorm(4, 0, 1),
    betaa1INBR = rnorm(23, 0, 1),
    betaa2INBR = rnorm(23, 0, 1),
    betab1INBR = rnorm(23, 0, 1),
    betab2INBR = rnorm(23, 0, 1),
    betacINBR = rnorm(23, 0, 1),
    z = rep(0, times = 23)
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
            targetNodes = c("betaa1INBR", "betaa2INBR", "betab1INBR", "betab2INBR", "betacINBR"),    ## coefficients for selection
            indicatorNodes = c("z", "z", "z", "z", "z"),   ## indicators paired with coefficients
            control = list(mean = 0, scale = 2))

config$addMonitors("betaa1INBR", "betaa2INBR", "betab1INBR", "betab2INBR", "betacINBR", "z")

#config$printSamplers(c("beta", "z", "a1", "a2", "b1", "b2", "c1", "mean.p", "betaINFCUB", "betaINFADULT"))

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
saveRDS(run, "outputs/m1_DropOneInbr_AllParameters_runsamples.rds")
#run <- readRDS("outputs/Sex3InfectionInbr_AllParameters_InbrB2runsamples.rds")

samples <- as.matrix(run$samples)
#saveRDS(samples, "outputs/Sex3Infection_AllParameters_samples.rds")

MCMCsummary(run$samples)
mcmcplot(run$samples)
MCMCtrace(run$samples)

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
saveRDS(res, "outputs/m1_DropOneInbr_AllParameters_PosteriorModelProbs.rds")
#res <- readRDS("outputs/Sex3Infection_AllParameters_PosteriorModelProbs.rds")
#samples <- as.matrix(samples)
samples <- as.data.frame(samples)

z_indicators <- samples %>%
  select(c(52:66)) %>%
  colSums()

z_indicators <- data.frame(z_indicators/sum(res$N))
z_indicators$parameter <- c("Inbr_a1", "Inbr:Sex_a1", "Inbr:Inf_a1",
                            "Inbr_a2", "Inbr:Sex_a2", "Inbr:Inf_a2",
                            "Inbr_b1", "Inbr:Sex_b1", "Inbr:Inf_b1",
                            "Inbr_b2", "Inbr:Sex_b2", "Inbr:Inf_b2",
                            "Inbr_c1", "Inbr:Sex_c1", "Inbr:Inf_c1")
colnames(z_indicators) <- c("Inclusion_Prob", "parameter")

z_indicators

ggplot(z_indicators, aes(x = parameter, y = Inclusion_Prob)) +
  geom_point()
