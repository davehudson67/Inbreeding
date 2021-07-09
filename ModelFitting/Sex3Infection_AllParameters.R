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
    tD[i] ~ dsilerNim(a1mult[i], a2mult[i], b1mult[i], b2mult[i], c1mult[i])
    
    log(a1mult[i]) <- log(a1) + betaSEX[1] * sex[i] * zSEX[1] + 
      (betaINFCUB[1] * infection[i, 3] + betaINFADULT[1] * infection[i, 2]) * zINF[1] + 
      sex[i] * (betaSEXINFCUB[1] * infection[i, 3] + betaSEXINFADULT[1] * infection[i, 2]) * zSEXINF[1]
    
    log(a2mult[i]) <- log(a2) + betaSEX[2] * sex[i] * zSEX[2] + 
      (betaINFCUB[2] * infection[i, 3] + betaINFADULT[2] * infection[i, 2]) * zINF[2] + 
      sex[i] * (betaSEXINFCUB[2] * infection[i, 3] + betaSEXINFADULT[2] * infection[i, 2]) * zSEXINF[2]
    
    log(b1mult[i]) <- log(b1) + betaSEX[3] * sex[i] * zSEX[3] + 
      (betaINFCUB[3] * infection[i, 3] + betaINFADULT[3] * infection[i, 2]) * zINF[3] + 
      sex[i] * (betaSEXINFCUB[3] * infection[i, 3] + betaSEXINFADULT[3] * infection[i, 2]) * zSEXINF[3]
    
    log(b2mult[i]) <- log(b2) + betaSEX[4] * sex[i] * zSEX[4] + 
      (betaINFCUB[4] * infection[i, 3] + betaINFADULT[4] * infection[i, 2]) * zINF[4] + 
      sex[i] * (betaSEXINFCUB[4] * infection[i, 3] + betaSEXINFADULT[4] * infection[i, 2]) * zSEXINF[4]
    
    log(c1mult[i]) <- log(c1) + betaSEX[5] * sex[i] * zSEX[5] + 
      (betaINFCUB[5] * infection[i, 3] + betaINFADULT[5] * infection[i, 2]) * zINF[5] + 
      sex[i] * (betaSEXINFCUB[5] * infection[i, 3] + betaSEXINFADULT[5] * infection[i, 2]) * zSEXINF[5]
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  for (k in 1:5) {
    betaSEX[k] ~ dnorm(0, sd = 1)
    betaINFCUB[k] ~ dnorm(0, sd = 1)
    betaINFADULT[k] ~ dnorm(0, sd = 1)
    betaSEXINFCUB[k] ~ dnorm(0, sd = 1)
    betaSEXINFADULT[k] ~ dnorm(0, sd = 1)
    zSEX[k] ~ dbern(0.5)
    zINF[k] ~ dbern(0.5)
    zSEXINF[k] ~ dbern(0.5)
    constraint_dataSEXINF[k] ~ dconstraint(zSEXINF[k] <= zSEX[k] * zINF[k])
  }  
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  b2 ~ dexp(1)
  c1 ~ dexp(1)
  mean.p ~ dunif(0, 1)
  
})

## set up data
consts <- list(nind = nind, tM = tM, sex = sex, infection = infection)

data <- list(
    y = y, cint = cint, censored = censored, tD = tD, dind = dind,
    constraint_dataSEXINF = rep(1, 5))

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
initFn <- function(cint, censored, sex, infection3, model) {
    ## get ML estimates as initial values
    optFn <- function(pars, t, sex, infection3) {
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
            pars <- optim(rexp(5, 10), optFn, t = tD, sex = sex, infection3 = infection3, control = list(fnscale = -1))
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
            betaINFCUB = rnorm(5, 0, 1),
            betaINFADULT = rnorm(5, 0, 1),
            betaSEXINFCUB = rnorm(5, 0, 1),
            betaSEXINFADULT = rnorm(5, 0, 1),
            zSEX = rep(0, 5),
            zINF = rep(0, 5),
            zSEXINF = rep(0, 5)
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
    inits[[k]] <- initFn(cint, censored, sex, infection3, cIndicatorModel)
}

## configure MCMC
config <- configureMCMC(model)

## load in custom RJ-MCMC samplers
source("ModelFitting/MCMC_RJ_multi.R")

## Add reversible jump
configureRJ_multi(conf = config,   ## model configuration
    targetNodes = c("betaSEX", "betaINFCUB", "betaINFADULT", "betaSEXINFCUB", "betaSEXINFADULT"),
    indicatorNodes = c("zSEX", "zINF", "zINF", "zSEXINF", "zSEXINF"),
    control = list(mean = 0, scale = 1))
                      
config$addMonitors("betaSEX", "betaINFCUB", "betaINFADULT", "betaSEXINFCUB", "betaSEXINFADULT")
config

rIndicatorMCMC <- buildMCMC(config)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

system.time(run <- runMCMC(cIndicatorMCMC, 
    niter = 20000, 
    nburnin = 5000, 
    nchains = 1, 
    inits = inits[[1]],
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
