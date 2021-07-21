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

## set inbreeding category 
median(inbr)

inbr_level0 <- quantile(inbr, probs = 0.25)
inbr_level1 <- quantile(inbr, probs = 0.5)
inbr_level2 <- quantile(inbr, probs = 0.75)

MLH30 <- if_else(inbr < inbr_level0, 1, 0)
MLH31 <- if_else(inbr >= inbr_level0 & inbr < inbr_level1, 1, 0)
MLH32 <- if_else(inbr >= inbr_level1 & inbr < inbr_level2, 1, 0)
MLH33 <- if_else(inbr >= inbr_level2, 1, 0)
MLH4 <- data.frame(MLH30, MLH31, MLH32, MLH33)

mlh3 <- NA
for(i in 1:nrow(MLH3)) {
  mlh3[i] <- which(MLH3[i, ] > 0) 
}
mlh3 <- as.factor(mlh3)
summary(mlh3)

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
      #sex[i] * (betaSEXINFCUB[1] * infection[i, 3] + betaSEXINFADULT[1] * infection[i, 2]) * zSEXINF[1] +
      (betaINBRmidlow[1] * inbr[i, 2] + betaINBRmidhigh[1] * inbr[i, 3] + betaINBRhigh[1] * inbr[i, 4]) * zINBR[1]
      #betaSEXINBR[1] * sex[i] * inbr[i] * zSEXINBR[1] +
      #(betaINFINBRCUB[1] * infection[i, 3] + betaINFINBRADULT[1] * infection[i, 2]) * inbr[i] * zINFINBR[1]
    
    log(a2mult[i]) <- log(a2) + betaSEX[2] * sex[i] * zSEX[2] + 
      (betaINFCUB[2] * infection[i, 3] + betaINFADULT[2] * infection[i, 2]) * zINF[2] + 
      #sex[i] * (betaSEXINFCUB[2] * infection[i, 3] + betaSEXINFADULT[2] * infection[i, 2]) * zSEXINF[2] +
      (betaINBRmidlow[2] * inbr[i, 2] + betaINBRhigh[2] * inbr[i, 3] + betaINBRhigh[2] * inbr[i, 4]) * zINBR[2]
      #betaSEXINBR[2] * sex[i] * inbr[i] * zSEXINBR[2] +
      #(betaINFINBRCUB[2] * infection[i, 3] + betaINFINBRADULT[2] * infection[i, 2]) * inbr[i] * zINFINBR[2]
    
    log(b1mult[i]) <- log(b1) + betaSEX[3] * sex[i] * zSEX[3] + 
      (betaINFCUB[3] * infection[i, 3] + betaINFADULT[3] * infection[i, 2]) * zINF[3] + 
      #sex[i] * (betaSEXINFCUB[3] * infection[i, 3] + betaSEXINFADULT[3] * infection[i, 2]) * zSEXINF[3] +
      (betaINBRmidlow[3] * inbr[i, 2] + betaINBRhigh[3] * inbr[i, 3] + betaINBRhigh[3] * inbr[i, 4]) * zINBR[3]
      #betaSEXINBR[3] * sex[i] * inbr[i] * zSEXINBR[3] +
      #(betaINFINBRCUB[3] * infection[i, 3] + betaINFINBRADULT[3] * infection[i, 2]) * inbr[i] * zINFINBR[3]
    
    log(b2mult[i]) <- log(b2) + betaSEX[4] * sex[i] * zSEX[4] + 
      (betaINFCUB[4] * infection[i, 3] + betaINFADULT[4] * infection[i, 2]) * zINF[4] + 
      #sex[i] * (betaSEXINFCUB[4] * infection[i, 3] + betaSEXINFADULT[4] * infection[i, 2]) * zSEXINF[4] +
      (betaINBRmidlow[4] * inbr[i, 2] + betaINBRhigh[4] * inbr[i, 3] + betaINBRhigh[4] * inbr[i, 4]) * zINBR[4]
      #betaSEXINBR[4] * sex[i] * inbr[i] * zSEXINBR[4] +
      #(betaINFINBRCUB[4] * infection[i, 3] + betaINFINBRADULT[4] * infection[i, 2]) * inbr[i] * zINFINBR[4]
    
    log(c1mult[i]) <- log(c1) + betaSEX[5] * sex[i] * zSEX[5] + 
      (betaINFCUB[5] * infection[i, 3] + betaINFADULT[5] * infection[i, 2]) * zINF[5] + 
      #sex[i] * (betaSEXINFCUB[5] * infection[i, 3] + betaSEXINFADULT[5] * infection[i, 2]) * zSEXINF[5] +
      (betaINBRmidlow[5] * inbr[i, 2] + betaINBRhigh[5] * inbr[i, 3] + betaINBRhigh[5] * inbr[i, 4]) * zINBR[5]
      #betaSEXINBR[5] * sex[i] * inbr[i] * zSEXINBR[5] +
      #(betaINFINBRCUB[5] * infection[i, 3] + betaINFINBRADULT[5] * infection[i, 2]) * inbr[i] * zINFINBR[5]
    
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  for (k in 1:5) {
    betaSEX[k] ~ dnorm(0, sd = 1)
    betaINFCUB[k] ~ dnorm(0, sd = 1)
    betaINFADULT[k] ~ dnorm(0, sd = 1)
    betaINBRmidlow[k] ~ dnorm(0, sd = 1)
    betaINBRmidhigh[k] ~ dnorm(0, sd = 1)
    betaINBRhigh[k] ~ dnorm(0, sd = 1)
    #betaSEXINFCUB[k] ~ dnorm(0, sd = 1)
    #betaSEXINFADULT[k] ~ dnorm(0, sd = 1)
    #betaSEXINBR[k] ~ dnorm(0, sd = 1)
    #betaINFINBRCUB[k] ~ dnorm(0, sd = 1)
    #betaINFINBRADULT[k] ~ dnorm(0, sd = 1)
    zSEX[k] ~ dbern(0.5)
    zINF[k] ~ dbern(0.5)
    zINBR[k] ~ dbern(0.5)
    #zSEXINF[k] ~ dbern(0.5)
    #zSEXINBR[k] ~ dbern(0.5)
    #zINFINBR[k] ~ dbern(0.5)
    #constraint_dataSEXINF[k] ~ dconstraint(zSEXINF[k] <= zSEX[k] * zINF[k])
    #constraint_dataSEXINBR[k] ~ dconstraint(zSEXINBR[k] <= zSEX[k] * zINBR[k])
    #constraint_dataINFINBR[k] ~ dconstraint(zINFINBR[k] <= zINF[k] * zINBR[k])
    
  }  
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  b2 ~ dexp(1)
  c1 ~ dexp(1)
  mean.p ~ dunif(0, 1)
  
})

## set up data
consts <- list(nind = nind, tM = tM, sex = sex, infection = infection, inbr = MLH4)

data <- list(
    y = y, cint = cint, censored = censored, tD = tD, dind = dind)
#    constraint_dataSEXINF = rep(1, 5), constraint_dataSEXINBR = rep(1, 5),
#    constraint_dataINFINBR = rep(1, 5))

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
            #betaSEXINFCUB = rnorm(5, 0, 1),
            #betaSEXINFADULT = rnorm(5, 0, 1),
            betaINBRmidlow  = rnorm(5, 0, 1),
            betaINBRmidhigh  = rnorm(5, 0, 1),
            betaINBRhigh  = rnorm(5, 0, 1),
            #betaSEXINBR = rnorm(5, 0, 1),
            #betaINFINBRADULT = rnorm(5, 0, 1),
            #betaINFINBRCUB = rnorm(5, 0, 1),
            zSEX = rep(0, 5),
            zINF = rep(0, 5),
            #zSEXINF = rep(0, 5),
            zINBR = rep(0, 5)
            #zSEXINBR = rep(0, 5),
            #zINFINBR = rep(0, 5)
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
config$removeSamplers(c("a1", "a2", "b1", "b2", "c1"))
config$addSampler(target = c("a1"), type = 'slice', control = list(sliceWidth = 0.5, adaptInterval = 50))
config$addSampler(target = c("a2"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 20))
config$addSampler(target = c("b1"), type = 'slice', control = list(sliceWidth = 0.5, adaptInterval = 50))
config$addSampler(target = c("b2"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 20))
config$addSampler(target = c("c1"), type = 'slice', control = list(sliceWidth = 0.5, adaptInterval = 50))

## load in custom RJ-MCMC samplers
source("ModelFitting/MCMC_RJ_multi.R")

## Add reversible jump
configureRJ_multi(conf = config,   ## model configuration
    targetNodes = c("betaSEX", "betaINFCUB", "betaINFADULT", "betaINBRmidlow", "betaINBRmidhigh", "betaINBRhigh"),
    indicatorNodes = c("zSEX", "zINF", "zINF", "zINBR", "zINBR", "zINBR"),
    control = list(mean = 0, scale = 1))
                      
config$addMonitors("betaSEX", "betaINFCUB", "betaINFADULT", "betaINBRmidlow", "betaINBRmidhigh", "betaINBRhigh")
config

rIndicatorMCMC <- buildMCMC(config)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

system.time(run <- runMCMC(cIndicatorMCMC, 
    niter = 100000, 
    nburnin = 10000, 
    nchains = 2, 
    inits = inits[[1]],
    progressBar = TRUE, 
    summary = TRUE, 
    samplesAsCodaMCMC = TRUE, 
    thin = 1))

## save mcmc ouput
saveRDS(run, "outputs/FullModel_NoInteractions_z_all_InbrQuart_runsamples.rds")
#run <- readRDS("outputs/FullModel_NoInteractions_z_all_runsamples.rds")

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
zNames <- model$expandNodeNames(c('z', 'zSEX', 'zINF', 'zINBR'))
zCols <- which(colnames(samples) %in% zNames)
binary <- as.data.table((samples[, zCols] != 0) + 0)
res <- binary[ , .N, by=names(binary)]
res <- res[order(N, decreasing = T)]
res <- res[, prob := N/dim(samples)[1]]
res
res
saveRDS(res, "outputs/FullModel_InbrCont_NoInteractions_z_all_PosteriorModelProbs.rds")
#res <- readRDS("outputs/Sex3Infection_AllParameters_PosteriorModelProbs.rds")

samples <- as.data.frame(samples)

z_indicators <- samples %>%
  select(c(27:41)) %>%
  colSums()

z_indicators <- data.frame(z_indicators/sum(res$N))
z_indicators$parameter <- c("a1_Inf", "a2_Inf", "b1_Inf", "b2_Inf", "c_Inf",
                            "a1_Sex", "a2_Sex", "b1_Sex", "b2_Sex", "c_Sex",
                            "a1_SexInf", "a2SexInf", "b1SexInf", "b2SexInf", "cSexInf")
                            
colnames(z_indicators) <- c("Inclusion_Prob", "parameter")

z_indicators

ggplot(z_indicators, aes(x = parameter, y = Inclusion_Prob)) +
  geom_point()
