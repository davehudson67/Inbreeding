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
load("Data/badgerSexInb_FullLifeInfvUninf.RData")
#load("Data/badgerSexInb_AdultInfvCubInfvUninf.RData")
load("Data/badgerSexInb_ICubvUninf.RData")

## load distributions
source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_Siler.R")
source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_SilerNim.R")
source("../SimulationStudy/FirstPaperFiles/ModelComparison_FUNCTIONS.R")

## set seed
set.seed(seeds[11])

## set up plot output file
#pdf("outputs/Sex3Infection_AllParameters.pdf")

code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood for interval-truncated siler
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1mult[i], a2mult[i], b1mult[i], b2mult[i], c1mult[i])
    
    log(a1mult[i]) <- log(a1) + betaSEX[1] * sex[i] * zSEX[1] + 
      betaINFCUB[1] * infection[i]* zINF[1] + 
      betaINBR[1] * inbr[i] * zINBR[1]
    
    log(a2mult[i]) <- log(a2) + betaSEX[2] * sex[i] * zSEX[2] + 
      betaINFCUB[2] * infection[i] * zINF[2] + 
      betaINBR[2] * inbr[i] * zINBR[2]
    
    log(b1mult[i]) <- log(b1) + betaSEX[3] * sex[i] * zSEX[3] + 
      betaINFCUB[3] * infection[i] * zINF[3] + 
      betaINBR[3] * inbr[i] * zINBR[3]
    
    log(b2mult[i]) <- log(b2) + betaSEX[4] * sex[i] * zSEX[4] + 
      betaINFCUB[4] * infection[i]  * zINF[4] + 
      betaINBR[4] * inbr[i] * zINBR[4]
    
    log(c1mult[i]) <- log(c1) + betaSEX[5] * sex[i] * zSEX[5] + 
      betaINFCUB[5] * infection[i] * zINF[5] + 
      betaINBR[5] * inbr[i] * zINBR[5]
    
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  for (k in 1:5) {
    betaSEX[k] ~ dnorm(0, sd = 1)
    betaINFCUB[k] ~ dnorm(0, sd = 1)
    betaINBR[k] ~dnorm(0, sd = 1)
    zSEX[k] ~ dbern(0.5)
    zINF[k] ~ dbern(0.5)
    zINBR[k] ~ dbern(0.5)
    
  }  
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  b2 ~ dexp(1)
  c1 ~ dexp(1)
  mean.p ~ dunif(0, 1)
  
})

## set up data
consts <- list(nind = nind, tM = tM, sex = sex, infection = infection, inbr = inbrCAT)

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
            #betaINFADULT = rnorm(5, 0, 1),
            #betaSEXINFCUB = rnorm(5, 0, 1),
            #betaSEXINFADULT = rnorm(5, 0, 1),
            betaINBR  = rnorm(5, 0, 1),
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
config$addSampler(target = c("a1"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 30))
config$addSampler(target = c("a2"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 30))
config$addSampler(target = c("b1"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 30))
config$addSampler(target = c("b2"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 30))
config$addSampler(target = c("c1"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 30))

## load in custom RJ-MCMC samplers
source("ModelFitting/MCMC_RJ_multi.R")

## Add reversible jump
configureRJ_multi(conf = config,   ## model configuration
    targetNodes = c("betaSEX", "betaINFCUB", "betaINBR"),
    indicatorNodes = c("zSEX", "zINF", "zINBR"),
    control = list(mean = 0, scale = 1))
                      
config$addMonitors("betaSEX", "betaINFCUB", "betaINBR")
config

rIndicatorMCMC <- buildMCMC(config)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

system.time(run <- runMCMC(cIndicatorMCMC, 
    niter = 300000, 
    nburnin = 24000, 
    nchains = 2, 
    inits = inits[[1]],
    progressBar = TRUE, 
    summary = TRUE, 
    samplesAsCodaMCMC = TRUE, 
    thin = 1))

## save mcmc ouput
saveRDS(run, "outputs/FullModel_InbrCont_NoInteractions_z_all_runsamples.rds")
run <- readRDS("outputs/FullModel_InbrCont_NoInteractions_z_all_runsamples.rds")

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
saveRDS(res, "outputs/FullModel_NoInteractions_z_all_PosteriorModelProbs1.rds")
res <- readRDS("outputs/FullModel_NoInteractions_z_all_PosteriorModelProbs1.rds")

samples <- as.data.frame(samples)

z_indicators <- samples %>%
  select(c(27:41)) %>%
  colSums()

z_indicators <- data.frame(z_indicators/sum(res$N))
#z_indicators$parameter <- c("a1_Inf", "a2_Inf", "b1_Inf", "b2_Inf", "c_Inf",
#                            "a1_Sex", "a2_Sex", "b1_Sex", "b2_Sex", "c_Sex",
#                            "a1_SexInf", "a2SexInf", "b1SexInf", "b2SexInf", "cSexInf")
                            
z_indicators$z <- rownames(z_indicators)
colnames(z_indicators) <- c("Inclusion_Prob", "z")
z_indicators$variable <- rep(c("Inbreeding", "Infection", "Sex"), each = 5)

incl <- ggplot(z_indicators, aes(x = z, y = Inclusion_Prob)) +
  geom_point(aes(colour = variable)) + 
  geom_hline(yintercept = 0.5, colour = "red") +
  scale_y_continuous(limits = c(0,1)) + 
  labs(title = "Main effects only model | Continuous inbreeding", subtitle = "[1] = a1, [2] = a2, [3] = b1, [4] = b2, [5] = c")
incl
