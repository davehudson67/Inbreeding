## load libraries
library(data.table)
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
library(MCMCvis)

rm(list=ls())

## set seed
set.seed(44)

## source necessary R functions
source("../../../FirstPaperFiles/Distributions/Dist_GompertzMakeham.R")
source("../../../FirstPaperFiles/Distributions/Dist_GompertzMakehamNim.R")
source("../../../FirstPaperFiles/ModelComparison_FUNCTIONS.R")

## read in data
CH <- readRDS("Badgers_foranalysis_NoSG3323.rds")
sex <- as.numeric(CH$sex) - 1

tKD <- CH$death
tB <- CH$birth

## extract max possible capture time
tM <- CH$max_captures

## extract last alive time
tL <- CH$last_seen

## some checks
stopifnot(all(tL >= tB))
stopifnot(all(tKD[!is.na(tKD)] >= tL[!is.na(tKD)]))

## normalise to survival times
## (necessary at the moment due to censoring
## constraints)
tKD <- tKD - tB
tL <- tL - tB

## define censoring matrices
cint <- cbind(tL, tKD)
cint[is.na(tKD), 2] <- cint[is.na(tKD), 1]
cint[is.na(tKD), 1] <- 0
colnames(cint) <- NULL
censored <- ifelse(!is.na(tKD), 1, 2)
CH$censored <- censored
tD <- rep(NA, length(tKD))
dind <- rep(1, length(tKD))

## extract number of captures
y <- CH$captures
names(y) <- NULL

## some checks
stopifnot(all(tM >= y))

## set up nind
nind <- length(y)

## set up zL/zU input matrix
zL <- NULL
for (i in 1:length(censored)){
  if (censored[i] == 2){
    zL[i] <- cint[i, 2]
  }else{
    zL[i] <- cint[i, 1]
  }
}
zU <- NULL
for (i in 1:length(censored)){
  if (censored[i] == 2){
    zU[i] <- NA
  }else{
    zU[i] <- cint[i, 2]
  }
}
zU[is.na(zU)] <- 0

save.image("badgerSex.RData")
load("../badgerSex.RData")

## set up plot output file
#pdf("RJoutputs/SexDiffFULLRJMCMC_GM.pdf")

code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood for interval-truncated gompertz
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dgompzMakeNim(a * amult[i], b * bmult[i], c * cmult[i])
    
    log(amult[i]) <- beta[1] * sex[i] * z[1]
    log(bmult[i]) <- beta[2] * sex[i] * z[2]
    log(cmult[i]) <- beta[3] * sex[i] * z[3]

        ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  for (j in 1:3) {
    beta[j] ~ dnorm(0, sd = 1)
    z[j] ~ dbern(0.5)
    #    psi[j] <- dunif(0, 1)
  }
  
  a ~ dexp(1)
  b ~ dexp(1)
  c ~ dexp(1)
  mean.p ~ dunif(0, 1)
  
})

## set up data
consts <- list(nind = nind, tM = tM, sex = sex)
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

initFn <- function(cint, censored, sex) {
  ## get ML estimates as initial values
  optFn <- function(pars, t, sex) {
    if(any(pars[c(1:3)] < 0)) {
      return(NA)
    }
    llM <- sum(dGompzMake(t[sex == 1], a = pars[1] * exp(pars[4]), b = pars[2] * exp(pars[5]), c = pars[3] * exp(pars[6]), log = TRUE))
    llF <- sum(dGompzMake(t[sex == 0], a = pars[1], b = pars[2], c = pars[3], log = TRUE))
    llM + llF
  }
  pars <- list(convergence = 1)
  k <- 0
  while(pars$convergence != 0 & k < 50) {
    ## sample missing values
    tD <- tinitFn(cint, censored)
    pars <- optim(c(rexp(3, 10), rnorm(3, 0, 1)), optFn, t = tD, sex = sex, control = list(fnscale = -1))
    k <- k + 1
  }
  if(k == 50) {
    stop("Can't sample initial values")
  }
  pars <- pars$par
  list(
    tD = tD,
    a = pars[1],
    b = pars[2],
    c = pars[3],
    mean.p = runif(1, 0, 1),
    beta = c(pars[4], pars[5], pars[6]),
    z = rep(0, times = 3)
  )
}


## build the model
rIndicatorModel <- nimbleModel(code, constants = consts,
                               data = data, 
                               inits = initFn(cint, censored, sex))

#indicatorModelConf <- configureMCMC(rIndicatorModel)
config <- configureMCMC(rIndicatorModel, monitors = c("a", "b", "c", "mean.p"), thin = 1)
config$removeSamplers(c("a", "b", "c"))
#config$addSampler(target = c("a1", "a2", "b1", "b2", "c1"), type = 'AF_slice', control = 50)
config$addSampler(target = c("a", "c"), type = 'AF_slice', control = 20)
config$addSampler(target = c("b"), type = 'slice')

## Add reversible jump
configureRJ(conf = config,   ## model configuration
            targetNodes = c("beta"),    ## coefficients for selection
            indicatorNodes = c("z"),    ## indicators paired with coefficients
            control = list(mean = 0, scale = 2))

config$addMonitors("beta", "z")

rIndicatorMCMC <- buildMCMC(config)
cIndicatorModel <- compileNimble(rIndicatorModel)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = rIndicatorModel)

system.time(run <- runMCMC(cIndicatorMCMC, 
                           niter = 200000, 
                           nburnin = 24000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))


MCMCsummary(run$samples)
plot(run$samples)
MCMCtrace(run$samples)
samples <- as.matrix(run$samples)

dev.off()
#MCMCplot(run$samples)
#MCMCtrace(run$samples, pdf = F)

## save MCMC
saveRDS(run, "RJoutputs/GM_samples_RJMCMC_Sex.rds")
samples <- readRDS("RJoutputs/samples_RJMCMC_Sex.rds")


## Marginal probabilities of inclusion for each variable
zNames <- rIndicatorModel$expandNodeNames('z')
zCols <- which(colnames(samples) %in% zNames)
binary <- as.data.table((samples[, zCols] != 0) + 0)
res <- binary[ , .N, by=names(binary)]
res <- res[order(N, decreasing = T)]
res <- res[, prob := N/dim(samples)[1]]

res

saveRDS(res, "RJoutputs/posteriorModelProbs_100Kc05.rds")

dev.off()


