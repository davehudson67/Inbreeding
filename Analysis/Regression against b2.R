## load libraries
library(nimble)
library(coda)
library(mcmcplots)
library(lamW)
library(data.table)
library(tidyverse)

source("Distributions/Siler/Dist_SilerNim.R")
source("Distributions/Siler/Dist_Siler.R")

## read data:
 

#####################################################################
##                                                                 ##
##            Never Positive badgers                               ##
##                                                                 ##
#####################################################################

## read data:
CH.all <- as.data.table(readRDS("Data/CHCVall.rds"))
#CH <- CH.all

## select group
levels(CH.all$status)
levels(CH.all$status.cub)
CH <- CH.all[CH.all$status == 'NPF' | CH.all$status == 'NPM'] #never positive badgers
#CH <- CH.all[CH.all$status.cub == 'CPF' | CH.all$status.cub == 'CPM'] #cub positive
#CH <- CH.all[CH.all$status == 'PF' | CH.all$status == 'PM'] # positive
inbreed <- CH$f_inbreed

## extract death and birth times
tKD <- CH$death 
tB <- CH$birth

## remove extra variable from CH
colnames(CH)
CH <- CH[,31:183]

## extract max possible death time
tM <- ifelse(is.na(tKD), ncol(CH), tKD)

## extract last alive time
tL <- as.numeric(apply(CH, 1, function(x) max(which(x == 1))))

## some checks
stopifnot(all(tL > tB))
stopifnot(all(tKD[!is.na(tKD)] > tL[!is.na(tKD)]))
stopifnot(all(tM >= tL))

## normalise to survival times
## (necessary at the moment due to censoring
## constraints)
tM <- tM - tB
tKD <- tKD - tB
tL <- tL - tB

## define censoring matrices
cint <- cbind(tL, tKD)
cint[is.na(tKD), 2] <- cint[is.na(tKD), 1]
cint[is.na(tKD), 1] <- 0
colnames(cint) <- NULL
censored <- ifelse(!is.na(tKD), 1, 2)
tD <- rep(NA, length(tKD))
dind <- rep(1, length(tKD))

## extract number of captures
y <- rowSums(CH)
names(y) <- NULL

## some checks
stopifnot(all(tM >= y))

## set up nind
nind <- length(y)

## code for NIMBLE model with censoring
code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood for interval-truncated Siler
    b2[i] <- beta0.b2 + (beta1.b2 * inbreed[i])
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1, a2, b1, b2[i], c)
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  beta0.b2 ~ dexp(1)
  beta1.b2 ~ dnorm(0, 100)
  c ~ dexp(1)
  mean.p ~ dunif(0, 1)
  
})


## set up other components of model
consts <- list(nind = nind, tM = tM)
data <- list(y = y, cint = cint, 
             censored = censored, tD = tD, dind = dind, inbreed = inbreed)

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


initFn <- function(cint, censored, inbreeding) {
  ## get ML estimates as initial values
  optFn <- function(pars, t) {
    if(any(pars < 0)) {
      return(NA)
    }
    sum(dSiler(t, a1 = pars[1], a2 = pars[2], b1 = pars[3], b2 = (pars[4] + (pars[5] * inbreeding)), c = pars[6], log = TRUE))
  }
  pars <- list(convergence = 1)
  k <- 0
  while(pars$convergence != 0 & k < 20) {
    ## sample missing values
    tD <- tinitFn(cint, censored)
    pars <- optim(rexp(6, 10), optFn, t = tD, control = list(fnscale = -1))
    k <- k + 1
  }
  if(k == 20) {
    stop("Can't sample initial values")
  }
  pars <- pars$par
  list(
    tD = tD,
    a1 = pars[1], 
    a2 = pars[2], 
    b1 = pars[3],
    beta0.b2 = pars[4],
    beta1.b2 = pars[5],
    c = pars[6],
    mean.p = runif(1, 0, 1)
  )
}

inits <- initFn(cint, censored, inbreeding = inbreed)

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = inits)

## compile the model
cModel <- compileNimble(model, showCompilerOutput = TRUE)

## try with adaptive slice sampler
config <- configureMCMC(cModel, monitors = c("a1", "a2", "b1", "beta0.b2", "beta1.b2", "c", "mean.p"), thin = 1)
config$removeSamplers(c("a1", "b1", "c", "beta0.b2", "beta1.b2", "a2"))
config$addSampler(target = c("a1", "b1", "c"), type = 'AF_slice')
config$addSampler(target = c("a2", "beta0.b2", "beta1.b2"), type = 'AF_slice')

#Check monitors and samplers
config$printMonitors()
config$printSamplers(c("a1", "a2", "b1", "b2", "c"))

#Build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

#Run the model
system.time(runAF <- runMCMC(cbuilt,  
                             niter = 100000, 
                             nburnin = 4000, 
                             nchains = 2, 
                             progressBar = TRUE, 
                             summary = TRUE, 
                             samplesAsCodaMCMC = TRUE, 
                             thin = 1))

runAF$summary
plot(runAF$samples)
NPrun <- runAF
saveRDS(NPrun, file = "Data/NPrun_regression.rds")

#####################################################################
##                                                                 ##
##                  Positive badgers                               ##
##                                                                 ##
#####################################################################

## read data:
CH.all <- as.data.table(readRDS("Data/CHCVall.rds"))
#CH <- CH.all

## select group
levels(CH.all$status)
levels(CH.all$status.cub)
#CH <- CH.all[CH.all$status == 'NPF' | CH.all$status == 'NPM'] #never positive badgers
#CH <- CH.all[CH.all$status.cub == 'CPF' | CH.all$status.cub == 'CPM'] #cub positive
CH <- CH.all[CH.all$status == 'PF' | CH.all$status == 'PM'] # positive
inbreed <- CH$f_inbreed

## extract death and birth times
tKD <- CH$death 
tB <- CH$birth

## remove extra variable from CH
colnames(CH)
CH <- CH[,31:183]

## extract max possible death time
tM <- ifelse(is.na(tKD), ncol(CH), tKD)

## extract last alive time
tL <- as.numeric(apply(CH, 1, function(x) max(which(x == 1))))

## some checks
stopifnot(all(tL > tB))
stopifnot(all(tKD[!is.na(tKD)] > tL[!is.na(tKD)]))
stopifnot(all(tM >= tL))

## normalise to survival times
## (necessary at the moment due to censoring
## constraints)
tM <- tM - tB
tKD <- tKD - tB
tL <- tL - tB

## define censoring matrices
cint <- cbind(tL, tKD)
cint[is.na(tKD), 2] <- cint[is.na(tKD), 1]
cint[is.na(tKD), 1] <- 0
colnames(cint) <- NULL
censored <- ifelse(!is.na(tKD), 1, 2)
tD <- rep(NA, length(tKD))
dind <- rep(1, length(tKD))

## extract number of captures
y <- rowSums(CH)
names(y) <- NULL

## some checks
stopifnot(all(tM >= y))

## set up nind
nind <- length(y)

## code for NIMBLE model with censoring
code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood for interval-truncated Siler
    b2[i] <- beta0.b2 + (beta1.b2 * inbreed[i])
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1, a2, b1, b2[i], c)
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  beta0.b2 ~ dexp(1)
  beta1.b2 ~ dnorm(0, 100)
  c ~ dexp(1)
  mean.p ~ dunif(0, 1)
  
})


## set up other components of model
consts <- list(nind = nind, tM = tM)
data <- list(y = y, cint = cint, 
             censored = censored, tD = tD, dind = dind, inbreed = inbreed)

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


initFn <- function(cint, censored, inbreeding) {
  ## get ML estimates as initial values
  optFn <- function(pars, t) {
    if(any(pars < 0)) {
      return(NA)
    }
    sum(dSiler(t, a1 = pars[1], a2 = pars[2], b1 = pars[3], b2 = (pars[4] + (pars[5] * inbreeding)), c = pars[6], log = TRUE))
  }
  pars <- list(convergence = 1)
  k <- 0
  while(pars$convergence != 0 & k < 20) {
    ## sample missing values
    tD <- tinitFn(cint, censored)
    pars <- optim(rexp(6, 10), optFn, t = tD, control = list(fnscale = -1))
    k <- k + 1
  }
  if(k == 20) {
    stop("Can't sample initial values")
  }
  pars <- pars$par
  list(
    tD = tD,
    a1 = pars[1], 
    a2 = pars[2], 
    b1 = pars[3],
    beta0.b2 = pars[4],
    beta1.b2 = pars[5],
    c = pars[6],
    mean.p = runif(1, 0, 1)
  )
}

inits <- initFn(cint, censored, inbreeding = inbreed)

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = inits)

## compile the model
cModel <- compileNimble(model, showCompilerOutput = TRUE)

## try with adaptive slice sampler
config <- configureMCMC(cModel, monitors = c("a1", "a2", "b1", "beta0.b2", "beta1.b2", "c", "mean.p"), thin = 1)
config$removeSamplers(c("a1", "b1", "c", "beta0.b2", "beta1.b2", "a2"))
config$addSampler(target = c("a1", "b1", "c"), type = 'AF_slice')
config$addSampler(target = c("a2", "beta0.b2", "beta1.b2"), type = 'AF_slice')

#Check monitors and samplers
config$printMonitors()
config$printSamplers(c("a1", "a2", "b1", "b2", "c"))

#Build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

#Run the model
system.time(runAF <- runMCMC(cbuilt,  
                             niter = 100000, 
                             nburnin = 4000, 
                             nchains = 2, 
                             progressBar = TRUE, 
                             summary = TRUE, 
                             samplesAsCodaMCMC = TRUE, 
                             thin = 1))

runAF$summary
plot(runAF$samples)
Prun <- runAF
saveRDS(Prun, file = "Data/Prun_regression.rds")

#####################################################################
##                                                                 ##
##                Cub positive badgers                             ##
##                                                                 ##
#####################################################################

## read data:
CH.all <- as.data.table(readRDS("Data/CHCVall.rds"))
#CH <- CH.all

## select group
levels(CH.all$status)
levels(CH.all$status.cub)
#CH <- CH.all[CH.all$status == 'NPF' | CH.all$status == 'NPM'] #never positive badgers
CH <- CH.all[CH.all$status.cub == 'CPF' | CH.all$status.cub == 'CPM'] #cub positive
#CH <- CH.all[CH.all$status == 'PF' | CH.all$status == 'PM'] # positive
inbreed <- CH$f_inbreed

## extract death and birth times
tKD <- CH$death 
tB <- CH$birth

## remove extra variable from CH
colnames(CH)
CH <- CH[,31:183]

## extract max possible death time
tM <- ifelse(is.na(tKD), ncol(CH), tKD)

## extract last alive time
tL <- as.numeric(apply(CH, 1, function(x) max(which(x == 1))))

## some checks
stopifnot(all(tL > tB))
stopifnot(all(tKD[!is.na(tKD)] > tL[!is.na(tKD)]))
stopifnot(all(tM >= tL))

## normalise to survival times
## (necessary at the moment due to censoring
## constraints)
tM <- tM - tB
tKD <- tKD - tB
tL <- tL - tB

## define censoring matrices
cint <- cbind(tL, tKD)
cint[is.na(tKD), 2] <- cint[is.na(tKD), 1]
cint[is.na(tKD), 1] <- 0
colnames(cint) <- NULL
censored <- ifelse(!is.na(tKD), 1, 2)
tD <- rep(NA, length(tKD))
dind <- rep(1, length(tKD))

## extract number of captures
y <- rowSums(CH)
names(y) <- NULL

## some checks
stopifnot(all(tM >= y))

## set up nind
nind <- length(y)

## code for NIMBLE model with censoring
code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood for interval-truncated Siler
    b2[i] <- beta0.b2 + (beta1.b2 * inbreed[i])
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1, a2, b1, b2[i], c)
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  beta0.b2 ~ dexp(1)
  beta1.b2 ~ dnorm(0, 100)
  c ~ dexp(1)
  mean.p ~ dunif(0, 1)
  
})


## set up other components of model
consts <- list(nind = nind, tM = tM)
data <- list(y = y, cint = cint, 
             censored = censored, tD = tD, dind = dind, inbreed = inbreed)

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


initFn <- function(cint, censored, inbreeding) {
  ## get ML estimates as initial values
  optFn <- function(pars, t) {
    if(any(pars < 0)) {
      return(NA)
    }
    sum(dSiler(t, a1 = pars[1], a2 = pars[2], b1 = pars[3], b2 = (pars[4] + (pars[5] * inbreeding)), c = pars[6], log = TRUE))
  }
  pars <- list(convergence = 1)
  k <- 0
  while(pars$convergence != 0 & k < 20) {
    ## sample missing values
    tD <- tinitFn(cint, censored)
    pars <- optim(rexp(6, 10), optFn, t = tD, control = list(fnscale = -1))
    k <- k + 1
  }
  if(k == 20) {
    stop("Can't sample initial values")
  }
  pars <- pars$par
  list(
    tD = tD,
    a1 = pars[1], 
    a2 = pars[2], 
    b1 = pars[3],
    beta0.b2 = pars[4],
    beta1.b2 = pars[5],
    c = pars[6],
    mean.p = runif(1, 0, 1)
  )
}

inits <- initFn(cint, censored, inbreeding = inbreed)

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = inits)

## compile the model
cModel <- compileNimble(model, showCompilerOutput = TRUE)

## try with adaptive slice sampler
config <- configureMCMC(cModel, monitors = c("a1", "a2", "b1", "beta0.b2", "beta1.b2", "c", "mean.p"), thin = 1)
config$removeSamplers(c("a1", "b1", "c", "beta0.b2", "beta1.b2", "a2"))
config$addSampler(target = c("a1", "b1", "c"), type = 'AF_slice')
config$addSampler(target = c("a2", "beta0.b2", "beta1.b2"), type = 'AF_slice')

#Check monitors and samplers
config$printMonitors()
config$printSamplers(c("a1", "a2", "b1", "b2", "c"))

#Build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

#Run the model
system.time(runAF <- runMCMC(cbuilt,  
                             niter = 100000, 
                             nburnin = 4000, 
                             nchains = 2, 
                             progressBar = TRUE, 
                             summary = TRUE, 
                             samplesAsCodaMCMC = TRUE, 
                             thin = 1))

runAF$summary
plot(runAF$samples)
Prun <- runAF
saveRDS(Prun, file = "Data/Prun_regression.rds")

## Comparisons of NP and P badgers

NPrun$summary
Prun$summary
