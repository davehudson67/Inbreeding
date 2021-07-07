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
CH.all <- as.data.table(readRDS("Data/CHCVall.rds"))

#####################################################################
##                                                                 ##
##            Inbred badgers                                       ##
##                                                                 ##
#####################################################################

## select badgers 
#CH <- CH.all[CH.all$status == 'NPF' | CH.all$status == 'NPM'] #never positive badgers
#CH <- CH.all[CH.all$status.cub == 'CPF' | CH.all$status.cub == 'CPM'] #cub positive
CH <- CH.all[CH.all$status == 'PF' | CH.all$status == 'PM'] # positive

## select group
levels(CH$Inb) # 1 = homozygous; 2 = heterozygous

CH <- CH[CH$Inb == 2] #more inbred badgers
#CH <- CH.all[CH.all$Inb == 1] #more outbred badgers

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
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1, a2, b1, b2, c)
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  b2 ~ dexp(1)
  c ~ dexp(1)
  mean.p ~ dunif(0, 1)
  
})


## set up other components of model
consts <- list(nind = nind, tM = tM)
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
    if(any(pars < 0)) {
      return(NA)
    }
    sum(dSiler(t, a1 = pars[1], a2 = pars[2], b1 = pars[3], b2 = pars[4], c = pars[5], log = TRUE))
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
  list(
    tD = tD,
    a1 = pars[1], 
    a2 = pars[2], 
    b1 = pars[3],
    b2 = pars[4],
    c = pars[5],
    mean.p = runif(1, 0, 1)
  )
}

inits <- initFn(cint, censored)

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = inits)

## compile the model
cModel <- compileNimble(model, showCompilerOutput = TRUE)

## try with adaptive slice sampler
config <- configureMCMC(cModel, monitors = c("a1", "a2", "b1", "b2", "c", "mean.p"), thin = 1)
config$removeSamplers(c("a1", "b1", "c", "b2", "a2"))
config$addSampler(target = c("a1", "b1", "c"), type = 'AF_slice')
config$addSampler(target = c("a2", "b2"), type = 'AF_slice')

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
Inb2_run <- runAF
#saveRDS(Inb2_run, file = "Data/Inb2_cat.rds")

##########################################################################################################################

## read data:
CH.all <- as.data.table(readRDS("Data/CHCVall.rds"))

#####################################################################
##                                                                 ##
##            Outbred badgers                                      ##
##                                                                 ##
#####################################################################

## select badgers 
#CH <- CH.all[CH.all$status == 'NPF' | CH.all$status == 'NPM'] #never positive badgers
#CH <- CH.all[CH.all$status.cub == 'CPF' | CH.all$status.cub == 'CPM'] #cub positive
CH <- CH.all[CH.all$status == 'PF' | CH.all$status == 'PM'] # positive

## select group
levels(CH.all$Inb) # 1 = homozygous; 2 = heterozygous

#CH <- CH.all[CH.all$Inb == 2] #more inbred badgers
CH <- CH.all[CH.all$Inb == 1] #more outbred badgers

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
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1, a2, b1, b2, c)
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  b2 ~ dexp(1)
  c ~ dexp(1)
  mean.p ~ dunif(0, 1)
  
})


## set up other components of model
consts <- list(nind = nind, tM = tM)
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
    if(any(pars < 0)) {
      return(NA)
    }
    sum(dSiler(t, a1 = pars[1], a2 = pars[2], b1 = pars[3], b2 = pars[4], c = pars[5], log = TRUE))
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
  list(
    tD = tD,
    a1 = pars[1], 
    a2 = pars[2], 
    b1 = pars[3],
    b2 = pars[4],
    c = pars[5],
    mean.p = runif(1, 0, 1)
  )
}

inits <- initFn(cint, censored)

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = inits)

## compile the model
cModel <- compileNimble(model, showCompilerOutput = TRUE)

## try with adaptive slice sampler
config <- configureMCMC(cModel, monitors = c("a1", "a2", "b1", "b2", "c", "mean.p"), thin = 1)
config$removeSamplers(c("a1", "b1", "c", "b2", "a2"))
config$addSampler(target = c("a1", "b1", "c"), type = 'AF_slice')
config$addSampler(target = c("a2", "b2"), type = 'AF_slice')

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
Inb1_run <- runAF
#saveRDS(Inb1_run, file = "Data/Inb1_cat.rds")


###########################################################################################################################

## Plot survival and mortality curves

hom.samples <- Inb2_run$samples
het.samples <- Inb1_run$samples

## Set age variable (quarter years)
x <- 0:60

## extract samples
het.samples <- as.matrix(het.samples)[, 1:5]
hom.samples <- as.matrix(hom.samples)[, 1:5]

#Siler Mortality rate
het.mort <- apply(het.samples, 1, function(pars, x) {
  ## extract pars
  a1 <- pars[1]
  a2 <- pars[2]
  b1 <- pars[3]
  b2 <- pars[4]
  c <- pars[5]
  
  ## return predictions
  exp(a1-(b1*x)) + c + exp(a2+(b2*x))
}, x = x)

## extract mean and 95% intervals
het.mort <- apply(het.mort, 1, function(x) {
  c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

hom.mort <- apply(hom.samples, 1, function(pars, x) {
  ## extract pars
  a1 <- pars[1]
  a2 <- pars[2]
  b1 <- pars[3]
  b2 <- pars[4]
  c <- pars[5]
  
  ## return predictions
  exp(a1-(b1*x)) + c + exp(a2+(b2*x))
}, x = x)

## extract mean and 95% intervals
hom.mort <- apply(hom.mort, 1, function(x) {
  c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

het.mort<-as.data.frame(t(het.mort))
het.mort$age<-seq(0:60)
het.mort$inb<-"heterozygous"
hom.mort<-as.data.frame(t(hom.mort))
hom.mort$age<-seq(0:60)
hom.mort$inb<-"homozygous"

inb.samples<-bind_rows(het.mort, hom.mort)
colnames(inb.samples)<-c("Mean", "lwr", "upr","Age", "Inbreeding_category")

saveRDS(inb.samples, file = "Data/inb_samplesPm.rds")

ggplot(inb.samples, aes(x=Age, y=Mean, col=Inbreeding_category)) +
  geom_line() +
geom_ribbon(data=inb.samples,aes(ymin=lwr,ymax=upr, fill=Inbreeding_category),alpha=0.3)


##########################################################################################################################

#Siler survival function
het.surv <- apply(het.samples, 1, function(pars, x) {
  ## extract pars
  a1 <- pars[1]
  a2 <- pars[2]
  b1 <- pars[3]
  b2 <- pars[4]
  c <- pars[5]
  
  ## return predictions
  pSiler(q=x, a1, a2, b1, b2, c, lower.tail = FALSE)
}, x = x)

## extract mean and 95% intervals
het.surv <- apply(het.surv, 1, function(x) {
  c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

hom.surv <- apply(hom.samples, 1, function(pars, x) {
  ## extract pars
  a1 <- pars[1]
  a2 <- pars[2]
  b1 <- pars[3]
  b2 <- pars[4]
  c <- pars[5]
  
  ## return predictions
  pSiler(q=x, a1, a2, b1, b2, c, lower.tail = FALSE)
}, x = x)

## extract mean and 95% intervals
hom.surv <- apply(hom.surv, 1, function(x) {
  c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

het.surv<-as.data.frame(t(het.surv))
het.surv$age<-seq(0:60)
het.surv$inb<-"heterozygous"
hom.surv<-as.data.frame(t(hom.surv))
hom.surv$age<-seq(0:60)
hom.surv$inb<-"homozygous"

inb.samples<-bind_rows(het.surv, hom.surv)
colnames(inb.samples)<-c("Survival_Probability", "lwr", "upr","Age", "Inbreeding_category")

saveRDS(inb.samples, file = "Data/inb_samplesPs.rds")

ggplot(inb.samples, aes(x=Age, y=Survival_Probability, col=Inbreeding_category)) +
  geom_line() +
  geom_ribbon(data=inb.samples,aes(ymin=lwr,ymax=upr, fill=Inbreeding_category),alpha=0.3)

##########################################################################################################################

#Siler pdf
het.pdf <- apply(het.samples, 1, function(pars, x) {
  ## extract pars
  a1 <- pars[1]
  a2 <- pars[2]
  b1 <- pars[3]
  b2 <- pars[4]
  c <- pars[5]
  
  ## return predictions
  dSiler(x, a1, a2, b1, b2, c)
}, x = x)

## extract mean and 95% intervals
het.pdf <- apply(het.pdf, 1, function(x) {
  c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

hom.pdf <- apply(hom.samples, 1, function(pars, x) {
  ## extract pars
  a1 <- pars[1]
  a2 <- pars[2]
  b1 <- pars[3]
  b2 <- pars[4]
  c <- pars[5]
  
  ## return predictions
  dSiler(x, a1, a2, b1, b2, c)
}, x = x)

## extract mean and 95% intervals
hom.pdf <- apply(hom.pdf, 1, function(x) {
  c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

het.pdf<-as.data.frame(t(het.pdf))
het.pdf$age<-seq(0:60)
het.pdf$inb<-"heterozygous"
hom.pdf<-as.data.frame(t(hom.pdf))
hom.pdf$age<-seq(0:60)
hom.pdf$inb<-"homozygous"

inb.samples<-bind_rows(het.pdf, hom.pdf)
colnames(inb.samples)<-c("Mean", "lwr", "upr","Age", "Inbreeding_category")

saveRDS(inb.samples, file = "Data/inb_samplesPPDF.rds")

ggplot(inb.samples, aes(x=Age, y=Mean, col=Inbreeding_category)) +
  geom_line() +
  geom_ribbon(data=inb.samples,aes(ymin=lwr, ymax=upr, fill=Inbreeding_category),alpha=0.3)

#########################################################################################################################

## KL discrepancies between pdf of low and hihg inbreeding estimates
#library(LaplacesDemon)

h.a1<-h.MCMC$summary$all.chains[1,1]
h.a2<-h.MCMC$summary$all.chains[2,1]
h.b1<-h.MCMC$summary$all.chains[3,1]
h.b2<-h.MCMC$summary$all.chains[4,1]
h.c<-h.MCMC$summary$all.chains[5,1]

l.a1<-l.MCMC$summary$all.chains[1,1]
l.a2<-l.MCMC$summary$all.chains[2,1]
l.b1<-l.MCMC$summary$all.chains[3,1]
l.b2<-l.MCMC$summary$all.chains[4,1]
l.c<-l.MCMC$summary$all.chains[5,1]


low.dist <- dSiler(1:80, l.a1, l.a2, l.b1, l.b2, l.c)
high.dist <- dSiler(1:80, h.a1, h.a2, h.b1, h.b2, h.c)
K <- KLD(low.dist, high.dist)
k <- K$intrinsic.discrepancy

#McCulloch adjustment
q <- (1 + (1 - exp(-2 * k)) ^ 0.5) / 2 
q

## KL discrepancies between parameters inbreeding estimates
## check a1
K <- KLD(l.samples[,1], h.samples[,1])
k <- K$intrinsic.discrepancy
#McCulloch adjustment
q <- (1 + (1 - exp(-2 * k)) ^ 0.5) / 2 
q

## check a2
K <- KLD(l.samples[,2], h.samples[,2])
k <- K$intrinsic.discrepancy
#McCulloch adjustment
q <- (1 + (1 - exp(-2 * k)) ^ 0.5) / 2 
q

## check b1
K <- KLD(l.samples[,3], h.samples[,3])
k <- K$intrinsic.discrepancy
#McCulloch adjustment
q <- (1 + (1 - exp(-2 * k)) ^ 0.5) / 2 
q

## check b2
K <- KLD(l.samples[,4], h.samples[,4])
k <- K$intrinsic.discrepancy
#McCulloch adjustment
q <- (1 + (1 - exp(-2 * k)) ^ 0.5) / 2 
q

## check c
K <- KLD(l.samples[,5], h.samples[,5])
k <- K$intrinsic.discrepancy
#McCulloch adjustment
q <- (1 + (1 - exp(-2 * k)) ^ 0.5) / 2 
q
