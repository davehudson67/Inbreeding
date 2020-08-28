## load libraries
library(nimble)
library(coda)
library(mcmcplots)
library(lamW)

source("Distributions/Siler/Dist_SilerNim.R")
source("Distributions/Siler/Dist_Siler.R")

## read data:
CH.all <- as.data.table(readRDS("Data/CHCVall.rds"))
#CH <- CH.all

## select group
levels(CH.all$status)
levels(CH.all$status.cub)
CH <- CH.all[CH.all$status == 'PF' | CH.all$status == 'PM'] #positive badgers
#CH <- CH.all[CH.all$status.cub == 'CPF' | CH.all$status.cub == 'CPM'] #cub positive
#CH <- CH.all[CH.all$status == 'NPF' | CH.all$status.cub == 'NPM'] #never positive

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
l.MCMC <- runAF

#Plot mcmcm
#het.samples <- runAF$samples
#mid.samples <- runAF$samples
hom.samples <- h.MCMC$samples
het.samples <- l.MCMC$samples

#mcmcplot(hom.samples)
#png("traceAF%d.png")
plot(het.samples)
#dev.off()

#png("pairsAF%d.png")
#pairs(samples)
#dev.off()

## save samples
#saveRDS(samples, "newSiler.rds")

#Set age variable (quarter years)
x <- 0:60

## extract samples
l.samples <- as.matrix(het.samples)[, 1:5]
#m.samples <- as.matrix(mid.samples)[, 1:5]
h.samples <- as.matrix(hom.samples)[, 1:5]


#Siler Mortality rate
l.mort <- apply(l.samples, 1, function(pars, x) {
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
l.mort <- apply(l.mort, 1, function(x) {
    c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

m.mort <- apply(m.samples, 1, function(pars, x, inbreed) {
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
m.mort <- apply(m.mort, 1, function(x) {
  c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

h.mort <- apply(h.samples, 1, function(pars, x) {
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
h.mort <- apply(h.mort, 1, function(x) {
  c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

l.mort<-as.data.frame(t(l.mort))
l.mort$age<-seq(0:60)
l.mort$inb<-"het"
m.mort<-as.data.frame(t(m.mort))
m.mort$age<-seq(0:60)
m.mort$inb<-"mid"
h.mort<-as.data.frame(t(h.mort))
h.mort$age<-seq(0:60)
h.mort$inb<-"hom"

inb.samples<-bind_rows(l.mort, h.mort)

colnames(inb.samples)<-c("mean", "lwr", "upr","age", "inb")

ggplot(inb.samples, aes(x=age, y=mean, col=inb)) +
  geom_line() 
  geom_ribbon(data=inb.samples,aes(ymin=lwr,ymax=upr, fill=inb),alpha=0.3)

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


#Siler survival function
surv <- apply(samples, 1, function(pars, x) {
    ## extract pars
    a1 <- pars[1]
    a2 <- pars[2]
    b1 <- pars[3]
    b2 <- pars[4]
    c <- pars[5]
    
    ## return predictions
    exp(((exp(a1))/b1)*(exp(-b1*x)-1) - c*x + ((exp(a2))/b2)*(1-exp(b2*x)))
}, x = x)

## extract mean and 95% intervals
surv <- apply(surv, 1, function(x) {
    c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

## produce plots
pdf("survcurves.pdf", width = 10, height = 5)

par(mfrow = c(1, 2))

#Draw mortality curve
plot(x, h.mort[1, ], type = 'l', main = "Hazard function")
lines(x, l.mort[2, ], lty = 2)
lines(x, l.mort[3, ], lty = 2)

#Draw survival curve
plot(x, surv[1, ], type = "l", main = "Survivor function")
lines(x, surv[2, ], lty = 2)
lines(x, surv[3, ], lty = 2)

dev.off()

