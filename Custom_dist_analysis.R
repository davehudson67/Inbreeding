## load libraries
library(nimble)
library(coda)
library(mcmcplots)
library(lamW)

source("Dist_Siler.R")

#read data:
CH<-low
tKD<-low.dead
tB <- low.birth

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
y <- apply(CH, 1, sum)
names(y) <- NULL

## some checks
stopifnot(all(tM >= y))

## set up nind
nind <- length(y)

## code for NIMBLE model with censoring
CJS.code <- nimbleCode({

    ## survival components for dead badgers
    for (i in 1:nind) {
      
        ## likelihood for interval-truncated Siler
        censored[i] ~ dinterval(tD[i], cint[i, ])
        tD[i] ~ dsiler(a1, a2, b1, b2, c)
        
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
CJS.Consts <- list(nind = nind, tM = tM)
CJS.data <- list(y = y, cint = cint, 
                 censored = censored, tD = tD, dind = dind)
tinit <- apply(cbind(cint, censored), 1, function(x) {
  if(x[3] == 2) {
    y <- x[2] + rexp(1, 0.1)
  } else {
    y <- runif(1, x[1], x[2])
  }
  y
})
CJS.inits <- list(
  tD = tinit,
  a1 = 0.1, 
  a2 = 0.1, 
  b1 = 0.1,
  b2 = 0.1,
  mean.p = 0.5,
  c = 0.05
)

## define the model, data, inits and constants
CJSModel <- nimbleModel(code = CJS.code, constants = CJS.Consts, data = CJS.data, inits = CJS.inits, name = "CJS")

## compile the model
cCJSModel <- compileNimble(CJSModel, showCompilerOutput = TRUE)

## try with adaptive slice sampler
CJSconfig <- configureMCMC(cCJSModel, monitors = c("a1", "a2", "b1", "b2", "c", "mean.p"), thin = 1)
CJSconfig$removeSamplers(c("a1", "b1", "c", "b2", "a2"))
CJSconfig$addSampler(target = c("a1", "b1", "c"), type = 'AF_slice')
CJSconfig$addSampler(target = c("a2", "b2"), type = 'AF_slice')


#Check monitors and samplers
CJSconfig$printMonitors()
CJSconfig$printSamplers(c("a1", "a2", "b1", "b2", "c"))

#Build the model
CJSbuilt <- buildMCMC(CJSconfig)
cCJSbuilt <- compileNimble(CJSbuilt)

#Run the model
system.time(runAF <- runMCMC(cCJSbuilt,  
    niter = 100000, 
    nburnin = 4000, 
    nchains = 2, 
    progressBar = TRUE, 
    summary = TRUE, 
    samplesAsCodaMCMC = TRUE, 
    thin = 1))
runAF$summary

#Plot mcmcm
hom.samples <- runAF$samples
mcmcplot(hom.samples)
#png("traceAF%d.png")
#plot(samples)
#dev.off()

#png("pairsAF%d.png")
#pairs(samples)
#dev.off()

## save samples
#saveRDS(samples, "newSiler.rds")

#Set age variable (quarter years)
x <- 0:80

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
l.mort$age<-seq(0:80)
l.mort$inb<-"het"
m.mort<-as.data.frame(t(m.mort))
m.mort$age<-seq(0:80)
m.mort$inb<-"mid"
h.mort<-as.data.frame(t(h.mort))
h.mort$age<-seq(0:80)
h.mort$inb<-"hom"

inb.samples<-bind_rows(l.mort, h.mort)

colnames(inb.samples)<-c("mean", "lwr", "upr","age", "inb")

ggplot(inb.samples, aes(x=age, y=mean, col=inb)) +
  geom_line() +
  geom_ribbon(data=inb.samples,aes(ymin=lwr,ymax=upr, fill=inb),alpha=0.3)

hist(gene$mlh)

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

