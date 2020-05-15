## load libraries
library(nimble)
library(coda)
library(mcmcplots)

source("/home2/ISAD/dh448/Desktop/Simulation study/Custom_Siler.R")

getwd()

#reset columns to start from 0
colnames(CP.CJS.array)<-seq(1:131)
colnames(CP.age.array)<-seq(1:131)

#inbreeding
inbreed<- CP.CVdata$Inbreed

## read in data
CH <- CP.CJS.array
tB <- CP.CVdata$birth_yr-9
names(tB) <- NULL
tKD <- CP.dead-9
age <- CP.age.array

#if death date is unknown change to NA
tKD[tKD<0]<-NA

## extract max possible death time
tM <- ifelse(is.na(tKD), ncol(CH), tKD)

## extract last alive time
tL <- apply(CH, 1, function(x) max(which(x == 1)))
names(tL) <- NULL

## normalise to survival times
## (necessary at the moment due to censoring
## constraints)
tM <- tM - tB
tL <- tL - tB
tKD <- tKD - tB

## define censoring matrices
cint <- cbind(tL, tM)
censored <- ifelse(!is.na(tKD), 1, 2)
tD <- rep(NA, length(tKD))
dind <- rep(1, length(tKD))

## extract number of captures
y <- apply(CH, 1, sum)
names(y) <- NULL

## set up nind
nind <- length(y)

## set up initial values
tinit <- apply(cint, 1, function(x) {
    runif(1, x[1], x[2])
})

## custom distribution

## probability density function
dsiler <- nimbleFunction(
    run = function(x = double(0), a1 = double(0),
        a2 = double(0), b1 = double(0), b2 = double(0), 
        c = double(0), 
        log = integer(0, default = 0)) {
        returnType(double(0))
        logS <- (exp(a1) / b1) * (exp(-b1 * x) - 1) - c * x + (exp(a2) / b2) * (1 - exp(b2 * x))
        logH <- log(exp(a1 - b1 * x) + c + exp(a2 + b2 * x))
        logProb <- logH + logS
        if(log) return(logProb)
        else return(exp(logProb)) 
    })

## function to produce random samples
rsiler <- nimbleFunction(
  run = function(n = integer(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c = double(0)) {
    returnType(double(0))
    if(n != 1) print("rsiler only allows n = 1; using n = 1.")
    print("No sampler implemented yet")
    return(NaN)
  })

## cumulative distribution function (and survivor function)
psiler <- nimbleFunction(
    run = function(q = double(0), a1 = double(0),
        a2 = double(0), b1 = double(0), b2 = double(0), 
        c = double(0), 
        lower.tail = integer(0, default = 1), 
        log.p = integer(0, default = 0)) {
        returnType(double(0))
        logS <- (exp(a1) / b1) * (exp(-b1 * q) - 1) - c * q + (exp(a2) / b2) * (1 - exp(b2 * q))
        if(!lower.tail) { 
            if(log.p) return(logS)
            else return(exp(logS))
        } else {
            p <- 1 - exp(logS)
            if(!log.p) return(p)
            else return(log(p))
        }
    })
    
### functions to optimise 
#qsilerObjective <- nimbleFunction(
#    run = function(par = double(1)) {
#        return(psiler(par, 0, 1))
#        returnType(double(0))
#    }
#)
#optimizer <- nimbleFunction(
#    run = function(method = character(0), fnscale = double(0)) {
#        control <- optimDefaultControl()
#        control$fnscale <- fnscale
#        par <- c(0.1, -0.1)
#        return(optim(par, objectiveFunction, method = method, control = control))
#        returnType(optimResultNimbleList())
#    }
#)
#cOptimizer <- compileNimble(optimizer)
#cOptimizer(method = 'BFGS', fnscale = -1)

## quantile function (not yet working)
qsiler <- nimbleFunction(
    run = function(p = double(0), a1 = double(0),
        a2 = double(0), b1 = double(0), b2 = double(0), 
        c = double(0),
        lower.tail = integer(0, default = 1), 
        log.p = integer(0, default = 0)) {
        returnType(double(0))
        if(log.p) p <- exp(p)
        if(!lower.tail) p <- 1 - p
        print("qsiler() not specified")
        return(NaN)
    })

## register distributions with NIMBLE
registerDistributions(list(
    dsiler = list(
        BUGSdist = "dsiler(a1, a2, b1, b2, c)",
        Rdist = "dsiler(a1, a2, b1, b2, c)",
        pqAvail = TRUE, 
        range = c(0, Inf)
        )
    ))

## code for NIMBLE model with censoring
CJS.code <- nimbleCode({

    ## survival components for dead badgers
    for (i in 1:nind) {
        ## likelihood for interval-truncated Siler
        censored[i] ~ dinterval(tD[i], cint[i, ])
        tD[i] ~ dsiler(a1, a2, b1, b2, c)
        
        ## sampling component
        pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), cint[i, 2]) - y[i]) * log(1 - mean.p))
        dind[i] ~ dbern(pd[i])
        
        b2 <- beta0.b2 + (beta1.b2 *  inbreed[i])
        
    }

    ## priors
    a1 ~ dnorm(-2, 0.01)
    a2 ~ dnorm(-3, 0.01)
    b1 ~ T(dnorm(0.01, 0.01), 0, )
    c ~ T(dnorm(0, 0.01), 0, )
    mean.p ~ dunif(0, 1)
    
    beta0.b2 ~ dnorm(0,0.001)
    beta1.b2 ~ dnorm(0,0.001)
})

## set up other components of model
CJS.Consts <- list(nind = nind)
CJS.data <- list(y = y, cint = cint, 
    censored = censored, tD = tD, dind = dind, inbreed=inbreed)
CJS.inits <- list(
    tD = tinit,
    a1 = -2, 
    a2 = -3, 
    b1 = 0.01,
    beta0.b2 = 0.01,
    beta1.b2 = 0.01,
    mean.p = runif(1, 0, 1),
    c = 0.05
)

## define the model, data, inits and constants
CJSModel <- nimbleModel(code = CJS.code, constants = CJS.Consts, data = CJS.data, inits = CJS.inits, name = "CJS")

## compile the model
cCJSModel <- compileNimble(CJSModel, showCompilerOutput = TRUE)

## try with adaptive slice sampler
CJSconfig <- configureMCMC(cCJSModel, monitors = c("a1", "a2", "b1", "beta0.b2", "beta1.b2", "c", "mean.p"), thin = 1)
CJSconfig$removeSamplers(c("a1", "a2", "b1", "c","beta0.b2", "beta1.b2"))
CJSconfig$addSampler(target = c("a1", "a2", "b1", "c", "beta0.b2", "beta1.b2"), type = 'AF_slice')

#Check monitors and samplers
CJSconfig$printMonitors()
CJSconfig$printSamplers(c("a1", "a2", "b1", "beta0.b2", "beta1.b2", "c"))

#Build the model
CJSbuilt <- buildMCMC(CJSconfig)
cCJSbuilt <- compileNimble(CJSbuilt)

#Run the model
system.time(runAF <- runMCMC(cCJSbuilt, 
    niter = 100000, 
    nburnin = 2000, 
    nchains = 2, 
    progressBar = TRUE, 
    summary = TRUE, 
    samplesAsCodaMCMC = TRUE, 
    thin = 1))
runAF$summary

#Plot mcmcm
samples.inbreed.uninf <- runAF$samples
mcmcplot(samples.inbreed.uninf)
png("traceAF%d.png")
plot(samples)
dev.off()

png("pairsAF%d.png")
pairs(samples)
dev.off()

## save samples
saveRDS(samples, "newSiler.rds")

#Set age variable (quarter years)
x <- 0:80

## extract samples
samples <- as.matrix(samples)[, 1:5]

#Siler Mortality rate
mort <- apply(samples, 1, function(pars, x) {
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
mort <- apply(mort, 1, function(x) {
    c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

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
plot(x, mort[1, ], type = 'l', main = "Hazard function")
lines(x, mort[2, ], lty = 2)
lines(x, mort[3, ], lty = 2)

#Draw survival curve
plot(x, surv[1, ], type = "l", main = "Survivor function")
lines(x, surv[2, ], lty = 2)
lines(x, surv[3, ], lty = 2)

dev.off()

