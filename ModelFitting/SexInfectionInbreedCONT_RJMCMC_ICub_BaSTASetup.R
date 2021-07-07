##########################################################
##                                                      ###
##       Linear regression of Inb against b2            ###
##                                                      ###
###########################################################

rm(list=ls())
#load("Data/badgerSexInb_FullLifeInfvUninf.RData")
load("Data/badgerSexInb_ICubvUninf.RData")
set.seed(seeds[12])

## set up plot output file
pdf("outputs/SexInfectionInbreed_CONT_ICub_BaSTASetup.pdf")

code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood for interval-truncated siler
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1 * a1mult[i], a2, b1 * b1mult[i], b2 * b2mult[i], c1 * c1mult[i])
    
    log(a1mult[i]) <- beta[1] * infection[i] 

#    log(a2mult[i]) <- beta[4] * sex[i] + 
#      beta[5] * infection[i]  + 
#      beta[6] * sex[i] * infection[i]
    
    log(b1mult[i]) <- beta[2] * sex[i] + 
      beta[3] * infection[i] + 
      beta[4] * sex[i] * infection[i]
    
    log(b2mult[i]) <- betab2[1] * inbr[i] * z[1] + 
      betab2[2] * infection[i] * z[2] + 
      betab2[3] * sex[i] * z[3] +
      betab2[4] * sex[i] * inbr[i] * z[4] +
      betab2[5] * infection[i] * inbr[i] * z[5]
    
    log(c1mult[i]) <- beta[5] * sex[i] + 
      beta[6] * infection[i] + 
      beta[7] * sex[i] * infection[i]
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  for (j in 1:7){
    beta[j] ~ dnorm(0, sd = 1)
  }
  
  for (k in 1:5){
    betab2[k] ~ dnorm(0, sd = 1)
    z[k] ~ dbern(0.5)
  }
  
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  b2 ~ dexp(1)
  c1 ~ dexp(1)
  mean.p ~ dunif(0, 1)
  constraint_dataSexInb4 ~ dconstraint(z[4] <= z[1] * z[3])
  constraint_dataInfInb4 ~ dconstraint(z[5] <= z[1] * z[2])
  
})


## set up data
consts <- list(nind = nind, tM = tM, sex = sex, infection = infection, inbr = inbr)

data <- list(y = y, cint = cint, 
             censored = censored, tD = tD, dind = dind, 
             constraint_dataSexInb4 = 1, constraint_dataInfInb4 = 1)

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
    if(any(pars[c(1:5)] < 0)) {
      return(NA)
    }
    
    sum(dSiler(t, a1 = pars[1], a2 = pars[2], b1 = pars[3], b2 = pars[4] * exp(pars[6]), c1 = pars[5], log = TRUE))
    
  }
  pars <- list(convergence = 1)
  k <- 0
  while(pars$convergence != 0 & k < 50) {
    ## sample missing values
    tD <- tinitFn(cint, censored)
    
    ## optimise to interval-censored only
    pars <- optim(c(rexp(5, 10), rnorm(1, 0, 1)), optFn, t = tD, control = list(fnscale = -1))
    k <- k + 1
  }
  if(k == 50) {
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
    beta = rep(pars[6], times = 7),
    betab2 = rep(pars[6], times = 5),
    z = rep(0, times = 5)
  )
}


## build the model
model <- nimbleModel(code, constants = consts,
                     data = data, 
                     inits = initFn(cint, censored))

## configure model
config <- configureMCMC(model)

config$removeSamplers(c("a1", "a2", "b1", "b2", "c1", "beta"))
config$addSampler(target = c("a1", "a2", "b1", "b2", "c1", "beta"), type = 'AF_slice', control = 50)
config$addSampler(target = c("a1", "b1", "c1"), type = 'AF_slice', control = 20)
config$addSampler(target = c("b2"), type = 'slice')

## Add reversible jump
configureRJ(conf = config,   ## model configuration
            targetNodes = c("betab2"),    ## coefficients for selection
            indicatorNodes = c("z"),   ## indicators paired with coefficients
            control = list(mean = 0.5, scale = .5))

config$addMonitors("betab2", "z")

config$printSamplers(c("betab2", "z", "a1", "a2", "b1", "b2", "c1", "beta", "mean.p"))

rIndicatorMCMC <- buildMCMC(config)
cIndicatorModel <- compileNimble(model)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

system.time(run <- runMCMC(cIndicatorMCMC, 
                           niter = 20000, 
                           nburnin = 5000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))


MCMCsummary(run$samples)
plot(run$samples)
#MCMCtrace(run$samples)
samples <- as.matrix(run$samples)

#samples <- samples[sample.int(nrow(samples), ceiling(nrow(samples) * 0.1)), ]
#samples %>%
#  as.data.frame() %>%
#  ggpairs()

dev.off()
#MCMCplot(run$samples)
#MCMCtrace(run$samples, pdf = F)

## save MCMC
saveRDS(samples, "outputs/SexInfectionInbreed_CONT_Samples_ICub_BaSTASetup.rds")
#samples <- readRDS("outputsRJ/samples_RJMCMC_SexInfectionFULLMODEL.rds")

## Marginal probabilities of inclusion for each variable
zNames <- model$expandNodeNames('z')
zCols <- which(colnames(samples) %in% zNames)
binary <- as.data.table((samples[, zCols] != 0) + 0)
res <- binary[ , .N, by=names(binary)]
res <- res[order(N, decreasing = T)]
res <- res[, prob := N/dim(samples)[1]]

res

saveRDS(res, "outputs/SexInfectionInbreed_CONT_PosteriorModelProbs_ICub_BaSTASetup.rds")
