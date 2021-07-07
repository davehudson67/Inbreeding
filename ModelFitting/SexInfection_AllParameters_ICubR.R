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
pdf("outputs/SexInfection_AllParameters_ICub.pdf")

code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood for interval-truncated siler
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1 * a1mult[i], a2 * a2mult[i], b1 * b1mult[i], b2 * b2mult[i], c1 * c1mult[i])
    
    log(a1mult[i]) <- beta[1] * sex[i] * z[1] + 
      beta[2] * infection[i] * z[2] + 
      beta[3] * sex[i] * infection[i] * z[3]
    
    log(a2mult[i]) <- beta[4] * sex[i] * z[4] + 
      beta[5] * infection[i] * z[5] + 
      beta[6] * sex[i] * infection[i] * z[6]
    
    log(b1mult[i]) <- beta[7] * sex[i] * z[7] + 
      beta[8] * infection[i] * z[8] + 
      beta[9] * sex[i] * infection[i] * z[9]
    
    log(b2mult[i]) <- beta[10] * sex[i] * z[10] + 
      beta[11] * infection[i] * z[11] + 
      beta[12] * sex[i] * infection[i] * z[12]
    
    log(c1mult[i]) <- beta[13] * sex[i] * z[13] + 
      beta[14] * infection[i] * z[14] + 
      beta[15] * sex[i] * infection[i] * z[15]
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  for (j in 1:15){
    beta[j] ~ dnorm(0, sd = 1)
    z[j] ~ dbern(0.5)
  }
  
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  b2 ~ dexp(1)
  c1 ~ dexp(1)
  mean.p ~ dunif(0, 1)
  constraint_dataSexInf1 ~ dconstraint(z[3] <= z[1] * z[2])
  constraint_dataSexInf2 ~ dconstraint(z[6] <= z[4] * z[5])
  constraint_dataSexInf3 ~ dconstraint(z[9] <= z[7] * z[8])
  constraint_dataSexInf4 ~ dconstraint(z[12] <= z[10] * z[11])
  constraint_dataSexInf5 ~ dconstraint(z[14] <= z[13] * z[14])
  
})


## set up data
consts <- list(nind = nind, tM = tM, sex = sex, infection = infection, inbr = inbr)

data <- list(y = y, cint = cint, 
             censored = censored, tD = tD, dind = dind, constraint_dataSexInf1 = 1,
             constraint_dataSexInf2 = 1,
             constraint_dataSexInf3 = 1,
             constraint_dataSexInf4 = 1,
             constraint_dataSexInf5 = 1)

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
    beta = rep(pars[6], times = 15),
    z = rep(0, times = 15)
  )
}


## build the model
model <- nimbleModel(code, constants = consts,
                     data = data, 
                     inits = initFn(cint, censored))

## configure model
config <- configureMCMC(model)

config$removeSamplers(c("a1", "a2", "b1", "b2", "c1"))
config$addSampler(target = c("a1", "a2", "b1", "b2", "c1"), type = 'AF_slice', control = 50)
config$addSampler(target = c("a1", "b1", "c1"), type = 'AF_slice', control = 20)
config$addSampler(target = c("b2"), type = 'slice')

## Add reversible jump
configureRJ(conf = config,   ## model configuration
            targetNodes = c("beta"),    ## coefficients for selection
            indicatorNodes = c("z"),   ## indicators paired with coefficients
            control = list(mean = 0.5, scale = .5))

config$addMonitors("beta", "z")

config$printSamplers(c("beta", "z", "a1", "a2", "b1", "b2", "c1", "mean.p"))

rIndicatorMCMC <- buildMCMC(config)
cIndicatorModel <- compileNimble(model)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

system.time(run <- runMCMC(cIndicatorMCMC, 
                           niter = 70000, 
                           nburnin = 10000, 
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
saveRDS(samples, "outputs/SexInfection_AllParameters_ICub.rds")
#samples <- readRDS("outputsRJ/samples_RJMCMC_SexInfectionFULLMODEL.rds")

## Marginal probabilities of inclusion for each variable
zNames <- model$expandNodeNames('z')
zCols <- which(colnames(samples) %in% zNames)
binary <- as.data.table((samples[, zCols] != 0) + 0)
res <- binary[ , .N, by=names(binary)]
res <- res[order(N, decreasing = T)]
res <- res[, prob := N/dim(samples)[1]]
res
res
saveRDS(res, "outputs/SexInfection_AllParameters_PosteriorModelProbs_ICub.rds")

samples <- as.data.frame(samples)

z_indicators <- samples %>%
  select(c(22:36)) %>%
  colSums()

z_indicators <- data.frame(z_indicators/120000)
z_indicators$parameter <- c("a1_Sex", "a1_Infection", "a1_SexInfection", 
                            "a2_Sex", "a2_Infection", "a2SexInfection",
                            "b1_Sex", "b1_Infection", "b1SexInfection",
                            "b2_Sex", "b2_Infection", "b2SexInfection",
                            "c_Sex", "c_Infection", "cSexInfection")
colnames(z_indicators) <- c("Inclusion_Prob", "parameter")

ggplot(z_indicators, aes(x = parameter, y = Inclusion_Prob)) +
  geom_point()
