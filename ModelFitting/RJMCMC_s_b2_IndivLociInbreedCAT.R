##########################################################
##                                                      ###
##       Linear regression of Inb against b2            ###
##                                                      ###
###########################################################

rm(list=ls())
load("badgerSexInb.RData")
set.seed(seeds[12])

## set up plot output file
pdf("outputs/Siler_b2_RJMCMC_IndivLociInbreedCAT.pdf")

## separate into vectors so each can be used as categorical variable in analysis
inb1 <- CH$V1
inb2 <- CH$V2
inb3 <- CH$V3
inb4 <- CH$V4
inb5 <- CH$V5
inb6 <- CH$V6
inb7 <- CH$V7
inb8 <- CH$V8
inb9 <- CH$V9
inb10 <- CH$V10
inb11 <- CH$V11
inb12 <- CH$V12
inb13 <- CH$V13
inb14 <- CH$V14
inb15 <- CH$V15
inb16 <- CH$V16
inb17 <- CH$V17
inb18 <- CH$V18
inb19 <- CH$V19
inb20 <- CH$V20
inb21 <- CH$V21
inb22 <- CH$V22



code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood for interval-truncated gompertz
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1, a2, b1, b2 * b2mult[i], c1)
    
    log(b2mult[i]) <- beta[1] * inb1[i] * z[1] +
      beta[2] * inb2[i] * z[2] +
      beta[3] * inb3[i] * z[3] +
      beta[4] * inb4[i] * z[4] +
      beta[5] * inb5[i] * z[5] +
      beta[6] * inb6[i] * z[6] +
      beta[7] * inb7[i] * z[7] +
      beta[8] * inb8[i] * z[8] +
      beta[9] * inb9[i] * z[9] +
      beta[10] * inb10[i] * z[10] +
      beta[11] * inb11[i] * z[11] +
      beta[12] * inb12[i] * z[12] +
      beta[13] * inb13[i] * z[13] +
      beta[14] * inb14[i] * z[14] +
      beta[15] * inb15[i] * z[15] +
      beta[16] * inb16[i] * z[16] +
      beta[17] * inb17[i] * z[17] +
      beta[18] * inb18[i] * z[18] +
      beta[19] * inb19[i] * z[19] +
      beta[20] * inb20[i] * z[20] +
      beta[21] * inb21[i] * z[21] +
      beta[22] * inb22[i] * z[22] +
      beta[23] * inbhetCAT[i] * z[23]
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  for (j in 1:23){
    beta[j] ~ dnorm(0, sd = 1)
    z[j] ~ dbern(0.5)
  }
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  b2 ~ dexp(1)
  c1 ~ dexp(1)
  mean.p ~ dunif(0, 1)
  constraint_data ~ dconstraint(sum(z[1:22]) <= 1)
  
})

## set up data
consts <- list(nind = nind, tM = tM, inb1 = inb1,
               inb2 = inb2, inb3 = inb3, inb4 = inb4, inb5 = inb5, inb6 = inb6, inb7 = inb7, inb8 = inb8,
               inb9 = inb9, inb10 = inb10, inb11 = inb11, inb12 = inb12, inb13 = inb13, inb14 = inb14,
               inb15 = inb15, inb16 = inb16, inb17 = inb17, inb18 = inb18, inb19 = inb19, inb20 = inb20,
               inb21 = inb21, inb22 = inb22, inbhetCAT = inbhetCAT)

data <- list(y = y, cint = cint, 
             censored = censored, tD = tD, dind = dind, constraint_data = 1)

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
    
#    llF <- sum(dSiler(t[sex == 0], a1 = pars[1], a2 = pars[2], b1 = pars[3], 
#                      b2 = pars[4], c1 = pars[5], log = TRUE))
#    llM + llF
  }
  pars <- list(convergence = 1)
  k <- 0
  while(pars$convergence != 0 & k < 50) {
    ## sample missing values
    tD <- tinitFn(cint, censored)
    ## sample sex proportion
    #    Pm <- runif(1, 0.4, 0.6)
    ## sample missing sex indicators
    #    sexI <- rbinom(length(censored), size = 1, prob = Pm)
    #    sexI[!is.na(sex)] <- sex[!is.na(sex)]
    ## optimise to interval-censored only
    pars <- optim(c(rexp(5, 10), rnorm(1, 0, 1)), optFn, t = tD, control = list(fnscale = -1))
    k <- k + 1
  }
  if(k == 50) {
    stop("Can't sample initial values")
  }
  pars <- pars$par
  ## check log-likelihoods
  #  ll <- sum(dGompertz(tD[sexI == 0], a = pars[1], b = pars[3], log = TRUE))
  #  ll <- ll + sum(dGompertz(tD[sexI == 1], a = pars[2], b = pars[3], log = TRUE))
  #  stopifnot(is.finite(ll))
  ## reformat sex initial conditions correctly
  #  sexI[!is.na(sex)] <- NA
  ## output initial values
  list(
    tD = tD,
    #    sex = sexI,
    #    Pm = Pm,
    a1 = pars[1],
    a2 = pars[2],
    b1 = pars[3],
    b2 = pars[4],
    c1 = pars[5],
    mean.p = runif(1, 0, 1),
    beta = rep(pars[6], times = 23),
    z = rep(0, times = 23)
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
config$addSampler(target = c("a1", "c1"), type = 'AF_slice', control = 20)
#config$addSampler(target = c("b2"), type = 'slice')

## Add reversible jump
configureRJ(conf = config,   ## model configuration
            targetNodes = c("beta[1]", "beta[2]","beta[3]", "beta[4]", "beta[5]",
                            "beta[6]", "beta[7]","beta[8]", "beta[9]", "beta[10]",
                            "beta[11]", "beta[12]","beta[13]", "beta[14]", "beta[15]",
                            "beta[16]", "beta[17]","beta[18]", "beta[19]", "beta[20]",
                            "beta[21]", "beta[22]", "beta[23]" ),    ## coefficients for selection
            indicatorNodes = c("z[1]", "z[2]", "z[3]", "z[4]", "z[5]", "z[6]", "z[7]", "z[8]", "z[9]", "z[10]",
                               "z[11]", "z[12]", "z[13]", "z[14]", "z[15]", "z[16]", "z[17]", "z[18]", "z[19]", "z[20]", 
                               "z[21]", "z[22]", "z[23]"),    ## indicators paired with coefficients
            control = list(mean = 0.5, scale = .5))

config$addMonitors("beta[1]", "beta[2]","beta[3]", "beta[4]", "beta[5]",
                   "beta[6]", "beta[7]","beta[8]", "beta[9]", "beta[10]",
                   "beta[11]", "beta[12]","beta[13]", "beta[14]", "beta[15]",
                   "beta[16]", "beta[17]","beta[18]", "beta[19]", "beta[20]",
                   "beta[21]", "beta[22]", "beta[23]", "z[1]", "z[2]", "z[3]", "z[4]", "z[5]", "z[6]", "z[7]", "z[8]", "z[9]", "z[10]",
                   "z[11]", "z[12]", "z[13]", "z[14]", "z[15]", "z[16]", "z[17]", "z[18]", "z[19]", "z[20]", 
                   "z[21]", "z[22]", "z[23]")

config$printSamplers(c("beta[1]", "beta[2]","beta[3]", "beta[4]", "beta[5]",
                       "beta[6]", "beta[7]","beta[8]", "beta[9]", "beta[10]",
                       "beta[11]", "beta[12]","beta[13]", "beta[14]", "beta[15]",
                       "beta[16]", "beta[17]","beta[18]", "beta[19]", "beta[20]",
                       "beta[21]", "beta[22]", "beta[23]", "z[1]", "z[2]", "z[3]", "z[4]", "z[5]", "z[6]", "z[7]", "z[8]", "z[9]", "z[10]",
                       "z[11]", "z[12]", "z[13]", "z[14]", "z[15]", "z[16]", "z[17]", "z[18]", "z[19]", "z[20]", 
                       "z[21]", "z[22]", "z[23]", "a1", "a2", "b1", "b2", "c1", "mean.p"))

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

dev.off()
#MCMCplot(run$samples)
#MCMCtrace(run$samples, pdf = F)

## save MCMC
saveRDS(samples, "outputs/samples_b2_RJMCMC_IndivLociInbreedCAT.rds")
#samples <- readRDS("outputsRJ/samples_RJMCMC_SexInfectionFULLMODEL.rds")

## Marginal probabilities of inclusion for each variable
zNames <- model$expandNodeNames('z')
zCols <- which(colnames(samples) %in% zNames)
binary <- as.data.table((samples[, zCols] != 0) + 0)
res <- binary[ , .N, by=names(binary)]
res <- res[order(N, decreasing = T)]
res <- res[, prob := N/dim(samples)[1]]

res

saveRDS(res, "outputs/posteriorModelProbs_RJMCMC_IndivLociInbreedCAT.rds")

dev.off()


