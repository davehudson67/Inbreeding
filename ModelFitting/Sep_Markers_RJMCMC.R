##########################################################
##                                                      ###
##       Linear regression of Inb against b2            ###
##                                                      ###
###########################################################

rm(list=ls())
load("badgerSexInb.RData")
inb <- readRDS("inbr_markerdropped.rds")
inb$ID <- rownames(inb)
CH <- inner_join(CH, inb, by = "ID")

set.seed(seeds[12])

## set up plot output file
pdf("outputs/Siler_b2_RJMCMC_Inbreed_markerdropped.pdf")

code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood for interval-truncated gompertz
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1, a2, b1, b2 * b2mult[i], c1)
    
    log(b2mult[i]) <- beta[1] * inbr1[i] * z[1] +
                      beta[2] * inbr2[i] * z[2] +
                      beta[3] * inbr3[i] * z[3] +
                      beta[4] * inbr4[i] * z[4] +
                      beta[5] * inbr5[i] * z[5] +
                      beta[6] * inbr6[i] * z[6] +
                      beta[7] * inbr7[i] * z[7] +
                      beta[8] * inbr8[i] * z[8] +
                      beta[9] * inbr9[i] * z[9] +
                      beta[10] * inbr10[i] * z[10] +
                      beta[11] * inbr11[i] * z[11] +
                      beta[12] * inbr12[i] * z[12] +
                      beta[13] * inbr13[i] * z[13] +
                      beta[14] * inbr14[i] * z[14] +
                      beta[15] * inbr15[i] * z[15] +
                      beta[16] * inbr16[i] * z[16] +
                      beta[17] * inbr17[i] * z[17] +
                      beta[18] * inbr18[i] * z[18] +
                      beta[19] * inbr19[i] * z[19] +
                      beta[20] * inbr20[i] * z[20] +
                      beta[21] * inbr21[i] * z[21] +
                      beta[22] * inbr22[i] * z[22] +
                      beta[23] * inbr[i] * z[23]
                      #beta[24] * infection[i] * z[24] + 
                      #beta[25] * sex[i] * z[25] +
                      #beta[26] * sex[i] * inbr[i] * z[26] + 
                      #beta[27] * infection[i] * inbr[i] * z[27]
    
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
  constraint_data ~ dconstraint(sum(z[c(1:22, 23)]) <= 1)
  
})

## set up data
consts <- list(nind = nind, tM = tM, inbr = inbr, inbr1 = CH$hom1, inbr2 = CH$hom2, inbr3 = CH$hom3, inbr4 = CH$hom4, inbr5 = CH$hom5, inbr6 = CH$hom6, inbr7 = CH$hom7,
               inbr8 = CH$hom8, inbr9 = CH$hom9, inbr10 = CH$hom10, inbr11 = CH$hom11, inbr12 = CH$hom12, inbr13 = CH$hom13, inbr14 = CH$hom14, inbr15 = CH$hom15,
               inbr16 = CH$hom16, inbr17 = CH$hom17, inbr18 = CH$hom18, inbr19 = CH$hom19, inbr20 = CH$hom20, inbr21 = CH$hom21, inbr22 = CH$hom22)

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

samples <- samples[sample.int(nrow(samples), ceiling(nrow(samples) * 0.1)), ]
samples %>%
  as.data.frame() %>%
  ggpairs()


dev.off()
#MCMCplot(run$samples)
#MCMCtrace(run$samples, pdf = F)

## save MCMC
saveRDS(samples, "outputs/Samples_b2_RJMCMC_SexInfectionInbreed.rds")
#samples <- readRDS("outputsRJ/samples_RJMCMC_SexInfectionFULLMODEL.rds")

## Marginal probabilities of inclusion for each variable
library(data.table)
zNames <- model$expandNodeNames('z')
zCols <- which(colnames(samples) %in% zNames)
binary <- as.data.table((samples[, zCols] != 0) + 0)
res <- binary[ , .N, by=names(binary)]
res <- res[order(N, decreasing = T)]
res <- res[, prob := N/dim(samples)[1]]

res

saveRDS(res, "outputs/posteriorModelProbs_RJMCMC_SexInfectionInbreed.rds")

dev.off()
