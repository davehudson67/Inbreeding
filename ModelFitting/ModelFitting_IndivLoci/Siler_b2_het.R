###########################################################
##                                                      ###
##        Now fit Siler model - inb diff b2             ###
##                                                      ###
###########################################################

rm(list=ls())
load("badgerSexInb.RData")
set.seed(seeds[11])

names(gene)

## set up plot output file
pdf("IndivLociOutputs/Siler_b2het.pdf")

## set up inbreeding loci
inb <- CH$het
summary(inb)

code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood for interval-truncated gompertz
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dsilerNim(a1, a2, b1, b2 * b2mult[i], c1)
    
    log(b2mult[i]) <- betab2 * inb[i]
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  betab2 ~ dnorm(0, 1)
  a1 ~ dexp(1)
  a2 ~ dexp(1)
  b1 ~ dexp(1)
  b2 ~ dexp(1)
  c1 ~ dexp(1)
  mean.p ~ dunif(0, 1)
  
})

## set up data
consts <- list(nind = nind, tM = tM, inb = inb)
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
initFn <- function(cint, censored, inb) {
  ## get ML estimates as initial values
  optFn <- function(pars, t, inb) {
    if(any(pars[1:5] < 0)) {
      
      return(NA)
    }
    ll_ho <- sum(dSiler(t, a1 = pars[1], a2 = pars[2], b1 = pars[3], 
                        b2 = pars[4] * exp(pars[6]), c1 = pars[5], log = TRUE))
    
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
    pars <- optim(c(rexp(5, 10), rnorm(1, 0, 1)), optFn, t = tD, inb = inb, control = list(fnscale = -1))
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
    betab2 = pars[6]
  )
}

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = initFn(cint, censored, inb))

## compile the model
cModel <- compileNimble(model, showCompilerOutput = TRUE)

## try with adaptive slice sampler
config <- configureMCMC(cModel, monitors = c("a1", "a2", "b1", "b2", "betab2", "c1", "mean.p"), thin = 1)
config$removeSamplers(c("a1", "a2", "b1", "b2", "betab2", "c1"))
config$addSampler(target = c("a1", "a2", "b1", "b2", "betab2", "c1"), type = 'AF_slice', control = 50)
config$addSampler(target = c("a1", "c1", "b1"), type = 'AF_slice', control = 20)
#config$addSampler(target = c("b1", "c1"), type = 'AF_slice', control = 20)
#config$addSampler(target = c("b2"), type = 'slice')

## check monitors and samplers
config$printMonitors()
config$printSamplers(c("a1", "a2", "b1", "b2", "c1", "betab2"))

## build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

## run the model
system.time(run <- runMCMC(cbuilt, 
                           niter = 20000, 
                           nburnin = 5000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

run$summary
## plot mcmcm
plot(run$samples)
samples <- run$samples

## save MCMC
saveRDS(samples, "IndivLociOutputs/samples_s_b2het.rds")
#samples <- readRDS("ISoutputs/samples_s_SexDiffb2.rds")

## pairs plot
samples <- as.matrix(samples)
samples <- samples[sample.int(nrow(samples), ceiling(nrow(samples) * 0.1)), ]
samples %>%
  as.data.frame() %>%
  ggpairs()

## fit range of finite mixture models
mod <- densityMclust(samples)

## summary of finite mixture models
summary(mod)
plot(mod, what = "BIC")

## take random samples from mixture
nimp <- 10000
nmix <- rbinom(1, size = nimp, prob = 0.95)
props <- sim(mod$modelName, mod$parameters, nmix)
props <- props[, -1, drop = FALSE]
colnames(props) <- c("a1", "a2", "b1", "b2", "betab2", "c1", "mean.p")

## take random samples from prior (to create defense mixture)
dmp <- runif((nimp - nmix), 0, 1)
db <- rnorm((nimp - nmix), 0, 1)
defense <- as.data.frame(matrix(rexp(5 * (nimp - nmix), 1), ncol = 5))
defense <- defense %>%
  mutate(mean.p = dmp) %>%
  mutate(beta = db) %>%
  select("V1", "V2", "V3", "V4", "beta", "V5", "mean.p")
colnames(defense) <- c("a1", "a2", "b1", "b2", "betab2", "c1", "mean.p")

## adjust column names of samples to match IS
#colnames(samples) <- c("Pm", "a1f", "a1m", "a2", "b1", "b2", "c1", "mean.p")

## check IS distribution against posterior samples
as.data.frame(props) %>%
  mutate(type = "IS") %>%
  rbind(as.data.frame(samples) %>%
          mutate(type = "Post")) %>%
  ggpairs(mapping = aes(colour = type, alpha = 0.5), upper = list(continuous = "density"), columns = 1:7)

## combine defense and importance samples
props <- rbind(props, defense)

## this next function takes vectors for y, zL, zU, and censored
## and scalar tM and parameters and calculates log-likelihood
loglike <- function(y, zL, zU, censored, tM, a1, a2, b1, b2, betab2, c1, p, inb) {
  ## loop over individuals
  fsurv1 <- pmap_dbl(list(y, zL, zU, censored, tM, inb), function(y, zL, zU, censored, tM, a1, a2, b1, b2, betab2, c1, p, inb) {
    
    ll <- y * log(p) + (zL - y) * log(1 - p)
    
    if(censored == 1){                             #interval
      t <- (zL + 1):zU
        temp <- pSiler(t, a1, a2, b1, (b2 * exp(betab2)), c1) - pSiler(t - 1, a1, a2, b1, (b2 * exp(betab2)), c1)
        temp1 <- log(temp) + (t - 1 - zL) * log(1 - p)
        ll <- ll + log_sum_exp_marg(temp1, mn = FALSE)
    } else {                                     #right censored
      if(zL < tM) {
        t <- (zL + 1):tM
          temp <- pSiler(t, a1, a2, b1, (b2 * exp(betab2)), c1) - pSiler(t - 1, a1, a2, b1, (b2 * exp(betab2)), c1)
          temp1 <- log(temp) + (t - 1 - zL) * log(1 - p)
          ## tail component (after last capture time)
          temp <- pSiler(tM, a1, a2, b1, (b2 * exp(betab2)), c1, lower.tail = FALSE)
          temp <- log(temp) + (tM - zL) * log(1 - p)
          temp1 <- c(temp1, temp)
          ll <- ll + log_sum_exp_marg(temp1, mn = FALSE)
      } else {                                      #zL = tM
          temp <- pSiler(tM, a1, a2, b1, (b2 * exp(betab2)), c1, lower.tail = FALSE)
          ll <- ll + log(temp) + (tM - zL) * log(1 - p)
      }
      ll}
  }, p = p, a1 = a1, a2 = a2, b1 = b1, b2 = b2, betab2 = betab2, c1 = c1)
  sum(fsurv1)
}

## calculate log-likelihoods in parallel
logimpweight <- apply(props, 1, list)
logimpweight <- purrr::map(logimpweight, 1)
logimpweight <- mclapply(logimpweight,
                         function(pars, y, zL, zU, censored, tM, p, inb) {
                           loglike(y, zL, zU, censored, tM, pars[1], pars[2], pars[3], pars[4], pars[5], pars[6], pars[7], inb)
                         }, y = y, zL = zL, zU = zU, censored = censored, tM = tM, inb = inb, mc.cores = 24)
logimpweight <- reduce(logimpweight, base::c)
# ## calculate log-likelihoods in serial
# logimpweight.s <- pmap_dbl(list(mvn.s[,1], mvn.s[,2], mvn.s[,3], mvn.s[,4], mvn.s[,5], mvn.s[,6]), 
#     function(a1, a2, b1, b2, c, p, y, zL, zU, censored, tM) {
#         loglike.s(y, zL, zU, censored, tM, a1, a2, b1, b2, c, p)
#     }, y = y, zL = zL, zU = zU, censored = censored, tM = tM)

## add prior densities
logimpweight <- logimpweight + 
  dexp(props[, 1], 1, log = TRUE) + dexp(props[, 2], 1, log = TRUE) +
  dexp(props[, 3], 1, log = TRUE) + dexp(props[, 4], 1, log = TRUE) + 
  dnorm(props[, 5], 0, 1, log = TRUE) + dexp(props[, 6], 1, log = TRUE) + 
  dunif(props[, 7], 0, 1, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(mod$modelName, props, FALSE, mod$parameters) + 0.05 * exp(dexp(props[, 1], 1, log = TRUE) +
                                                                              dexp(props[, 2], 1, log = TRUE) +
                                                                              dexp(props[, 3], 1, log = TRUE) +
                                                                              dexp(props[, 4], 1, log = TRUE) +
                                                                              dnorm(props[, 5], 0, 1, log = TRUE) +
                                                                              dexp(props[, 6], 1, log = TRUE) +
                                                                              dunif(props[, 7], 0, 1, log = TRUE)))

saveRDS(logimpweight, "IndivLociOutputs/logimpweight_s_b2het.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## bootstrap the importance weights to create 95% intervals
imp.boot <- BootsPlot(logimpweight, 5000, TRUE)

dev.off()