###########################################################
##                                                      ###
##        Now fit Gompertz Makeham model                ###
##                                                      ###
###########################################################
library(tidyverse)
library(nimble)
library(purrr)
library(mclust)
library(GGally)
library(boot)
library(parallel)
library(MCMCvis)
dev.off()
#load("badgerInitialModelComparison.RData")
## set up plot output file
pdf("InitialModelComparison_AllBadgers/outputs/GompertzMakeham.pdf")

code <- nimbleCode({
  ## survival components for dead badgers
  for (i in 1:nind) {
    ## likelihood for interval-truncated gompertz makeham
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dgompzMakeNim(a, b, c1)
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  a ~ dexp(1)
  b ~ dexp(1)
  c1 ~ dexp(1)
  mean.p ~ dunif(0,1)
})

## set up other components of model
consts <- list(nind = nind, tM = tM)
data <- list(y = y, cint = cint, 
             censored = censored, tD = tD, dind = dind)

## find overdispersed initial values
tinitFn <- function(cint, censored) {
  apply(cbind(cint, censored), 1, function(x) {
    if(x[3] == 2) {
      y <- x[2] + 1
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
    sum(dGompzMake(t, a = pars[1], b = pars[2], c1 = pars[3], log = TRUE))
  }
  pars <- list(convergence = 1)
  k <- 0
  while(pars$convergence != 0 & k < 20) {
    ## sample missing values
    tD <- tinitFn(cint, censored)
    pars <- optim(rexp(3, 10), optFn, t = tD, control = list(fnscale = -1))
    k <- k + 1
  }
  if(k == 20) {
    stop("Can't sample initial values")
  }
  pars <- pars$par
  list(
    tD = tD,
    a = pars[1], 
    b = pars[2],
    c1 = pars[3],
    mean.p = runif(1, 0, 1)
  )
}

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = initFn(cint, censored))

## compile the model
cmodel <- compileNimble(model)

## set monitor
#help(samplers)
config <- configureMCMC(cmodel, monitors = c("a", "b", "c1", "mean.p"), thin = 1)
config$removeSamplers(c("a", "b", "c1"))
config$addSampler(target = c("a", "c1"), type = 'AF_slice', control = list(sliceAdaptFactorInterval = 250))
#config$addSampler(target = c("a"), type = 'slice')
config$addSampler(target = c("b"), type = 'slice')
#config$addSampler(target = c("c1"), type = 'slice')

## check monitors and samplers
config$printMonitors()
config$printSamplers(c("a", "b", "c1"))

## build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

## run the model
system.time(run <- runMCMC(cbuilt, 
                              niter = 500000, 
                              nburnin = 19000, 
                              nchains = 2, 
                              progressBar = TRUE, 
                              summary = TRUE, 
                              samplesAsCodaMCMC = TRUE, 
                              thin = 1))

## plot mcmcm
plot(run$samples)
samples <- run$samples
MCMCsummary(samples)
## save MCMC
saveRDS(run, "InitialModelComparison_AllBadgers/outputs/samples_gm.rds")
#samples <- readRDS("outputs/samples_s.rds")

## predictive plot against data
t_pred <- seq(0, max(cint), length.out = 100)
preds <- as.matrix(samples)[1:2000, ] %>%
  apply(1, function(pars, t) {
    pGompzMake(t, a = pars[1], b = pars[2], c1 = pars[3], lower.tail = FALSE)
  }, t = t_pred) %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, )

## plot ignoring sex differences
pred_plot <- km_plot[[1]] +
  geom_line(aes(x = t, y = Median), data = preds) +
  geom_ribbon(aes(x = t, y = Median, ymin = LCI, ymax = UCI), data = preds, alpha = 0.3)
pred_plot

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
props <- mclust::sim(mod$modelName, mod$parameters, nmix)
props <- props[, -1, drop = FALSE]
colnames(props) <- c("a", "b", "c1", "mean.p")

## take random samples from prior (to create defense mixture)
dmp <- runif((nimp - nmix), 0, 1)
defense <- as.data.frame(matrix(rexp(3 * (nimp - nmix), 1), ncol = 3))
defense <- defense %>%
  mutate(mean.p = dmp)
colnames(defense) <- c("a", "b", "c1", "mean.p")

## check IS distribution against posterior samples
as.data.frame(props) %>%
  mutate(type = "IS") %>%
  rbind(as.data.frame(samples) %>%
          mutate(type = "Post")) %>%
  ggpairs(mapping = aes(colour = type, alpha = 0.5), upper = list(continuous = "density"), columns = 1:4)

## combine defense and importance samples
props <- rbind(props, defense)

## calculate log likelihood
loglike <- function(y, zL, zU, censored, tM, a, b, c1, p) {
  ## loop over individuals
  fsurv1 <- pmap_dbl(list(y, zL, zU, censored, tM), function(y, zL, zU, censored, tM, a, b, c1, p) {
    if(censored == 1){                             #interval
      ll <- y * log(p) + (zL - y) * log(1 - p)
      t <- (zL + 1):zU
      temp <- pGompzMake(t, a, b, c1) - pGompzMake(t - 1, a, b, c1)
      temp1 <- log(temp) + (t - 1 - zL) * log(1 - p)
      ll <- ll + log_sum_exp_marg(temp1, mn = FALSE)
    } else {                                     #right censored
      ll <- y * log(p) + (zL - y) * log(1 - p)
      if(zL < tM) {
        t <- (zL + 1):tM
        temp <- pGompzMake(t, a, b, c1) - pGompzMake(t - 1, a, b, c1)
        temp1 <- log(temp) + (t - 1 - zL) * log(1 - p)
        ## tail component (after last capture time)
        temp <- pGompzMake(tM, a, b, c1, lower.tail = FALSE)
        temp <- log(temp) + (tM - zL) * log(1 - p)
        temp1 <- c(temp1, temp)
        ll <- ll + log_sum_exp_marg(temp1, mn = FALSE)
      } else {
        temp <- pGompzMake(tM, a, b, c1, lower.tail = FALSE)
        ll <- ll + log(temp) + (tM - zL) * log(1 - p)
      }
    }
    ll
  }, a = a, b = b, c1 = c1, p = p)
  sum(fsurv1)
}

## calculate log-likelihoods in parallel
logimpweight <- apply(props, 1, list)
logimpweight <- purrr::map(logimpweight, 1)
logimpweight <- mclapply(logimpweight,
                            function(pars, y, zL, zU, censored, tM) {
                              loglike(y, zL, zU, censored, tM, pars[1], pars[2], pars[3], pars[4])
                            }, y = y, zL = zL, zU = zU, censored = censored, tM = tM, mc.cores = 24)
logimpweight <- reduce(logimpweight, base::c)
# ## calculate log-likelihoods in serial
# logimpweight.gm <- pmap_dbl(list(mvn.gm[,1], mvn.gm[,2], mvn.gm[,3], mvn.gm[,4]),
#     function(a, b, c, p, y, zL, zU, censored, tM) {
#         loglike.gm(y, zL, zU, censored, tM, a, b, c, p)
#     }, y = y, zL = zL, zU = zU, censored = censored, tM = tM)

## add prior densities
logimpweight <- logimpweight + dexp(props[, 1], 1, log = TRUE) + 
  dexp(props[, 2], 1, log = TRUE) + dexp(props[, 3], 1, log = TRUE) +
  dunif(props[, 4], 0, 1, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(props, mod$modelName, mod$parameters, FALSE) + 0.05 * exp(dexp(props[, 1], 1, log = TRUE) +
                                                                              dexp(props[, 2], 1, log = TRUE) +
                                                                              dexp(props[, 3], 1, log = TRUE) +
                                                                              dunif(props[, 4], 0, 1, log = TRUE)))

saveRDS(logimpweight, "InitialModelComparison_AllBadgers/outputs/logimpweight_gm.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## bootstrap the importance weights to create 95% intervals
imp.boot <- BootsPlot(logimpweight, 5000, TRUE)

dev.off()
