###################################################
##                                               ##
##   refit the data using an exponential model   ##
##                                               ##
###################################################
library(tidyverse)
library(nimble)
library(purrr)
library(mclust)
library(GGally)
library(boot)
library(parallel)
library(MCMCvis)

## set up plot output file
pdf("InitialModelComparison_AllBadgers/outputs/Exponential.pdf")

code <- nimbleCode({
  ## survival components for dead badgers
  for (i in 1:nind) {
    ## likelihood for interval-truncated exponential
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dexp(r)
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  r ~ dexp(1)
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
    if(pars < 0) {
      return(NA)
    }
    sum(dExpo(t, r = pars, log = TRUE))
  }
  ## sample missing values
  tD <- tinitFn(cint, censored)
  pars <- optimise(optFn, c(0, 100), t = tD, maximum = TRUE)
  pars <- pars$maximum
  list(
    tD = tD,
    r = pars,
    mean.p = runif(1, 0, 1)
  )
}

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = initFn(cint, censored))

## compile the model
cmodel <- compileNimble(model)

## set monitor
config <- configureMCMC(cmodel, monitors = c("r", "mean.p"), thin = 1)

## check monitors and samplers
config$printMonitors()
config$printSamplers("r")

## build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

## run the model
system.time(run_e <- runMCMC(cbuilt, 
                             niter = 500000, 
                             nburnin = 19000, 
                             nchains = 2, 
                             progressBar = TRUE, 
                             summary = TRUE, 
                             samplesAsCodaMCMC = TRUE, 
                             thin = 1))
run_e$summary

## plot mcmcm
plot(run_e$samples)
samples <- run_e$samples

## save MCMC
saveRDS(run_e, "InitialModelComparison_AllBadgers/outputs/samples_e.rds")
#samples <- readRDS("outputs/samples_e.rds")

## predictive plot against data
t_pred <- seq(0, max(cint), length.out = 100)
preds <- as.matrix(samples)%>%
  apply(1, function(pars, t) {
    pexp(t, r = pars[1], lower.tail = FALSE)
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
props <- sim(mod$modelName, mod$parameters, nmix)
props <- props[, -1, drop = FALSE]
colnames(props) <- c("mean.p", "r")

## take random samples from prior (to create defense mixture)
dmp <- runif((nimp - nmix), 0, 1)
defense <- rexp(1 * (nimp - nmix), 1)
defense <- cbind(dmp, defense)
colnames(defense) <- c("mean.p", "r")

## check IS distribution against posterior samples
as.data.frame(props) %>%
  mutate(type = "IS") %>%
  rbind(as.data.frame(samples) %>%
          mutate(type = "Post")) %>%
  ggplot(aes(x = r, colour = type)) + geom_density()

## combine defense and importance samples
props <- rbind(props, defense)

## this next function takes vectors for y, zL, zU and censored
## and scalar tM and parameters and calculates log-likelihood
loglike <- function(y, zL, zU, censored, tM, r, p) {
  ## loop over individuals
  fsurv1 <- pmap_dbl(list(y, zL, zU, censored, tM), function(y, zL, zU, censored, tM, r, p) {
    if(censored == 1) {                             #interval
      ll <- y * log(p) + (zL - y) * log(1 - p)
      t <- (zL + 1):zU
      temp <- pexp(t, r) - pexp(t - 1, r)
      temp1 <- log(temp) + (t - 1 - zL) * log(1 - p)
      ll <- ll + log_sum_exp_marg(temp1, mn = FALSE)
    } else {                                     #right censored
      ll <- y * log(p) + (zL - y) * log(1 - p)
      if(zL < tM) {
        t <- (zL + 1):tM
        temp <- pexp(t, r) - pexp(t - 1, r)
        temp1 <- log(temp) + (t - 1 - zL) * log(1 - p)
        ## tail component (after last capture time)
        temp <- pexp(tM, r, lower.tail = FALSE)
        temp <- log(temp) + (tM - zL) * log(1 - p)
        temp1 <- c(temp1, temp)
        ll <- ll + log_sum_exp_marg(temp1, mn = FALSE)
      } else {
        temp <- pexp(tM, r, lower.tail = FALSE)
        ll <- ll + log(temp) + (tM - zL) * log(1 - p)
      }
    }
    ll
  }, r = r, p = p)
  sum(fsurv1)
}

## calculate log-likelihoods in parallel
logimpweight <- apply(props, 1, list)
logimpweight <- purrr::map(logimpweight, 1)
logimpweight <- mclapply(logimpweight,
                         function(pars, y, zL, zU, censored, tM) {
                           loglike(y, zL, zU, censored, tM, pars[2], pars[1])
                         }, y = y, zL = zL, zU = zU, censored = censored, tM = tM, mc.cores = 24)
logimpweight <- reduce(logimpweight, base::c)

# ## calculate log-likelihoods in serial
## this now loops over the importance samples
#logimpweight <- map2_dbl(props[, 1], props[, 2], function(p, r, y, zL, zU, censored, tM) {
#  loglike(y, zL, zU, censored, tM, r, p)
#}, y = y, zL = zL, zU = zU, censored = censored, tM = tM)

## add prior densities
logimpweight <- logimpweight + dexp(props[, 2], 1, log = TRUE) + dunif(props[, 1], 0, 1, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(props, mod$modelName, mod$parameters, FALSE) + 0.05 * exp(dexp(props[, 2], 1, log = TRUE) +
                                                                              dunif(props[, 1], 0, 1, log = TRUE)))

saveRDS(logimpweight, "InitialModelComparison_AllBadgers/outputs/logimpweight_e.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## calculate log-marginal likelihood
logmarg.e <- log_sum_exp_marg(logimpweight)

## bootstrap the importance weights to create 95% intervals
imp.boot.e <- BootsPlot(logimpweight, 5000, TRUE)

dev.off()
