## load libraries
library(tidyverse)
library(coda)

## source Siler functions
source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_Siler.R")

## load posterior samples
post <- readRDS("outputs/FullModel_z_all_runsamples.rds")
post <- as.matrix(post$samples)
post <- sample_n(as.data.frame(post), 50000)
post <- as.matrix(post)




## Need to produce different trajectories for...

# 1 - Outbred Uninfected male 
# 2 - Outbred Infected as cub male
# 3 - Outbred Infected as adult male
# 4 - Outbred Uninfected female 
# 5 - Outbred Infected as cub female
# 6 - Outbred Infected as adult female
# 7 - Inbred Uninfected male 
# 8 - Inbred Infected as cub male
# 9 - Inbred Infected as adult male
# 10 - Inbred Uninfected female 
# 11 - Inbred Infected as cub female
# 12 - Inbred Infected as adult female

####################################################################################################################################################################################
## create points to predict to for different groups
newdata <- expand.grid(t = 0:50, sex = 1, infection1 = 0, infection2 = 0, inbr = 0)
newdata$sexinb <- newdata$sex * newdata$inbr
newdata$infinb1 <- newdata$infection1 * newdata$inbr
newdata$infinb2 <- newdata$infection2 * newdata$inbr
newdata$sexinf1 <- newdata$sex * newdata$infection1
newdata$sexinf2 <- newdata$sex * newdata$infection2

## loop over Siler parameters
parnms <- c("a1", "a2", "b1", "b2", "c1")
pars <- list()
for(i in 1:5) {
  ## extract samples for linear predictor for e.g. a1
  ## making sure they line up with newdata
  ## (IN PARTICULAR INFECTION[, 2] AND [, 3] SEEM TO BE IN A DIFFERENT ORDER TO THE BETAS
  ##  WHICH IS FINE BUT THEY NEED TO MATCH UP - I'M CALLING INFECTION1 CUBS AND INFECTION2 ADULTS)
  parPost <- post[, match(paste0("beta", c("SEX", "INFCUB", "INFADULT", "INBR", "SEXINBR", "INFINBRCUB", "INFINBRADULT", "SEXINFCUB", "SEXINFADULT"), "[", i, "]"), colnames(post))]
  
  ## extract inclusion samples
  zpost <- post[, match(paste0("z", c("INBR", "INF", "SEX", "INFINBR", "SEXINBR", "SEXINF"), "[", i, "]"), colnames(post))]
  
  ## now multiply by relevant z indexes
  ## (expanding zpost matrix to make elementwise
  ## multiplication straightforward)
  zpost <- cbind(
    zpost[, paste0("zSEX[", i, "]")], ## Sex
    zpost[, paste0("zINF[", i, "]")], ## infection Cub
    zpost[, paste0("zINF[", i, "]")], ## infection Adult
    zpost[, paste0("zINBR[", i, "]")], ## inbreeding
    zpost[, paste0("zSEXINBR[", i, "]")], ## sex/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf1/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf2/inbreeding
    zpost[, paste0("zSEXINF[", i, "]")], ## inf1/sex
    zpost[, paste0("zSEXINF[", i, "]")] ## inf2/sex
  )
  
  ## posterior samples for a1 linear predictor
  ## using matrix multiplication because it's faster
  ## (first brackets uses ELEMENTWISE multiplication (*)
  ## and then we use MATRIX multiplication e.g. %*%)
  parLinpred <- (parPost * zpost) %*% t(newdata[, -1])
  
  ## add intercept
  parLinpred <- parLinpred + log(post[, parnms[i]])
  
  ## transformed parameters
  pars[[i]] <- exp(parLinpred)
}

## for efficient sampling of the Siler collapse down to large matrix
pars <- map(pars, t) %>%
  map(as.vector) %>%
  reduce(cbind) %>%
  {cbind(rep(newdata$t, nrow(zpost)), .)}

## Siler sampling
postSurv <- pSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6], lower.tail = FALSE)
postDens <- dSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6])
postMort <- postDens/postSurv


postSurv1 <- postSurv %>%    
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x),
      LCI = quantile(x, probs = 0.975, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Outbred") %>%
  mutate(Infection = "Uninfected") %>%
  mutate(Sex = "Male")

postMort1 <- postMort %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      LCI = quantile(x, probs = 0.975,na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Outbred") %>%
  mutate(Infection = "Uninfected") %>%
  mutate(Sex = "Male")

##########################################################################################################################################################################################

## create points to predict to for different groups
newdata <- expand.grid(t = 0:50, sex = 1, infection1 = 1, infection2 = 0, inbr = 0)
newdata$sexinb <- newdata$sex * newdata$inbr
newdata$infinb1 <- newdata$infection1 * newdata$inbr
newdata$infinb2 <- newdata$infection2 * newdata$inbr
newdata$sexinf1 <- newdata$sex * newdata$infection1
newdata$sexinf2 <- newdata$sex * newdata$infection2

## loop over Siler parameters
parnms <- c("a1", "a2", "b1", "b2", "c1")
pars <- list()
for(i in 1:5) {
  ## extract samples for linear predictor for e.g. a1
  ## making sure they line up with newdata
  ## (IN PARTICULAR INFECTION[, 2] AND [, 3] SEEM TO BE IN A DIFFERENT ORDER TO THE BETAS
  ##  WHICH IS FINE BUT THEY NEED TO MATCH UP - I'M CALLING INFECTION1 CUBS AND INFECTION2 ADULTS)
  parPost <- post[, match(paste0("beta", c("SEX", "INFCUB", "INFADULT", "INBR", "SEXINBR", "INFINBRCUB", "INFINBRADULT", "SEXINFCUB", "SEXINFADULT"), "[", i, "]"), colnames(post))]  
  ## extract inclusion samples
  zpost <- post[, match(paste0("z", c("INBR", "INF", "SEX", "INFINBR", "SEXINBR", "SEXINF"), "[", i, "]"), colnames(post))]
  
  ## now multiply by relevant z indexes
  ## (expanding zpost matrix to make elementwise
  ## multiplication straightforward)
  zpost <- cbind(
    zpost[, paste0("zSEX[", i, "]")], ## Sex
    zpost[, paste0("zINF[", i, "]")], ## infection Cub
    zpost[, paste0("zINF[", i, "]")], ## infection Adult
    zpost[, paste0("zINBR[", i, "]")], ## inbreeding
    zpost[, paste0("zSEXINBR[", i, "]")], ## sex/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf1/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf2/inbreeding
    zpost[, paste0("zSEXINF[", i, "]")], ## inf1/sex
    zpost[, paste0("zSEXINF[", i, "]")] ## inf2/sex
  )
  
  ## posterior samples for a1 linear predictor
  ## using matrix multiplication because it's faster
  ## (first brackets uses ELEMENTWISE multiplication (*)
  ## and then we use MATRIX multiplication e.g. %*%)
  parLinpred <- (parPost * zpost) %*% t(newdata[, -1])
  
  ## add intercept
  parLinpred <- parLinpred + log(post[, parnms[i]])
  
  ## transformed parameters
  pars[[i]] <- exp(parLinpred)
}

## for efficient sampling of the Siler collapse down to large matrix
pars <- map(pars, t) %>%
  map(as.vector) %>%
  reduce(cbind) %>%
  {cbind(rep(newdata$t, nrow(zpost)), .)}

## Siler sampling
postSurv <- pSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6], lower.tail = FALSE)
postDens <- dSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6])
postMort <- postDens/postSurv


postSurv2 <- postSurv %>%    
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x),
      LCI = quantile(x, probs = 0.975, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Outbred") %>%
  mutate(Infection = "Infected Cub") %>%
  mutate(Sex = "Male")

postMort2 <- postMort %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      LCI = quantile(x, probs = 0.975,na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Outbred") %>%
  mutate(Infection = "Infected Cub") %>%
  mutate(Sex = "Male")

###################################################################################################################################################################################

## create points to predict to for different groups
newdata <- expand.grid(t = 0:50, sex = 1, infection1 = 0, infection2 = 1, inbr = 0)
newdata$sexinb <- newdata$sex * newdata$inbr
newdata$infinb1 <- newdata$infection1 * newdata$inbr
newdata$infinb2 <- newdata$infection2 * newdata$inbr
newdata$sexinf1 <- newdata$sex * newdata$infection1
newdata$sexinf2 <- newdata$sex * newdata$infection2

## loop over Siler parameters
parnms <- c("a1", "a2", "b1", "b2", "c1")
pars <- list()
for(i in 1:5) {
  ## extract samples for linear predictor for e.g. a1
  ## making sure they line up with newdata
  ## (IN PARTICULAR INFECTION[, 2] AND [, 3] SEEM TO BE IN A DIFFERENT ORDER TO THE BETAS
  ##  WHICH IS FINE BUT THEY NEED TO MATCH UP - I'M CALLING INFECTION1 CUBS AND INFECTION2 ADULTS)
  parPost <- post[, match(paste0("beta", c("SEX", "INFCUB", "INFADULT", "INBR", "SEXINBR", "INFINBRCUB", "INFINBRADULT", "SEXINFCUB", "SEXINFADULT"), "[", i, "]"), colnames(post))]  
  ## extract inclusion samples
  zpost <- post[, match(paste0("z", c("INBR", "INF", "SEX", "INFINBR", "SEXINBR", "SEXINF"), "[", i, "]"), colnames(post))]
  
  ## now multiply by relevant z indexes
  ## (expanding zpost matrix to make elementwise
  ## multiplication straightforward)
  zpost <- cbind(
    zpost[, paste0("zSEX[", i, "]")], ## Sex
    zpost[, paste0("zINF[", i, "]")], ## infection Cub
    zpost[, paste0("zINF[", i, "]")], ## infection Adult
    zpost[, paste0("zINBR[", i, "]")], ## inbreeding
    zpost[, paste0("zSEXINBR[", i, "]")], ## sex/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf1/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf2/inbreeding
    zpost[, paste0("zSEXINF[", i, "]")], ## inf1/sex
    zpost[, paste0("zSEXINF[", i, "]")] ## inf2/sex
  )
  
  ## posterior samples for a1 linear predictor
  ## using matrix multiplication because it's faster
  ## (first brackets uses ELEMENTWISE multiplication (*)
  ## and then we use MATRIX multiplication e.g. %*%)
  parLinpred <- (parPost * zpost) %*% t(newdata[, -1])
  
  ## add intercept
  parLinpred <- parLinpred + log(post[, parnms[i]])
  
  ## transformed parameters
  pars[[i]] <- exp(parLinpred)
}

## for efficient sampling of the Siler collapse down to large matrix
pars <- map(pars, t) %>%
  map(as.vector) %>%
  reduce(cbind) %>%
  {cbind(rep(newdata$t, nrow(zpost)), .)}

## Siler sampling
postSurv <- pSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6], lower.tail = FALSE)
postDens <- dSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6])
postMort <- postDens/postSurv


postSurv3 <- postSurv %>%    
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x),
      LCI = quantile(x, probs = 0.975, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Outbred") %>%
  mutate(Infection = "Infected Adult") %>%
  mutate(Sex = "Male")

postMort3 <- postMort %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      LCI = quantile(x, probs = 0.975,na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Outbred") %>%
  mutate(Infection = "Infected Adult") %>%
  mutate(Sex = "Male")

##########################################################################################################################################################################################

## create points to predict to for different groups
newdata <- expand.grid(t = 0:50, sex = 0, infection1 = 0, infection2 = 0, inbr = 0)
newdata$sexinb <- newdata$sex * newdata$inbr
newdata$infinb1 <- newdata$infection1 * newdata$inbr
newdata$infinb2 <- newdata$infection2 * newdata$inbr
newdata$sexinf1 <- newdata$sex * newdata$infection1
newdata$sexinf2 <- newdata$sex * newdata$infection2

## loop over Siler parameters
parnms <- c("a1", "a2", "b1", "b2", "c1")
pars <- list()
for(i in 1:5) {
  ## extract samples for linear predictor for e.g. a1
  ## making sure they line up with newdata
  ## (IN PARTICULAR INFECTION[, 2] AND [, 3] SEEM TO BE IN A DIFFERENT ORDER TO THE BETAS
  ##  WHICH IS FINE BUT THEY NEED TO MATCH UP - I'M CALLING INFECTION1 CUBS AND INFECTION2 ADULTS)
  parPost <- post[, match(paste0("beta", c("SEX", "INFCUB", "INFADULT", "INBR", "SEXINBR", "INFINBRCUB", "INFINBRADULT", "SEXINFCUB", "SEXINFADULT"), "[", i, "]"), colnames(post))]  
  ## extract inclusion samples
  zpost <- post[, match(paste0("z", c("INBR", "INF", "SEX", "INFINBR", "SEXINBR", "SEXINF"), "[", i, "]"), colnames(post))]
  
  ## now multiply by relevant z indexes
  ## (expanding zpost matrix to make elementwise
  ## multiplication straightforward)
  zpost <- cbind(
    zpost[, paste0("zSEX[", i, "]")], ## Sex
    zpost[, paste0("zINF[", i, "]")], ## infection Cub
    zpost[, paste0("zINF[", i, "]")], ## infection Adult
    zpost[, paste0("zINBR[", i, "]")], ## inbreeding
    zpost[, paste0("zSEXINBR[", i, "]")], ## sex/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf1/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf2/inbreeding
    zpost[, paste0("zSEXINF[", i, "]")], ## inf1/sex
    zpost[, paste0("zSEXINF[", i, "]")] ## inf2/sex
  )
  
  ## posterior samples for a1 linear predictor
  ## using matrix multiplication because it's faster
  ## (first brackets uses ELEMENTWISE multiplication (*)
  ## and then we use MATRIX multiplication e.g. %*%)
  parLinpred <- (parPost * zpost) %*% t(newdata[, -1])
  
  ## add intercept
  parLinpred <- parLinpred + log(post[, parnms[i]])
  
  ## transformed parameters
  pars[[i]] <- exp(parLinpred)
}

## for efficient sampling of the Siler collapse down to large matrix
pars <- map(pars, t) %>%
  map(as.vector) %>%
  reduce(cbind) %>%
  {cbind(rep(newdata$t, nrow(zpost)), .)}

## Siler sampling
postSurv <- pSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6], lower.tail = FALSE)
postDens <- dSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6])
postMort <- postDens/postSurv


postSurv4 <- postSurv %>%    
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x),
      LCI = quantile(x, probs = 0.975, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Outbred") %>%
  mutate(Infection = "Uninfected") %>%
  mutate(Sex = "Female")

postMort4 <- postMort %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      LCI = quantile(x, probs = 0.975,na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Outbred") %>%
  mutate(Infection = "Uninfected") %>%
  mutate(Sex = "Female")

###################################################################################################################################################################################

## create points to predict to for different groups
newdata <- expand.grid(t = 0:50, sex = 0, infection1 = 1, infection2 = 0, inbr = 0)
newdata$sexinb <- newdata$sex * newdata$inbr
newdata$infinb1 <- newdata$infection1 * newdata$inbr
newdata$infinb2 <- newdata$infection2 * newdata$inbr
newdata$sexinf1 <- newdata$sex * newdata$infection1
newdata$sexinf2 <- newdata$sex * newdata$infection2

## loop over Siler parameters
parnms <- c("a1", "a2", "b1", "b2", "c1")
pars <- list()
for(i in 1:5) {
  ## extract samples for linear predictor for e.g. a1
  ## making sure they line up with newdata
  ## (IN PARTICULAR INFECTION[, 2] AND [, 3] SEEM TO BE IN A DIFFERENT ORDER TO THE BETAS
  ##  WHICH IS FINE BUT THEY NEED TO MATCH UP - I'M CALLING INFECTION1 CUBS AND INFECTION2 ADULTS)
  parPost <- post[, match(paste0("beta", c("SEX", "INFCUB", "INFADULT", "INBR", "SEXINBR", "INFINBRCUB", "INFINBRADULT", "SEXINFCUB", "SEXINFADULT"), "[", i, "]"), colnames(post))]  
  ## extract inclusion samples
  zpost <- post[, match(paste0("z", c("INBR", "INF", "SEX", "INFINBR", "SEXINBR", "SEXINF"), "[", i, "]"), colnames(post))]
  
  ## now multiply by relevant z indexes
  ## (expanding zpost matrix to make elementwise
  ## multiplication straightforward)
  zpost <- cbind(
    zpost[, paste0("zSEX[", i, "]")], ## Sex
    zpost[, paste0("zINF[", i, "]")], ## infection Cub
    zpost[, paste0("zINF[", i, "]")], ## infection Adult
    zpost[, paste0("zINBR[", i, "]")], ## inbreeding
    zpost[, paste0("zSEXINBR[", i, "]")], ## sex/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf1/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf2/inbreeding
    zpost[, paste0("zSEXINF[", i, "]")], ## inf1/sex
    zpost[, paste0("zSEXINF[", i, "]")] ## inf2/sex
  )
  
  ## posterior samples for a1 linear predictor
  ## using matrix multiplication because it's faster
  ## (first brackets uses ELEMENTWISE multiplication (*)
  ## and then we use MATRIX multiplication e.g. %*%)
  parLinpred <- (parPost * zpost) %*% t(newdata[, -1])
  
  ## add intercept
  parLinpred <- parLinpred + log(post[, parnms[i]])
  
  ## transformed parameters
  pars[[i]] <- exp(parLinpred)
}

## for efficient sampling of the Siler collapse down to large matrix
pars <- map(pars, t) %>%
  map(as.vector) %>%
  reduce(cbind) %>%
  {cbind(rep(newdata$t, nrow(zpost)), .)}

## Siler sampling
postSurv <- pSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6], lower.tail = FALSE)
postDens <- dSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6])
postMort <- postDens/postSurv


postSurv5 <- postSurv %>%    
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x),
      LCI = quantile(x, probs = 0.975, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Outbred") %>%
  mutate(Infection = "Infected Cub") %>%
  mutate(Sex = "Female")

postMort5 <- postMort %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      LCI = quantile(x, probs = 0.975,na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Outbred") %>%
  mutate(Infection = "Infected Cub") %>%
  mutate(Sex = "Female")

###################################################################################################################################################################################

## create points to predict to for different groups
newdata <- expand.grid(t = 0:50, sex = 0, infection1 = 0, infection2 = 1, inbr = 0)
newdata$sexinb <- newdata$sex * newdata$inbr
newdata$infinb1 <- newdata$infection1 * newdata$inbr
newdata$infinb2 <- newdata$infection2 * newdata$inbr
newdata$sexinf1 <- newdata$sex * newdata$infection1
newdata$sexinf2 <- newdata$sex * newdata$infection2

## loop over Siler parameters
parnms <- c("a1", "a2", "b1", "b2", "c1")
pars <- list()
for(i in 1:5) {
  ## extract samples for linear predictor for e.g. a1
  ## making sure they line up with newdata
  ## (IN PARTICULAR INFECTION[, 2] AND [, 3] SEEM TO BE IN A DIFFERENT ORDER TO THE BETAS
  ##  WHICH IS FINE BUT THEY NEED TO MATCH UP - I'M CALLING INFECTION1 CUBS AND INFECTION2 ADULTS)
  parPost <- post[, match(paste0("beta", c("SEX", "INFCUB", "INFADULT", "INBR", "SEXINBR", "INFINBRCUB", "INFINBRADULT", "SEXINFCUB", "SEXINFADULT"), "[", i, "]"), colnames(post))]  
  ## extract inclusion samples
  zpost <- post[, match(paste0("z", c("INBR", "INF", "SEX", "INFINBR", "SEXINBR", "SEXINF"), "[", i, "]"), colnames(post))]
  
  ## now multiply by relevant z indexes
  ## (expanding zpost matrix to make elementwise
  ## multiplication straightforward)
  zpost <- cbind(
    zpost[, paste0("zSEX[", i, "]")], ## Sex
    zpost[, paste0("zINF[", i, "]")], ## infection Cub
    zpost[, paste0("zINF[", i, "]")], ## infection Adult
    zpost[, paste0("zINBR[", i, "]")], ## inbreeding
    zpost[, paste0("zSEXINBR[", i, "]")], ## sex/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf1/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf2/inbreeding
    zpost[, paste0("zSEXINF[", i, "]")], ## inf1/sex
    zpost[, paste0("zSEXINF[", i, "]")] ## inf2/sex
  )
  
  ## posterior samples for a1 linear predictor
  ## using matrix multiplication because it's faster
  ## (first brackets uses ELEMENTWISE multiplication (*)
  ## and then we use MATRIX multiplication e.g. %*%)
  parLinpred <- (parPost * zpost) %*% t(newdata[, -1])
  
  ## add intercept
  parLinpred <- parLinpred + log(post[, parnms[i]])
  
  ## transformed parameters
  pars[[i]] <- exp(parLinpred)
}

## for efficient sampling of the Siler collapse down to large matrix
pars <- map(pars, t) %>%
  map(as.vector) %>%
  reduce(cbind) %>%
  {cbind(rep(newdata$t, nrow(zpost)), .)}

## Siler sampling
postSurv <- pSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6], lower.tail = FALSE)
postDens <- dSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6])
postMort <- postDens/postSurv


postSurv6 <- postSurv %>%    
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x),
      LCI = quantile(x, probs = 0.975, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Outbred") %>%
  mutate(Infection = "Infected Adult") %>%
  mutate(Sex = "Female")

postMort6 <- postMort %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      LCI = quantile(x, probs = 0.975,na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Outbred") %>%
  mutate(Infection = "Infected Adult") %>%
  mutate(Sex = "Female")

##########################################################################################################################################################################################
##########################################################################################################################################################################################

## create points to predict to for different groups
newdata <- expand.grid(t = 0:50, sex = 1, infection1 = 0, infection2 = 0, inbr = 1)
newdata$sexinb <- newdata$sex * newdata$inbr
newdata$infinb1 <- newdata$infection1 * newdata$inbr
newdata$infinb2 <- newdata$infection2 * newdata$inbr
newdata$sexinf1 <- newdata$sex * newdata$infection1
newdata$sexinf2 <- newdata$sex * newdata$infection2

## loop over Siler parameters
parnms <- c("a1", "a2", "b1", "b2", "c1")
pars <- list()
for(i in 1:5) {
  ## extract samples for linear predictor for e.g. a1
  ## making sure they line up with newdata
  ## (IN PARTICULAR INFECTION[, 2] AND [, 3] SEEM TO BE IN A DIFFERENT ORDER TO THE BETAS
  ##  WHICH IS FINE BUT THEY NEED TO MATCH UP - I'M CALLING INFECTION1 CUBS AND INFECTION2 ADULTS)
  parPost <- post[, match(paste0("beta", c("SEX", "INFCUB", "INFADULT", "INBR", "SEXINBR", "INFINBRCUB", "INFINBRADULT", "SEXINFCUB", "SEXINFADULT"), "[", i, "]"), colnames(post))]  
  ## extract inclusion samples
  zpost <- post[, match(paste0("z", c("INBR", "INF", "SEX", "INFINBR", "SEXINBR", "SEXINF"), "[", i, "]"), colnames(post))]
  
  ## now multiply by relevant z indexes
  ## (expanding zpost matrix to make elementwise
  ## multiplication straightforward)
  zpost <- cbind(
    zpost[, paste0("zSEX[", i, "]")], ## Sex
    zpost[, paste0("zINF[", i, "]")], ## infection Cub
    zpost[, paste0("zINF[", i, "]")], ## infection Adult
    zpost[, paste0("zINBR[", i, "]")], ## inbreeding
    zpost[, paste0("zSEXINBR[", i, "]")], ## sex/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf1/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf2/inbreeding
    zpost[, paste0("zSEXINF[", i, "]")], ## inf1/sex
    zpost[, paste0("zSEXINF[", i, "]")] ## inf2/sex
  )
  
  ## posterior samples for a1 linear predictor
  ## using matrix multiplication because it's faster
  ## (first brackets uses ELEMENTWISE multiplication (*)
  ## and then we use MATRIX multiplication e.g. %*%)
  parLinpred <- (parPost * zpost) %*% t(newdata[, -1])
  
  ## add intercept
  parLinpred <- parLinpred + log(post[, parnms[i]])
  
  ## transformed parameters
  pars[[i]] <- exp(parLinpred)
}

## for efficient sampling of the Siler collapse down to large matrix
pars <- map(pars, t) %>%
  map(as.vector) %>%
  reduce(cbind) %>%
  {cbind(rep(newdata$t, nrow(zpost)), .)}

## Siler sampling
postSurv <- pSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6], lower.tail = FALSE)
postDens <- dSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6])
postMort <- postDens/postSurv


postSurv7 <- postSurv %>%    
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      LCI = quantile(x, probs = 0.975, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Inbred") %>%
  mutate(Infection = "Uninfected") %>%
  mutate(Sex = "Male")

postMort7 <- postMort %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      LCI = quantile(x, probs = 0.975,na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Inbred") %>%
  mutate(Infection = "Uninfected") %>%
  mutate(Sex = "Male")

##########################################################################################################################################################################################

## create points to predict to for different groups
newdata <- expand.grid(t = 0:50, sex = 1, infection1 = 1, infection2 = 0, inbr = 1)
newdata$sexinb <- newdata$sex * newdata$inbr
newdata$infinb1 <- newdata$infection1 * newdata$inbr
newdata$infinb2 <- newdata$infection2 * newdata$inbr
newdata$sexinf1 <- newdata$sex * newdata$infection1
newdata$sexinf2 <- newdata$sex * newdata$infection2

## loop over Siler parameters
parnms <- c("a1", "a2", "b1", "b2", "c1")
pars <- list()
for(i in 1:5) {
  ## extract samples for linear predictor for e.g. a1
  ## making sure they line up with newdata
  ## (IN PARTICULAR INFECTION[, 2] AND [, 3] SEEM TO BE IN A DIFFERENT ORDER TO THE BETAS
  ##  WHICH IS FINE BUT THEY NEED TO MATCH UP - I'M CALLING INFECTION1 CUBS AND INFECTION2 ADULTS)
  parPost <- post[, match(paste0("beta", c("SEX", "INFCUB", "INFADULT", "INBR", "SEXINBR", "INFINBRCUB", "INFINBRADULT", "SEXINFCUB", "SEXINFADULT"), "[", i, "]"), colnames(post))]  
  ## extract inclusion samples
  zpost <- post[, match(paste0("z", c("INBR", "INF", "SEX", "INFINBR", "SEXINBR", "SEXINF"), "[", i, "]"), colnames(post))]
  
  ## now multiply by relevant z indexes
  ## (expanding zpost matrix to make elementwise
  ## multiplication straightforward)
  zpost <- cbind(
    zpost[, paste0("zSEX[", i, "]")], ## Sex
    zpost[, paste0("zINF[", i, "]")], ## infection Cub
    zpost[, paste0("zINF[", i, "]")], ## infection Adult
    zpost[, paste0("zINBR[", i, "]")], ## inbreeding
    zpost[, paste0("zSEXINBR[", i, "]")], ## sex/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf1/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf2/inbreeding
    zpost[, paste0("zSEXINF[", i, "]")], ## inf1/sex
    zpost[, paste0("zSEXINF[", i, "]")] ## inf2/sex
  )
  
  ## posterior samples for a1 linear predictor
  ## using matrix multiplication because it's faster
  ## (first brackets uses ELEMENTWISE multiplication (*)
  ## and then we use MATRIX multiplication e.g. %*%)
  parLinpred <- (parPost * zpost) %*% t(newdata[, -1])
  
  ## add intercept
  parLinpred <- parLinpred + log(post[, parnms[i]])
  
  ## transformed parameters
  pars[[i]] <- exp(parLinpred)
}

## for efficient sampling of the Siler collapse down to large matrix
pars <- map(pars, t) %>%
  map(as.vector) %>%
  reduce(cbind) %>%
  {cbind(rep(newdata$t, nrow(zpost)), .)}

## Siler sampling
postSurv <- pSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6], lower.tail = FALSE)
postDens <- dSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6])
postMort <- postDens/postSurv


postSurv8 <- postSurv %>%    
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x),
      LCI = quantile(x, probs = 0.975, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Inbred") %>%
  mutate(Infection = "Infected Cub") %>%
  mutate(Sex = "Male")

postMort8 <- postMort %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      LCI = quantile(x, probs = 0.975,na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Inbred") %>%
  mutate(Infection = "Infected Cub") %>%
  mutate(Sex = "Male")

###################################################################################################################################################################################

## create points to predict to for different groups
newdata <- expand.grid(t = 0:50, sex = 1, infection1 = 0, infection2 = 1, inbr = 1)
newdata$sexinb <- newdata$sex * newdata$inbr
newdata$infinb1 <- newdata$infection1 * newdata$inbr
newdata$infinb2 <- newdata$infection2 * newdata$inbr
newdata$sexinf1 <- newdata$sex * newdata$infection1
newdata$sexinf2 <- newdata$sex * newdata$infection2

## loop over Siler parameters
parnms <- c("a1", "a2", "b1", "b2", "c1")
pars <- list()
for(i in 1:5) {
  ## extract samples for linear predictor for e.g. a1
  ## making sure they line up with newdata
  ## (IN PARTICULAR INFECTION[, 2] AND [, 3] SEEM TO BE IN A DIFFERENT ORDER TO THE BETAS
  ##  WHICH IS FINE BUT THEY NEED TO MATCH UP - I'M CALLING INFECTION1 CUBS AND INFECTION2 ADULTS)
  parPost <- post[, match(paste0("beta", c("SEX", "INFCUB", "INFADULT", "INBR", "SEXINBR", "INFINBRCUB", "INFINBRADULT", "SEXINFCUB", "SEXINFADULT"), "[", i, "]"), colnames(post))]  
  ## extract inclusion samples
  zpost <- post[, match(paste0("z", c("INBR", "INF", "SEX", "INFINBR", "SEXINBR", "SEXINF"), "[", i, "]"), colnames(post))]
  
  ## now multiply by relevant z indexes
  ## (expanding zpost matrix to make elementwise
  ## multiplication straightforward)
  zpost <- cbind(
    zpost[, paste0("zSEX[", i, "]")], ## Sex
    zpost[, paste0("zINF[", i, "]")], ## infection Cub
    zpost[, paste0("zINF[", i, "]")], ## infection Adult
    zpost[, paste0("zINBR[", i, "]")], ## inbreeding
    zpost[, paste0("zSEXINBR[", i, "]")], ## sex/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf1/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf2/inbreeding
    zpost[, paste0("zSEXINF[", i, "]")], ## inf1/sex
    zpost[, paste0("zSEXINF[", i, "]")] ## inf2/sex
  )
  
  ## posterior samples for a1 linear predictor
  ## using matrix multiplication because it's faster
  ## (first brackets uses ELEMENTWISE multiplication (*)
  ## and then we use MATRIX multiplication e.g. %*%)
  parLinpred <- (parPost * zpost) %*% t(newdata[, -1])
  
  ## add intercept
  parLinpred <- parLinpred + log(post[, parnms[i]])
  
  ## transformed parameters
  pars[[i]] <- exp(parLinpred)
}

## for efficient sampling of the Siler collapse down to large matrix
pars <- map(pars, t) %>%
  map(as.vector) %>%
  reduce(cbind) %>%
  {cbind(rep(newdata$t, nrow(zpost)), .)}

## Siler sampling
postSurv <- pSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6], lower.tail = FALSE)
postDens <- dSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6])
postMort <- postDens/postSurv


postSurv9 <- postSurv %>%    
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x),
      LCI = quantile(x, probs = 0.975, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Inbred") %>%
  mutate(Infection = "Infected Adult") %>%
  mutate(Sex = "Male")

postMort9 <- postMort %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      LCI = quantile(x, probs = 0.975,na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Inbred") %>%
  mutate(Infection = "Infected Adult") %>%
  mutate(Sex = "Male")

##########################################################################################################################################################################################

## create points to predict to for different groups
newdata <- expand.grid(t = 0:50, sex = 0, infection1 = 0, infection2 = 0, inbr = 1)
newdata$sexinb <- newdata$sex * newdata$inbr
newdata$infinb1 <- newdata$infection1 * newdata$inbr
newdata$infinb2 <- newdata$infection2 * newdata$inbr
newdata$sexinf1 <- newdata$sex * newdata$infection1
newdata$sexinf2 <- newdata$sex * newdata$infection2

## loop over Siler parameters
parnms <- c("a1", "a2", "b1", "b2", "c1")
pars <- list()
for(i in 1:5) {
  ## extract samples for linear predictor for e.g. a1
  ## making sure they line up with newdata
  ## (IN PARTICULAR INFECTION[, 2] AND [, 3] SEEM TO BE IN A DIFFERENT ORDER TO THE BETAS
  ##  WHICH IS FINE BUT THEY NEED TO MATCH UP - I'M CALLING INFECTION1 CUBS AND INFECTION2 ADULTS)
  parPost <- post[, match(paste0("beta", c("SEX", "INFCUB", "INFADULT", "INBR", "SEXINBR", "INFINBRCUB", "INFINBRADULT", "SEXINFCUB", "SEXINFADULT"), "[", i, "]"), colnames(post))]  
  ## extract inclusion samples
  zpost <- post[, match(paste0("z", c("INBR", "INF", "SEX", "INFINBR", "SEXINBR", "SEXINF"), "[", i, "]"), colnames(post))]
  
  ## now multiply by relevant z indexes
  ## (expanding zpost matrix to make elementwise
  ## multiplication straightforward)
  zpost <- cbind(
    zpost[, paste0("zSEX[", i, "]")], ## Sex
    zpost[, paste0("zINF[", i, "]")], ## infection Cub
    zpost[, paste0("zINF[", i, "]")], ## infection Adult
    zpost[, paste0("zINBR[", i, "]")], ## inbreeding
    zpost[, paste0("zSEXINBR[", i, "]")], ## sex/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf1/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf2/inbreeding
    zpost[, paste0("zSEXINF[", i, "]")], ## inf1/sex
    zpost[, paste0("zSEXINF[", i, "]")] ## inf2/sex
  )
  
  ## posterior samples for a1 linear predictor
  ## using matrix multiplication because it's faster
  ## (first brackets uses ELEMENTWISE multiplication (*)
  ## and then we use MATRIX multiplication e.g. %*%)
  parLinpred <- (parPost * zpost) %*% t(newdata[, -1])
  
  ## add intercept
  parLinpred <- parLinpred + log(post[, parnms[i]])
  
  ## transformed parameters
  pars[[i]] <- exp(parLinpred)
}

## for efficient sampling of the Siler collapse down to large matrix
pars <- map(pars, t) %>%
  map(as.vector) %>%
  reduce(cbind) %>%
  {cbind(rep(newdata$t, nrow(zpost)), .)}

## Siler sampling
postSurv <- pSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6], lower.tail = FALSE)
postDens <- dSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6])
postMort <- postDens/postSurv


postSurv10 <- postSurv %>%    
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x),
      LCI = quantile(x, probs = 0.975, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Inbred") %>%
  mutate(Infection = "Uninfected") %>%
  mutate(Sex = "Female")

postMort10 <- postMort %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      LCI = quantile(x, probs = 0.975,na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Inbred") %>%
  mutate(Infection = "Uninfected") %>%
  mutate(Sex = "Female")

###################################################################################################################################################################################

## create points to predict to for different groups
newdata <- expand.grid(t = 0:50, sex = 0, infection1 = 1, infection2 = 0, inbr = 1)
newdata$sexinb <- newdata$sex * newdata$inbr
newdata$infinb1 <- newdata$infection1 * newdata$inbr
newdata$infinb2 <- newdata$infection2 * newdata$inbr
newdata$sexinf1 <- newdata$sex * newdata$infection1
newdata$sexinf2 <- newdata$sex * newdata$infection2

## loop over Siler parameters
parnms <- c("a1", "a2", "b1", "b2", "c1")
pars <- list()
for(i in 1:5) {
  ## extract samples for linear predictor for e.g. a1
  ## making sure they line up with newdata
  ## (IN PARTICULAR INFECTION[, 2] AND [, 3] SEEM TO BE IN A DIFFERENT ORDER TO THE BETAS
  ##  WHICH IS FINE BUT THEY NEED TO MATCH UP - I'M CALLING INFECTION1 CUBS AND INFECTION2 ADULTS)
  parPost <- post[, match(paste0("beta", c("SEX", "INFCUB", "INFADULT", "INBR", "SEXINBR", "INFINBRCUB", "INFINBRADULT", "SEXINFCUB", "SEXINFADULT"), "[", i, "]"), colnames(post))]  
  ## extract inclusion samples
  zpost <- post[, match(paste0("z", c("INBR", "INF", "SEX", "INFINBR", "SEXINBR", "SEXINF"), "[", i, "]"), colnames(post))]
  
  ## now multiply by relevant z indexes
  ## (expanding zpost matrix to make elementwise
  ## multiplication straightforward)
  zpost <- cbind(
    zpost[, paste0("zSEX[", i, "]")], ## Sex
    zpost[, paste0("zINF[", i, "]")], ## infection Cub
    zpost[, paste0("zINF[", i, "]")], ## infection Adult
    zpost[, paste0("zINBR[", i, "]")], ## inbreeding
    zpost[, paste0("zSEXINBR[", i, "]")], ## sex/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf1/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf2/inbreeding
    zpost[, paste0("zSEXINF[", i, "]")], ## inf1/sex
    zpost[, paste0("zSEXINF[", i, "]")] ## inf2/sex
  )
  
  ## posterior samples for a1 linear predictor
  ## using matrix multiplication because it's faster
  ## (first brackets uses ELEMENTWISE multiplication (*)
  ## and then we use MATRIX multiplication e.g. %*%)
  parLinpred <- (parPost * zpost) %*% t(newdata[, -1])
  
  ## add intercept
  parLinpred <- parLinpred + log(post[, parnms[i]])
  
  ## transformed parameters
  pars[[i]] <- exp(parLinpred)
}

## for efficient sampling of the Siler collapse down to large matrix
pars <- map(pars, t) %>%
  map(as.vector) %>%
  reduce(cbind) %>%
  {cbind(rep(newdata$t, nrow(zpost)), .)}

## Siler sampling
postSurv <- pSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6], lower.tail = FALSE)
postDens <- dSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6])
postMort <- postDens/postSurv


postSurv11 <- postSurv %>%    
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x),
      LCI = quantile(x, probs = 0.975, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Inbred") %>%
  mutate(Infection = "Infected Cub") %>%
  mutate(Sex = "Female")

postMort11 <- postMort %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      LCI = quantile(x, probs = 0.975,na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Inbred") %>%
  mutate(Infection = "Infected Cub") %>%
  mutate(Sex = "Female")

###################################################################################################################################################################################

## create points to predict to for different groups
newdata <- expand.grid(t = 0:50, sex = 0, infection1 = 0, infection2 = 1, inbr = 1)
newdata$sexinb <- newdata$sex * newdata$inbr
newdata$infinb1 <- newdata$infection1 * newdata$inbr
newdata$infinb2 <- newdata$infection2 * newdata$inbr
newdata$sexinf1 <- newdata$sex * newdata$infection1
newdata$sexinf2 <- newdata$sex * newdata$infection2

## loop over Siler parameters
parnms <- c("a1", "a2", "b1", "b2", "c1")
pars <- list()
for(i in 1:5) {
  ## extract samples for linear predictor for e.g. a1
  ## making sure they line up with newdata
  ## (IN PARTICULAR INFECTION[, 2] AND [, 3] SEEM TO BE IN A DIFFERENT ORDER TO THE BETAS
  ##  WHICH IS FINE BUT THEY NEED TO MATCH UP - I'M CALLING INFECTION1 CUBS AND INFECTION2 ADULTS)
  parPost <- post[, match(paste0("beta", c("SEX", "INFCUB", "INFADULT", "INBR", "SEXINBR", "INFINBRCUB", "INFINBRADULT", "SEXINFCUB", "SEXINFADULT"), "[", i, "]"), colnames(post))]  
  ## extract inclusion samples
  zpost <- post[, match(paste0("z", c("INBR", "INF", "SEX", "INFINBR", "SEXINBR", "SEXINF"), "[", i, "]"), colnames(post))]
  
  ## now multiply by relevant z indexes
  ## (expanding zpost matrix to make elementwise
  ## multiplication straightforward)
  zpost <- cbind(
    zpost[, paste0("zSEX[", i, "]")], ## Sex
    zpost[, paste0("zINF[", i, "]")], ## infection Cub
    zpost[, paste0("zINF[", i, "]")], ## infection Adult
    zpost[, paste0("zINBR[", i, "]")], ## inbreeding
    zpost[, paste0("zSEXINBR[", i, "]")], ## sex/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf1/inbreeding
    zpost[, paste0("zINFINBR[", i, "]")], ## inf2/inbreeding
    zpost[, paste0("zSEXINF[", i, "]")], ## inf1/sex
    zpost[, paste0("zSEXINF[", i, "]")] ## inf2/sex
  )
  
  ## posterior samples for a1 linear predictor
  ## using matrix multiplication because it's faster
  ## (first brackets uses ELEMENTWISE multiplication (*)
  ## and then we use MATRIX multiplication e.g. %*%)
  parLinpred <- (parPost * zpost) %*% t(newdata[, -1])
  
  ## add intercept
  parLinpred <- parLinpred + log(post[, parnms[i]])
  
  ## transformed parameters
  pars[[i]] <- exp(parLinpred)
}

## for efficient sampling of the Siler collapse down to large matrix
pars <- map(pars, t) %>%
  map(as.vector) %>%
  reduce(cbind) %>%
  {cbind(rep(newdata$t, nrow(zpost)), .)}

## Siler sampling
postSurv <- pSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6], lower.tail = FALSE)
postDens <- dSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6])
postMort <- postDens/postSurv


postSurv12 <- postSurv %>%    
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x),
      LCI = quantile(x, probs = 0.975, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Inbred") %>%
  mutate(Infection = "Infected Adult") %>%
  mutate(Sex = "Female")

postMort12 <- postMort %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      LCI = quantile(x, probs = 0.975, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = newdata$t) %>%
  mutate(MLH = "Inbred") %>%
  mutate(Infection = "Infected Adult") %>%
  mutate(Sex = "Female")

###################################################################################################################################################################################

postSurv <- bind_rows(postSurv1, postSurv2, postSurv3, postSurv4, postSurv5, postSurv6, postSurv7, postSurv8,
                      postSurv9, postSurv10, postSurv11, postSurv12)

postMort <- bind_rows(postMort1, postMort2, postMort3, postMort4, postMort5, postMort6, postMort7, postMort8,
                      postMort9, postMort10, postMort11, postMort12)

saveRDS(postSurv, "outputs/PostSurv_FullModel_inbrCONT_zAll.rds")
saveRDS(postMort, "outputs/PostMort_FullModel_inbrCONT_zAll.rds")

postMort <- readRDS("outputs/PostMort_FullModel_zAll.rds")
saveRDS(postMort, "outputs/PostMort_FullModel_inbrCONT_zAll.rds")



ggplot(postMort, aes(x = t)) +
  geom_line(aes(y = Median, colour = MLH)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = MLH), alpha = 0.2) +
  facet_wrap(Sex ~ Infection) +
  labs(title = "Mortality") +
  scale_y_continuous(limits = c(0, 0.5))



ggplot(postSurv, aes(x = t)) +
  geom_line(aes(y = Median, colour = MLH)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = MLH), alpha = 0.2) +
  facet_wrap(Sex ~ Infection) +
  labs(title = "Survival")

dev.off()
scale_x_continuous(limits = c(0, 50)) +
  scale_y_continuous(limits = c(0, 0.5))
