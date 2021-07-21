## load libraries
library(tidyverse)
library(coda)

## source Siler functions
source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_Siler.R")

## load posterior samples
post <- readRDS("outputs/Sex3InfectionInbrCATmed_AllParameters_runsamples.rds")
post <- as.matrix(post$samples)

## take a chain to test
post <- post[[1]]#[1:10, ]

## create points to predict to for uninfected outbred males
newdata <- expand.grid(t = 0:80, sex = 1, infection1 = 0, infection2 = 0, inbr = 0)
newdata$sexinb <- newdata$sex * newdata$inbr
newdata$infinb1 <- newdata$infection1 * newdata$inbr
newdata$infinb2 <- newdata$infection2 * newdata$inbr

## extract inclusion samples
zpost <- post[, match(paste0("z", c("INBR", "SEXINBR", "INFINBR"), "[1]"), colnames(post))]

## now multiply by relevant z indexes
## (expanding zpost matrix to make elementwise
## multiplication straightforward)
zpost <- cbind(
    rep(1, nrow(zpost)), ## sex
    rep(1, nrow(zpost)), ## infection1
    rep(1, nrow(zpost)), ## infection2
    zpost[, "zINBR[1]"], ## inbreeding
    zpost[, "zSEXINBR[1]"], ## sex/inbreeding
    zpost[, "zINFINBR[1]"], ## inf1/inbreeding
    zpost[, "zINFINBR[1]"] ## inf2/inbreeding
)

## loop over Siler parameters
parnms <- c("a1", "a2", "b1", "b2", "c1")
pars <- list()
for(i in 1:5) {
    ## extract samples for linear predictor for e.g. a1
    ## making sure they line up with newdata
    ## (IN PARTICULAR INFECTION[, 2] AND [, 3] SEEM TO BE IN A DIFFERENT ORDER TO THE BETAS
    ##  WHICH IS FINE BUT THEY NEED TO MATCH UP - I'M CALLING INFECTION1 CUBS AND INFECTION2 ADULTS)
    parPost[[i]] <- post[, match(paste0("beta", c("SEX", "INFCUB", "INFADULT", "INBR", "SEXINBR", "INFINBRCUB", "INFINBRADULT"), "[", i, "]"), colnames(post))]
    zpost[[i]] <- post[, match(paste0("z", c("INBR", "SEXINBR", "INFINBR"), "[", i, "]"), colnames(post))]
    
    zpost[[i]] <- cbind(
        rep(1, nrow(zpost)), ## sex
        rep(1, nrow(zpost)), ## infection1
        rep(1, nrow(zpost)), ## infection2
        zpost[, paste0("zINBR", "[", i, "]")], ## inbreeding
        zpost[, paste0("zSEXINBR", "[", i, "]")], ## sex/inbreeding
        zpost[, paste0("zINFINBRCUB", "[", i, "]")], ## inf1/inbreeding
        zpost[, paste0("zINFINBRADULT", "[", i, "]")] ## inf2/inbreeding
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
postSurv <- pSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6], lower.tail = FALSE) %>%
    matrix(nrow = nrow(newdata)) %>%
    apply(1, function(x) {
        c(LCI = quantile(x, probs = 0.025),
        median = median(x),
        LCI = quantile(x, probs = 0.975))
    }) %>%
    t() %>%
    as_tibble() %>%
    set_names(c("LCI", "Median", "UCI")) %>%
    mutate(t = newdata$t)
ggplot(postSurv, aes(x = t)) +
    geom_line(aes(y = Median)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5)
