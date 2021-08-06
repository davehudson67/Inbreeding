## load libraries
library(tidyverse)
library(coda)
library(patchwork)

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
newdata <- expand.grid(t = 0:50, sex = 0:1, infection1 = 0:1, infection2 = 0:1, inbr = 0:1) %>%
    filter(!(infection1 == 1 & infection2 == 1))
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

## Siler sampling (calculate on the log-scale to try and avoid rounding errors)
postSurv <- pSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6], lower.tail = FALSE, log = TRUE)
postDens <- dSiler(pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], pars[, 6], log = TRUE)
postMort <- exp(postDens - postSurv)
postSurv <- exp(postSurv)

## posterior predictive summaries
postSurv <- postSurv %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025),#, na.rm = TRUE),
      Median = median(x),
      UCI = quantile(x, probs = 0.975))#, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  cbind(newdata[, 1:5]) %>%
  mutate(MLH = ifelse(inbr == 0, "Outbred", "Inbred")) %>%
  mutate(Infection = ifelse(infection1 == 1, "Infected (cub)", "Uninfected")) %>%
  mutate(Infection = ifelse(infection2 == 1, "Infected (adult)", Infection)) %>%
  mutate(Sex = ifelse(sex == 1, "Male", "Female")) %>%
  select(-sex, -infection1, -infection2, -inbr)

postMort <- postMort %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025),#, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      UCI = quantile(x, probs = 0.975))#, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  cbind(newdata[, 1:5]) %>%
  mutate(MLH = ifelse(inbr == 0, "Outbred", "Inbred")) %>%
  mutate(Infection = ifelse(infection1 == 1, "Infected (cub)", "Uninfected")) %>%
  mutate(Infection = ifelse(infection2 == 1, "Infected (adult)", Infection)) %>%
  mutate(Sex = ifelse(sex == 1, "Male", "Female")) %>%
  select(-sex, -infection1, -infection2, -inbr)

## save outputs
saveRDS(postSurv, "outputs/PostSurv_FullModel_inbrCONT_zAll.rds")
saveRDS(postMort, "outputs/PostMort_FullModel_inbrCONT_zAll.rds")

## mortality plot
p1 <- mutate(postMort, Infection = factor(Infection)) %>%
  mutate(Infection = factor(Infection, levels = rev(levels(Infection)))) %>%
  ggplot(aes(x = t)) +
    geom_line(aes(y = Median, colour = MLH)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = MLH), alpha = 0.2) +
    facet_grid(Sex ~ Infection) +
    labs(title = "Mortality") +
    # scale_y_continuous(limits = c(0, 0.5)) +
    ylab("Posterior predictive mortality rate") +
    xlab("Time")

## survival plot
p2 <- mutate(postSurv, Infection = factor(Infection)) %>%
  mutate(Infection = factor(Infection, levels = rev(levels(Infection)))) %>%
  ggplot(aes(x = t)) +
    geom_line(aes(y = Median, colour = MLH)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = MLH), alpha = 0.2) +
    facet_grid(Sex ~ Infection) +
    labs(title = "Survival") +
    scale_y_continuous(limits = c(0, 1)) +
    ylab("Posterior predictive probability") +
    xlab("Time")

p <- p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("predictivePosteriors.pdf", p, width = 10, height = 5)
