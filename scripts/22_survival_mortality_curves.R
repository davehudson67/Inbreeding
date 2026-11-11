## load libraries
library(tidyverse)
library(coda)
library(patchwork)
library(MCMCvis)
library(bayestestR)

## source Siler functions
source("../../SimulationStudy/FirstPaperFiles/Distributions/Dist_Siler.R")

## load posterior samples
run <- readRDS("Outputs/SilerFullModel_InbrCONT_runsamples.rds")
post <- as.matrix(run$samples)
post <- sample_n(as.data.frame(post), 50000)
post <- as.data.frame(post)

#post$'betaINBR[1]' <- 0
#post$'betaINBR[2]' <- 0
#post$'betaINBR[4]' <- 0
#post$'betaINBR[3]' <- post$betaINBR
#post$'betaINBR[5]' <- 0
#post$'zINBR[1]' <- 0
#post$'zINBR[2]' <- 0
#post$'zINBR[4]' <- 0
#post$'zINBR[3]' <- post$zINBR
#post$'zINBR[5]' <- 0

#post$'betaINBRINF[1]' <- 0
#post$'betaINBRINF[2]' <- 0
#post$'betaINBRINF[4]' <- 0
#post$'betaINBRINF[3]' <- post$betaINBRINF
#post$'betaINBRINF[5]' <- 0
#post$'zINBRINF[1]' <- 0
#post$'zINBRINF[2]' <- 0
#post$'zINBRINF[4]' <- 0
#post$'zINBRINF[3]' <- post$zINBRINF
#post$'zINBRINF[5]' <- 0

#post$'betaINBRSEX[1]' <- 0
#post$'betaINBRSEX[2]' <- 0
#post$'betaINBRSEX[4]' <- 0
#post$'betaINBRSEX[3]' <- post$betaINBRSEX
#post$'betaINBRSEX[5]' <- 0
#post$'zINBRSEX[1]' <- 0
#post$'zINBRSEX[2]' <- 0
#post$'zINBRSEX[4]' <- 0
#post$'zINBRSEX[3]' <- post$zINBRSEX
#post$'zINBRSEX[5]' <- 0

#post$'betaINBRINFSEX[1]' <- 0
#post$'betaINBRINFSEX[2]' <- 0
#post$'betaINBRINFSEX[4]' <- 0
#post$'betaINBRINFSEX[3]' <- post$betaINBRINFSEX
#post$'betaINBRINFSEX[5]' <- 0
#post$'zINBRINFSEX[1]' <- 0
#post$'zINBRINFSEX[2]' <- 0
#post$'zINBRINFSEX[4]' <- 0
#post$'zINBRINFSEX[3]' <- post$zINBRINFSEX
#post$'zINBRINFSEX[5]' <- 0

post$'zSEX[1]' <- 1
post$'zSEX[2]' <- 1
post$'zSEX[3]' <- 1
post$'zSEX[4]' <- 1
post$'zSEX[5]' <- 1

post$'zINF[1]' <- 1
post$'zINF[2]' <- 1
post$'zINF[3]' <- 1
post$'zINF[4]' <- 1
post$'zINF[5]' <- 1

#post <- select(post, -c(betaINBR, zINBR, betaINBRINF, zINBRINF, betaINBRSEX, zINBRSEX, betaINBRINFSEX, zINBRINFSEX))

post <- as.matrix(post)

## Need to produce different trajectories for...

# 1 - Outbred Uninfected male 
# 2 - Outbred Infected as cub male
# 3 - Outbred Uninfected female 
# 4 - Outbred Infected as cub female
# 5 - Inbred Uninfected male
# 6 - Inbred Infected as cub male
# 7 - Inbred Uninfected female 
# 8 - Inbred Infected as cub female

####################################################################################################################################################################################
## create points to predict to for different groups
newdata <- expand.grid(t = 0:50, sex = 0:1, inf = 0:1, inbr = seq(0.2, 0.8, 0.2))

#filter(!(infection1 == 1 & infection2 == 1))
newdata$inbrinf <- newdata$inbr * newdata$inf
newdata$inbrsex <- newdata$inbr * newdata$sex
newdata$inbrinfsex <- newdata$inbr * newdata$inf * newdata$sex
#newdata$sexinf2 <- newdata$sex * newdata$infection2

## loop over Siler parameters
parnms <- c("a1", "a2", "b1", "b2", "c1")
pars <- list()

for(i in 1:5) {
  ## extract samples for linear predictor for e.g. a1
  ## making sure they line up with newdata
  ## (IN PARTICULAR INFECTION[, 2] AND [, 3] SEEM TO BE IN A DIFFERENT ORDER TO THE BETAS
  ##  WHICH IS FINE BUT THEY NEED TO MATCH UP - I'M CALLING INFECTION1 CUBS AND INFECTION2 ADULTS)
  parPost <- post[, match(paste0("beta", c("SEX", "INF", "INBR", "INBRINF", "INBRSEX", "INBRINFSEX"), "[", i, "]"), colnames(post))]
  
  ## extract inclusion samples
  zpost <- post[, match(paste0("z", c("SEX", "INF", "INBR", "INBRINF", "INBRSEX", "INBRINFSEX"), "[", i, "]"), colnames(post))]
  
  ## now multiply by relevant z indexes
  ## (expanding zpost matrix to make elementwise
  ## multiplication straightforward)
  zpost <- cbind(
    zpost[, paste0("zSEX[", i, "]")], ## Sex
    zpost[, paste0("zINF[", i, "]")], ## infection Cub
    zpost[, paste0("zINBR[", i, "]")], ## inbreeding
    zpost[, paste0("zINBRINF[", i, "]")], ## inbreeding*infectionCUB
    zpost[, paste0("zINBRSEX[", i, "]")], ## inbreeding
    zpost[, paste0("zINBRINFSEX[", i, "]")] ## inbreeding*infectionCUB
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
pars <- purrr::map(pars, t) %>%
  purrr::map(as.vector) %>%
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
  cbind(newdata[, 1:4]) %>%
  mutate(MLH = as.factor(inbr)) %>%
  mutate(Infection = ifelse(inf == 1, "Infected as a cub", "Uninfected life")) %>%
  #mutate(Infection = ifelse(infection2 == 1, "Infected (adult)", Infection)) %>%
  mutate(Sex = ifelse(sex == 1, "Male", "Female")) %>%
  select(-sex, -inf, -inbr)

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
  cbind(newdata[, 1:4]) %>%
  mutate(MLH = as.factor(inbr)) %>%
  mutate(Infection = ifelse(inf == 1, "Infected as a cub", "Uninfected life")) %>%
  #mutate(Infection = ifelse(infection2 == 1, "Infected (adult)", Infection)) %>%
  mutate(Sex = ifelse(sex == 1, "Male", "Female")) %>%
  select(-sex, -inf, -inbr)

## save outputs
saveRDS(postSurv, "outputs/PostSurv_FullModel_inbrCONT_levels_zAll.rds")
saveRDS(postMort, "outputs/PostMort_FullModel_inbrCONT_levels_zAll.rds")

## mortality plot
p1 <- mutate(postMort, Infection = factor(Infection)) %>%
  mutate(Infection = factor(Infection, levels = rev(levels(Infection)))) %>%
  ggplot(aes(x = t)) +
  geom_line(aes(y = Median, colour = MLH)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = MLH), alpha = 0.2) +
  facet_grid(Sex ~ Infection) +
  labs(title = "Mortality") + #, subtitle = "tD ~ Sex + Infection on all parameters (Inbreeding and inbr:infection on b2) SG in PH framework") +
  #scale_y_continuous(limits = c(0, 1.25)) +
  ylab("Posterior predictive mortality rate") +
  xlab("Time") +
  coord_cartesian(ylim = c(0, 1.25), xlim = c(0,50))

p1

## survival plot
p2 <- mutate(postSurv, Infection = factor(Infection)) %>%
  mutate(Infection = factor(Infection, levels = rev(levels(Infection)))) %>%
  ggplot(aes(x = t)) +
  geom_line(aes(y = Median, colour = MLH)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = MLH), alpha = 0.2) +
  facet_grid(Sex ~ Infection) +
  labs(title = "Survival") + # subtitle = "tD ~ Sex + Infection on all parameters (Inbreeding and inbr:infection on b2) SG in PH framework") +
  scale_y_continuous(limits = c(0, 1)) +
  ylab("Posterior predictive probability") +
  xlab("Time")

p <- p1 / p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
p
ggsave("predictivePosteriors_SexInfectionInbrCONTlevels.pdf", p, width = 10, height = 5)

p1 <- mutate(postMort, Infection = factor(Infection)) %>%
  mutate(Infection = factor(Infection, levels = rev(levels(Infection)))) %>%
  ggplot(aes(x = t)) +
  geom_line(aes(y = Median, colour = MLH)) +
  #geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = MLH), alpha = 0.2) +
  facet_grid(Sex ~ Infection) +
  labs(title = "Mortality", subtitle = "tD ~ Sex + Infection on all parameters (Inbreeding and inbr:infection on b2) SG in PH framework") +
  scale_y_continuous(limits = c(0, 0.1)) +
  scale_x_continuous(limits = c(0, 5)) +
  ylab("Posterior predictive mortality rate") +
  xlab("Time")
p1

#####  Posterior plots  ###############################################################################################

## reload samples
post <- readRDS("Outputs/SexInfall_Inbrb2_SGa1a2c1_3wayandsexinfremoved_runsamples.rds")

post <- as.matrix(run$samples)
postdf <- as.data.frame(post)

#### INBR
betaINBR <- postdf %>%
  filter(zINBR == 1) %>%
  select(betaINBR)

ci_hdi <- ci(betaINBR, method = "HDI", ci = 0.9)
ci_eti <- ci(betaINBR, method = "ETI", ci = 0.9)

i <- ggplot(betaINBR, aes(x = betaINBR)) +
      geom_density(colour = "blue") +
      geom_vline(xintercept = ci_hdi$CI_low, colour = "red") +
      geom_vline(xintercept = ci_hdi$CI_high, colour = "red")
      #geom_vline(xintercept = ci_eti$CI_low, colour = "green") +
      #geom_vline(xintercept = ci_eti$CI_high, colour = "green")
i
#### INBRINF
betaINBRINF <- postdf %>%
  filter(zINBRINF == 1) %>%
  select(betaINBRINF)

ci_hdi <- ci(betaINBRINF, method = "HDI", ci = 0.9)
ci_eti <- ci(betaINBRINF, method = "ETI", ci = 0.9)

ii <- ggplot(betaINBRINF, aes(x = betaINBRINF)) +
  geom_density(colour = "blue") +
  geom_vline(xintercept = ci_hdi$CI_low, colour = "red") +
  geom_vline(xintercept = ci_hdi$CI_high, colour = "red")
  #geom_vline(xintercept = ci_eti$CI_low, colour = "green") +
  #geom_vline(xintercept = ci_eti$CI_high, colour = "green")

#### INBRSEX
betaINBRSEX <- postdf %>%
  filter(zINBRSEX == 1) %>%
  select(betaINBRSEX)

ci_hdi <- ci(betaINBRSEX, method = "HDI", ci = 0.9)
ci_eti <- ci(betaINBRSEX, method = "ETI", ci = 0.9)

is <- ggplot(betaINBRSEX, aes(x = betaINBRSEX)) +
  geom_density(colour = "blue") +
  geom_vline(xintercept = ci_hdi$CI_low, colour = "red") +
  geom_vline(xintercept = ci_hdi$CI_high, colour = "red")
  #geom_vline(xintercept = ci_eti$CI_low, colour = "green")
  #geom_vline(xintercept = ci_eti$CI_high, colour = "green")

#### INBRINFSEX
betaINBRINFSEX <- postdf %>%
  filter(zINBRINFSEX == 1) %>%
  select(betaINBRINFSEX)

ci_hdi <- ci(betaINBRINFSEX, method = "HDI", ci = 0.9)
ci_eti <- ci(betaINBRINFSEX, method = "ETI", ci = 0.9)

iis <- ggplot(betaINBRINFSEX, aes(x = betaINBRINFSEX)) +
  geom_density(colour = "blue") +
  geom_vline(xintercept = ci_hdi$CI_low, colour = "red") +
  geom_vline(xintercept = ci_hdi$CI_high, colour = "red")
  #geom_vline(xintercept = ci_eti$CI_low, colour = "green")
  #geom_vline(xintercept = ci_eti$CI_high, colour = "green")

iis

(i|ii)/(is|iis)
