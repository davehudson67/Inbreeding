## load libraries
library(tidyverse)
library(coda)
library(patchwork)
library(ggplot2)

## source Siler functions
source("../Distributions/Dist_Siler.R")

## load data
load("../Data/badgerSexInb_ICubvUninf_SG.RData")

## load posterior samples
run <- readRDS("Outputs/SilerFullModel_NULLtD_runsamples.rds")

post <- as.matrix(run$samples)
#post <- sample_n(as.data.frame(post), 10000)
post <- as.data.frame(post)

## extract samples for inferred lifespan
tD <- post %>%
  select(starts_with("tD"))

## take median
tD <- apply(tD, 2, median)

## combine data
CH$tD <- tD

inferredAgeatDeath <- CH %>%
  select(hom, tD) %>%
  filter(tD > 24)

#plot(x = inferredAgeatDeath$hom, y = inferredAgeatDeath$tD, )
ggplot(inferredAgeatDeath, aes(x = hom, y = tD)) +
  geom_point()

CH$Sex <- CH$sex
CH$Infection <- CH$infected_as_cub
CH$Sex <- ifelse(CH$Sex == "FEMALE", "Female", "Male")
CH$Infection <- ifelse(CH$Infection == 1, "Test positive as a cub", "Never positive")

## plot
p <- CH %>%
  filter(tD > 24) %>%
  ggplot(aes(x = hom, y = tD)) +
  geom_point() +
  facet_grid(Sex ~ Infection) +
  scale_y_continuous(limits = c(24, 56), breaks = c(24, 32, 40, 48, 56))

p

#saveRDS(CH, "CH_inferredDeath.rds")


##
## generate median lifespans for different category of badger
##

## load posterior samples
run <- readRDS("Outputs/SilerFullModel_InbrCONT_runsamples.rds")
post <- as.matrix(run$samples)
post <- sample_n(as.data.frame(post), 2000)
post <- as.data.frame(post)

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

post <- as.matrix(post)

####################################################################################################################################################################################
## create points to predict to for different groups
newdata <- expand.grid(t = 0:50, sex = 0:1, infection = 0:1, inbr = seq(0.02, 0.7, 0.2))  #
newdata$inbrsex <- newdata$sex * newdata$inbr
newdata$inbrinf <- newdata$infection * newdata$inbr
newdata$inbrinfsex <- newdata$sex * newdata$infection * newdata$inbr

## loop over Siler parameters
parnms <- c("a1", "a2", "b1", "b2", "c1")
pars <- list()
for(i in 1:5) {
  ## extract samples for linear predictor for e.g. a1
  ## making sure they line up with newdata
  ## (IN PARTICULAR INFECTION[, 2] AND [, 3] SEEM TO BE IN A DIFFERENT ORDER TO THE BETAS
  ##  WHICH IS FINE BUT THEY NEED TO MATCH UP - I'M CALLING INFECTION1 CUBS AND INFECTION2 ADULTS)
  parPost <- post[, match(paste0("beta", c("SEX", "INF", "INBR", "INBRSEX", "INBRINF", "INBRINFSEX"), "[", i, "]"), colnames(post))]
  
  ## extract inclusion samples
  zpost <- post[, match(paste0("z", c("SEX", "INF", "INBR", "INBRSEX", "INBRINF", "INBRINFSEX"), "[", i, "]"), colnames(post))]
  
  ## now multiply by relevant z indexes
  ## (expanding zpost matrix to make elementwise
  ## multiplication straightforward)
  zpost <- cbind(
    zpost[, paste0("zSEX[", i, "]")], ## Sex
    zpost[, paste0("zINF[", i, "]")], ## infection
    zpost[, paste0("zINBR[", i, "]")], ## inbreeding
    zpost[, paste0("zINBRSEX[", i, "]")], ## inf1/sex
    zpost[, paste0("zINBRINF[", i, "]")], ## inf1/inbreeding
    zpost[, paste0("zINBRINFSEX[", i, "]")] ## sex/inbreeding
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

pars <- as.data.frame(pars)
colnames(pars) <- c("t", "a1", "a2", "b1", "b2", "c")  
pars <- pars %>%
  filter(t == 1) %>%
  select(-t)

## Siler sampling (calculate on the log-scale to try and avoid rounding errors)
postLifespan <- qSiler(rep(0.5, nrow(pars)), pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], zL = 24)

newdata <- newdata %>%
  filter(t == 1) %>%
  #slice(rep(1:n(), each = 8000)) %>%
  select(-t)

## posterior predictive summaries
postLifespan1 <- postLifespan %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025),#, na.rm = TRUE),
      Median = median(x),
      UCI = quantile(x, probs = 0.975))#, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  cbind(newdata[, 1:3]) %>%
  mutate(MLH = inbr) %>%
  mutate(Infection = ifelse(infection == 1, "Test positive as a cub", "Never positive")) %>%
  #  mutate(Infection = ifelse(infection2 == 1, "Infected (adult)", Infection)) %>%
  mutate(Sex = ifelse(sex == 1, "Male", "Female")) %>%
  select(-sex, -infection, -inbr)

median(CH$hom)

## plot
p3_data <- mutate(postLifespan1, Infection = factor(Infection)) %>%
  mutate(Infection = factor(Infection, levels = rev(levels(Infection)))) %>%
  unite(group, c(Infection, Sex), sep = " ", remove = FALSE)

p1 <- p +  
  geom_line(data = p3_data, aes(x = MLH, y = Median, colour = "redbrick")) +
  ylab("Lifespan (yrs)") +
  xlab("Multi-locus Homozygosity") + 
  theme(text = element_text(size = 16)) +
  facet_grid(Infection ~ Sex) +
  theme(legend.position="none")
  
##################################################################################

## Add categorical

run <- readRDS("Outputs/SilerFullModel_InbrCAT_runsamples.rds")
post <- as.matrix(run$samples)
post <- sample_n(as.data.frame(post), 2000)

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

post <- as.matrix(post)

####################################################################################################################################################################################
## create points to predict to for different groups
newdata <- expand.grid(t = 0:50, sex = 0:1, infection = 0:1, inbr = 0:1)  #seq(0.2, 0.8, 0.2)
newdata$inbrsex <- newdata$sex * newdata$inbr
newdata$inbrinf <- newdata$infection * newdata$inbr
newdata$inbrinfsex <- newdata$sex * newdata$infection * newdata$inbr

## loop over Siler parameters
parnms <- c("a1", "a2", "b1", "b2", "c1")
pars <- list()
for(i in 1:5) {
  ## extract samples for linear predictor for e.g. a1
  ## making sure they line up with newdata
  ## (IN PARTICULAR INFECTION[, 2] AND [, 3] SEEM TO BE IN A DIFFERENT ORDER TO THE BETAS
  ##  WHICH IS FINE BUT THEY NEED TO MATCH UP - I'M CALLING INFECTION1 CUBS AND INFECTION2 ADULTS)
  parPost <- post[, match(paste0("beta", c("SEX", "INF", "INBR", "INBRSEX", "INBRINF", "INBRINFSEX"), "[", i, "]"), colnames(post))]
  
  ## extract inclusion samples
  zpost <- post[, match(paste0("z", c("SEX", "INF", "INBR", "INBRSEX", "INBRINF", "INBRINFSEX"), "[", i, "]"), colnames(post))]
  
  ## now multiply by relevant z indexes
  ## (expanding zpost matrix to make elementwise
  ## multiplication straightforward)
  zpost <- cbind(
    zpost[, paste0("zSEX[", i, "]")], ## Sex
    zpost[, paste0("zINF[", i, "]")], ## infection
    zpost[, paste0("zINBR[", i, "]")], ## inbreeding
    zpost[, paste0("zINBRSEX[", i, "]")], ## inf1/sex
    zpost[, paste0("zINBRINF[", i, "]")], ## inf1/inbreeding
    zpost[, paste0("zINBRINFSEX[", i, "]")] ## sex/inbreeding
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

pars <- as.data.frame(pars)
colnames(pars) <- c("t", "a1", "a2", "b1", "b2", "c")  
pars <- pars %>%
  filter(t == 1) %>%
  select(-t)

## Siler sampling (calculate on the log-scale to try and avoid rounding errors)
postLifespan <- qSiler(rep(0.5, nrow(pars)), pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5], zL = 24)

newdata <- newdata %>%
  filter(t == 1) %>%
  select(-t)

## posterior predictive summaries
postLifespan1 <- postLifespan %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025),#, na.rm = TRUE),
      Median = median(x),
      UCI = quantile(x, probs = 0.975))#, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  cbind(newdata[, 1:3]) %>%
  mutate(Inbreeding = ifelse(inbr == 0, "Outbred", "Inbred")) %>%
  mutate(Infection = ifelse(infection == 1, "Test positive as a cub", "Never positive")) %>%
  #  mutate(Infection = ifelse(infection2 == 1, "Infected (adult)", Infection)) %>%
  mutate(Sex = ifelse(sex == 1, "Male", "Female")) %>%
  select(-sex, -infection, -inbr)

## plot
p2_data <- mutate(postLifespan1, Infection = factor(Infection)) %>%
  mutate(Infection = factor(Infection, levels = rev(levels(Infection)))) %>%
  unite(group, c(Inbreeding, Infection, Sex), sep = " ", remove = FALSE)

p2 <- p +
  geom_hline(data = p2_data, aes(yintercept = Median, colour = Inbreeding)) +
  #  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = MLH), alpha = 0.2) +
  #  facet_grid(Sex ~ Infection) +
  #  labs(title = "Survival", subtitle = "tD ~ Sex * Infection + Inbreeding:Sex + Inbreeding:Infection") +
  #  scale_y_continuous(limits = c(0, 1)) +
  ylab("Lifespan (1/4 yrs)") +
  xlab("Multi-locus Homozygosity") +
  facet_grid(Infection ~ Sex)

## split p2 data
p2_inbred <- p2_data %>%
  filter(Inbreeding == "Inbred")
p2_outbred <- p2_data %>%
  filter(Inbreeding == "Outbred")

p3 <- p1 +
  #geom_hline(data = p2_data, aes(yintercept = Median, colour = Inbreeding)) +
  geom_segment(data = p2_outbred, aes(x = 0.02, xend = 0.318, y = Median, yend = Median, colour = Inbreeding)) +
  geom_segment(data = p2_inbred, aes(x = 0.318, xend = 0.6, y = Median, yend = Median, colour = Inbreeding)) +
  #  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = MLH), alpha = 0.2) +
  #  facet_grid(Sex ~ Infection) +
  #  labs(title = "Survival", subtitle = "tD ~ Sex * Infection + Inbreeding:Sex + Inbreeding:Infection") +
  #  scale_y_continuous(limits = c(0, 1)) +
  ylab("Lifespan (yrs)") +
  xlab("Multi-locus Homozygosity") +
  facet_grid(Infection ~ Sex) +
  scale_colour_manual(values = c("#F8766D", "#619CFF", "#00BA38"), drop = FALSE)
p3
p1
median(CH$hom)
