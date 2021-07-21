library(tidyverse)
library(nimble)

load("Data/badgerSexInb_AdultInfvCubInfvUninf.RData")
samples <- readRDS("outputs/Sex3InfectionInbrCATmed_AllParameters_runsamples.rds")
samples <- samples$samples[1:20000,]

saveRDS(samples, "posteriorSamples.rds")


names(samples)

colnames(samples) <- c("a1", "a2", "b1", "b2", "betaINBR1", "betaINBR2", "betaINBR3", "betaINBR4", "betaINBR5",  "betaINFADULT1", "betaINFADULT2", 
                       "betaINFADULT3", "betaINFADULT4", "betaINFADULT5", "betaINFCUB1", "betaINFCUB2", "betaINFCUB3", "betaINFCUB4", "betaINFCUB5", 
                       "betaINFINBRADULT1", "betaINFINBRADULT2", "betaINFINBRADULT3", "betaINFINBRADULT4", "betaINFINBRADULT5", "betaINFINBRCUB1", "betaINFINBRCUB2", 
                       "betaINFINBRCUB3", "betaINFINBRCUB4", "betaINFINBRCUB5",  "betaSex1", "betaSex2", "betaSex3", "betaSex4", "betaSex5",
                       "betaSexINBR1", "betaSexINBR2", "betaSexINBR3", "betaSexINBR4", "betaSexINBR5", "c1", "mean.p", "zINBR1", "zINBR2", "zINBR3", "zINBR4", "zINBR5",
                       "zINFINBR1", "zINFINBR2", "zINFINBR3", "zINFINBR4", "zINFINBR5", "zSEXINBR1", "zSEXINBR2", "zSEXINBR3", "zSEXINBR4", "zSEXINBR5") 

source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_Siler.R")
source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_SilerNim.R")
source("../SimulationStudy/FirstPaperFiles/ModelComparison_FUNCTIONS.R")

t_pred <- seq(from = 0, to = max(cint), length.out = 50)


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


## group 1 - outbred uninfected males
dens1 <- samples %>%
  filter(zINBR1 == 0, zINBR2 == 0, zINBR3 == 0, zINBR4 == 0, zINBR5 == 0,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
    apply(1, function(pars, t) {
    dSiler(t, a1 = pars[1] * exp(pars[30]), 
              a2 = pars[2] * exp(pars[31]),
              b1 = pars[3] * exp(pars[32]),
              b2 = pars[4] * exp(pars[33]),
              c1 = pars[40] * exp(pars[34]))
           }, t = t_pred)

surv1 <- samples %>%
  filter(zINBR1 == 0, zINBR2 == 0, zINBR3 == 0, zINBR4 == 0, zINBR5 == 0,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[1] * exp(pars[30]), 
           a2 = pars[2] * exp(pars[31]),
           b1 = pars[3] * exp(pars[32]),
           b2 = pars[4] * exp(pars[33]),
           c1 = pars[40] * exp(pars[34]), lower.tail = FALSE)
  }, t = t_pred)

mort1 <- dens1/surv1

survival1 <- surv1 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Outbred, Uninfected, Male")

mortality1 <- mort1 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Outbred, Uninfected, Male")

#################
## Group 2 - Outbred infected as cub male

dens2 <- samples %>%
  filter(zINBR1 == 0, zINBR2 == 0, zINBR3 == 0, zINBR4 == 0, zINBR5 == 0,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    dSiler(t, a1 = pars[1] * exp(pars[30] * pars[15]), 
           a2 = pars[2] * exp(pars[31] * pars[16]),
           b1 = pars[3] * exp(pars[32] * pars[17]),
           b2 = pars[4] * exp(pars[33] * pars[18]),
           c1 = pars[40] * exp(pars[34] * pars[19]))
  }, t = t_pred)

surv2 <- samples %>%
  filter(zINBR1 == 0, zINBR2 == 0, zINBR3 == 0, zINBR4 == 0, zINBR5 == 0,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[1] * exp(pars[30] * pars[15]), 
           a2 = pars[2] * exp(pars[31] * pars[16]),
           b1 = pars[3] * exp(pars[32] * pars[17]),
           b2 = pars[4] * exp(pars[33] * pars[18]),
           c1 = pars[40] * exp(pars[34] * pars[19]), lower.tail = FALSE)
  }, t = t_pred)

mort2 <- dens2/surv2

survival2 <- surv2 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Outbred, Infected_as_cub, Male")

mortality2 <- mort2 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Outbred, Infected_as_cub, Male")

#################
## Group 3 - Outbred infected as adult male

dens3 <- samples %>%
  filter(zINBR1 == 0, zINBR2 == 0, zINBR3 == 0, zINBR4 == 0, zINBR5 == 0,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    dSiler(t, a1 = pars[1] * exp(pars[30] * pars[10]), 
           a2 = pars[2] * exp(pars[31] * pars[11]),
           b1 = pars[3] * exp(pars[32] * pars[12]),
           b2 = pars[4] * exp(pars[33] * pars[13]),
           c1 = pars[40] * exp(pars[34] * pars[14]))
  }, t = t_pred)

surv3 <- samples %>%
  filter(zINBR1 == 0, zINBR2 == 0, zINBR3 == 0, zINBR4 == 0, zINBR5 == 0,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[1] * exp(pars[30] * pars[10]), 
           a2 = pars[2] * exp(pars[31] * pars[11]),
           b1 = pars[3] * exp(pars[32] * pars[12]),
           b2 = pars[4] * exp(pars[33] * pars[13]),
           c1 = pars[40] * exp(pars[34] * pars[14]), lower.tail = FALSE)
  }, t = t_pred)

which(surv3 == 0)


mort3 <- dens3/surv3

survival3 <- surv3 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Outbred, Infected_as_adult, Male")

mortality3 <- mort3 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Outbred, Infected_as_adult, Male")

#################
## Group 4 - Outbred uninfected female

dens4 <- samples %>%
  filter(zINBR1 == 0, zINBR2 == 0, zINBR3 == 0, zINBR4 == 0, zINBR5 == 0,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    dSiler(t, a1 = pars[1], 
           a2 = pars[2],
           b1 = pars[3],
           b2 = pars[4],
           c1 = pars[40])
  }, t = t_pred)

surv4 <- samples %>%
  filter(zINBR1 == 0, zINBR2 == 0, zINBR3 == 0, zINBR4 == 0, zINBR5 == 0,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[1], 
           a2 = pars[2],
           b1 = pars[3],
           b2 = pars[4],
           c1 = pars[40], lower.tail = FALSE)
  }, t = t_pred)

mort4 <- dens4/surv4

survival4 <- surv4 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Outbred, Uninfected, Female")

mortality4 <- mort4 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Outbred, Uninfected, Female")

#################
## Group 5 - Outbred infected as cub female

dens5 <- samples %>%
  filter(zINBR1 == 0, zINBR2 == 0, zINBR3 == 0, zINBR4 == 0, zINBR5 == 0,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    dSiler(t, a1 = pars[1] * exp(pars[15]), 
           a2 = pars[2] * exp(pars[16]),
           b1 = pars[3] * exp(pars[17]),
           b2 = pars[4] * exp(pars[18]),
           c1 = pars[40] * exp(pars[19]))
  }, t = t_pred)

surv5 <- samples %>%
  filter(zINBR1 == 0, zINBR2 == 0, zINBR3 == 0, zINBR4 == 0, zINBR5 == 0,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[1] * exp(pars[15]), 
           a2 = pars[2] * exp(pars[16]),
           b1 = pars[3] * exp(pars[17]),
           b2 = pars[4] * exp(pars[18]),
           c1 = pars[40] * exp(pars[19]), lower.tail = FALSE)
  }, t = t_pred)

mort5 <- dens5/surv5

survival5 <- surv5 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Outbred, Infected_as_cub, Female")

mortality5 <- mort5 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Outbred, Infected_as_cub, Female")

#################
## Group 6 - Outbred infected as adult female

dens6 <- samples %>%
  filter(zINBR1 == 0, zINBR2 == 0, zINBR3 == 0, zINBR4 == 0, zINBR5 == 0,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    dSiler(t, a1 = pars[1] * exp(pars[10]), 
           a2 = pars[2] * exp(pars[11]),
           b1 = pars[3] * exp(pars[12]),
           b2 = pars[4] * exp(pars[13]),
           c1 = pars[40] * exp(pars[14]))
  }, t = t_pred)

surv6 <- samples %>%
  filter(zINBR1 == 0, zINBR2 == 0, zINBR3 == 0, zINBR4 == 0, zINBR5 == 0,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[1] * exp(pars[10]), 
           a2 = pars[2] * exp(pars[11]),
           b1 = pars[3] * exp(pars[12]),
           b2 = pars[4] * exp(pars[13]),
           c1 = pars[40] * exp(pars[14]), lower.tail = FALSE)
  }, t = t_pred)

mort6 <- dens6/surv6

survival6 <- surv6 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Outbred, Infected_as_adult, Female")

mortality6 <- mort6 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Outbred, Infected_as_adult, Female")

## Group 7 - Inbred uninfected male

dens7 <- samples %>%
  filter(zINBR1 == 1, zINBR2 == 1, zINBR3 == 1, zINBR4 == 1, zINBR5 == 1,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    dSiler(t, a1 = pars[1] * exp(pars[5] + pars[30]), 
           a2 = pars[2] * exp(pars[6] + pars[31]),
           b1 = pars[3] * exp(pars[7] + pars[32]),
           b2 = pars[4] * exp(pars[8] + pars[33]),
           c1 = pars[40] * exp(pars[9] + pars[34]))
  }, t = t_pred)

surv7 <- samples %>%
  filter(zINBR1 == 1, zINBR2 == 1, zINBR3 == 1, zINBR4 == 1, zINBR5 == 1,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[1] * exp(pars[5] + pars[30]), 
           a2 = pars[2] * exp(pars[6] + pars[31]),
           b1 = pars[3] * exp(pars[7] + pars[32]),
           b2 = pars[4] * exp(pars[8] + pars[33]),
           c1 = pars[40] * exp(pars[9] + pars[34]), lower.tail = FALSE)
  }, t = t_pred)

mort7 <- dens7/surv7

survival7 <- surv7 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Inbred, Uninfected, Male")

mortality7 <- mort7 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Outbred, Uninfected, Male")

##########################
## Group 8 - Inbred infected as cub male

dens8 <- samples %>%
  filter(zINBR1 == 1, zINBR2 == 1, zINBR3 == 1, zINBR4 == 1, zINBR5 == 1,
         zSEXINBR1 == 1, zSEXINBR2 == 1, zSEXINBR3 == 1, zSEXINBR4 == 1, zSEXINBR5 == 1,
         zINFINBR1 == 1, zINFINBR2 == 1, zINFINBR3 == 1, zINFINBR4 == 1, zINFINBR5 == 1) %>%
  apply(1, function(pars, t) {
    dSiler(t, a1 = pars[1] * exp(pars[5] + pars[30] + pars[25] + pars[15] + pars[35]), 
           a2 = pars[2] * exp(pars[6] + pars[31] + pars[26]+ pars[16] + pars[36]),
           b1 = pars[3] * exp(pars[7] + pars[32] + pars[27]+ pars[17] + pars[37]),
           b2 = pars[4] * exp(pars[8] + pars[33] + pars[28]+ pars[18] + pars[38]),
           c1 = pars[40] * exp(pars[9] + pars[34] + pars[29]+ pars[19] + pars[39]))
  }, t = t_pred)

surv8 <- samples %>%
  filter(zINBR1 == 1, zINBR2 == 1, zINBR3 == 1, zINBR4 == 1, zINBR5 == 1,
         zSEXINBR1 == 1, zSEXINBR2 == 1, zSEXINBR3 == 1, zSEXINBR4 == 1, zSEXINBR5 == 1,
         zINFINBR1 == 1, zINFINBR2 == 1, zINFINBR3 == 1, zINFINBR4 == 1, zINFINBR5 == 1) %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[1] * exp(pars[5] + pars[30] + pars[25] + pars[15] + pars[35]), 
           a2 = pars[2] * exp(pars[6] + pars[31] + pars[26]+ pars[16] + pars[36]),
           b1 = pars[3] * exp(pars[7] + pars[32] + pars[27]+ pars[17] + pars[37]),
           b2 = pars[4] * exp(pars[8] + pars[33] + pars[28]+ pars[18] + pars[38]),
           c1 = pars[40] * exp(pars[9] + pars[34] + pars[29]+ pars[19] + pars[39]), lower.tail = FALSE)
  }, t = t_pred)

mort8 <- dens8/surv8

survival8 <- surv8 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Inbred, infected as cub, Male")

mortality8 <- mort8 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Inbred, infected as cub, Male")

##########################
## Group 9 - Inbred infected as adult male

dens9 <- samples %>%
  filter(zINBR1 == 1, zINBR2 == 1, zINBR3 == 1, zINBR4 == 1, zINBR5 == 1,
         zSEXINBR1 == 1, zSEXINBR2 == 1, zSEXINBR3 == 1, zSEXINBR4 == 1, zSEXINBR5 == 1,
         zINFINBR1 == 1, zINFINBR2 == 1, zINFINBR3 == 1, zINFINBR4 == 1, zINFINBR5 == 1) %>%
  apply(1, function(pars, t) {
    dSiler(t, a1 = pars[1] * exp(pars[5] + pars[30] + pars[20] + pars[10] + pars[35]), 
           a2 = pars[2] * exp(pars[6] + pars[31] + pars[21]+ pars[11] + pars[36]),
           b1 = pars[3] * exp(pars[7] + pars[32] + pars[22]+ pars[12] + pars[37]),
           b2 = pars[4] * exp(pars[8] + pars[33] + pars[23]+ pars[13] + pars[38]),
           c1 = pars[40] * exp(pars[9] + pars[34] + pars[24]+ pars[14] + pars[39]))
  }, t = t_pred)

surv9 <- samples %>%
  filter(zINBR1 == 1, zINBR2 == 1, zINBR3 == 1, zINBR4 == 1, zINBR5 == 1,
         zSEXINBR1 == 1, zSEXINBR2 == 1, zSEXINBR3 == 1, zSEXINBR4 == 1, zSEXINBR5 == 1,
         zINFINBR1 == 1, zINFINBR2 == 1, zINFINBR3 == 1, zINFINBR4 == 1, zINFINBR5 == 1) %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[1] * exp(pars[5] + pars[30] + pars[20] + pars[10] + pars[35]), 
           a2 = pars[2] * exp(pars[6] + pars[31] + pars[21]+ pars[11] + pars[36]),
           b1 = pars[3] * exp(pars[7] + pars[32] + pars[22]+ pars[12] + pars[37]),
           b2 = pars[4] * exp(pars[8] + pars[33] + pars[23]+ pars[13] + pars[38]),
           c1 = pars[40] * exp(pars[9] + pars[34] + pars[24]+ pars[14] + pars[39]), lower.tail = FALSE)
  }, t = t_pred)

mort9 <- dens9/surv9

survival9 <- surv9 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Inbred, infected as adult, Male")

mortality9 <- mort9 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Inbred, infected as adult, Male")

##########################
## Group 10 - Inbred uninfected female

dens10 <- samples %>%
  filter(zINBR1 == 1, zINBR2 == 1, zINBR3 == 1, zINBR4 == 1, zINBR5 == 1,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    dSiler(t, a1 = pars[1] * exp(pars[5]), 
           a2 = pars[2] * exp(pars[6]),
           b1 = pars[3] * exp(pars[7]),
           b2 = pars[4] * exp(pars[8]),
           c1 = pars[40] * exp(pars[9]))
  }, t = t_pred)

surv10 <- samples %>%
  filter(zINBR1 == 1, zINBR2 == 1, zINBR3 == 1, zINBR4 == 1, zINBR5 == 1,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[1] * exp(pars[5]), 
           a2 = pars[2] * exp(pars[6]),
           b1 = pars[3] * exp(pars[7]),
           b2 = pars[4] * exp(pars[8]),
           c1 = pars[40] * exp(pars[9]), lower.tail = FALSE)
  }, t = t_pred)

mort10 <- dens10/surv10

survival10 <- surv10 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Inbred, uninfected, female")

mortality10 <- mort10 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Inbred, uninfected, female")

##########################
## Group 11 - Inbred infected as cub female

dens11 <- samples %>%
  filter(zINBR1 == 1, zINBR2 == 1, zINBR3 == 1, zINBR4 == 1, zINBR5 == 1,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 1, zINFINBR2 == 1, zINFINBR3 == 1, zINFINBR4 == 1, zINFINBR5 == 1) %>%
  apply(1, function(pars, t) {
    dSiler(t, a1 = pars[1] * exp(pars[5] + pars[25] + pars[15]), 
           a2 = pars[2] * exp(pars[6] + pars[26]+ pars[16]),
           b1 = pars[3] * exp(pars[7] + pars[27]+ pars[17]),
           b2 = pars[4] * exp(pars[8] + pars[28]+ pars[18]),
           c1 = pars[40] * exp(pars[9] + pars[29]+ pars[19]))
  }, t = t_pred)

surv11 <- samples %>%
  filter(zINBR1 == 1, zINBR2 == 1, zINBR3 == 1, zINBR4 == 1, zINBR5 == 1,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 1, zINFINBR2 == 1, zINFINBR3 == 1, zINFINBR4 == 1, zINFINBR5 == 1) %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[1] * exp(pars[5] + pars[25] + pars[15]), 
           a2 = pars[2] * exp(pars[6] + pars[26]+ pars[16]),
           b1 = pars[3] * exp(pars[7] + pars[27]+ pars[17]),
           b2 = pars[4] * exp(pars[8] + pars[28]+ pars[18]),
           c1 = pars[40] * exp(pars[9] + pars[29]+ pars[19]), lower.tail = FALSE)
  }, t = t_pred)

mort11 <- dens11/surv11

survival11 <- surv11 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Inbred, infected as cub, female")

mortality11 <- mort11 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Inbred, infected as cub, feale")

#####################################################
## Group 12 - Inbred infected as adult female

dens12 <- samples %>%
  filter(zINBR1 == 1, zINBR2 == 1, zINBR3 == 1, zINBR4 == 1, zINBR5 == 1,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 1, zINFINBR2 == 1, zINFINBR3 == 1, zINFINBR4 == 1, zINFINBR5 == 1) %>%
  apply(1, function(pars, t) {
    dSiler(t, a1 = pars[1] * exp(pars[5] + pars[20] + pars[10]), 
           a2 = pars[2] * exp(pars[6] + pars[21]+ pars[11]),
           b1 = pars[3] * exp(pars[7] + pars[22]+ pars[12]),
           b2 = pars[4] * exp(pars[8] + pars[23]+ pars[13]),
           c1 = pars[40] * exp(pars[9] + pars[24]+ pars[14]))
  }, t = t_pred)

surv12 <- samples %>%
  filter(zINBR1 == 1, zINBR2 == 1, zINBR3 == 1, zINBR4 == 1, zINBR5 == 1,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 1, zINFINBR2 == 1, zINFINBR3 == 1, zINFINBR4 == 1, zINFINBR5 == 1) %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[1] * exp(pars[5] + pars[20] + pars[10]), 
           a2 = pars[2] * exp(pars[6] + pars[21]+ pars[11]),
           b1 = pars[3] * exp(pars[7] + pars[22]+ pars[12]),
           b2 = pars[4] * exp(pars[8] + pars[23]+ pars[13]),
           c1 = pars[40] * exp(pars[9] + pars[24]+ pars[14]), lower.tail = FALSE)
  }, t = t_pred)

mort12 <- dens12/surv12

survival12 <- surv12 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Inbred, infected as adult, female")

mortality12 <- mort12 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Inbred, infected as adult, female")








ggplot(survival3) +
  geom_line(aes(x=t, y=Median))
     
ggplot(mortality2) +
  geom_line(aes(x=t, y=Median))


















dens2 <- samples %>%
  filter(zINBR1 == 0, zINBR2 == 0, zINBR3 == 0, zINBR4 == 0, zINBR5 == 0,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    dSiler(t, a1 = pars[1] + exp(pars[30] + pars[15]), 
           a2 = pars[2] + exp(pars[31] + pars[16]),
           b1 = pars[3] + exp(pars[32] + pars[17]),
           b2 = pars[4] + exp(pars[33] + pars[18]),
           c1 = pars[40] + exp(pars[34] + pars[19]))
  }, t = t_pred)

surv1 <- samples %>%
  filter(zINBR1 == 0, zINBR2 == 0, zINBR3 == 0, zINBR4 == 0, zINBR5 == 0,
         zSEXINBR1 == 0, zSEXINBR2 == 0, zSEXINBR3 == 0, zSEXINBR4 == 0, zSEXINBR5 == 0,
         zINFINBR1 == 0, zINFINBR2 == 0, zINFINBR3 == 0, zINFINBR4 == 0, zINFINBR5 == 0) %>%
  apply(1, function(pars, t) {
    pSiler(t, a1 = pars[1] + exp(pars[30]), 
           a2 = pars[2] + exp(pars[31]),
           b1 = pars[3] + exp(pars[32]),
           b2 = pars[4] + exp(pars[33]),
           c1 = pars[40] + exp(pars[34]), lower.tail = FALSE)
  }, t = t_pred)

survival1 <- surv1 %>%
  apply(1, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(t = t_pred, ) %>%
  mutate(category = "Outbred, Uninfected, Male")

ggplot(survival1) +
  geom_line(aes(x=t, y=Median))

