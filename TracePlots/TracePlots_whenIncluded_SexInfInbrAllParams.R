library(tidyverse)
library(patchwork)
library(MCMCvis)
library(mcmcplots)

samples <- readRDS("outputs/Sex3InfectionInbrCATmed_AllParameters_runsamples.rds")

mcmcplot(samples$samples, parms = c("a1", "a2", "b1", "b2", "c1"))

samples <- as.matrix(samples$samples, chains = TRUE)
samples <- as.data.frame(samples)

names(samples)
colnames(samples) <- c("CHAIN", "a1", "a2", "b1", "b2", "betaINBR1", "betaINBR2", "betaINBR3", "betaINBR4", "betaINBR5",  "betaINFADULT1", "betaINFADULT2", 
                       "betaINFADULT3", "betaINFADULT4", "betaINFADULT5", "betaINFCUB1", "betaINFCUB2", "betaINFCUB3", "betaINFCUB4", "betaINFCUB5", 
                       "betaINFINBRADULT1", "betaINFINBRADULT2", "betaINFINBRADULT3", "betaINFINBRADULT4", "betaINFINBRADULT5", "betaINFINBRCUB1", "betaINFINBRCUB2", 
                       "betaINFINBRCUB3", "betaINFINBRCUB4", "betaINFINBRCUB5",  "betaSex1", "betaSex2", "betaSex3", "betaSex4", "betaSex5",
                       "betaSexINBR1", "betaSexINBR2", "betaSexINBR3", "betaSexINBR4", "betaSexINBR5", "c1", "mean.p", "zINBR1", "zINBR2", "zINBR3", "zINBR4", "zINBR5",
                       "zINFINBR1", "zINFINBR2", "zINFINBR3", "zINFINBR4", "zINFINBR5", "zSEXINBR1", "zSEXINBR2", "zSEXINBR3", "zSEXINBR4", "zSEXINBR5") 

samples$CHAIN <- as.factor(samples$CHAIN)

##  beta INBR ## 

betaINBR1 <- samples %>%
  filter(zINBR1 == 1) %>%
  select(CHAIN, betaINBR1) %>%
  mutate(iteration = 1:length(betaINBR1))

betaINBR2 <- samples %>%
  filter(zINBR2 == 1) %>%
  select(CHAIN, betaINBR2) %>%
  mutate(iteration = 1:length(betaINBR2))

betaINBR3 <- samples %>%
  filter(zINBR3 == 1) %>%
  select(CHAIN, betaINBR3) %>%
  mutate(iteration = 1:length(betaINBR3))

betaINBR4 <- samples %>%
  filter(zINBR4 == 1) %>%
  select(CHAIN, betaINBR4) %>%
  mutate(iteration = 1:length(betaINBR4))

betaINBR5 <- samples %>%
  filter(zINBR5 == 1) %>%
  select(CHAIN, betaINBR5) %>%
  mutate(iteration = 1:length(betaINBR5))

betaINBR1 <- as.data.frame(betaINBR1)
b1 <- ggplot(betaINBR1) +
  geom_line(aes(y = betaINBR1, x = iteration, colour = CHAIN))
b1d <- ggplot(betaINBR1) +
  geom_density(aes(x = betaINBR1, colour = CHAIN))

betaINBR2 <- as.data.frame(betaINBR2)
betaINBR2$iteration <- seq(1:nrow(betaINBR2))
b2 <- ggplot(betaINBR2) +
  geom_line(aes(y = betaINBR2, x = iteration, colour = CHAIN))
b2d <- ggplot(betaINBR2) +
  geom_density(aes(x = betaINBR2, colour = CHAIN))

betaINBR3 <- as.data.frame(betaINBR3)
betaINBR3$iteration <- seq(1:nrow(betaINBR3))
b3 <- ggplot(betaINBR3) +
  geom_line(aes(y = betaINBR3, x = iteration, colour = CHAIN))
b3d <- ggplot(betaINBR3) +
  geom_density(aes(x = betaINBR3, colour = CHAIN))

betaINBR4 <- as.data.frame(betaINBR4)
betaINBR4$iteration <- seq(1:nrow(betaINBR4))
b4 <- ggplot(betaINBR4) +
  geom_line(aes(y = betaINBR4, x = iteration, colour = CHAIN))
b4d <- ggplot(betaINBR4) +
  geom_density(aes(x = betaINBR4, colour = CHAIN))


betaINBR5 <- as.data.frame(betaINBR5)
betaINBR5$iteration <- seq(1:nrow(betaINBR5))
b5 <- ggplot(betaINBR5) +
  geom_line(aes(y = betaINBR5, x = iteration, colour = CHAIN))
b5d <- ggplot(betaINBR5) +
  geom_density(aes(x = betaINBR5, colour = CHAIN))

## plot

b1/b2/b3|b1d/b2d/b3d
b4/b5|b4d/b5d

###############################################################################
## beta Inf:Inbr

betaINFINBRADULT1 <- samples %>%
  filter(zINFINBR1 == 1) %>%
  select(CHAIN, betaINFINBRADULT1) %>%
  mutate(iteration = 1:length(betaINFINBRADULT1))

betaINFINBRADULT2 <- samples %>%
  filter(zINFINBR2 == 1) %>%
  select(CHAIN, betaINFINBRADULT2) %>%
  mutate(iteration = 1:length(betaINFINBRADULT2))

betaINFINBRADULT3 <- samples %>%
  filter(zINFINBR3 == 1) %>%
  select(CHAIN, betaINFINBRADULT3) %>%
  mutate(iteration = 1:length(betaINFINBRADULT3))

betaINFINBRADULT4 <- samples %>%
  filter(zINFINBR4 == 1) %>%
  select(CHAIN, betaINFINBRADULT4) %>%
  mutate(iteration = 1:length(betaINFINBRADULT4))

betaINFINBRADULT5 <- samples %>%
  filter(zINFINBR5 == 1) %>%
  select(CHAIN, betaINFINBRADULT5) %>%
  mutate(iteration = 1:length(betaINFINBRADULT5))

betaINFINBRCUB1 <- samples %>%
  filter(zINFINBR1 == 1) %>%
  select(CHAIN, betaINFINBRCUB1) %>%
  mutate(iteration = 1:length(betaINFINBRCUB1))

betaINFINBRCUB2 <- samples %>%
  filter(zINFINBR2 == 1) %>%
  select(CHAIN, betaINFINBRCUB2) %>%
  mutate(iteration = 1:length(betaINFINBRCUB2))

betaINFINBRCUB3 <- samples %>%
  filter(zINFINBR3 == 1) %>%
  select(CHAIN, betaINFINBRCUB3) %>%
  mutate(iteration = 1:length(betaINFINBRCUB3))

betaINFINBRCUB4 <- samples %>%
  filter(zINFINBR4 == 1) %>%
  select(CHAIN, betaINFINBRCUB4) %>%
  mutate(iteration = 1:length(betaINFINBRCUB4))

betaINFINBRCUB5 <- samples %>%
  filter(zINFINBR5 == 1) %>%
  select(CHAIN, betaINFINBRCUB5) %>%
  mutate(iteration = 1:length(betaINFINBRCUB5))

##

betaINFINBRADULT1 <- as.data.frame(betaINFINBRADULT1)
BIIA1 <- ggplot(betaINFINBRADULT1) +
  geom_line(aes(y = betaINFINBRADULT1, x = iteration, colour = CHAIN))
BIIA1d <- ggplot(betaINFINBRADULT1) +
  geom_density(aes(x = betaINFINBRADULT1, colour = CHAIN))

betaINFINBRADULT2 <- as.data.frame(betaINFINBRADULT2)
betaINFINBRADULT2$iteration <- seq(1:nrow(betaINFINBRADULT2))
BIIA2 <- ggplot(betaINFINBRADULT2) +
  geom_line(aes(y = betaINFINBRADULT2, x = iteration, colour = CHAIN))
BII2Ad <- ggplot(betaINFINBRADULT2) +
  geom_density(aes(x = betaINFINBRADULT2, colour = CHAIN))

betaINFINBRADULT3 <- as.data.frame(betaINFINBRADULT3)
betaINFINBRADULT3$iteration <- seq(1:nrow(betaINFINBRADULT3))
BIIA3 <- ggplot(betaINFINBRADULT3) +
  geom_line(aes(y = betaINFINBRADULT3, x = iteration, colour = CHAIN))
BII3Ad <- ggplot(betaINFINBRADULT3) +
  geom_density(aes(x = betaINFINBRADULT3, colour = CHAIN))

betaINFINBRADULT4 <- as.data.frame(betaINFINBRADULT4)
betaINFINBRADULT4$iteration <- seq(1:nrow(betaINFINBRADULT4))
BIIA4 <- ggplot(betaINFINBRADULT4) +
  geom_line(aes(y = betaINFINBRADULT4, x = iteration, colour = CHAIN))
BII4Ad <- ggplot(betaINFINBRADULT4) +
  geom_density(aes(x = betaINFINBRADULT4, colour = CHAIN))

betaINFINBRADULT5 <- as.data.frame(betaINFINBRADULT5)
betaINFINBRADULT5$iteration <- seq(1:nrow(betaINFINBRADULT5))
BIIA5 <- ggplot(betaINFINBRADULT5) +
  geom_line(aes(y = betaINFINBRADULT5, x = iteration, colour = CHAIN))
BII5Ad <- ggplot(betaINFINBRADULT5) +
  geom_density(aes(x = betaINFINBRADULT5, colour = CHAIN))

betaINFINBRCUB1 <- as.data.frame(betaINFINBRCUB1)
BII1C <- ggplot(betaINFINBRCUB1) +
  geom_line(aes(y = betaINFINBRCUB1, x = iteration, colour = CHAIN))
BII1Cd <- ggplot(betaINFINBRCUB1) +
  geom_density(aes(x = betaINFINBRCUB1, colour = CHAIN))

betaINFINBRCUB2 <- as.data.frame(betaINFINBRCUB2)
betaINFINBRCUB2$iteration <- seq(1:nrow(betaINFINBRCUB2))
BIIC2 <- ggplot(betaINFINBRCUB2) +
  geom_line(aes(y = betaINFINBRCUB2, x = iteration, colour = CHAIN))
BII2Cd <- ggplot(betaINFINBRCUB2) +
  geom_density(aes(x = betaINFINBRCUB2, colour = CHAIN))

betaINFINBRCUB3 <- as.data.frame(betaINFINBRCUB3)
betaINFINBRCUB3$iteration <- seq(1:nrow(betaINFINBRCUB3))
BIIC3 <- ggplot(betaINFINBRCUB3) +
  geom_line(aes(y = betaINFINBRCUB3, x = iteration, colour = CHAIN))
BII3Cd <- ggplot(betaINFINBRCUB3) +
  geom_density(aes(x = betaINFINBRCUB3, colour = CHAIN))

betaINFINBRCUB4 <- as.data.frame(betaINFINBRCUB4)
betaINFINBRCUB4$iteration <- seq(1:nrow(betaINFINBRCUB4))
BIIC4 <- ggplot(betaINFINBRCUB4) +
  geom_line(aes(y = betaINFINBRCUB4, x = iteration, colour = CHAIN))
BII4Cd <- ggplot(betaINFINBRCUB4) +
  geom_density(aes(x = betaINFINBRCUB4, colour = CHAIN))

betaINFINBRCUB5 <- as.data.frame(betaINFINBRCUB5)
betaINFINBRCUB5$iteration <- seq(1:nrow(betaINFINBRCUB5))
BIIC5 <- ggplot(betaINFINBRCUB5) +
  geom_line(aes(y = betaINFINBRCUB5, x = iteration, colour = CHAIN))
BII5Cd <- ggplot(betaINFINBRCUB5) +
  geom_density(aes(x = betaINFINBRCUB5, colour = CHAIN))


#############################################################

## beta Sex:Inbr

betaSEXINBR1 <- samples %>%
  filter(zSEXINBR1 == 1) %>%
  select(CHAIN, betaSexINBR1) %>%
  mutate(iteration = 1:length(betaSexINBR1))

betaSEXINBR2 <- samples %>%
  filter(zSEXINBR2 == 1) %>%
  select(CHAIN, betaSexINBR2) %>%
  mutate(iteration = 1:length(betaSexINBR2))

betaSEXINBR3 <- samples %>%
  filter(zSEXINBR3 == 1) %>%
  select(CHAIN, betaSexINBR3) %>%
  mutate(iteration = 1:length(betaSexINBR3))

betaSEXINBR4 <- samples %>%
  filter(zSEXINBR4 == 1) %>%
  select(CHAIN, betaSexINBR4) %>%
  mutate(iteration = 1:length(betaSexINBR4))

betaSEXINBR5 <- samples %>%
  filter(zSEXINBR5 == 1) %>%
  select(CHAIN, betaSexINBR5) %>%
  mutate(iteration = 1:length(betaSexINBR5))

betaSEXINBR1 <- as.data.frame(betaSEXINBR1)
bSI1 <- ggplot(betaSEXINBR1) +
  geom_line(aes(y = betaSexINBR1, x = iteration, colour = CHAIN))
bSI1d <- ggplot(betaSEXINBR1) +
  geom_density(aes(x = betaSexINBR1, colour = CHAIN))

betaSEXINBR2 <- as.data.frame(betaSEXINBR2)
betaSEXINBR2$iteration <- seq(1:nrow(betaSEXINBR2))
bSI2 <- ggplot(betaSEXINBR2) +
  geom_line(aes(y = betaSexINBR2, x = iteration, colour = CHAIN))
bSI2d <- ggplot(betaSEXINBR2) +
  geom_density(aes(x = betaSexINBR2, colour = CHAIN))

betaSEXINBR3 <- as.data.frame(betaSEXINBR3)
betaSEXINBR3$iteration <- seq(1:nrow(betaSEXINBR3))
bSI3 <- ggplot(betaSEXINBR3) +
  geom_line(aes(y = betaSexINBR3, x = iteration, colour = CHAIN))
bSI3d <- ggplot(betaSEXINBR3) +
  geom_density(aes(x = betaSexINBR3, colour = CHAIN))

betaSEXINBR4 <- as.data.frame(betaSEXINBR4)
betaSEXINBR4$iteration <- seq(1:nrow(betaSEXINBR4))
bSI4 <- ggplot(betaSEXINBR4) +
  geom_line(aes(y = betaSexINBR4, x = iteration, colour = CHAIN))
bSI4d <- ggplot(betaSEXINBR4) +
  geom_density(aes(x = betaSexINBR4, colour = CHAIN))

betaSEXINBR5 <- as.data.frame(betaSEXINBR5)
betaSEXINBR5$iteration <- seq(1:nrow(betaSEXINBR5))
bSI5 <- ggplot(betaSEXINBR5) +
  geom_line(aes(y = betaSexINBR5, x = iteration, colour = CHAIN))
bSI5d <- ggplot(betaSEXINBR5) +
  geom_density(aes(x = betaSexINBR5, colour = CHAIN))


##########################################################################

## beta INF Adult

betaINFADULT1 <- samples %>%
  select(CHAIN, betaINFADULT1) %>%
  mutate(iteration = 1:length(betaINFADULT1))

betaINFADULT2 <- samples %>%
  select(CHAIN, betaINFADULT2) %>%
  mutate(iteration = 1:length(betaINFADULT2))

betaINFADULT3 <- samples %>%
  select(CHAIN, betaINFADULT3) %>%
  mutate(iteration = 1:length(betaINFADULT3))

betaINFADULT4 <- samples %>%
  select(CHAIN, betaINFADULT4) %>%
  mutate(iteration = 1:length(betaINFADULT4))

betaINFADULT5 <- samples %>%
  select(CHAIN, betaINFADULT5) %>%
  mutate(iteration = 1:length(betaINFADULT5))

## beta INF Cub

betaINFCUB1 <- samples %>%
  select(CHAIN, betaINFCUB1) %>%
  mutate(iteration = 1:length(betaINFCUB1))

betaINFCUB2 <- samples %>%
  select(CHAIN, betaINFCUB2) %>%
  mutate(iteration = 1:length(betaINFCUB2))

betaINFCUB3 <- samples %>%
  select(CHAIN, betaINFCUB3) %>%
  mutate(iteration = 1:length(betaINFCUB3))

betaINFCUB4 <- samples %>%
  select(CHAIN, betaINFCUB4) %>%
  mutate(iteration = 1:length(betaINFCUB4))

betaINFCUB5 <- samples %>%
  select(CHAIN, betaINFCUB5) %>%
  mutate(iteration = 1:length(betaINFCUB5))

##

betaSex1 <- samples %>%
  select(CHAIN, betaSex1) %>%
  mutate(iteration = 1:length(betaSex1))

betaSex2 <- samples %>%
  select(CHAIN, betaSex2) %>%
  mutate(iteration = 1:length(betaSex2))

betaSex3 <- samples %>%
  select(CHAIN, betaSex3) %>%
  mutate(iteration = 1:length(betaSex3))

betaSex4 <- samples %>%
  select(CHAIN, betaSex4) %>%
  mutate(iteration = 1:length(betaSex4))

betaSex5 <- samples %>%
  select(CHAIN, betaSex5) %>%
  mutate(iteration = 1:length(betaSex5))


betaINFADULT1 <- as.data.frame(betaINFADULT1)
betaINFADULT1$iteration <- seq(1:nrow(betaINFADULT1))
biAA1 <- ggplot(betaINFADULT1) +
  geom_line(aes(x = iteration, y = betaINFADULT1, colour = CHAIN))
biAA1d <- ggplot(betaINFADULT1) + 
  geom_density(aes(x = betaINFADULT1, colour = CHAIN))

betaINFADULT2 <- as.data.frame(betaINFADULT2)
betaINFADULT2$iteration <- seq(1:nrow(betaINFADULT2))
biAA2 <- ggplot(betaINFADULT2) +
  geom_line(aes(x = iteration, y = betaINFADULT2, colour = CHAIN))
biAA2d <- ggplot(betaINFADULT2) + 
  geom_density(aes(x = betaINFADULT2, colour = CHAIN))

betaINFADULT3 <- as.data.frame(betaINFADULT3)
betaINFADULT3$iteration <- seq(1:nrow(betaINFADULT3))
biAB1 <- ggplot(betaINFADULT3) +
  geom_line(aes(x = iteration, y = betaINFADULT3, colour = CHAIN))
biAB1d <- ggplot(betaINFADULT3) + 
  geom_density(aes(x = betaINFADULT3, colour = CHAIN))

betaINFADULT4 <- as.data.frame(betaINFADULT4)
betaINFADULT4$iteration <- seq(1:nrow(betaINFADULT4))
biAB2 <- ggplot(betaINFADULT4) +
  geom_line(aes(x = iteration, y = betaINFADULT4, colour = CHAIN))
biAB2d <- ggplot(betaINFADULT4) + 
  geom_density(aes(x = betaINFADULT4, colour = CHAIN))

betaINFADULT5 <- as.data.frame(betaINFADULT5)
betaINFADULT5$iteration <- seq(1:nrow(betaINFADULT5))
biAC1 <- ggplot(betaINFADULT5) +
  geom_line(aes(x = iteration, y = betaINFADULT5, colour = CHAIN))
biAC1d <- ggplot(betaINFADULT5) + 
  geom_density(aes(x = betaINFADULT5, colour = CHAIN))

##

betaINFCUB1 <- as.data.frame(betaINFCUB1)
betaINFCUB1$iteration <- seq(1:nrow(betaINFCUB1))
biCA1 <- ggplot(betaINFCUB1) +
  geom_line(aes(x = iteration, y = betaINFCUB1, colour = CHAIN))
biCA1d <- ggplot(betaINFCUB1) + 
  geom_density(aes(x = betaINFCUB1, colour = CHAIN))

betaINFCUB2 <- as.data.frame(betaINFCUB2)
betaINFCUB2$iteration <- seq(1:nrow(betaINFCUB2))
biCA2 <- ggplot(betaINFCUB2) +
  geom_line(aes(x = iteration, y = betaINFCUB2, colour = CHAIN))
biCA2d <- ggplot(betaINFCUB2) + 
  geom_density(aes(x = betaINFCUB2, colour = CHAIN))

betaINFCUB3 <- as.data.frame(betaINFCUB3)
betaINFCUB3$iteration <- seq(1:nrow(betaINFCUB3))
biCB1 <- ggplot(betaINFCUB3) +
  geom_line(aes(x = iteration, y = betaINFCUB3, colour = CHAIN))
biCB1d <- ggplot(betaINFCUB3) + 
  geom_density(aes(x = betaINFCUB3, colour = CHAIN))

betaINFCUB4 <- as.data.frame(betaINFCUB4)
betaINFCUB4$iteration <- seq(1:nrow(betaINFCUB4))
biCB2 <- ggplot(betaINFCUB4) +
  geom_line(aes(x = iteration, y = betaINFCUB4, colour = CHAIN))
biCB2d <- ggplot(betaINFCUB4) + 
  geom_density(aes(x = betaINFCUB4, colour = CHAIN))

betaINFCUB5 <- as.data.frame(betaINFCUB5)
betaINFCUB5$iteration <- seq(1:nrow(betaINFCUB5))
biCC1 <- ggplot(betaINFCUB5) +
  geom_line(aes(x = iteration, y = betaINFCUB5, colour = CHAIN))
biCC1d <- ggplot(betaINFCUB5) + 
  geom_density(aes(x = betaINFCUB5, colour = CHAIN))

##

betaSEX1 <- as.data.frame(betaSex1)
betaSEX1$iteration <- seq(1:nrow(betaSEX1))
bsA1 <- ggplot(betaSEX1) +
  geom_line(aes(x = iteration, y = betaSex1, colour = CHAIN))
bsA1d <- ggplot(betaSEX1) + 
  geom_density(aes(x = betaSex1, colour = CHAIN))

betaSEX2 <- as.data.frame(betaSex2)
betaSEX2$iteration <- seq(1:nrow(betaSEX2))
bsA2 <- ggplot(betaSEX2) +
  geom_line(aes(x = iteration, y = betaSex2, colour = CHAIN))
bsA2d <- ggplot(betaSEX2) + 
  geom_density(aes(x = betaSex2, colour = CHAIN))

betaSEX3 <- as.data.frame(betaSex3)
betaSEX3$iteration <- seq(1:nrow(betaSEX3))
bsB1 <- ggplot(betaSEX3) +
  geom_line(aes(x = iteration, y = betaSex3, colour = CHAIN))
bsB1d <- ggplot(betaSEX3) + 
  geom_density(aes(x = betaSex3, colour = CHAIN))

betaSEX4 <- as.data.frame(betaSex4)
betaSEX4$iteration <- seq(1:nrow(betaSEX4))
bsB2 <- ggplot(betaSEX4) +
  geom_line(aes(x = iteration, y = betaSex4, colour = CHAIN))
bsB2d <- ggplot(betaSEX4) + 
  geom_density(aes(x = betaSex4, colour = CHAIN))

betaSEX5 <- as.data.frame(betaSex5)
betaSEX5$iteration <- seq(1:nrow(betaSEX5))
bsC1 <- ggplot(betaSEX5) +
  geom_line(aes(x = iteration, y = betaSex5, colour = CHAIN))
bsC1d <- ggplot(betaSEX5) + 
  geom_density(aes(x = betaSex5, colour = CHAIN))

## Plot with patchwork

BIIA1/BIIA2/BIIA3|BIIA1d/BII2Ad/BII3Ad
BIIA4/BIIA5/BII1C|BII4Ad/BII5Ad/BII1Cd
BIIC2/BIIC3/BIIC4|BII1Cd/BII3Cd/BII4Cd
BIIC5/bSI1/bSI2|BII5Cd/bSI1d/bSI2d
bSI3/bSI4/bSI5|bSI3d/bSI4d/bSI5d
biAA1/biAA2/biAB1|biAA1d/biAA2d/biAB1d
biAB2/biAC1/biCA1|biAB2d/biAC1d/biCA1d
biCA2/biCB1/biCB2|biCA2d/biCB1d/biCB2d
biCC1/bsA1/bsA2|biCC1d/bsA1d/bsA2d
bsB1/bsB2/bsC1|bsB1d/bsB2d/bsC1d

##

