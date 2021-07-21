library(tidyverse)
library(patchwork)
library(MCMCvis)
library(mcmcplots)

#samples <- readRDS("outputs/FullModel_NoInteractions_z_all_runsamples.rds")
samples <- readRDS("outputs/FullModel_InbrCont_NoInteractions_z_all_runsamples.rds")

mcmcplot(samples$samples, parms = c("a1", "a2", "b1", "b2", "c1"))

#MCMCtrace(samples$samples, params = c("a1", "a2", "b1", "b2", "c1"), pdf = F)

samples <- as.matrix(samples$samples, chains = TRUE)
samples <- as.data.frame(samples)
#samples <- sample_n(samples, 50000)

names(samples)
colnames(samples) <- c("CHAIN", "a1", "a2", "b1", "b2", 
                       "betaINBR1", "betaINBR2", "betaINBR3", "betaINBR4", "betaINBR5",  
                       "betaINFADULT1", "betaINFADULT2", "betaINFADULT3", "betaINFADULT4", "betaINFADULT5", 
                       "betaINFCUB1", "betaINFCUB2", "betaINFCUB3", "betaINFCUB4", "betaINFCUB5", 
                       #"betaINFINBRADULT1", "betaINFINBRADULT2", "betaINFINBRADULT3", "betaINFINBRADULT4", "betaINFINBRADULT5", 
                       #"betaINFINBRCUB1", "betaINFINBRCUB2", "betaINFINBRCUB3", "betaINFINBRCUB4", "betaINFINBRCUB5",  
                       "betaSex1", "betaSex2", "betaSex3", "betaSex4", "betaSex5",
                       #"betaSexINBR1", "betaSexINBR2", "betaSexINBR3", "betaSexINBR4", "betaSexINBR5", 
                       #"betaSEXINFADULT1", "betaSEXINFADULT2", "betaSEXINFADULT3", "betaSEXINFADULT4", "betaSEXINFADULT5", 
                       #"betaSEXINFCUB1", "betaSEXINFCUB2", "betaSEXINFCUB3", "betaSEXINFCUB4", "betaSEXINFCUB5", 
                       "c1", "mean.p", "zINBR1", "zINBR2", "zINBR3", "zINBR4", "zINBR5",
                       "zINF1", "zINF2", "zINF3", "zINF4", "zINF5",
                       #"zINFINBR1", "zINFINBR2", "zINFINBR3", "zINFINBR4", "zINFINBR5",
                       "zSEX1", "zSEX2", "zSEX3", "zSEX4", "zSEX5")
                       #"zSEXINBR1", "zSEXINBR2", "zSEXINBR3", "zSEXINBR4", "zSEXINBR5",
                       #"zSEXINF1", "zSEXINF2", "zSEXINF3", "zSEXINF4", "zSEXINF5") 

samples$CHAIN <- as.factor(samples$CHAIN)
pdf("Posteriors_FullModel_NoInteractions_InbrCONT.pdf")


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
#b1/b2/b3|b1d/b2d/b3d
#b4/b5|b4d/b5d

##########################################################################

## beta INF Adult

betaINFADULT1 <- samples %>%
  filter(zINF1 == 1) %>%
  select(CHAIN, betaINFADULT1) %>%
  mutate(iteration = 1:length(betaINFADULT1))

betaINFADULT2 <- samples %>%
  filter(zINF2 == 1) %>%
  select(CHAIN, betaINFADULT2) %>%
  mutate(iteration = 1:length(betaINFADULT2))

betaINFADULT3 <- samples %>%
  filter(zINF3 == 1) %>%
  select(CHAIN, betaINFADULT3) %>%
  mutate(iteration = 1:length(betaINFADULT3))

betaINFADULT4 <- samples %>%
  filter(zINF4 == 1) %>%
  select(CHAIN, betaINFADULT4) %>%
  mutate(iteration = 1:length(betaINFADULT4))

betaINFADULT5 <- samples %>%
  filter(zINF5 == 1) %>%
  select(CHAIN, betaINFADULT5) %>%
  mutate(iteration = 1:length(betaINFADULT5))

## beta INF Cub

betaINFCUB1 <- samples %>%
  filter(zINF1 == 1) %>%
  select(CHAIN, betaINFCUB1) %>%
  mutate(iteration = 1:length(betaINFCUB1))

betaINFCUB2 <- samples %>%
  filter(zINF2 == 1) %>%
  select(CHAIN, betaINFCUB2) %>%
  mutate(iteration = 1:length(betaINFCUB2))

betaINFCUB3 <- samples %>%
  filter(zINF3 == 1) %>%
  select(CHAIN, betaINFCUB3) %>%
  mutate(iteration = 1:length(betaINFCUB3))

betaINFCUB4 <- samples %>%
  filter(zINF4 == 1) %>%
  select(CHAIN, betaINFCUB4) %>%
  mutate(iteration = 1:length(betaINFCUB4))

betaINFCUB5 <- samples %>%
  filter(zINF5 == 1) %>%
  select(CHAIN, betaINFCUB5) %>%
  mutate(iteration = 1:length(betaINFCUB5))

##

betaSex1 <- samples %>%
  filter(zSEX1 == 1) %>%
  select(CHAIN, betaSex1) %>%
  mutate(iteration = 1:length(betaSex1))

betaSex2 <- samples %>%
  filter(zSEX2 == 1) %>%
  select(CHAIN, betaSex2) %>%
  mutate(iteration = 1:length(betaSex2))

betaSex3 <- samples %>%
  filter(zSEX3 == 1) %>%
  select(CHAIN, betaSex3) %>%
  mutate(iteration = 1:length(betaSex3))

betaSex4 <- samples %>%
  filter(zSEX4 == 1) %>%
  select(CHAIN, betaSex4) %>%
  mutate(iteration = 1:length(betaSex4))

betaSex5 <- samples %>%
  filter(zSEX5 == 1) %>%
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

biAA1/biAA2/biAB1|biAA1d/biAA2d/biAB1d
biAB2/biAC1/biCA1|biAB2d/biAC1d/biCA1d
biCA2/biCB1/biCB2|biCA2d/biCB1d/biCB2d
biCC1/bsA1/bsA2|biCC1d/bsA1d/bsA2d
bsB1/bsB2/bsC1|bsB1d/bsB2d/bsC1d

dev.off()
##

## density only

b1d/b3d/b5d|b2d/b4d/biAA1d

biAA2d/biAB2d/biCA1d|biAB1d/biAC1d/biCA2d

biCB1d/biCC1d/bsA2d|biCB2d/bsA1d/bsB1d

bsB2d/a1/b1|bsC1d/a2/b2

c1

## main params

a1 <- samples %>%
  select(a1, CHAIN) %>%
  ggplot() +
  geom_density(aes(x = a1, colour = CHAIN))

a2 <- samples %>%
  select(a2, CHAIN) %>%
  ggplot() +
  geom_density(aes(x = a2, colour = CHAIN))

b1 <- samples %>%
  select(b1, CHAIN) %>%
  ggplot() +
  geom_density(aes(x = b1, colour = CHAIN))

b2 <- samples %>%
  select(b2, CHAIN) %>%
  ggplot() +
  geom_density(aes(x = b2, colour = CHAIN))

c1 <- samples %>%
  select(c1, CHAIN) %>%
  ggplot() +
  geom_density(aes(x = c1, colour = CHAIN))
