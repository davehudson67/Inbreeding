library(tidyverse)
library(patchwork)
library(MCMCvis)

samples <- readRDS("outputs/Sex3Infection_AllParameters_runsamples.rds")

mcmcplot(samples$samples, parms = c("a1", "a2", "b1", "b2", "c1"))

samples <- as.matrix(samples$samples, chains = TRUE)
samples <- as.data.frame(samples)

names(samples)
colnames(samples) <- c("CHAIN", "a1", "a2", "b1", "b2", "betaINFADULT1", "betaINFADULT2", "betaINFADULT3", "betaINFADULT4", "betaINFADULT5", "betaINFCUB1", "betaINFCUB2", "betaINFCUB3", 
              "betaINFCUB4", "betaINFCUB5", "betaSEX1", "betaSEX2", "betaSEX3", "betaSEX4", "betaSEX5", "betaSEXINFADULT1", "betaSEXINFADULT2", "betaSEXINFADULT3", "betaSEXINFADULT4", 
              "betaSEXINFADULT5", "betaSEXINFCUB1", "betaSEXINFCUB2", "betaSEXINFCUB3", "betaSEXINFCUB4", "betaSEXINFCUB5", "c1", "mean.p", "zINF1", "zINF2", "zINF3", "zINF4", "zINF5",
              "zSEX1", "zSEX2", "zSEX3", "zSEX4", "zSEX5", "zSEXINF1", "zSEXINF2", "zSEXINF3", "zSEXINF4", "zSEXINF5") 

samples$CHAIN <- as.factor(samples$CHAIN)

##

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

##

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

betaSEXINFCUB1 <- samples %>%
  filter(zSEXINF1 == 1) %>%
  select(CHAIN, betaSEXINFCUB1) %>%
  mutate(iteration = 1:length(betaSEXINFCUB1))

betaSEXINFCUB2 <- samples %>%
  filter(zSEXINF2 == 1) %>%
  select(CHAIN, betaSEXINFCUB2) %>%
  mutate(iteration = 1:length(betaSEXINFCUB2))

betaSEXINFCUB3 <- samples %>%
  filter(zSEXINF3 == 1) %>%
  select(CHAIN, betaSEXINFCUB3) %>%
  mutate(iteration = 1:length(betaSEXINFCUB3))

betaSEXINFCUB4 <- samples %>%
  filter(zSEXINF4 == 1) %>%
  select(CHAIN, betaSEXINFCUB4) %>%
  mutate(iteration = 1:length(betaSEXINFCUB4))

betaSEXINFCUB5 <- samples %>%
  filter(zSEXINF5 == 1) %>%
  select(CHAIN, betaSEXINFCUB5) %>%
  mutate(iteration = 1:length(betaSEXINFCUB5))

##

betaSEXINFADULT1 <- samples %>%
  filter(zSEXINF1 == 1) %>%
  select(CHAIN, betaSEXINFADULT1) %>%
  mutate(iteration = 1:length(betaSEXINFADULT1))

betaSEXINFADULT2 <- samples %>%
  filter(zSEXINF2 == 1) %>%
  select(CHAIN, betaSEXINFADULT2) %>%
  mutate(iteration = 1:length(betaSEXINFADULT2))

betaSEXINFADULT3 <- samples %>%
  filter(zSEXINF3 == 1) %>%
  select(CHAIN, betaSEXINFADULT3) %>%
  mutate(iteration = 1:length(betaSEXINFADULT3))

betaSEXINFADULT4 <- samples %>%
  filter(zSEXINF4 == 1) %>%
  select(CHAIN, betaSEXINFADULT4) %>%
  mutate(iteration = 1:length(betaSEXINFADULT4))

betaSEXINFADULT5 <- samples %>%
  filter(zSEXINF5 == 1) %>%
  select(CHAIN, betaSEXINFADULT5) %>%
  mutate(iteration = 1:length(betaSEXINFADULT5))

##

betaSEXINFADULT1 <- as.data.frame(betaSEXINFADULT1)
betaSEXINFADULT1$iteration <- seq(1:nrow(betaSEXINFADULT1))
bsiAA1 <- ggplot(betaSEXINFADULT1) +
  geom_line(aes(x = iteration, y = betaSEXINFADULT1, colour = CHAIN))
bsiAA1d <- ggplot(betaSEXINFADULT1) + 
  geom_density(aes(x = betaSEXINFADULT1, colour = CHAIN))

betaSEXINFADULT2 <- as.data.frame(betaSEXINFADULT2)
betaSEXINFADULT2$iteration <- seq(1:nrow(betaSEXINFADULT2))
bsiAA2 <- ggplot(betaSEXINFADULT2) +
  geom_line(aes(x = iteration, y = betaSEXINFADULT2, colour = CHAIN))
bsiAA2d <- ggplot(betaSEXINFADULT2) + 
  geom_density(aes(x = betaSEXINFADULT2, colour = CHAIN))

betaSEXINFADULT3 <- as.data.frame(betaSEXINFADULT3)
betaSEXINFADULT3$iteration <- seq(1:nrow(betaSEXINFADULT3))
bsiAB1 <- ggplot(betaSEXINFADULT3) +
  geom_line(aes(x = iteration, y = betaSEXINFADULT3, colour = CHAIN))
bsiAB1d <- ggplot(betaSEXINFADULT3) + 
  geom_density(aes(x = betaSEXINFADULT3, colour = CHAIN))

betaSEXINFADULT4 <- as.data.frame(betaSEXINFADULT4)
betaSEXINFADULT4$iteration <- seq(1:nrow(betaSEXINFADULT4))
bsiAB2 <- ggplot(betaSEXINFADULT4) +
  geom_line(aes(x = iteration, y = betaSEXINFADULT4, colour = CHAIN))
bsiAB2d <- ggplot(betaSEXINFADULT4) + 
  geom_density(aes(x = betaSEXINFADULT4, colour = CHAIN))

betaSEXINFADULT5 <- as.data.frame(betaSEXINFADULT5)
betaSEXINFADULT5$iteration <- seq(1:nrow(betaSEXINFADULT5))
bsiAC1 <- ggplot(betaSEXINFADULT5) +
  geom_line(aes(x = iteration, y = betaSEXINFADULT5, colour = CHAIN))
bsiAC1d <- ggplot(betaSEXINFADULT5) + 
  geom_density(aes(x = betaSEXINFADULT5, colour = CHAIN))

##

betaSEXINFCUB1 <- as.data.frame(betaSEXINFCUB1)
betaSEXINFCUB1$iteration <- seq(1:nrow(betaSEXINFCUB1))
bsiCA1 <- ggplot(betaSEXINFCUB1) +
  geom_line(aes(x = iteration, y = betaSEXINFCUB1, colour = CHAIN))
bsiCA1d <- ggplot(betaSEXINFCUB1) + 
  geom_density(aes(x = betaSEXINFCUB1, colour = CHAIN))

betaSEXINFCUB2 <- as.data.frame(betaSEXINFCUB2)
betaSEXINFCUB2$iteration <- seq(1:nrow(betaSEXINFCUB2))
bsiCA2 <- ggplot(betaSEXINFCUB2) +
  geom_line(aes(x = iteration, y = betaSEXINFCUB2, colour = CHAIN))
bsiCA2d <- ggplot(betaSEXINFCUB2) + 
  geom_density(aes(x = betaSEXINFCUB2, colour = CHAIN))

betaSEXINFCUB3 <- as.data.frame(betaSEXINFCUB3)
betaSEXINFCUB3$iteration <- seq(1:nrow(betaSEXINFCUB3))
bsiCB1 <- ggplot(betaSEXINFCUB3) +
  geom_line(aes(x = iteration, y = betaSEXINFCUB3, colour = CHAIN))
bsiCB1d <- ggplot(betaSEXINFCUB3) + 
  geom_density(aes(x = betaSEXINFCUB3, colour = CHAIN))

betaSEXINFCUB4 <- as.data.frame(betaSEXINFCUB4)
betaSEXINFCUB4$iteration <- seq(1:nrow(betaSEXINFCUB4))
bsiCB2 <- ggplot(betaSEXINFCUB4) +
  geom_line(aes(x = iteration, y = betaSEXINFCUB4, colour = CHAIN))
bsiCB2d <- ggplot(betaSEXINFCUB4) + 
  geom_density(aes(x = betaSEXINFCUB4, colour = CHAIN))

betaSEXINFCUB5 <- as.data.frame(betaSEXINFCUB5)
betaSEXINFCUB5$iteration <- seq(1:nrow(betaSEXINFCUB5))
bsiCC1 <- ggplot(betaSEXINFCUB5) +
  geom_line(aes(x = iteration, y = betaSEXINFCUB5, colour = CHAIN))
bsiCC1d <- ggplot(betaSEXINFCUB5) + 
  geom_density(aes(x = betaSEXINFCUB5, colour = CHAIN))

##

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

## Plot with patchwork

bsiAA1/bsiAA2/bsiAB1|bsiAA1d/bsiAA2d/bsiAB1d
bsiAB2/bsiAC1/bsiCA1|bsiAB2d/bsiAC1d/bsiCA1d
bsiCA2/bsiCB1/bsiCB2|bsiCA2d/bsiCB1d/bsiCB2d
bsiCC1/biAA1/biAA2|bsiCC1d/biAA1d/biAA2d
biAB1/biAB2/biAC1|biAB1d/biAB2d/biAC1d
biCA1/biCA2/biCB1|biCA1d/biCA2d/biCB1d
biCB2/biCC1|biCB2d/biCC1d

