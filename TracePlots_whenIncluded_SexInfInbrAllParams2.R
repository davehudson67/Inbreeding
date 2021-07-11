library(tidyverse)
library(patchwork)
library(MCMCvis)

samples <- readRDS("outputs/Sex3InfectionInbr_AllParameters_runsamples.rds")

mcmcplot(samples$samples, parms = c("a1", "a2", "b1", "b2", "c1"))

samples <- as.matrix(samples$samples, chains = TRUE)
samples <- as.data.frame(samples)

names(samples)
colnames(samples) <- c("CHAIN", "a1", "a2", "b1", "b2", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "beta8", "beta9", "beta10",
              "betaINBR1", "betaINBR2", "betaINBR3", "betaINBR4", "betaINBR5", "betaINBR6", "betaINBR7", "betaINBR8", "betaINBR9", "betaINBR10",
              "betaINBR11", "betaINBR12", "betaINBR13", "betaINBR14", "betaINBR15",  "betaINFADULT1", "betaINFADULT2", "betaINFADULT3", "betaINFADULT4", 
              "betaINFADULT5", "betaINFCUB1", "betaINFCUB2", "betaINFCUB3", "betaINFCUB4", "betaINFCUB5", "betaInbrInfADULT1", "betaInbrInfADULT2", 
              "betaInbrInfADULT3", "betaInbrInfADULT4", "betaInbrInfADULT5", "betaInbrInfCUB1", "betaInbrInfCUB2", "betaInbrInfCUB3", "betaInbrInfCUB4",
              "betaInbrInfCUB5", "c1", "mean.p", "z1", "z2", "z3", "z4", "z5", "z6", "z7", "z8", "z9", "z10", "z11",
              "z12", "z13", "z14", "z15") 

samples$CHAIN <- as.factor(samples$CHAIN)

beta1 <- samples %>%
  select(CHAIN, beta1) %>%
  mutate(iteration = 1:length(beta1))

beta2 <- samples %>%
  select(CHAIN, beta2) %>%
  mutate(iteration = 1:length(beta2))

beta3 <- samples %>%
  select(CHAIN, beta3) %>%
  mutate(iteration = 1:length(beta3))

beta4 <- samples %>%
  select(CHAIN, beta4) %>%
  mutate(iteration = 1:length(beta4))

beta5 <- samples %>%
  select(CHAIN, beta5) %>%
  mutate(iteration = 1:length(beta5))

beta6 <- samples %>%
  select(CHAIN, beta6) %>%
  mutate(iteration = 1:length(beta6))

beta7 <- samples %>%
  select(CHAIN, beta7) %>%
  mutate(iteration = 1:length(beta7))

beta8 <- samples %>%
  select(CHAIN, beta8) %>%
  mutate(iteration = 1:length(beta8))

beta9 <- samples %>%
  select(CHAIN, beta9) %>%
  mutate(iteration = 1:length(beta9))

beta10 <- samples %>%
  select(CHAIN, beta10) %>%
  mutate(iteration = 1:length(beta10))

##

betaINBR1 <- samples %>%
  filter(z1 == 1) %>%
  select(CHAIN, betaINBR1) %>%
  mutate(iteration = 1:length(betaINBR1))

betaINBR2 <- samples %>%
  filter(z2 == 1) %>%
  select(CHAIN, betaINBR2) %>%
  mutate(iteration = 1:length(betaINBR2))

betaINBR3 <- samples %>%
  filter(z3 == 1) %>%
  select(CHAIN, betaINBR3) %>%
  mutate(iteration = 1:length(betaINBR3))

betaINBR4 <- samples %>%
  filter(z4 == 1) %>%
  select(CHAIN, betaINBR4) %>%
  mutate(iteration = 1:length(betaINBR4))

betaINBR5 <- samples %>%
  filter(z5 == 1) %>%
  select(CHAIN, betaINBR5) %>%
  mutate(iteration = 1:length(betaINBR5))

betaINBR6 <- samples %>%
  filter(z6 == 1) %>%
  select(CHAIN, betaINBR6) %>%
  mutate(iteration = 1:length(betaINBR6))

betaINBR7 <- samples %>%
  filter(z7 == 1) %>%
  select(CHAIN, betaINBR7) %>%
  mutate(iteration = 1:length(betaINBR7))

betaINBR8 <- samples %>%
  filter(z8 == 1) %>%
  select(CHAIN, betaINBR8) %>%
  mutate(iteration = 1:length(betaINBR8))

betaINBR9 <- samples %>%
  filter(z9 == 1) %>%
  select(CHAIN, betaINBR9) %>%
  mutate(iteration = 1:length(betaINBR9))

betaINBR10 <- samples %>%
  filter(z10 == 1) %>%
  select(CHAIN, betaINBR10) %>%
  mutate(iteration = 1:length(betaINBR10))

betaINBR11 <- samples %>%
  filter(z11 == 1) %>%
  select(CHAIN, betaINBR11) %>%
  mutate(iteration = 1:length(betaINBR11))

betaINBR12 <- samples %>%
  filter(z12 == 1) %>%
  select(CHAIN, betaINBR12) %>%
  mutate(iteration = 1:length(betaINBR12))

betaINBR13 <- samples %>%
  filter(z13 == 1) %>%
  select(CHAIN, betaINBR13) %>%
  mutate(iteration = 1:length(betaINBR13))

betaINBR14 <- samples %>%
  filter(z14 == 1) %>%
  select(CHAIN, betaINBR14) %>%
  mutate(iteration = 1:length(betaINBR14))

betaINBR15 <- samples %>%
  filter(z15 == 1) %>%
  select(CHAIN, betaINBR15) %>%
  mutate(iteration = 1:length(betaINBR15))

##

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

##

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

betaInbrInfCUB1 <- samples %>%
  filter(z3 == 1) %>%
  select(CHAIN, betaInbrInfCUB1) %>%
  mutate(iteration = 1:length(betaInbrInfCUB1))

betaInbrInfCUB2 <- samples %>%
  filter(z6 == 1) %>%
  select(CHAIN, betaInbrInfCUB2) %>%
  mutate(iteration = 1:length(betaInbrInfCUB2))

betaInbrInfCUB3 <- samples %>%
  filter(z9 == 1) %>%
  select(CHAIN, betaInbrInfCUB3) %>%
  mutate(iteration = 1:length(betaInbrInfCUB3))

betaInbrInfCUB4 <- samples %>%
  filter(z12 == 1) %>%
  select(CHAIN, betaInbrInfCUB4) %>%
  mutate(iteration = 1:length(betaInbrInfCUB4))

betaInbrInfCUB5 <- samples %>%
  filter(z15 == 1) %>%
  select(CHAIN, betaInbrInfCUB5) %>%
  mutate(iteration = 1:length(betaInbrInfCUB5))

##

betaInbrInfADULT1 <- samples %>%
  filter(z3 == 1) %>%
  select(CHAIN, betaInbrInfADULT1) %>%
  mutate(iteration = 1:length(betaInbrInfADULT1))

betaInbrInfADULT2 <- samples %>%
  filter(z6 == 1) %>%
  select(CHAIN, betaInbrInfADULT2) %>%
  mutate(iteration = 1:length(betaInbrInfADULT2))

betaInbrInfADULT3 <- samples %>%
  filter(z9 == 1) %>%
  select(CHAIN, betaInbrInfADULT3) %>%
  mutate(iteration = 1:length(betaInbrInfADULT3))

betaInbrInfADULT4 <- samples %>%
  filter(z12 == 1) %>%
  select(CHAIN, betaInbrInfADULT4) %>%
  mutate(iteration = 1:length(betaInbrInfADULT4))

betaInbrInfADULT5 <- samples %>%
  filter(z15 == 1) %>%
  select(CHAIN, betaInbrInfADULT5) %>%
  mutate(iteration = 1:length(betaInbrInfADULT5))

##

beta1 <- as.data.frame(beta1)
b1 <- ggplot(beta1) +
  geom_line(aes(y = beta1, x = iteration, colour = CHAIN))
b1d <- ggplot(beta1) +
  geom_density(aes(x = beta1, colour = CHAIN))

beta2 <- as.data.frame(beta2)
beta2$iteration <- seq(1:nrow(beta2))
b2 <- ggplot(beta2) +
  geom_line(aes(y = beta2, x = iteration, colour = CHAIN))
b2d <- ggplot(beta2) +
  geom_density(aes(x = beta2, colour = CHAIN))

beta3 <- as.data.frame(beta3)
beta3$iteration <- seq(1:nrow(beta3))
b3 <- ggplot(beta3) +
  geom_line(aes(y = beta3, x = iteration, colour = CHAIN))
b3d <- ggplot(beta3) +
  geom_density(aes(x = beta3, colour = CHAIN))

beta4 <- as.data.frame(beta4)
beta4$iteration <- seq(1:nrow(beta4))
b4 <- ggplot(beta4) +
  geom_line(aes(y = beta4, x = iteration, colour = CHAIN))
b4d <- ggplot(beta4) +
  geom_density(aes(x = beta4, colour = CHAIN))


beta5 <- as.data.frame(beta5)
beta5$iteration <- seq(1:nrow(beta5))
b5 <- ggplot(beta5) +
  geom_line(aes(y = beta5, x = iteration, colour = CHAIN))
b5d <- ggplot(beta5) +
  geom_density(aes(x = beta5, colour = CHAIN))

beta6 <- as.data.frame(beta6)
beta6$iteration <- seq(1:nrow(beta6))
b6 <- ggplot(beta6) +
  geom_line(aes(y = beta6, x = iteration, colour = CHAIN))
b6d <- ggplot(beta6) +
  geom_density(aes(x = beta6, colour = CHAIN))

beta7 <- as.data.frame(beta7)
beta7$iteration <- seq(1:nrow(beta7))
b7 <- ggplot(beta7) +
  geom_line(aes(y = beta7, x = iteration, colour = CHAIN))
b7d <- ggplot(beta7) +
  geom_density(aes(x = beta7, colour = CHAIN))

beta8 <- as.data.frame(beta8)
beta8$iteration <- seq(1:nrow(beta8))
b8 <- ggplot(beta8) +
  geom_line(aes(y = beta8, x = iteration, colour = CHAIN))
b8d <- ggplot(beta8) +
  geom_density(aes(x = beta8, colour = CHAIN))

beta9 <- as.data.frame(beta9)
beta9$iteration <- seq(1:nrow(beta9))
b9 <- ggplot(beta9) +
  geom_line(aes(y = beta9, x = iteration, colour = CHAIN))
b9d <- ggplot(beta9) +
  geom_density(aes(x = beta9, colour = CHAIN))

beta10 <- as.data.frame(beta10)
beta10$iteration <- seq(1:nrow(beta10))
b10 <- ggplot(beta10) +
  geom_line(aes(y = beta10, x = iteration, colour = CHAIN))
b10d <- ggplot(beta10) +
  geom_density(aes(x = beta10, colour = CHAIN))

b1/b2/b3|b1d/b2d/b3d
b4/b5/b6|b4d/b5d/b6d
b7/b8/b9|b7d/b8d/b9d
b10|b10d

##

betaInbrInfADULT1 <- as.data.frame(betaInbrInfADULT1)
betaInbrInfADULT1$iteration <- seq(1:nrow(betaInbrInfADULT1))
bsiAA1 <- ggplot(betaInbrInfADULT1) +
  geom_line(aes(x = iteration, y = betaInbrInfADULT1, colour = CHAIN))
bsiAA1d <- ggplot(betaInbrInfADULT1) + 
  geom_density(aes(x = betaInbrInfADULT1, colour = CHAIN))

betaInbrInfADULT2 <- as.data.frame(betaInbrInfADULT2)
betaInbrInfADULT2$iteration <- seq(1:nrow(betaInbrInfADULT2))
bsiAA2 <- ggplot(betaInbrInfADULT2) +
  geom_line(aes(x = iteration, y = betaInbrInfADULT2, colour = CHAIN))
bsiAA2d <- ggplot(betaInbrInfADULT2) + 
  geom_density(aes(x = betaInbrInfADULT2, colour = CHAIN))

betaInbrInfADULT3 <- as.data.frame(betaInbrInfADULT3)
betaInbrInfADULT3$iteration <- seq(1:nrow(betaInbrInfADULT3))
bsiAB1 <- ggplot(betaInbrInfADULT3) +
  geom_line(aes(x = iteration, y = betaInbrInfADULT3, colour = CHAIN))
bsiAB1d <- ggplot(betaInbrInfADULT3) + 
  geom_density(aes(x = betaInbrInfADULT3, colour = CHAIN))

betaInbrInfADULT4 <- as.data.frame(betaInbrInfADULT4)
betaInbrInfADULT4$iteration <- seq(1:nrow(betaInbrInfADULT4))
bsiAB2 <- ggplot(betaInbrInfADULT4) +
  geom_line(aes(x = iteration, y = betaInbrInfADULT4, colour = CHAIN))
bsiAB2d <- ggplot(betaInbrInfADULT4) + 
  geom_density(aes(x = betaInbrInfADULT4, colour = CHAIN))

betaInbrInfADULT5 <- as.data.frame(betaInbrInfADULT5)
betaInbrInfADULT5$iteration <- seq(1:nrow(betaInbrInfADULT5))
bsiAC1 <- ggplot(betaInbrInfADULT5) +
  geom_line(aes(x = iteration, y = betaInbrInfADULT5, colour = CHAIN))
bsiAC1d <- ggplot(betaInbrInfADULT5) + 
  geom_density(aes(x = betaInbrInfADULT5, colour = CHAIN))

##

betaInbrInfCUB1 <- as.data.frame(betaInbrInfCUB1)
betaInbrInfCUB1$iteration <- seq(1:nrow(betaInbrInfCUB1))
bsiCA1 <- ggplot(betaInbrInfCUB1) +
  geom_line(aes(x = iteration, y = betaInbrInfCUB1, colour = CHAIN))
bsiCA1d <- ggplot(betaInbrInfCUB1) + 
  geom_density(aes(x = betaInbrInfCUB1, colour = CHAIN))

betaInbrInfCUB2 <- as.data.frame(betaInbrInfCUB2)
betaInbrInfCUB2$iteration <- seq(1:nrow(betaInbrInfCUB2))
bsiCA2 <- ggplot(betaInbrInfCUB2) +
  geom_line(aes(x = iteration, y = betaInbrInfCUB2, colour = CHAIN))
bsiCA2d <- ggplot(betaInbrInfCUB2) + 
  geom_density(aes(x = betaInbrInfCUB2, colour = CHAIN))

betaInbrInfCUB3 <- as.data.frame(betaInbrInfCUB3)
betaInbrInfCUB3$iteration <- seq(1:nrow(betaInbrInfCUB3))
bsiCB1 <- ggplot(betaInbrInfCUB3) +
  geom_line(aes(x = iteration, y = betaInbrInfCUB3, colour = CHAIN))
bsiCB1d <- ggplot(betaInbrInfCUB3) + 
  geom_density(aes(x = betaInbrInfCUB3, colour = CHAIN))

betaInbrInfCUB4 <- as.data.frame(betaInbrInfCUB4)
betaInbrInfCUB4$iteration <- seq(1:nrow(betaInbrInfCUB4))
bsiCB2 <- ggplot(betaInbrInfCUB4) +
  geom_line(aes(x = iteration, y = betaInbrInfCUB4, colour = CHAIN))
bsiCB2d <- ggplot(betaInbrInfCUB4) + 
  geom_density(aes(x = betaInbrInfCUB4, colour = CHAIN))

betaInbrInfCUB5 <- as.data.frame(betaInbrInfCUB5)
betaInbrInfCUB5$iteration <- seq(1:nrow(betaInbrInfCUB5))
bsiCC1 <- ggplot(betaInbrInfCUB5) +
  geom_line(aes(x = iteration, y = betaInbrInfCUB5, colour = CHAIN))
bsiCC1d <- ggplot(betaInbrInfCUB5) + 
  geom_density(aes(x = betaInbrInfCUB5, colour = CHAIN))

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

##

betaINBR1 <- as.data.frame(betaINBR1)
bI1 <- ggplot(betaINBR1) +
  geom_line(aes(y = betaINBR1, x = iteration, colour = CHAIN))
bI1d <- ggplot(betaINBR1) +
  geom_density(aes(x = betaINBR1, colour = CHAIN))

betaINBR2 <- as.data.frame(betaINBR2)
betaINBR2$iteration <- seq(1:nrow(betaINBR2))
bI2 <- ggplot(betaINBR2) +
  geom_line(aes(y = betaINBR2, x = iteration, colour = CHAIN))
bI2d <- ggplot(betaINBR2) +
  geom_density(aes(x = betaINBR2, colour = CHAIN))

betaINBR3 <- as.data.frame(betaINBR3)
betaINBR3$iteration <- seq(1:nrow(betaINBR3))
bI3 <- ggplot(betaINBR3) +
  geom_line(aes(y = betaINBR3, x = iteration, colour = CHAIN))
bI3d <- ggplot(betaINBR3) +
  geom_density(aes(x = betaINBR3, colour = CHAIN))

betaINBR4 <- as.data.frame(betaINBR4)
betaINBR4$iteration <- seq(1:nrow(betaINBR4))
bI4 <- ggplot(betaINBR4) +
  geom_line(aes(y = betaINBR4, x = iteration, colour = CHAIN))
bI4d <- ggplot(betaINBR4) +
  geom_density(aes(x = betaINBR4, colour = CHAIN))


betaINBR5 <- as.data.frame(betaINBR5)
betaINBR5$iteration <- seq(1:nrow(betaINBR5))
bI5 <- ggplot(betaINBR5) +
  geom_line(aes(y = betaINBR5, x = iteration, colour = CHAIN))
bI5d <- ggplot(betaINBR5) +
  geom_density(aes(x = betaINBR5, colour = CHAIN))

betaINBR6 <- as.data.frame(betaINBR6)
betaINBR6$iteration <- seq(1:nrow(betaINBR6))
bI6 <- ggplot(betaINBR6) +
  geom_line(aes(y = betaINBR6, x = iteration, colour = CHAIN))
bI6d <- ggplot(betaINBR6) +
  geom_density(aes(x = betaINBR6, colour = CHAIN))

betaINBR7 <- as.data.frame(betaINBR7)
betaINBR7$iteration <- seq(1:nrow(betaINBR7))
bI7 <- ggplot(betaINBR7) +
  geom_line(aes(y = betaINBR7, x = iteration, colour = CHAIN))
bI7d <- ggplot(betaINBR7) +
  geom_density(aes(x = betaINBR7, colour = CHAIN))

betaINBR8 <- as.data.frame(betaINBR8)
betaINBR8$iteration <- seq(1:nrow(betaINBR8))
bI8 <- ggplot(betaINBR8) +
  geom_line(aes(y = betaINBR8, x = iteration, colour = CHAIN))
bI8d <- ggplot(betaINBR8) +
  geom_density(aes(x = betaINBR8, colour = CHAIN))

betaINBR9 <- as.data.frame(betaINBR9)
betaINBR9$iteration <- seq(1:nrow(betaINBR9))
bI9 <- ggplot(betaINBR9) +
  geom_line(aes(y = betaINBR9, x = iteration, colour = CHAIN))
bI9d <- ggplot(betaINBR9) +
  geom_density(aes(x = betaINBR9, colour = CHAIN))

betaINBR10 <- as.data.frame(betaINBR10)
betaINBR10$iteration <- seq(1:nrow(betaINBR10))
bI10 <- ggplot(betaINBR10) +
  geom_line(aes(y = betaINBR10, x = iteration, colour = CHAIN))
bI10d <- ggplot(betaINBR10) +
  geom_density(aes(x = betaINBR10, colour = CHAIN))

betaINBR11 <- as.data.frame(betaINBR11)
bI11 <- ggplot(betaINBR11) +
  geom_line(aes(y = betaINBR11, x = iteration, colour = CHAIN))
bI11d <- ggplot(betaINBR11) +
  geom_density(aes(x = betaINBR11, colour = CHAIN))

betaINBR12 <- as.data.frame(betaINBR12)
betaINBR12$iteration <- seq(1:nrow(betaINBR12))
bI12 <- ggplot(betaINBR12) +
  geom_line(aes(y = betaINBR12, x = iteration, colour = CHAIN))
bI12d <- ggplot(betaINBR12) +
  geom_density(aes(x = betaINBR12, colour = CHAIN))

betaINBR13 <- as.data.frame(betaINBR13)
betaINBR13$iteration <- seq(1:nrow(betaINBR13))
bI13 <- ggplot(betaINBR13) +
  geom_line(aes(y = betaINBR13, x = iteration, colour = CHAIN))
bI13d <- ggplot(betaINBR13) +
  geom_density(aes(x = betaINBR13, colour = CHAIN))

betaINBR14 <- as.data.frame(betaINBR14)
betaINBR14$iteration <- seq(1:nrow(betaINBR14))
bI14 <- ggplot(betaINBR14) +
  geom_line(aes(y = betaINBR14, x = iteration, colour = CHAIN))
bI14d <- ggplot(betaINBR14) +
  geom_density(aes(x = betaINBR14, colour = CHAIN))

betaINBR15 <- as.data.frame(betaINBR15)
betaINBR15$iteration <- seq(1:nrow(betaINBR15))
bI15 <- ggplot(betaINBR15) +
  geom_line(aes(y = betaINBR15, x = iteration, colour = CHAIN))
bI15d <- ggplot(betaINBR15) +
  geom_density(aes(x = betaINBR15, colour = CHAIN))


bI1/bI2/bI3|bI1d/bI2d/bI3d
bI4/bI5/bI6|bI4d/bI5d/bI6d
bI7/bI8/bI9|bI7d/bI8d/bI9d
bI10/bI11/bI12|bI10d/bI11d/bI12d
bI13/bI14/bI15|bI13d/bI14d/bI15d


## betaINBR main effects
bI1/bI4/bI7|bI1d/bI4d/bI7d
bI10/bI13|bI10d/bI13d

bI13/bI14/bI15|bI13d/bI14d/bI15d


## correlation plots...
samples.corr <- samples %>%
  select(betaINBR1, betaINBR4, betaINBR7, betaINBR10, betaINBR13)

ggpairs(samples.corr)


b1_corr <- betaINBR1$betaINBR1
b4_corr <- betaINBR4$betaINBR4
b7_corr <- betaINBR7$betaINBR7
b10_corr <- betaINBR10$betaINBR10
b13_corr <- betaINBR13$betaINBR13


c_samples <- samples %>%
  filter(z1 == 1 & z4 == 1 & z7 == 1 & z10 == 1 & z13 ==1) %>%
  select(betaINBR1, betaINBR4, betaINBR7, betaINBR10, betaINBR13)


ggpairs(c_samples)


bind_cols()

bi1_p <- (length(which(b1_corr < 0)))/length(b1_corr)
bi4_p <- (length(which(b4_corr > 0)))/length(b4_corr)
bi7_p <- (length(which(b7_corr > 0)))/length(b7_corr)
bi10_p <- (length(which(b10_corr > 0)))/length(b10_corr)
bi13_p <- (length(which(b13_corr > 0)))/length(b13_corr)





length(c(b1_corr, b4_corr, b7_corr, b10_corr, b13_corr)))

length(b10_corr)

b1_corr <- sample(b1_corr, length(b10_corr))
b4_corr <- sample(b1_corr, length(b10_corr))
b7_corr <- sample(b1_corr, length(b10_corr))
b13_corr <- sample(b1_corr, length(b10_corr))

c_table <- data.frame(b1_corr, b4_corr, b7_corr, b10_corr, b13_corr)

ggpairs(c_table)





