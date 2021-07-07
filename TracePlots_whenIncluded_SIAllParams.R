library(tidyverse)
library(patchwork)
library(MCMCvis)

samples <- readRDS("outputs/Sex3Infection_AllParameters_runsamples.rds")

mcmcplot(samples$samples, parms = c("a1", "a2", "b1", "b2", "c1"))

samples <- as.matrix(samples$samples, chains = TRUE)
samples <- as.data.frame(samples)

names(samples)
colnames(samples) <- c("CHAIN", "a1", "a2", "b1", "b2", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "beta8", "beta9", "beta10", "beta11", "beta12", "beta13",
              "beta14", "beta15", "betaINFADULT1", "betaINFADULT2", "betaINFADULT3", "betaINFADULT4", "betaINFADULT5", "betaINFCUB1", "betaINFCUB2", "betaINFCUB3", 
              "betaINFCUB4", "betaINFCUB5", "betaSexINFADULT1", "betaSexINFADULT2", "betaSexINFADULT3", "betaSexINFADULT4", "betaSexINFADULT5", "betaSexINFCUB1",
              "betaSexINFCUB2", "betaSexINFCUB3", "betaSexINFCUB4", "betaSexINFCUB5", "c1", "mean.p", "z1", "z2", "z3", "z4", "z5", "z6", "z7", "z8", "z9", "z10", "z11",
              "z12", "z13", "z14", "z15") 

samples$CHAIN <- as.factor(samples$CHAIN)

beta1 <- samples %>%
  filter(z1 == 1) %>%
  select(CHAIN, beta1) %>%
  mutate(iteration = 1:length(beta1))

beta2 <- samples %>%
  filter(z2 == 1) %>%
  select(CHAIN, beta2) %>%
  mutate(iteration = 1:length(beta2))

beta3 <- samples %>%
  filter(z3 == 1) %>%
  select(CHAIN, beta3) %>%
  mutate(iteration = 1:length(beta3))

beta4 <- samples %>%
  filter(z4 == 1) %>%
  select(CHAIN, beta4) %>%
  mutate(iteration = 1:length(beta4))

beta5 <- samples %>%
  filter(z5 == 1) %>%
  select(CHAIN, beta5) %>%
  mutate(iteration = 1:length(beta5))

beta6 <- samples %>%
  filter(z6 == 1) %>%
  select(CHAIN, beta6) %>%
  mutate(iteration = 1:length(beta6))

beta7 <- samples %>%
  filter(z7 == 1) %>%
  select(CHAIN, beta7) %>%
  mutate(iteration = 1:length(beta7))

beta8 <- samples %>%
  filter(z8 == 1) %>%
  select(CHAIN, beta8) %>%
  mutate(iteration = 1:length(beta8))

beta9 <- samples %>%
  filter(z9 == 1) %>%
  select(CHAIN, beta9) %>%
  mutate(iteration = 1:length(beta9))

beta10 <- samples %>%
  filter(z10 == 1) %>%
  select(CHAIN, beta10) %>%
  mutate(iteration = 1:length(beta10))

beta11 <- samples %>%
  filter(z11 == 1) %>%
  select(CHAIN, beta11) %>%
  mutate(iteration = 1:length(beta11))

beta12 <- samples %>%
  filter(z12 == 1) %>%
  select(CHAIN, beta12) %>%
  mutate(iteration = 1:length(beta12))

beta13 <- samples %>%
  filter(z13 == 1) %>%
  select(CHAIN, beta13) %>%
  mutate(iteration = 1:length(beta13))

beta14 <- samples %>%
  filter(z14 == 1) %>%
  select(CHAIN, beta14) %>%
  mutate(iteration = 1:length(beta14))

beta15 <- samples %>%
  filter(z15 == 1) %>%
  select(CHAIN, beta15) %>%
  mutate(iteration = 1:length(beta15))

##

betaINFADULT1 <- samples %>%
  filter(z2 == 1) %>%
  select(CHAIN, betaINFADULT1) %>%
  mutate(iteration = 1:length(betaINFADULT1))

betaINFADULT2 <- samples %>%
  filter(z5 == 1) %>%
  select(CHAIN, betaINFADULT2) %>%
  mutate(iteration = 1:length(betaINFADULT2))

betaINFADULT3 <- samples %>%
  filter(z8 == 1) %>%
  select(CHAIN, betaINFADULT3) %>%
  mutate(iteration = 1:length(betaINFADULT3))

betaINFADULT4 <- samples %>%
  filter(z11 == 1) %>%
  select(CHAIN, betaINFADULT4) %>%
  mutate(iteration = 1:length(betaINFADULT4))

betaINFADULT5 <- samples %>%
  filter(z14 == 1) %>%
  select(CHAIN, betaINFADULT5) %>%
  mutate(iteration = 1:length(betaINFADULT5))

##

betaINFCUB1 <- samples %>%
  filter(z2 == 1) %>%
  select(CHAIN, betaINFCUB1) %>%
  mutate(iteration = 1:length(betaINFCUB1))

betaINFCUB2 <- samples %>%
  filter(z5 == 1) %>%
  select(CHAIN, betaINFCUB2) %>%
  mutate(iteration = 1:length(betaINFCUB2))

betaINFCUB3 <- samples %>%
  filter(z8 == 1) %>%
  select(CHAIN, betaINFCUB3) %>%
  mutate(iteration = 1:length(betaINFCUB3))

betaINFCUB4 <- samples %>%
  filter(z11 == 1) %>%
  select(CHAIN, betaINFCUB4) %>%
  mutate(iteration = 1:length(betaINFCUB4))

betaINFCUB5 <- samples %>%
  filter(z14 == 1) %>%
  select(CHAIN, betaINFCUB5) %>%
  mutate(iteration = 1:length(betaINFCUB5))

##

betaSexINFCUB1 <- samples %>%
  filter(z3 == 1) %>%
  select(CHAIN, betaSexINFCUB1) %>%
  mutate(iteration = 1:length(betaSexINFCUB1))

betaSexINFCUB2 <- samples %>%
  filter(z6 == 1) %>%
  select(CHAIN, betaSexINFCUB2) %>%
  mutate(iteration = 1:length(betaSexINFCUB2))

betaSexINFCUB3 <- samples %>%
  filter(z9 == 1) %>%
  select(CHAIN, betaSexINFCUB3) %>%
  mutate(iteration = 1:length(betaSexINFCUB3))

betaSexINFCUB4 <- samples %>%
  filter(z12 == 1) %>%
  select(CHAIN, betaSexINFCUB4) %>%
  mutate(iteration = 1:length(betaSexINFCUB4))

betaSexINFCUB5 <- samples %>%
  filter(z15 == 1) %>%
  select(CHAIN, betaSexINFCUB5) %>%
  mutate(iteration = 1:length(betaSexINFCUB5))

##

betaSexINFADULT1 <- samples %>%
  filter(z3 == 1) %>%
  select(CHAIN, betaSexINFADULT1) %>%
  mutate(iteration = 1:length(betaSexINFADULT1))

betaSexINFADULT2 <- samples %>%
  filter(z6 == 1) %>%
  select(CHAIN, betaSexINFADULT2) %>%
  mutate(iteration = 1:length(betaSexINFADULT2))

betaSexINFADULT3 <- samples %>%
  filter(z9 == 1) %>%
  select(CHAIN, betaSexINFADULT3) %>%
  mutate(iteration = 1:length(betaSexINFADULT3))

betaSexINFADULT4 <- samples %>%
  filter(z12 == 1) %>%
  select(CHAIN, betaSexINFADULT4) %>%
  mutate(iteration = 1:length(betaSexINFADULT4))

betaSexINFADULT5 <- samples %>%
  filter(z15 == 1) %>%
  select(CHAIN, betaSexINFADULT5) %>%
  mutate(iteration = 1:length(betaSexINFADULT5))

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

beta11 <- as.data.frame(beta11)
beta11$iteration <- seq(1:nrow(beta11))
b11 <- ggplot(beta11) +
  geom_line(aes(y = beta11, x = iteration, colour = CHAIN))
b11d <- ggplot(beta11) +
  geom_density(aes(x = beta11, colour = CHAIN))

beta12 <- as.data.frame(beta12)
beta12$iteration <- seq(1:nrow(beta12))
b12 <- ggplot(beta12) +
  geom_line(aes(y = beta12, x = iteration, colour = CHAIN))
b12d <- ggplot(beta12) +
  geom_density(aes(x = beta12, colour = CHAIN))

beta13 <- as.data.frame(beta13)
beta13$iteration <- seq(1:nrow(beta13))
b13 <- ggplot(beta13) +
  geom_line(aes(y = beta13, x = iteration, colour = CHAIN))
b13d <- ggplot(beta13) +
  geom_density(aes(x = beta13, colour = CHAIN))

beta14 <- as.data.frame(beta14)
beta14$iteration <- seq(1:nrow(beta14))
b14 <- ggplot(beta14) +
  geom_line(aes(y = beta14, x = iteration, colour = CHAIN))
b14d <- ggplot(beta14) +
  geom_density(aes(x = beta14, colour = CHAIN))

beta15 <- as.data.frame(beta15)
beta15$iteration <- seq(1:nrow(beta15))
b15 <- ggplot(beta15) +
  geom_line(aes(y = beta15, x = iteration, colour = CHAIN))
b15d <- ggplot(beta15) +
  geom_density(aes(x = beta15, colour = CHAIN))


b1/b2/b3|b1d/b2d/b3d
b4/b5/b6|b4d/b5d/b6d
b7/b8/b9|b7d/b8d/b9d
b10/b11/b12|b10d/b11d/b12d
b13/b14/b15|b13d/b14d/b15d

##

betaSexINFADULT1 <- as.data.frame(betaSexINFADULT1)
betaSexINFADULT1$iteration <- seq(1:nrow(betaSexINFADULT1))
bsiAA1 <- ggplot(betaSexINFADULT1) +
  geom_line(aes(x = iteration, y = betaSexINFADULT1, colour = CHAIN))
bsiAA1d <- ggplot(betaSexINFADULT1) + 
  geom_density(aes(x = betaSexINFADULT1, colour = CHAIN))

betaSexINFADULT2 <- as.data.frame(betaSexINFADULT2)
betaSexINFADULT2$iteration <- seq(1:nrow(betaSexINFADULT2))
bsiAA2 <- ggplot(betaSexINFADULT2) +
  geom_line(aes(x = iteration, y = betaSexINFADULT2, colour = CHAIN))
bsiAA2d <- ggplot(betaSexINFADULT2) + 
  geom_density(aes(x = betaSexINFADULT2, colour = CHAIN))

betaSexINFADULT3 <- as.data.frame(betaSexINFADULT3)
betaSexINFADULT3$iteration <- seq(1:nrow(betaSexINFADULT3))
bsiAB1 <- ggplot(betaSexINFADULT3) +
  geom_line(aes(x = iteration, y = betaSexINFADULT3, colour = CHAIN))
bsiAB1d <- ggplot(betaSexINFADULT3) + 
  geom_density(aes(x = betaSexINFADULT3, colour = CHAIN))

betaSexINFADULT4 <- as.data.frame(betaSexINFADULT4)
betaSexINFADULT4$iteration <- seq(1:nrow(betaSexINFADULT4))
bsiAB2 <- ggplot(betaSexINFADULT4) +
  geom_line(aes(x = iteration, y = betaSexINFADULT4, colour = CHAIN))
bsiAB2d <- ggplot(betaSexINFADULT4) + 
  geom_density(aes(x = betaSexINFADULT4, colour = CHAIN))

betaSexINFADULT5 <- as.data.frame(betaSexINFADULT5)
betaSexINFADULT5$iteration <- seq(1:nrow(betaSexINFADULT5))
bsiAC1 <- ggplot(betaSexINFADULT5) +
  geom_line(aes(x = iteration, y = betaSexINFADULT5, colour = CHAIN))
bsiAC1d <- ggplot(betaSexINFADULT5) + 
  geom_density(aes(x = betaSexINFADULT5, colour = CHAIN))

##

betaSexINFCUB1 <- as.data.frame(betaSexINFCUB1)
betaSexINFCUB1$iteration <- seq(1:nrow(betaSexINFCUB1))
bsiCA1 <- ggplot(betaSexINFCUB1) +
  geom_line(aes(x = iteration, y = betaSexINFCUB1, colour = CHAIN))
bsiCA1d <- ggplot(betaSexINFCUB1) + 
  geom_density(aes(x = betaSexINFCUB1, colour = CHAIN))

betaSexINFCUB2 <- as.data.frame(betaSexINFCUB2)
betaSexINFCUB2$iteration <- seq(1:nrow(betaSexINFCUB2))
bsiCA2 <- ggplot(betaSexINFCUB2) +
  geom_line(aes(x = iteration, y = betaSexINFCUB2, colour = CHAIN))
bsiCA2d <- ggplot(betaSexINFCUB2) + 
  geom_density(aes(x = betaSexINFCUB2, colour = CHAIN))

betaSexINFCUB3 <- as.data.frame(betaSexINFCUB3)
betaSexINFCUB3$iteration <- seq(1:nrow(betaSexINFCUB3))
bsiCB1 <- ggplot(betaSexINFCUB3) +
  geom_line(aes(x = iteration, y = betaSexINFCUB3, colour = CHAIN))
bsiCB1d <- ggplot(betaSexINFCUB3) + 
  geom_density(aes(x = betaSexINFCUB3, colour = CHAIN))

betaSexINFCUB4 <- as.data.frame(betaSexINFCUB4)
betaSexINFCUB4$iteration <- seq(1:nrow(betaSexINFCUB4))
bsiCB2 <- ggplot(betaSexINFCUB4) +
  geom_line(aes(x = iteration, y = betaSexINFCUB4, colour = CHAIN))
bsiCB2d <- ggplot(betaSexINFCUB4) + 
  geom_density(aes(x = betaSexINFCUB4, colour = CHAIN))

betaSexINFCUB5 <- as.data.frame(betaSexINFCUB5)
betaSexINFCUB5$iteration <- seq(1:nrow(betaSexINFCUB5))
bsiCC1 <- ggplot(betaSexINFCUB5) +
  geom_line(aes(x = iteration, y = betaSexINFCUB5, colour = CHAIN))
bsiCC1d <- ggplot(betaSexINFCUB5) + 
  geom_density(aes(x = betaSexINFCUB5, colour = CHAIN))

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

