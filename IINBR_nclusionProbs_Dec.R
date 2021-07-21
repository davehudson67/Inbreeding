r09 <- readRDS("outputs/FullModel_NoInteractions_z_all_InbrD09_runsamples.rds")
r08 <- readRDS("outputs/FullModel_NoInteractions_z_all_InbrD08_runsamples.rds")
r06 <- readRDS("outputs/FullModel_NoInteractions_z_all_InbrD06_runsamples.rds")
r04 <- readRDS("outputs/FullModel_NoInteractions_z_all_InbrD04_runsamples.rds")
r02 <- readRDS("outputs/FullModel_NoInteractions_z_all_InbrD02_runsamples.rds")
r05 <- readRDS("outputs/FullModel_NoInteractions_z_all_runsamples.rds")

#run <- readRDS("outputs/FullModel_NoInteractions_z_all_runsamples.rds")

r09 <- as.matrix(r09$samples)
r08 <- as.matrix(r08$samples)
r06 <- as.matrix(r06$samples)
r04 <- as.matrix(r04$samples)
r02 <- as.matrix(r02$samples)

#MCMCsummary(run$samples)

## Marginal probabilities of inclusion for each variable
zNames <- model$expandNodeNames(c('z', 'zSEX', 'zINF', 'zINBR'))
zCols <- which(colnames(r09) %in% zNames)


binary <- as.data.table((r08[, zCols] != 0) + 0)
res <- binary[ , .N, by=names(binary)]
res <- res[order(N, decreasing = T)]
res <- res[, prob := N/dim(r08)[1]]
res
res
sum(res$N)

saveRDS(res, "outputs/FullModel_InbrCont_NoInteractions_z_all_PosteriorModelProbs.rds")
#res <- readRDS("outputs/Sex3Infection_AllParameters_PosteriorModelProbs.rds")

r <- as.data.frame(r02)
z_indicators <- r %>%
  select(c(27:31)) %>%
  colSums()
z_indicators <- data.frame(z_indicators/sum(res$N))
z_indicators <- as.data.frame(t(z_indicators))
z_2 <- z_indicators %>%
  mutate(dec = 0.2)

r <- as.data.frame(r04)
z_indicators <- r %>%
  select(c(27:31)) %>%
  colSums()
z_indicators <- data.frame(z_indicators/sum(res$N))
z_indicators <- as.data.frame(t(z_indicators))
z_4 <- z_indicators %>%
  mutate(dec = 0.4)

r <- as.data.frame(r06)
z_indicators <- r %>%
  select(c(27:31)) %>%
  colSums()
z_indicators <- data.frame(z_indicators/sum(res$N))
z_indicators <- as.data.frame(t(z_indicators))
z_6 <- z_indicators %>%
  mutate(dec = 0.6)

r <- as.data.frame(r08)
z_indicators <- r %>%
  select(c(27:31)) %>%
  colSums()
z_indicators <- data.frame(z_indicators/sum(res$N))
z_indicators <- as.data.frame(t(z_indicators))
z_8 <- z_indicators %>%
  mutate(dec = 0.8)

r <- as.data.frame(r09)
z_indicators <- r %>%
  select(c(27:31)) %>%
  colSums()
z_indicators <- data.frame(z_indicators/sum(res$N))
z_indicators <- as.data.frame(t(z_indicators))
z_9 <- z_indicators %>%
  mutate(dec = 0.9)

z <- bind_rows(z_2, z_4, z_6, z_8, z_9)
colnames(z) <- c("a1", "a2", "b1", "b2", "c", "dec")

zl <- pivot_longer(z, cols = c("a1", "a2", "b1", "b2", "c"), names_to = "Inbr_on")

#zl$dec <- as.factor(zl$dec)
ggplot(zl, aes(x = dec, y = value)) +
  geom_line(aes(colour = Inbr_on)) +
  ylab("Inclusion_prob")
