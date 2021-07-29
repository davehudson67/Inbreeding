library(haven)
library(survival)

rm(list=ls())

## load data
load("Data/badgerSexInb_AdultInfvCubInfvUninf.RData")


CH$death[is.na(CH$death)] <- CH$last_seen[is.na(CH$death)] + rnorm(1, 2, 1)
CH$dur <- CH$death - CH$birth
CH$censored[CH$censored == 2] <- 0

badgers <- CH %>%
  filter(censored > 0)


badgers <- CH

sex <- badgers$sex
#inbreeding <- ifelse(badgers$hom > median(badgers$hom), 1, 0)
infection <- rep(0, times = nrow(badgers))
infection[badgers$infected_as_cub == 1] <- 1
infection[badgers$infected_lifetime > badgers$infected_as_cub] <- 2


inbr_level0 <- quantile(badgers$hom, probs = 0.333)
inbr_level1 <- quantile(badgers$hom, probs = 0.666)
inbr <- rep(0, times = nrow(badgers))

inbr[badgers$hom >= inbr_level0] <- 1
inbr[badgers$hom >= inbr_level1] <- 2

  
badgersS <- Surv(badgers$dur, badgers$censored)


fit <- survfit(badgersS ~ inbr + sex)

ggsurvplot(fit, data = badgersS, risk.table = TRUE,
           conf.int = FALSE, pval = TRUE, title = "Kaplan-Meier plot of survival probability",
           legend = "right")
           
           
library(muhaz)
plot(muhaz(badgers$dur[inbr == 1], badgers$censored[inbr == 1]))
plot(muhaz(badgers$dur[inbr == 0], badgers$censored[inbr == 0]))

library(KMsurv)
library(biostat3)

lt <- lifetab2(badgersS ~ 1, badgersS)

plot(lt$`inbr=0`, type = "b")
           
           , legend.title = "Sex", legend.labs = c("Female","Male", "Pup"),
           pval.coord = c(1000,0.75), surv.median.line = "v",
           censor.shape=124, ggtheme = theme_bw(), font.main = c(16, "darkblue"), 
           subtitle = "Censored individuals shown by vertical lines, dashed line indicates median age.
           (p-value corresponds to log-rank comparison of survival curves)")
