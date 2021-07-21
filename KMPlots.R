library(haven)
library(survival)

rm(list=ls())

## load data
load("Data/badgerSexInb_AdultInfvCubInfvUninf.RData")


CH$death[is.na(CH$death)] <- CH$birth[is.na(CH$death)] + 12
CH$dur <- CH$death - CH$birth
CH$censored[CH$censored == 2] <- 0

badgers <- CH %>%
  filter(censored > 0)

badgers <- CH

sex <- badgers$sex
inbreeding <- ifelse(badgers$hom > median(badgers$hom), 1, 0)
infection <- rep(0, times = nrow(badgers))
infection[badgers$infected_as_cub == 1] <- 1
infection[badgers$infected_lifetime > badgers$infected_as_cub] <- 2

  
badgers <- Surv(badgers$dur, badgers$censored)


fit <- survfit(badgers ~ inbreeding+ sex)

ggsurvplot(fit, data = badgers, risk.table = FALSE,
           conf.int = FALSE, pval = TRUE, title = "Kaplan-Meier plot of survival probability",
           legend = "right")
           
           
           
           , legend.title = "Sex", legend.labs = c("Female","Male", "Pup"),
           pval.coord = c(1000,0.75), surv.median.line = "v",
           censor.shape=124, ggtheme = theme_bw(), font.main = c(16, "darkblue"), 
           subtitle = "Censored individuals shown by vertical lines, dashed line indicates median age.
           (p-value corresponds to log-rank comparison of survival curves)")
