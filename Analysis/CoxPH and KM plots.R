library(tidyverse)
library(data.table)
library(survival)
library(haven)
library(survminer)

source("TCH/TCH for Nimble CP_Inbreeding.R")

## sort data
PosBad <- as.data.table(P.CVdata[,c(1,9,16,17,19,20,21,25,29)])

## create duration variable
PosBad$death[is.na(PosBad$death)] <- 147
PosBad$dur <- PosBad$death - PosBad$birth_yr

## cox PH
PosBad.s <- Surv(PosBad$dur, PosBad$knowndeath)
PosBad.sfit <- coxph(PosBad.s ~ PosBad$f_inbreed * PosBad$sex, ties = "breslow")
summary(PosBad.sfit)


## select only cub positive badgers
PosBad.CP <- PosBad[PosBad$infected.as.cub==1]

## cox PH
PosBad.CP.s <- Surv(PosBad.CP$dur, PosBad.CP$knowndeath)
PosBad.CP.sfit <- coxph(PosBad.CP.s ~ PosBad.CP$f_inbreed * PosBad.CP$sex, ties = "breslow")
summary(PosBad.CP.sfit)

##KM plots
fit1 <- survfit(Surv(PosBad$dur, PosBad$knowndeath) ~ PosBad$Inb)
#Survival curves
ggsurvplot(fit1, data = PosBad, risk.table = TRUE,
           conf.int = TRUE, pval = TRUE, title = "Kaplan-Meier plot of survival probability",
           legend = "right", legend.title = "Inbreeding", legend.labs = c("Low","High"),
           pval.coord = c(2,0.75), surv.median.line = "v",
           censor.shape=124, ggtheme = theme_bw(), font.main = c(16, "darkblue"), 
           subtitle = "Censored individuals shown by vertical lines, dashed line indicates median age.
           (p-value corresponds to log-rank comparison of survival curves)")

fit2 <- survfit(Surv(PosBad.CP$dur, PosBad.CP$knowndeath) ~ PosBad.CP$Inb)
#Survival curves
ggsurvplot(fit2, data = PosBad.CP, risk.table = TRUE,
           conf.int = TRUE, pval = TRUE, title = "Kaplan-Meier plot of survival probability",
           legend = "right", legend.title = "Inbreeding", legend.labs = c("Low","High"),
           pval.coord = c(2,0.75), surv.median.line = "v",
           censor.shape=124, ggtheme = theme_bw(), font.main = c(16, "darkblue"), 
           subtitle = "Censored individuals shown by vertical lines, dashed line indicates median age.
           (p-value corresponds to log-rank comparison of survival curves)")
