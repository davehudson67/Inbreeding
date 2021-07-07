## load libraries
library(tidyverse)
library(boot)
library(lamW)
library(nimble)
rm(list = ls())

## source necessary R distributions
source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_Gompertz.R")
source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_GompertzNim.R")
source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_GompertzMakeham.R")
source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_GompertzMakehamNim.R")
source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_Siler.R")
source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_SilerNim.R")
source("../SimulationStudy/FirstPaperFiles/Distributions/Dist_Expo.R")
source("../SimulationStudy/FirstPaperFiles/ModelComparison_FUNCTIONS.R")

## source additional R functions
source("../SimulationStudy/FirstPaperFiles/ModelComparison_FUNCTIONS.R")

## load data
load("badgerSexInb.RData")

## set seed according to model
set.seed(seeds[16])

## set up plot output file
pdf("outputs/ModelComparisons_stage1.pdf")

###########################################################
##                                                      ###
##          Now conduct model comparisons               ###
##                                                      ###
###########################################################

## load IS samples

logimpweight_s <- readRDS("outputs/logimpweight_s_NoInbDiff.rds")
#logimpweight_s_SexDiffa1 <- readRDS("ISoutputs/logimpweight_s_SexDiffa1.rds")
#logimpweight_s_SexDiffa2 <- readRDS("ISoutputs/logimpweight_s_SexDiffa2.rds")
#logimpweight_s_SexDiffb1 <- readRDS("ISoutputs/logimpweight_s_SexDiffb1.rds")
logimpweight_s_InbDiffb2 <- readRDS("outputs/logimpweight_s_InbDiffb2.rds")
#logimpweight_s_SexDiffc1 <- readRDS("ISoutputs/logimpweight_s_SexDiffc1.rds")
#logimpweight_s_SexDiffFULL <- readRDS("ISoutputs/logimpweight_s_SexDiffFULLMODEL.rds")

## generate log marginal likelihoods
logmarg_s <- log_sum_exp_marg(logimpweight_s)
#logmarg_s_InbDiffa1 <- log_sum_exp_marg(logimpweight_s_InbDiffa1)
#logmarg_s_InbDiffa2 <- log_sum_exp_marg(logimpweight_s_InbDiffa2)
#logmarg_s_InbDiffb1 <- log_sum_exp_marg(logimpweight_s_InbDiffb1)
logmarg_s_InbDiffb2 <- log_sum_exp_marg(logimpweight_s_InbDiffb2)
#logmarg_s_InbDiffc1 <- log_sum_exp_marg(logimpweight_s_InbDiffc1)
#logmarg_s_InbDiffFULL <- log_sum_exp_marg(logimpweight_s_InbDiffFULL)

## bootstrap samples
imp_boot_s <- BootsPlot(logimpweight_s, 5000)
#imp_boot_s_InbDiffa1 <- BootsPlot(logimpweight_s_InbDiffa1, 5000)
#imp_boot_s_InbDiffa2 <- BootsPlot(logimpweight_s_InbDiffa2, 5000)
#imp_boot_s_InbDiffb1 <- BootsPlot(logimpweight_s_InbDiffb1, 5000)
imp_boot_s_InbDiffb2 <- BootsPlot(logimpweight_s_InbDiffb2, 5000)
#imp_boot_s_InbDiffc1 <- BootsPlot(logimpweight_s_InbDiffc1, 5000)
#imp_boot_s_InbDiffFULL <- BootsPlot(logimpweight_s_InbDiffFULL, 5000)

## add prior model weights
priorp <- 1/2
p_s <- logmarg_s + log(priorp)
#p_s_InbDiffa1 <- logmarg_s_InbDiffa1 + log(priorp)
#p_s_InbDiffa2 <- logmarg_s_InbDiffa2 + log(priorp)
#p_s_InbDiffb1 <- logmarg_s_InbDiffb1 + log(priorp)
p_s_InbDiffb2 <- logmarg_s_InbDiffb2 + log(priorp)
#p_s_InbDiffc1 <- logmarg_s_InbDiffc1 + log(priorp)
#p_s_InbDiffFULL <- logmarg_s_InbDiffFULL + log(priorp)


p <- c(p_s, p_s_InbDiffb2)

pd <- log_sum_exp_marg(p, mn = FALSE)

## normalise
p <- p - pd
p <- exp(p)
p

## plot marginal likelihoods
mods <- list(
  S = imp_boot_s,
#  S_sda1 = imp_boot_s_InbDiffa1,
#  S_sda2 = imp_boot_s_InbDiffa2,
#  S_sdb1 = imp_boot_s_InbDiffb1,
  S_sdb2 = imp_boot_s_InbDiffb2
#  S_sdc1 = imp_boot_s_InbDiffc1,
#  S_sdFULL = imp_boot_s_InbDiffFULL
)

MargLike.plot(mods)

## which models within log(20) of best
logmarg <- map_dbl(mods, "logmarg")
bestind <- which(logmarg == max(logmarg))
logmargLCI <- mods[[bestind]]$LCI
logmarg <- map_dbl(mods, "UCI")
logmarg <- logmarg[map_lgl(logmarg, ~ . >= logmargLCI - log(20))]
logmarg

## turn graphics device off
dev.off()
