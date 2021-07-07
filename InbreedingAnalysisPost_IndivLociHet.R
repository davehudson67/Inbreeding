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
pdf("outputs/ModelComparisons_IndivLociHetIR.pdf")

###########################################################
##                                                      ###
##          Now conduct model comparisons               ###
##                                                      ###
###########################################################

## load IS samples
logimpweight_l1 <- readRDS("IndivLociOutputs/logimpweight_s_b2L1.rds")
logimpweight_l2 <- readRDS("IndivLociOutputs/logimpweight_s_b2L2.rds")
logimpweight_l3 <- readRDS("IndivLociOutputs/logimpweight_s_b2L3.rds")
logimpweight_l4 <- readRDS("IndivLociOutputs/logimpweight_s_b2L4.rds")
logimpweight_l5 <- readRDS("IndivLociOutputs/logimpweight_s_b2L5.rds")
logimpweight_l6 <- readRDS("IndivLociOutputs/logimpweight_s_b2L6.rds")
logimpweight_l7 <- readRDS("IndivLociOutputs/logimpweight_s_b2L7.rds")
logimpweight_l8 <- readRDS("IndivLociOutputs/logimpweight_s_b2L8.rds")
logimpweight_l9 <- readRDS("IndivLociOutputs/logimpweight_s_b2L9.rds")
logimpweight_l10 <- readRDS("IndivLociOutputs/logimpweight_s_b2L10.rds")
logimpweight_l11 <- readRDS("IndivLociOutputs/logimpweight_s_b2L11.rds")
logimpweight_l12 <- readRDS("IndivLociOutputs/logimpweight_s_b2L12.rds")
logimpweight_l13 <- readRDS("IndivLociOutputs/logimpweight_s_b2L13.rds")
logimpweight_l14 <- readRDS("IndivLociOutputs/logimpweight_s_b2L14.rds")
logimpweight_l15 <- readRDS("IndivLociOutputs/logimpweight_s_b2L15.rds")
logimpweight_l16 <- readRDS("IndivLociOutputs/logimpweight_s_b2L16.rds")
logimpweight_l17 <- readRDS("IndivLociOutputs/logimpweight_s_b2L17.rds")
logimpweight_l18 <- readRDS("IndivLociOutputs/logimpweight_s_b2L18.rds")
logimpweight_l19 <- readRDS("IndivLociOutputs/logimpweight_s_b2L19.rds")
logimpweight_l20 <- readRDS("IndivLociOutputs/logimpweight_s_b2L20.rds")
logimpweight_l21 <- readRDS("IndivLociOutputs/logimpweight_s_b2L21.rds")
logimpweight_l22 <- readRDS("IndivLociOutputs/logimpweight_s_b2L22.rds")
logimpweight_b2het <- readRDS("IndivLociOutputs/logimpweight_s_b2het.rds")
logimpweight_b2IR <- readRDS("IndivLociOutputs/logimpweight_s_b2IR.rds")
logimpweight_b2hetCAT <- readRDS("IndivLociOutputs/logimpweight_s_b2hetCAT.rds")

## generate log marginal likelihoods
logmarg_l1 <- log_sum_exp_marg(logimpweight_l1)
logmarg_l2 <- log_sum_exp_marg(logimpweight_l2)
logmarg_l3 <- log_sum_exp_marg(logimpweight_l3)
logmarg_l4 <- log_sum_exp_marg(logimpweight_l4)
logmarg_l5 <- log_sum_exp_marg(logimpweight_l5)
logmarg_l6 <- log_sum_exp_marg(logimpweight_l6)
logmarg_l7 <- log_sum_exp_marg(logimpweight_l7)
logmarg_l8 <- log_sum_exp_marg(logimpweight_l8)
logmarg_l9 <- log_sum_exp_marg(logimpweight_l9)
logmarg_l10 <- log_sum_exp_marg(logimpweight_l10)
logmarg_l11 <- log_sum_exp_marg(logimpweight_l11)
logmarg_l12 <- log_sum_exp_marg(logimpweight_l12)
logmarg_l13 <- log_sum_exp_marg(logimpweight_l13)
logmarg_l14 <- log_sum_exp_marg(logimpweight_l14)
logmarg_l15 <- log_sum_exp_marg(logimpweight_l15)
logmarg_l16 <- log_sum_exp_marg(logimpweight_l16)
logmarg_l17 <- log_sum_exp_marg(logimpweight_l17)
logmarg_l18 <- log_sum_exp_marg(logimpweight_l18)
logmarg_l19 <- log_sum_exp_marg(logimpweight_l19)
logmarg_l20 <- log_sum_exp_marg(logimpweight_l20)
logmarg_l21 <- log_sum_exp_marg(logimpweight_l21)
logmarg_l22 <- log_sum_exp_marg(logimpweight_l22)
logmarg_het <- log_sum_exp_marg(logimpweight_b2het)
logmarg_IR <- log_sum_exp_marg(logimpweight_b2IR)
logmarg_hetCAT <- log_sum_exp_marg(logimpweight_b2hetCAT)

## bootstrap samples
imp_boot_l1 <- BootsPlot(logimpweight_l1, 5000)
imp_boot_l2 <- BootsPlot(logimpweight_l2, 5000)
imp_boot_l3 <- BootsPlot(logimpweight_l3, 5000)
imp_boot_l4 <- BootsPlot(logimpweight_l4, 5000)
imp_boot_l5 <- BootsPlot(logimpweight_l5, 5000)
imp_boot_l6 <- BootsPlot(logimpweight_l6, 5000)
imp_boot_l7 <- BootsPlot(logimpweight_l7, 5000)
imp_boot_l8 <- BootsPlot(logimpweight_l8, 5000)
imp_boot_l9 <- BootsPlot(logimpweight_l9, 5000)
imp_boot_l10 <- BootsPlot(logimpweight_l10, 5000)
imp_boot_l11 <- BootsPlot(logimpweight_l11, 5000)
imp_boot_l12 <- BootsPlot(logimpweight_l12, 5000)
imp_boot_l13 <- BootsPlot(logimpweight_l13, 5000)
imp_boot_l14 <- BootsPlot(logimpweight_l14, 5000)
imp_boot_l15 <- BootsPlot(logimpweight_l15, 5000)
imp_boot_l16 <- BootsPlot(logimpweight_l16, 5000)
imp_boot_l17 <- BootsPlot(logimpweight_l17, 5000)
imp_boot_l18 <- BootsPlot(logimpweight_l18, 5000)
imp_boot_l19 <- BootsPlot(logimpweight_l19, 5000)
imp_boot_l20 <- BootsPlot(logimpweight_l20, 5000)
imp_boot_l21 <- BootsPlot(logimpweight_l21, 5000)
imp_boot_l22 <- BootsPlot(logimpweight_l22, 5000)
imp_boot_het <- BootsPlot(logimpweight_b2het, 5000)
imp_boot_IR <- BootsPlot(logimpweight_b2IR, 5000)
imp_boot_hetCAT <- BootsPlot(logimpweight_b2hetCAT, 5000)

## add prior model weights
priorp <- 1/25
p_l1 <- logmarg_l1 + log(priorp)
p_l2 <- logmarg_l2 + log(priorp)
p_l3 <- logmarg_l3 + log(priorp)
p_l4 <- logmarg_l4 + log(priorp)
p_l5 <- logmarg_l5 + log(priorp)
p_l6 <- logmarg_l6 + log(priorp)
p_l7 <- logmarg_l7 + log(priorp)
p_l8 <- logmarg_l8 + log(priorp)
p_l9 <- logmarg_l9 + log(priorp)
p_l10 <- logmarg_l10 + log(priorp)
p_l11 <- logmarg_l11 + log(priorp)
p_l12 <- logmarg_l12 + log(priorp)
p_l13 <- logmarg_l13 + log(priorp)
p_l14 <- logmarg_l14 + log(priorp)
p_l15 <- logmarg_l15 + log(priorp)
p_l16 <- logmarg_l16 + log(priorp)
p_l17 <- logmarg_l17 + log(priorp)
p_l18 <- logmarg_l18 + log(priorp)
p_l19 <- logmarg_l19 + log(priorp)
p_l20 <- logmarg_l20 + log(priorp)
p_l21 <- logmarg_l21 + log(priorp)
p_l22 <- logmarg_l22 + log(priorp)
p_het <- logmarg_het + log(priorp)
p_IR <- logmarg_IR + log(priorp)
p_hetCAT <- logmarg_hetCAT + log(priorp)

p <- c(p_l1, p_l2, p_l3, p_l4, p_l5, p_l6, p_l7, p_l8, p_l9, p_l10, p_l11,
       p_l12, p_l13, p_l14, p_l15, p_l16, p_l17, p_l18, p_l19, p_l20, p_l21,
       p_l22, p_het, p_IR, p_hetCAT)

pd <- log_sum_exp_marg(p, mn = FALSE)

## normalise
p <- p - pd
p <- exp(p)
p

## plot marginal likelihoods
mods <- list(
  p_l1 = imp_boot_l1,
  p_l2 = imp_boot_l2,
  p_l3 = imp_boot_l3,
  p_l4 = imp_boot_l4,
  p_l5 = imp_boot_l5,
  p_l6 = imp_boot_l6,
  p_l7 = imp_boot_l7,
  p_l8 = imp_boot_l8,
  p_l9 = imp_boot_l9,
  p_l10 = imp_boot_l10,
  p_l11 = imp_boot_l11,
  p_l12 = imp_boot_l12,
  p_l13 = imp_boot_l13,
  p_l14 = imp_boot_l14,
  p_l15 = imp_boot_l15,
  p_l16 = imp_boot_l16,
  p_l17 = imp_boot_l17,
  p_l18 = imp_boot_l18,
  p_l19 = imp_boot_l19,
  p_l20 = imp_boot_l20,
  p_l21 = imp_boot_l21,
  p_l22 = imp_boot_l22,
  p_het = imp_boot_het,
  p_IR = imp_boot_IR,
  p_hetCAT = imp_boot_hetCAT
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
