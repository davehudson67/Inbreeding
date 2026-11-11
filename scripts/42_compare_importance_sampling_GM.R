## load libraries
library(tidyverse)
library(boot)
library(lamW)
library(nimble)
rm(list = ls())

## source necessary R distributions
source("../../../FirstPaperFiles/Distributions/Dist_Gompertz.R")
source("../../../FirstPaperFiles/Distributions/Dist_GompertzNim.R")
source("../../../FirstPaperFiles/Distributions/Dist_GompertzMakeham.R")
source("../../../FirstPaperFiles/Distributions/Dist_GompertzMakehamNim.R")
source("../../../FirstPaperFiles/Distributions/Dist_Siler.R")
source("../../../FirstPaperFiles/Distributions/Dist_SilerNim.R")
source("../../../FirstPaperFiles/Distributions/Dist_Expo.R")

## source additional R functions
source("../../../FirstPaperFiles/ModelComparison_FUNCTIONS.R")

## load data
load("../badgerSex.RData")

## set seed according to model
#set.seed(seeds[16])

## set up plot output file
#pdf("ISoutputs/ModelComparisons_StageAllModels.pdf")

###########################################################
##                                                      ###
##          Now conduct model comparisons               ###
##                                                      ###
###########################################################

## load IS samples
logimpweight_gm_SexDiffFULL <- readRDS("ISoutputs/GM/logimpweight_gm_SexDiffabc.rds")
logimpweight_gm_NoSexDiff <- readRDS("ISoutputs/GM/logimpweight_gm_NoSexDiff.rds")

logimpweight_gm_SexDiffa <- readRDS("ISoutputs/GM/logimpweight_gm_SexDiffa.rds")
logimpweight_gm_SexDiffb <- readRDS("ISoutputs/GM/logimpweight_gm_SexDiffb.rds")
logimpweight_gm_SexDiffc <- readRDS("ISoutputs/GM/logimpweight_gm_SexDiffc.rds")
logimpweight_gm_SexDiffab <- readRDS("ISoutputs/GM/logimpweight_gm_SexDiffab.rds")
logimpweight_gm_SexDiffac <- readRDS("ISoutputs/GM/logimpweight_gm_SexDiffac.rds")
logimpweight_gm_SexDiffbc <- readRDS("ISoutputs/GM/logimpweight_gm_SexDiffbc.rds")

## generate log marginal likelihoods
logmarg_gm_SexDiffFULL <- log_sum_exp_marg(logimpweight_gm_SexDiffFULL)
logmarg_gm_NoSexDiff <- log_sum_exp_marg(logimpweight_gm_NoSexDiff)
logmarg_gm_SexDiffa <- log_sum_exp_marg(logimpweight_gm_SexDiffa)
logmarg_gm_SexDiffb <- log_sum_exp_marg(logimpweight_gm_SexDiffb)
logmarg_gm_SexDiffc <- log_sum_exp_marg(logimpweight_gm_SexDiffc)
logmarg_gm_SexDiffab <- log_sum_exp_marg(logimpweight_gm_SexDiffab)
logmarg_gm_SexDiffac <- log_sum_exp_marg(logimpweight_gm_SexDiffac)
logmarg_gm_SexDiffbc <- log_sum_exp_marg(logimpweight_gm_SexDiffbc)

## bootstrap samples
imp_boot_gm_SexDiffFULL <- BootsPlot(logimpweight_gm_SexDiffFULL, 5000)
imp_boot_gm_NoSexDiff <- BootsPlot(logimpweight_gm_NoSexDiff, 5000)
imp_boot_gm_SexDiffa <- BootsPlot(logimpweight_gm_SexDiffa, 5000)
imp_boot_gm_SexDiffb <- BootsPlot(logimpweight_gm_SexDiffb, 5000)
imp_boot_gm_SexDiffc <- BootsPlot(logimpweight_gm_SexDiffc, 5000)
imp_boot_gm_SexDiffab <- BootsPlot(logimpweight_gm_SexDiffab, 5000)
imp_boot_gm_SexDiffac <- BootsPlot(logimpweight_gm_SexDiffac, 5000)
imp_boot_gm_SexDiffbc <- BootsPlot(logimpweight_gm_SexDiffbc, 5000)

## add prior model weights
priorp <- 1/8

p_gm_SexDiffFULL <- logmarg_gm_SexDiffFULL + log(priorp)
p_gm_NoSexDiff <- logmarg_gm_NoSexDiff + log(priorp)
p_gm_SexDiffa <- logmarg_gm_SexDiffa + log(priorp)
p_gm_SexDiffb <- logmarg_gm_SexDiffb + log(priorp)
p_gm_SexDiffc <- logmarg_gm_SexDiffc + log(priorp)
p_gm_SexDiffab <- logmarg_gm_SexDiffab + log(priorp)
p_gm_SexDiffac <- logmarg_gm_SexDiffac + log(priorp)
p_gm_SexDiffbc <- logmarg_gm_SexDiffbc + log(priorp)


p <- c(p_gm_SexDiffFULL,
       p_gm_NoSexDiff,
       p_gm_SexDiffa,
       p_gm_SexDiffb,
       p_gm_SexDiffc,
       p_gm_SexDiffab,
       p_gm_SexDiffac,
       p_gm_SexDiffbc)

pd <- log_sum_exp_marg(p, mn = FALSE)

## normalise
p <- p - pd
p <- exp(p)
p

## plot marginal likelihoods
mods <- list(
  GMSexDiffFULL = imp_boot_gm_SexDiffFULL,
  GMNoSexDiff = imp_boot_gm_NoSexDiff,
  GMSexDiffa = imp_boot_gm_SexDiffa,
  GMSexDiffb = imp_boot_gm_SexDiffb,
  GMSexDiffc = imp_boot_gm_SexDiffc,
  GMSexDiffab = imp_boot_gm_SexDiffab,
  GMSexDiffac = imp_boot_gm_SexDiffac,
  GMSexDiffbc = imp_boot_gm_SexDiffbc
)

marg <- MargLike.plot(mods)
marg
## which models within log(20) of best
logmarg <- map_dbl(mods, "logmarg")
bestind <- which(logmarg == max(logmarg))
logmargLCI <- mods[[bestind]]$LCI
logmarg <- map_dbl(mods, "UCI")
logmarg <- logmarg[map_lgl(logmarg, ~ . >= logmargLCI - log(20))]
logmarg

names <- names(logmarg)
names(logmarg) <- NULL
logmarg <- unlist(logmarg) 

logmarg_table <- as.data.frame(matrix(NA, nrow = 4, ncol = 3))
logmarg_table[, 1] <- names
logmarg_table[, 2] <- logmarg

pb <- sort(p, decreasing = TRUE)
sum(pb)

pb <- pb[1:4]
logmarg_table <- arrange(logmarg_table, desc(V2) )

logmarg_table[, 3] <- pb
colnames(logmarg_table) <- c("Model", "Logmarg", "Prob")
logmarg_table

saveRDS(logmarg_table, file = "ISoutputs/logmarg_table_GM.rds")

#RJout <- readRDS("RJoutputs/posteriorModelProbs.rds")
#RJout

## turn graphics device off
dev.off()

samples <- t(rbind(logimpweight_gm_NoSexDiff, logimpweight_gm_SexDiffa, logimpweight_gm_SexDiffb, logimpweight_gm_SexDiffc,
                   logimpweight_gm_SexDiffab, logimpweight_gm_SexDiffac, logimpweight_gm_SexDiffbc, logimpweight_gm_SexDiffFULL))

## set up bootsrap function
bootsIS <- function(samples, nboot, nmodel){
  lists <- list()
  priorp <- 1/nmodel
  
  for (i in 1:nboot) {
    samp_ints <- sample(1:10000, 10000, replace = TRUE)
    temp_post <- samples[samp_ints, ]
    logmarg <- apply(temp_post, 2, log_sum_exp_marg)
    logmarg <- logmarg + log(priorp)
    pd <- log_sum_exp_marg(logmarg, mn = FALSE)
    p <- logmarg - pd
    p <- exp(p)
    lists[[i]] <- p
    
  }
  t <- as.data.frame(do.call(rbind, lists))
  return(t)
}

## run function
out <- bootsIS(samples, 2000, 8)

## create CIs
IS_ModProbs <- out %>%
  apply(2, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(model = colnames(out))

## adjust model names
IS_ModProbs$model <- gsub("logimpweight_gm_", "", IS_ModProbs$model) 
#IS_ModProbs$model[1] <- "SilerNoSexDiff"

IS_ModProbs <- arrange(IS_ModProbs, desc(Median))

## save
saveRDS(IS_ModProbs, "ISoutputs/logmarg_tableFULL.rds")

