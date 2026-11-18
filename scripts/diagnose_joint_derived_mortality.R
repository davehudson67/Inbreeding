# =============================================================================
# Three-chain derived mortality diagnostics for ONE joint SEI-survival model.
# READ ONLY. Does not refit model or change samples/prior/likelihood.
#
# Strong h4/f4/sigmaD posterior correlations mean marginal parameter Rhat
# alone does not indicate whether an age-specific predicted hazard mixes well.
# Compute age-specific baseline hazard and total frailty SD per draw.
# Compare original THREE chains for Rhat/ESS and half-chain drift.
# The short blocked pilot is used solely as a descriptive visual comparator.
# =============================================================================
options(stringsAsFactors=FALSE)
stopifnot(requireNamespace("posterior",quietly=TRUE))
fail <- function(...)stop(paste0(...),call.=FALSE)
root <- "outputs/joint_sei_survival/full_joint_20261117"
p <- file.path(root,sprintf("chain_%02d",1:3),"posterior_samples.rds")
if(!all(file.exists(p)))fail("Original full chains missing")
chains <- lapply(p,readRDS)
if(!all(vapply(chains,function(z)z$N==1603,logical(1))))
  fail("Wrong cohort in source chains")
if(length(unique(vapply(chains,function(z)
     paste(unname(as.character(z$input_md5)),collapse=":"),
     character(1))))!=1L)
  fail("Different data across source chains")
if(length(unique(vapply(chains,function(z)
     paste(unname(as.character(z$model_code_md5)),collapse=":"),
     character(1))))!=1L)
  fail("Different model code across original chains")
M <- lapply(chains,function(z)as.matrix(z$samples))
if(!all(vapply(M,function(m)identical(colnames(m),colnames(M[[1]])),
               logical(1))))
  fail("Different monitored parameters")
need <- c("h4","f4","bGM","sigmaD","load31","load32",
          "betaD[1]","betaD[2]","betaD[3]","betaD[4]",
          "betaSE[4]","betaEI[4]")
if(any(!need%in%colnames(M[[1]])))
  fail("Missing needed parameters: ",
       paste(setdiff(need,colnames(M[[1]])),collapse=", "))
if(any(vapply(M,function(x)any(!is.finite(x)),logical(1))))
  fail("Nonfinite posterior samples")
if(length(unique(vapply(M,nrow,integer(1))))!=1L)
  fail("Unequal posterior chain lengths")

derive <- function(x){
  h <- function(age)
    x[,"h4"]*((1-x[,"f4"])*exp(x[,"bGM"]*(age-4))+
                      x[,"f4"])
  cbind(
    baseline_hazard_age_2=h(2),
    baseline_hazard_age_4=h(4),
    baseline_hazard_age_6=h(6),
    baseline_hazard_age_8=h(8),
    baseline_hazard_age_10=h(10),
    baseline_hazard_age_12=h(12),
    baseline_hazard_ratio_8_to_2=h(8)/h(2),
    baseline_makeham_component=x[,"h4"]*x[,"f4"],
    baseline_gompertz_at_age_4=x[,"h4"]*(1-x[,"f4"]),
    frailty_SD_death_total=sqrt(
      x[,"load31"]^2+x[,"load32"]^2+x[,"sigmaD"]^2),
    mortality_HR_homozygosity_1SD=exp(x[,"betaD[2]"]),
    mortality_HR_male=exp(x[,"betaD[1]"]),
    mortality_HR_exposed_vs_S=exp(x[,"betaD[3]"]),
    mortality_HR_infectious_vs_S=exp(x[,"betaD[4]"]),
    acquisition_OR_homozygosity_1SD=exp(x[,"betaSE[4]"]),
    progression_OR_homozygosity_1SD=exp(x[,"betaEI[4]"])
  )
}
D <- lapply(M,derive)
nm <- colnames(D[[1]])
arr <- array(NA_real_,
  dim=c(nrow(D[[1]]),3L,length(nm)),
  dimnames=list(NULL,paste0("chain",1:3),nm))
for(j in 1:3)arr[,j,] <- D[[j]]
ps <- as.data.frame(posterior::summarise_draws(
  posterior::as_draws_array(arr)))
idx <- match(nm,as.character(ps$variable))
if(anyNA(idx))fail("Missing derived draw summary")
all <- do.call(rbind,D)
sumtab <- data.frame(
 quantity=nm,
 median=unname(apply(all,2,median)),
 lower95=unname(apply(all,2,quantile,.025)),
 upper95=unname(apply(all,2,quantile,.975)),
 rhat=ps$rhat[idx],
 ess_bulk=ps$ess_bulk[idx],
 ess_tail=ps$ess_tail[idx],
 row.names=NULL)
bychain <- do.call(rbind,lapply(1:3,function(j){
  x <- D[[j]]
  half <- floor(nrow(x)/2)
  data.frame(
    chain=j,quantity=nm,
    median=unname(apply(x,2,median)),
    first_half=unname(apply(x[seq_len(half),,drop=FALSE],2,median)),
    second_half=unname(apply(
      x[seq.int(half+1,nrow(x)),,drop=FALSE],2,median)),
    row.names=NULL)
}))
cat("\n=== ONE JOINT MODEL: DERIVED MORTALITY CONVERGENCE ===\n")
cat("Original 3 chains; 2,000 retained draws each; n=1603.\n")
print(sumtab[,c("quantity","median","lower95","upper95",
               "rhat","ess_bulk","ess_tail")],
      row.names=FALSE,digits=4)
cat("\n=== SELECTED MORTALITY QUANTITIES BY CHAIN/HALF ===\n")
key <- c("baseline_hazard_age_4","baseline_hazard_age_8",
         "baseline_hazard_age_12","frailty_SD_death_total",
         "mortality_HR_homozygosity_1SD",
         "mortality_HR_exposed_vs_S",
         "mortality_HR_infectious_vs_S")
print(bychain[bychain$quantity%in%key,],
      row.names=FALSE,digits=4)

pp <- file.path(root,"chain_04","posterior_samples.rds")
pilotcomp <- NULL
if(file.exists(pp)){
  pilot <- readRDS(pp)
  if(identical(pilot$sampler_scheme,"blocked")&&pilot$N==1603L){
    x <- derive(as.matrix(pilot$samples))
    pilotcomp <- data.frame(
      quantity=colnames(x),
      pilot_median=unname(apply(x,2,median)),
      pilot_q025=unname(apply(x,2,quantile,.025)),
      pilot_q975=unname(apply(x,2,quantile,.975)),
      old_3chain_median=sumtab$median,
      row.names=NULL)
    cat("\n=== SHORT BLOCKED PILOT 4: DESCRIPTIVE ONLY ===\n")
    print(pilotcomp[pilotcomp$quantity%in%key,],
          row.names=FALSE,digits=4)
    cat("WARNING: One 600-draw pilot cannot establish between-chain",
        " convergence or improvement over original 3 chains.\n")
  } else warning("Chain 4 is not a valid blocked pilot; ignoring")
}

name <- paste0("derived_mortality_review_",
               format(Sys.time(),"%Y%m%dT%H%M%S"))
out <- file.path(root,name)
if(dir.exists(out))fail("Existing review folder: ",out)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(sumtab,file.path(out,"derived_mortality_convergence.csv"),
          row.names=FALSE)
write.csv(bychain,file.path(out,"derived_mortality_by_chain_half.csv"),
          row.names=FALSE)
if(!is.null(pilotcomp))
  write.csv(pilotcomp,file.path(out,"short_blocked_pilot_comparison.csv"),
            row.names=FALSE)
cat("\nSaved read-only derived review:",out,"\n")
cat("These are baseline hazards, not cause-specific mortality",
    " or posterior predictive validation.\n")
