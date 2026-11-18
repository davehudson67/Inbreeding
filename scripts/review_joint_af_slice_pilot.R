# ============================================================================
# Review completed AF_slice chain 5 against previous blocked pilot chain 4
# and the original three 16,000-iteration random-walk chains.
#
# READ ONLY. Original chain files and scientific model are never changed.
# One short AF_slice chain cannot establish convergence or unbiased mixing.
# coda spectral ESS and ESS per hour are DESCRIPTIVE, with different burn-in,
# thinning and chain starts; treat only large differences as useful clues.
# ============================================================================
options(stringsAsFactors=FALSE)
fail <- function(...)stop(paste0(...),call.=FALSE)
stopifnot(requireNamespace("coda",quietly=TRUE))
root <- "outputs/joint_sei_survival/full_joint_20261117"
paths <- file.path(root,sprintf("chain_%02d",1:5),"posterior_samples.rds")
if(!all(file.exists(paths)))
  fail("Missing chain files: ",paste(paths[!file.exists(paths)],collapse=", "))
runs <- lapply(paths,readRDS)
tag <- c("RW1","RW2","RW3","RW_block4","AF_slice5")
expected <- c("default","default","default","blocked","af_slice")
for(j in 1:5){
  z <- runs[[j]]
  scheme_matches <- if(j<=3L)
    (is.null(z$sampler_scheme) ||
       identical(z$sampler_scheme,"default")) else
       identical(z$sampler_scheme,expected[j])
  if(as.integer(z$chain)!=j || as.integer(z$N)!=1603L ||
     !identical(z$status,"UNVALIDATED_PROVISIONAL")||
     !scheme_matches)
    fail("Chain ",j," unexpected cohort, sampler or status")
  path <- dirname(paths[j])
  exit <- file.path(path,"process_exit_status.txt")
  marker <- file.path(path,"chain_finished.txt")
  if(!file.exists(exit)||!file.exists(marker)||
     !any(trimws(readLines(exit,warn=FALSE))=="EXIT_STATUS=0")||
     !any(grepl("^status UNVALIDATED_PROVISIONAL$",
                readLines(marker,warn=FALSE))))
    fail("Chain ",j," lacks successful process completion marker")
  if(!is.matrix(z$samples) ||
     nrow(z$samples)!=as.integer(z$saved_draws) ||
     any(!is.finite(z$samples)))
    fail("Nonfinite, absent or mismatched samples for chain ",j)
}
hashkey <- function(z,k)
  paste(unname(as.character(z[[k]])),collapse=":")
for(k in c("source_md5","input_md5"))
  if(length(unique(vapply(runs,hashkey,character(1),k=k)))!=1L)
    fail("Different full-data provenance across chains: ",k)
science <- c("R/joint_sei_survival_hierarchy_nimble.R",
             "R/joint_sei_survival_forward_nimble.R")
sci <- vapply(runs,function(z){
  h <- z$model_code_md5
  if(!all(science%in%names(h)))
    fail("Missing scientific model hashes")
  paste(unname(h[science]),collapse=":")
},character(1))
if(length(unique(sci))!=1L)
  fail("Scientific likelihood or priors changed; cannot compare samplers")
m <- lapply(runs,function(z)as.matrix(z$samples))
if(!all(vapply(m,function(x)identical(colnames(x),colnames(m[[1]])),
               logical(1))))
  fail("Different monitor columns across samplers")
pars <- c("h4","f4","bGM","sigmaD","load31","load32",
          "sigmaSE","sigmaEI","sigmaCapture",
          "betaSE[4]","betaEI[4]","betaD[2]",
          "betaD[3]","betaD[4]","alphaP","rho")
if(!all(pars%in%colnames(m[[1]])))
  fail("Required parameters not monitored: ",
       paste(setdiff(pars,colnames(m[[1]])),collapse=", "))
derive <- function(x){
  h <- function(age)
    x[,"h4"]*((1-x[,"f4"])*exp(x[,"bGM"]*(age-4))+
                  x[,"f4"])
  cbind(
    h4=x[,"h4"],f4=x[,"f4"],bGM=x[,"bGM"],
    sigmaD=x[,"sigmaD"],
    load31=x[,"load31"],load32=x[,"load32"],
    baseline_hazard_age4=h(4),
    baseline_hazard_age8=h(8),
    baseline_hazard_age12=h(12),
    combined_mortality_frailty_SD=sqrt(x[,"load31"]^2+
        x[,"load32"]^2+x[,"sigmaD"]^2),
    mortality_HR_exposed=exp(x[,"betaD[3]"]),
    mortality_HR_infectious=exp(x[,"betaD[4]"]),
    acquisition_OR_hom=exp(x[,"betaSE[4]"]),
    progression_OR_hom=exp(x[,"betaEI[4]"]),
    mortality_HR_hom=exp(x[,"betaD[2]"]),
    alphaP=x[,"alphaP"],rho=x[,"rho"]
  )
}
D <- lapply(m,derive)
nm <- colnames(D[[1]])
summarise_one <- function(x,index) {
  z <- runs[[index]]
  n <- nrow(x)
  mid <- n%/%2L
  ess <- suppressWarnings(coda::effectiveSize(coda::mcmc(x)))
  med <- apply(x,2,stats::median)
  qlo <- apply(x,2,stats::quantile,probs=.025)
  qhi <- apply(x,2,stats::quantile,probs=.975)
  h1 <- apply(x[seq_len(mid),,drop=FALSE],2,stats::median)
  h2 <- apply(x[seq.int(mid+1L,n),,drop=FALSE],2,stats::median)
  sdv <- apply(x,2,stats::sd)
  runtime <- as.numeric(z$elapsed_sampling_seconds)
  data.frame(
    chain=index,sampler=tag[index],
    parameter=nm,
    draws=n,seconds_per_iteration=runtime/as.integer(z$iterations),
    median=unname(med),
    q025=unname(qlo),q975=unname(qhi),
    half1=unname(h1),half2=unname(h2),
    half_drift_over_sd=unname(abs(h2-h1)/pmax(sdv,1e-9)),
    ESS_spectral=unname(as.numeric(ess[nm])),
    ESS_per_sampling_hour=unname(as.numeric(ess[nm]))/(runtime/3600),
    stringsAsFactors=FALSE
  )
}
tab <- do.call(rbind,lapply(1:5,function(j)
  summarise_one(D[[j]],j)))
cat("\n=== AF_SLICE SINGLE-MODEL 1603-BADGER PILOT REVIEW ===\n")
cat("All 5 chains passed exit, input and scientific source checks: PASS\n")
cat("IMPORTANT: three original full chains failed Rhat/ESS criteria.\n")
cat("One AF pilot and one RW_block pilot cannot demonstrate convergence.\n\n")
meta <- do.call(rbind,lapply(1:5,function(j) {
  x<-runs[[j]]
  data.frame(chain=j,sampler=tag[j],niter=as.integer(x$iterations),
             burnin=as.integer(x$burnin),thin=as.integer(x$thin),
             draws=nrow(x$samples),
             sampling_hours=as.numeric(x$elapsed_sampling_seconds)/3600,
             seconds_per_sweep=as.numeric(x$elapsed_sampling_seconds)/
                               as.integer(x$iterations),
             row.names=NULL)
}))
print(meta,row.names=FALSE,digits=4)

focus <- c("h4","f4","bGM","sigmaD","load31","load32",
           "baseline_hazard_age8","baseline_hazard_age12",
           "combined_mortality_frailty_SD",
           "acquisition_OR_hom","progression_OR_hom","mortality_HR_hom")
cat("\n=== AF PILOT: DRAWS, WITHIN-CHAIN DRIFT AND ESS ===\n")
print(tab[tab$chain==5 & tab$parameter%in%focus,
          c("parameter","median","q025","q975","half1","half2",
            "half_drift_over_sd","ESS_spectral","ESS_per_sampling_hour")],
      row.names=FALSE,digits=4)
cat("\n=== AF VS BLOCKED PILOT, DESCRIPTIVE EFFICIENCY ===\n")
print(tab[tab$chain%in%c(4,5) & tab$parameter%in%
            c("h4","f4","sigmaD","load31","load32",
              "baseline_hazard_age8",
              "combined_mortality_frailty_SD"),
          c("sampler","parameter","ESS_spectral",
            "ESS_per_sampling_hour","half_drift_over_sd")],
      row.names=FALSE,digits=4)
orig <- do.call(rbind,D[1:3])
af <- D[[5]]
cat("\n=== AF VS ORIGINAL THREE-CHAIN MEDIANS ===\n")
co <- data.frame(parameter=focus,
  original_three_chain_median=vapply(focus,function(k)
      median(orig[,k]),numeric(1)),
  AF_pilot_median=vapply(focus,function(k)
      median(af[,k]),numeric(1)),
  original_q025=vapply(focus,function(k)
      unname(quantile(orig[,k],.025)),numeric(1)),
  original_q975=vapply(focus,function(k)
      unname(quantile(orig[,k],.975)),numeric(1)))
co$AF_median_within_original_interval <-
  co$AF_pilot_median>=co$original_q025 &
  co$AF_pilot_median<=co$original_q975
print(co,row.names=FALSE,digits=4)

# Immutable timestamped outputs; never write to chain dirs.
out <- file.path(root,paste0("AF_slice_pilot_05_review_",
  format(Sys.time(),"%Y%m%dT%H%M%S")))
if(dir.exists(out))fail("Review path exists: ",out)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(meta,file.path(out,"sampler_runtime.csv"),row.names=FALSE)
write.csv(tab,file.path(out,"all_five_chain_descriptive.csv"),
          row.names=FALSE)
write.csv(co,file.path(out,"AF_vs_original_medians.csv"),
          row.names=FALSE)
writeLines(c(
  "AF_slice pilot-5 review; scientific model and data provenance matched.",
  "Single-chain spectral ESS and ESS/hour are descriptive and uncertain.",
  "Pilot 4 differs in burn-in/thinning and should not be treated as a controlled sampler comparison.",
  "AF_slice pilot has only 500 draws; three independent AF_slice chains needed for Rhat.",
  "Original 3 chains did not converge on mortality hazards and frailty.",
  "Single cause-independent death hazard includes 202 RTA and 11 management deaths.",
  "Biological genetic associations remain provisional; no causal effect established."
),file.path(out,"read_me_limitations.txt"))
cat("\nSaved read-only AF pilot review:",out,"\n")
cat("These figures are not a convergence verdict or validated effect estimates.\n")
