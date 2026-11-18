# =============================================================================
# Review the completed, separately saved full-1603 blocked-sampler chain 4.
# READ ONLY. Does not refit the scientific model, extend runs or overwrite chains.
#
# Compares pilot chain's trajectory with first three default-sampler chains
# for gross parameter drift and computational throughput. ONE PILOT CHAIN
# CANNOT demonstrate between-chain convergence or superior mixing.
# =============================================================================
options(stringsAsFactors=FALSE)
fail <- function(...)stop(paste0(...),call.=FALSE)
root <- "outputs/joint_sei_survival/full_joint_20261117"
pilot_dir <- file.path(root,"chain_04")
f <- file.path(pilot_dir,"posterior_samples.rds")
exitfile <- file.path(pilot_dir,"process_exit_status.txt")
finish <- file.path(pilot_dir,"chain_finished.txt")
logfile <- file.path(pilot_dir,"chain_console.log")
if(!all(file.exists(c(f,exitfile,finish)))){
  cat("Chain 4 output is not complete/verified.\n")
  if(file.exists(logfile)){
    cat("Last 25 log lines:\n")
    cat(paste(tail(readLines(logfile,warn=FALSE),25),
              collapse="\n"),"\n")
  }
  fail("Missing: ",paste(c(f,exitfile,finish)[
    !file.exists(c(f,exitfile,finish))],collapse=", "))
}
if(!any(trimws(readLines(exitfile,warn=FALSE))=="EXIT_STATUS=0"))
  fail("Chain 4 did not have a successful process exit.")
if(!any(grepl("^status UNVALIDATED_PROVISIONAL$",
              readLines(finish,warn=FALSE))))
  fail("Chain 4 lacks the saved-model completion marker.")
pilot <- readRDS(f)
if(!identical(pilot$sampler_scheme,"blocked")||
   as.integer(pilot$chain)!=4L ||
   as.integer(pilot$N)!=1603L ||
   !identical(pilot$status,"UNVALIDATED_PROVISIONAL"))
  fail("Chain 4 is not the expected blocked, provisional 1603-badger run")
m <- as.matrix(pilot$samples)
if(!is.numeric(m)||any(!is.finite(m))||
   nrow(m)!=as.integer(pilot$saved_draws)||
   nrow(m)<100L)
  fail("Pilot sample matrix absent, nonfinite or unexpectedly short")

orig_paths <- file.path(root,sprintf("chain_%02d",1:3),
                        "posterior_samples.rds")
if(!all(file.exists(orig_paths)))
  fail("Original 3 chains unavailable for read-only comparison")
old <- lapply(orig_paths,readRDS)
if(!all(vapply(old,function(z)identical(z$sampler_scheme,"default"),
               logical(1))) &&
   !all(vapply(old,function(z)is.null(z$sampler_scheme),
               logical(1))))
  fail("Original chain sampler provenance differs unexpectedly")
hashkey <- function(x,slot)
  paste(unname(as.character(x[[slot]])),collapse=":")
for(slot in c("input_md5","source_md5")){
  hashes <- vapply(c(old,list(pilot)),hashkey,character(1),slot=slot)
  if(length(unique(hashes))!=1L)
    fail("Pilot and previous chains differ on ",slot)
}
# Orchestration script legitimately changed when introducing block samplers.
# Require the hierarchy and HMM itself to be byte-identical in each chain.
scientific_files <- c("R/joint_sei_survival_hierarchy_nimble.R",
                       "R/joint_sei_survival_forward_nimble.R")
modelhash <- function(x) {
  h <- x$model_code_md5
  if(!all(scientific_files%in%names(h)))
    fail("Model hash record lacks hierarchy / likelihood sources")
  paste(unname(h[scientific_files]),collapse=":")
}
sci_hash <- vapply(c(old,list(pilot)),modelhash,character(1))
if(length(unique(sci_hash))!=1L)
  fail("Scientific hierarchy/likelihood code changed; cannot compare")

oldM <- lapply(old,function(z)as.matrix(z$samples))
if(!all(vapply(oldM,function(z)identical(colnames(z),colnames(m)),
               logical(1))))
  fail("Monitored variables differ between original and pilot chains")
param <- c("h4","f4","bGM","sigmaD","load31","load32",
           "sigmaSE","sigmaEI","sigmaCapture",
           "betaSE[4]","betaEI[4]","betaD[2]",
           "rho","alphaP")
if(!all(param%in%colnames(m)))
  fail("Pilot lacks required parameters: ",
       paste(setdiff(param,colnames(m)),collapse=", "))
old_combined <- do.call(rbind,oldM)
n <- nrow(m)
half <- floor(n/2)
pilot_table <- data.frame(
 parameter=param,
 pilot_median=vapply(param,function(k)median(m[,k]),numeric(1)),
 pilot_first_half=vapply(param,function(k)
   median(m[seq_len(half),k]),numeric(1)),
 pilot_second_half=vapply(param,function(k)
   median(m[seq.int(half+1L,n),k]),numeric(1)),
 pilot_q025=vapply(param,function(k)unname(
   quantile(m[,k],probs=.025)),numeric(1)),
 pilot_q975=vapply(param,function(k)unname(
   quantile(m[,k],probs=.975)),numeric(1)),
 original_combined_median=vapply(param,function(k)
   median(old_combined[,k]),numeric(1)),
 original_q025=vapply(param,function(k)
   unname(quantile(old_combined[,k],.025)),numeric(1)),
 original_q975=vapply(param,function(k)
   unname(quantile(old_combined[,k],.975)),numeric(1)),
 row.names=NULL
)
pilot_table$median_outside_original_95interval <- with(
  pilot_table,pilot_median<original_q025|pilot_median>original_q975)

# Qualitative warning signs only: 600 draws from ONE chain are far too few
# to judge whether the adaptive block proposal fixes cross-chain mixing.
pilot_table$half_median_diff <- pilot_table$pilot_second_half -
                                 pilot_table$pilot_first_half
pilot_table$pilot_sd <- vapply(param,function(k)sd(m[,k]),numeric(1))
pilot_table$abs_half_drift_over_sd <-
  abs(pilot_table$half_median_diff)/pmax(pilot_table$pilot_sd,1e-9)

old_speed <- vapply(old,function(z)
  as.numeric(z$elapsed_sampling_seconds)/as.numeric(z$iterations),
  numeric(1))
speed <- as.numeric(pilot$elapsed_sampling_seconds)/
         as.numeric(pilot$iterations)

cat("\n=== BLOCKED SAMPLER: FULL COHORT CHAIN 4 REVIEW ===\n")
cat("Original 3 chains and pilot share data and scientific code: PASS\n")
cat("Pilot: N =",pilot$N,
    "iterations =",pilot$iterations,
    "burnin =",pilot$burnin,"thin =",pilot$thin,
    "retained =",n,"\n")
cat("Sampling hours:",
    signif(as.numeric(pilot$elapsed_sampling_seconds)/3600,4),"\n")
cat("Seconds/sweep (pilot):",signif(speed,4),"\n")
cat("Seconds/sweep (3 original chains):",
    paste(signif(old_speed,4),collapse=", "),"\n")
cat("\n=== PROVISIONAL: PILOT MORTALITY CHAIN TRAJECTORY ===\n")
keep <- pilot_table$parameter%in%c("h4","f4","bGM",
                                   "sigmaD","load31","load32")
print(pilot_table[keep,c(
  "parameter","pilot_first_half","pilot_second_half",
  "pilot_median","original_combined_median",
  "abs_half_drift_over_sd",
  "median_outside_original_95interval")],
  row.names=FALSE,digits=4)
cat("\n=== PROVISIONAL: PILOT GENETIC SLOPES ===\n")
print(pilot_table[grepl("^beta",pilot_table$parameter),
  c("parameter","pilot_median","original_combined_median",
    "pilot_first_half","pilot_second_half")],
  row.names=FALSE,digits=4)

# Purely descriptive single-chain ESS: no 3-chain Rhat or cross-chain
# mixing improvement is computable from this one pilot.
if(requireNamespace("coda",quietly=TRUE)){
  ess <- coda::effectiveSize(coda::mcmc(m[,param,drop=FALSE]))
  cat("\nSingle-chain pilot spectral ESS (descriptive only):\n")
  print(data.frame(parameter=param,
                   ESS=round(as.numeric(ess[param]),1)),
        row.names=FALSE)
}

out <- file.path(root,
  paste0("blocked_pilot_04_review_",
         format(Sys.time(),"%Y%m%dT%H%M%S")))
if(dir.exists(out))fail("Output review folder already exists: ",out)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(pilot_table,file.path(out,"blocked_pilot_vs_original.csv"),
          row.names=FALSE)
writeLines(c(
  "Pilot full-cohort blocked sampler completed; scientific sources match.",
  paste("Pilot chain 4, retained draws",n),
  paste("Sampling seconds per sweep",signif(speed,6)),
  paste("Original chains mean seconds per sweep",
        signif(mean(old_speed),6)),
  "Do NOT infer convergence improvement from one short chain.",
  "Do NOT pool this pilot with original long chains as final posterior.",
  "All-cause mortality includes 202 RTA, 8 euthanised, 3 killed.",
  "Acquisition/progression frailty scales approach upper prior limits."
  ),file.path(out,"interpretation_warning.txt"))
cat("\nSaved non-destructive pilot review:",out,"\n")
cat("One blocked pilot can test execution and trajectory only;",
    " three independent blocked runs are required to compare Rhat/ESS.\n")
