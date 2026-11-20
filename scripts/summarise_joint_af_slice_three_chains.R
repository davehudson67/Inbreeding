# AF_slice-specific post-burn-in checks: sampler identical on all
# three chains and separate outputs from the original default sampler.
# Important: rank-normalized Rhat/ESS must be reassessed, independent
# of the short AF pilot-5 or blocked pilot-4.
# ============================================================================
# Three independent AF_slice chains (6-8) of SAME joint SEI-survival model.
#
# READ ONLY: does not rerun or modify MCMC, private inputs, or chain files.
# Requires three successfully finished, provenance-matched RDS chains.
# Writes an immutable timestamped diagnostics folder (all outputs gitignored).
# Modern rank-normalized Rhat / bulk and tail ESS via 'posterior' if installed;
# otherwise coda's classic Gelman-Rubin and spectral ESS, explicitly labelled.
#
# Important: this fit uses a proxy for trapping effort, original fixed assay
# accuracies and recorded deaths that may include management interventions.
# No unqualified natural-mortality interpretation is supported yet.
# ============================================================================
options(stringsAsFactors=FALSE)
root <- Sys.getenv("JOINT_RUN_DIR",
  "outputs/joint_sei_survival/full_joint_20261117")
chain_ids <- 6:8
sampler_expected <- "af_slice"
fail <- function(...) stop(paste0(...),call.=FALSE)
f <- file.path(root,sprintf("chain_%02d",chain_ids),
               "posterior_samples.rds")
done <- file.path(root,sprintf("chain_%02d",chain_ids),
                  "chain_finished.txt")
status <- file.path(root,sprintf("chain_%02d",chain_ids),
                    "process_exit_status.txt")
logs <- file.path(root,sprintf("chain_%02d",chain_ids),
                  "chain_console.log")
cat("=== SAME JOINT SEI-SURVIVAL: THREE INDEPENDENT AF_SLICE CHAINS 6-8 ===\n")
if(!all(file.exists(f))){
  for(i in which(!file.exists(f))) {
    cat("MISSING posterior chain",i,":",f[i],"\n")
    if(file.exists(logs[i])){
      cat("Last log lines:\n")
      cat(paste(tail(readLines(logs[i],warn=FALSE),12L),
                collapse="\n"),"\n")
    }
  }
  fail("All three completed chain RDS files are required")
}
if(!all(file.exists(done))||!all(file.exists(status)))
  fail("One or more chains lacks a completed-run or process-exit marker")
for(i in 1:3){
  ex <- readLines(status[i],warn=FALSE)
  if(!any(grepl("^EXIT_STATUS=0$",ex)))
    fail("Chain ",i," did not record successful exit: ",
         paste(ex,collapse="; "))
  finish <- readLines(done[i],warn=FALSE)
  if(!any(grepl("^status UNVALIDATED_PROVISIONAL$",finish)))
    fail("Chain ",i," lacks expected completed provisional status")
}
chains <- lapply(f,readRDS)
for(i in 1:3){
  z <- chains[[i]]
  needed <- c("chain","seed","N","iterations","burnin","thin",
              "saved_draws","samples","source_md5","input_md5",
              "model_code_md5","status")
  if(!all(needed%in%names(z)))
    fail("Chain ",i," lacks required provenance or samples: ",
         paste(setdiff(needed,names(z)),collapse=","))
  if(as.integer(z$chain)!=chain_ids[i]||as.integer(z$N)!=1603L||
     !identical(z$sampler_scheme,sampler_expected)||
     !identical(z$status,"UNVALIDATED_PROVISIONAL"))
    fail("Chain ",i," has wrong ID, cohort or provisional status")
  if(!is.matrix(z$samples)||!is.numeric(z$samples)||
     anyNA(z$samples)||any(!is.finite(z$samples))||
     nrow(z$samples)<100L||
     nrow(z$samples)!=as.integer(z$saved_draws))
    fail("Chain ",i," lacks valid finite posterior sample matrix")
}
keys <- c("N","iterations","burnin","thin")
for(k in keys)
  if(length(unique(vapply(chains,function(z)as.integer(z[[k]]),
                          integer(1))))!=1L)
    fail("Chains disagree on ",k,"; do NOT pool")
if(anyDuplicated(vapply(chains,function(z)as.integer(z$seed),
                        integer(1))))
  fail("Identical chain RNG seeds; independence not established")
check_hash <- function(name){
  hashes <- lapply(chains,function(z){
    h <- z[[name]]
    paste(unname(as.character(h)),collapse=":")
  })
  if(length(unique(unlist(hashes)))!=1L)
    fail("Chains disagree on ",name,"; do NOT pool")
}
check_hash("source_md5")
check_hash("input_md5")
check_hash("model_code_md5")
M <- lapply(chains,function(z)z$samples)
if(!all(vapply(M,function(m)
  identical(colnames(m),colnames(M[[1]])),logical(1))))
  fail("Monitored parameters differ between chains")
if(!all(vapply(M,nrow,integer(1))==nrow(M[[1]])))
  fail("Unequal retained draws, cannot construct parallel draw arrays")
if(anyDuplicated(colnames(M[[1]]))||
   is.null(colnames(M[[1]])))
  fail("Missing or duplicated posterior parameter labels")

cat("Chain IDs:",paste(chain_ids,collapse=", "),"; sampler:",sampler_expected,"\n")
cat("Genotyped badgers:",chains[[1]]$N,"\n")
cat("Iterations / burnin / thin:",
    chains[[1]]$iterations,chains[[1]]$burnin,
    chains[[1]]$thin,"\n")
cat("Retained draws by chain:",
    paste(vapply(M,nrow,integer(1)),collapse=", "),"\n")
cat("Seeds:",paste(vapply(chains,function(z)z$seed,integer(1)),
                    collapse=", "),"\n")
cat("Input and code hashes: identical across all 3 AF_slice chains: PASS\n")

required <- c("betaSE[4]","betaEI[4]","betaD[2]",
 "betaSE[3]","betaEI[3]","betaD[1]",
 "sigmaSE","sigmaEI","sigmaD",
 "load21","load31","load32","sigmaCapture",
 "h4","f4","bGM","entryRho","alphaP","rho",
 "sigmaYearSE","sigmaYearEI")
all_names <- colnames(M[[1]])
missing <- setdiff(required,all_names)
if(length(missing))
  fail("Missing required scientific parameters in saved MCMC: ",
       paste(missing,collapse=", "))
vary <- vapply(seq_along(all_names),function(j){
  v <- unlist(lapply(M,function(x)x[,j]))
  sd(v)>1e-12
},logical(1))
if(!any(vary))fail("All monitored parameters appear constant")
cat("Excluded constant reference terms:",
    paste(all_names[!vary],collapse=", "),"\n")
M <- lapply(M,function(x)x[,vary,drop=FALSE])
params <- colnames(M[[1]])
cat("Nonconstant monitored parameters:",length(params),"\n")

# Derived quantities are calculated at EACH posterior draw.
# The 3 biological random effects have triangular factor loadings
# SE = sigmaSE * z1
# EI = load21*z1 + sigmaEI*z2
#  D = load31*z1 + load32*z2 + sigmaD*z3
# Therefore SD(EI), SD(D) include the off-diagonal loadings.
derived <- lapply(M,function(x){
  se <- x[,"sigmaSE"]; ei <- x[,"sigmaEI"]
  dd <- x[,"sigmaD"]
  a21 <- x[,"load21"];a31 <- x[,"load31"]
  a32 <- x[,"load32"]
  sEI <- sqrt(a21^2+ei^2)
  sD <- sqrt(a31^2+a32^2+dd^2)
  # Same anchored Gompertz-Makeham hazard as hierarchy:
  # aGM = h4*(1-f4)*exp(-4*bGM), cGM = h4*f4.
  # These baseline hazards exclude sex, infection, genotype and frailty.
  hage <- function(age) x[,"h4"] * (
    (1-x[,"f4"])*exp(x[,"bGM"]*(age-4)) + x[,"f4"])
  r12 <- a21/sEI
  r13 <- a31/sD
  r23 <- (a21*a31+ei*a32)/(sEI*sD)
  cbind(
   acquisition_OR_hom_1SD=exp(x[,"betaSE[4]"]),
   progression_OR_hom_1SD=exp(x[,"betaEI[4]"]),
   mortality_HR_hom_1SD=exp(x[,"betaD[2]"]),
   acquisition_OR_male=exp(x[,"betaSE[3]"]),
   progression_OR_male=exp(x[,"betaEI[3]"]),
   mortality_HR_male=exp(x[,"betaD[1]"]),
   mortality_HR_exposed_vs_S=exp(x[,"betaD[3]"]),
   mortality_HR_infectious_vs_S=exp(x[,"betaD[4]"]),
   SD_frailty_SE=se,SD_frailty_EI=sEI,SD_frailty_D=sD,
   SD_frailty_capture=x[,"sigmaCapture"],
   baseline_hazard_age_2=hage(2),
   baseline_hazard_age_4=hage(4),
   baseline_hazard_age_6=hage(6),
   baseline_hazard_age_8=hage(8),
   baseline_hazard_age_10=hage(10),
   baseline_hazard_age_12=hage(12),
   baseline_hazard_ratio_8_to_2=hage(8)/hage(2),
   baseline_makeham_component=x[,"h4"]*x[,"f4"],
   baseline_gompertz_at_age_4=x[,"h4"]*(1-x[,"f4"]),
   cor_frailty_SE_EI=r12,
   cor_frailty_SE_D=r13,
   cor_frailty_EI_D=r23
  )
})

make_summary <- function(chains_m, label) {
  nm <- colnames(chains_m[[1]])
  if(any(vapply(chains_m,function(m) !identical(colnames(m),nm),
                logical(1))))
    fail("Derived chain columns differ for ",label)
  combined <- do.call(rbind,chains_m)
  stats <- data.frame(
    parameter=nm,
    mean=unname(colMeans(combined)),
    median=unname(apply(combined,2,median)),
    lower95=unname(apply(combined,2,quantile,probs=.025)),
    upper95=unname(apply(combined,2,quantile,probs=.975)),
    prob_positive=unname(colMeans(combined>0)),
    stringsAsFactors=FALSE
  )
  if(requireNamespace("posterior",quietly=TRUE)){
    arr <- array(NA_real_,
       dim=c(nrow(chains_m[[1]]),length(chains_m),length(nm)),
       dimnames=list(NULL,paste0("chain",seq_along(chains_m)),nm))
    for(j in seq_along(chains_m))arr[,j,] <- chains_m[[j]]
    post <- as.data.frame(
      posterior::summarise_draws(posterior::as_draws_array(arr)))
    match_idx <- match(nm,as.character(post$variable))
    if(anyNA(match_idx)||
       !all(c("rhat","ess_bulk","ess_tail")%in%names(post)))
      fail("Modern posterior diagnostic results lacked expected fields")
    stats$rhat <- as.numeric(post$rhat[match_idx])
    stats$ess_bulk <- as.numeric(post$ess_bulk[match_idx])
    stats$ess_tail <- as.numeric(post$ess_tail[match_idx])
    diag_type <- "rank-normalized split Rhat, bulk/tail ESS (posterior)"
  } else if(requireNamespace("coda",quietly=TRUE)){
    ml <- coda::mcmc.list(lapply(chains_m,coda::mcmc))
    gr <- tryCatch(coda::gelman.diag(
      ml,autoburnin=FALSE,multivariate=FALSE)$psrf[,1],
      error=function(e){
        warning("Classic Rhat could not be computed: ",
                conditionMessage(e))
        setNames(rep(NA_real_,length(nm)),nm)
      })
    ess <- coda::effectiveSize(ml)
    stats$rhat <- unname(as.numeric(gr[nm]))
    stats$ess_bulk <- unname(as.numeric(ess[nm]))
    stats$ess_tail <- NA_real_
    diag_type <- "CLASSIC Gelman-Rubin PSRF + coda spectral ESS (fallback)"
  } else {
    fail("Install posterior or coda to compute MCMC diagnostics")
  }
  attr(stats,"diag_type") <- diag_type
  stats
}
summary <- make_summary(M,"raw")
derived_summary <- make_summary(derived,"derived")
diag_type <- attr(summary,"diag_type")
cat("Diagnostics method:",diag_type,"\n")
if(any(!is.finite(summary$rhat))||
   any(!is.finite(summary$ess_bulk)))
  warning("One or more convergence diagnostics are missing/nonfinite")
cat("Max Rhat:",max(summary$rhat,na.rm=TRUE),"\n")
cat("Min effective sample size:",min(summary$ess_bulk,na.rm=TRUE),"\n")
poor <- summary[
  !is.finite(summary$rhat) | summary$rhat>1.01 |
  !is.finite(summary$ess_bulk) | summary$ess_bulk<400 |
  (!is.na(summary$ess_tail) & summary$ess_tail<400),]
poor <- poor[order(-poor$rhat),]
cat("Parameters failing nominal Rhat<=1.01 and ESS>=400:",
    nrow(poor),"/",nrow(summary),"\n")
if(nrow(poor))
  print(utils::head(
    poor[,c("parameter","median","rhat","ess_bulk","ess_tail")],
    25L),row.names=FALSE,digits=4)

# Immutable review output: retain every earlier diagnostic audit.
runstamp <- format(Sys.time(),"%Y%m%dT%H%M%S")
outdir <- file.path(root,paste0("AF_slice_three_chain_review_",runstamp))
if(dir.exists(outdir))fail("Diagnostic output directory already exists")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
utils::write.csv(summary,
  file.path(outdir,"joint_parameter_posterior.csv"),
  row.names=FALSE)
utils::write.csv(derived_summary,
  file.path(outdir,"joint_derived_effects.csv"),
  row.names=FALSE)
utils::write.csv(poor,
  file.path(outdir,"parameters_needing_convergence_review.csv"),
  row.names=FALSE)
metadata <- data.frame(
 chain=chain_ids,
 seed=vapply(chains,function(z)as.integer(z$seed),integer(1)),
 n=vapply(chains,function(z)as.integer(z$N),integer(1)),
 iterations=vapply(chains,function(z)as.integer(z$iterations),integer(1)),
 burnin=vapply(chains,function(z)as.integer(z$burnin),integer(1)),
 thin=vapply(chains,function(z)as.integer(z$thin),integer(1)),
 draws=vapply(chains,function(z)nrow(z$samples),integer(1)),
 sampling_seconds=vapply(chains,function(z)
    as.numeric(z$elapsed_sampling_seconds),numeric(1)),
 input_md5=vapply(chains,function(z)
    paste(unname(z$input_md5),collapse=":"),character(1)),
 model_md5=vapply(chains,function(z)
    paste(unname(z$model_code_md5),collapse=":"),character(1))
)
utils::write.csv(metadata,file.path(outdir,"chain_metadata.csv"),
                 row.names=FALSE)

# Trace plots: show within-chain persistence and disagreement.
trace_keep <- intersect(c(
  "betaSE[4]","betaEI[4]","betaD[2]",
  "sigmaSE","sigmaEI","sigmaD",
  "load21","load31","load32",
  "h4","f4","bGM","entryRho","rho"),params)
pdf(file.path(outdir,"core_parameter_traces.pdf"),width=11,height=8)
oldpar <- par(mfrow=c(3,2),mar=c(3.5,4.0,2.1,1.3))
colors <- c("#26547C","#EF476F","#06A77D")
for(p in trace_keep){
  yr <- range(unlist(lapply(M,function(m)m[,p])))
  plot(seq_len(nrow(M[[1]])),M[[1]][,p],type="l",
       xlab="Saved iteration",ylab=p,main=p,
       col=colors[1],ylim=yr)
  for(j in 2:3)
    lines(seq_len(nrow(M[[j]])),M[[j]][,p],col=colors[j])
  legend("topright",legend=paste("Chain",chain_ids),
         col=colors,lty=1,bty="n",cex=.7)
}
par(oldpar)
dev.off()

# By-chain medians on biological effects help identify disagreement.
derived_by_chain <- do.call(rbind,lapply(1:3,function(j){
  x <- derived[[j]]
  data.frame(
   chain=chain_ids[j],parameter=colnames(x),
   median=unname(apply(x,2,median)),
   lower95=unname(apply(x,2,quantile,probs=.025)),
   upper95=unname(apply(x,2,quantile,probs=.975)),
   row.names=NULL)
}))
utils::write.csv(derived_by_chain,
  file.path(outdir,"derived_effects_by_chain.csv"),
  row.names=FALSE)

# Retain all three chains' marginal distributions separately to distinguish
# strong between-chain disagreement from persistent within-chain movement.
parameter_by_chain <- do.call(rbind,lapply(1:3,function(j){
  x <- M[[j]]
  data.frame(
    chain=chain_ids[j],parameter=colnames(x),
    median=unname(apply(x,2,median)),
    lower95=unname(apply(x,2,quantile,probs=.025)),
    upper95=unname(apply(x,2,quantile,probs=.975)),
    first_half_median=unname(apply(
      x[seq_len(floor(nrow(x)/2)),,drop=FALSE],2,median)),
    second_half_median=unname(apply(
      x[seq.int(floor(nrow(x)/2)+1L,nrow(x)),,drop=FALSE],2,median)),
    row.names=NULL)
}))
utils::write.csv(parameter_by_chain,
  file.path(outdir,"parameters_by_chain_and_half.csv"),
  row.names=FALSE)
cat("\n=== KEY GENETIC SLOPES BY CHAIN (PROVISIONAL) ===\n")
key_slopes <- c("betaSE[4]","betaEI[4]","betaD[2]")
print(parameter_by_chain[
  parameter_by_chain$parameter%in%key_slopes,
  c("chain","parameter","median","first_half_median",
    "second_half_median")],row.names=FALSE,digits=4)

# Save a short machine-readable verdict, NOT a scientific conclusion.
method_modern <- grepl("rank-normalized",diag_type,fixed=TRUE)
passing <- method_modern && nrow(poor)==0 &&
  all(is.finite(derived_summary$rhat)) &&
  all(derived_summary$rhat<=1.01) &&
  all(is.finite(derived_summary$ess_bulk)) &&
  all(derived_summary$ess_bulk>=400) &&
  all(is.finite(derived_summary$ess_tail)) &&
  all(derived_summary$ess_tail>=400)
verdict <- if (passing) {
  "SCREEN_PASS: numeric mixing thresholds met, scientific validation pending"
} else {
  "REVIEW_REQUIRED: mixing/diagnostics or scientific validation pending"
}
writeLines(c(
  verdict,
  paste("Diagnostic method:",diag_type),
  paste("Number of monitored parameters:",nrow(summary)),
  paste("Number failing nominal thresholds:",nrow(poor)),
  "These checks cannot certify model identifiability, adequacy or causality.",
  "Mortality includes potentially deliberate removal/euthanasia.",
  "Sampling conditioned on genotype and first live encounter.",
  "Capture effort is a proxy; six assay properties were fixed.",
  "Parameter-recovery simulations and predictive checks remain outstanding."),
  file.path(outdir,"AF_slice_convergence_verdict.txt"))

cat("\n=== AF_SLICE: DERIVED BASELINE MORTALITY HAZARD CONVERGENCE ===\n")
hazard_names <- c(
  "baseline_hazard_age_2","baseline_hazard_age_4",
  "baseline_hazard_age_6","baseline_hazard_age_8",
  "baseline_hazard_age_10","baseline_hazard_age_12",
  "baseline_hazard_ratio_8_to_2",
  "baseline_makeham_component","baseline_gompertz_at_age_4",
  "SD_frailty_D",
  "mortality_HR_exposed_vs_S","mortality_HR_infectious_vs_S")
print(derived_summary[
  match(hazard_names,derived_summary$parameter),
  c("parameter","median","lower95","upper95",
    "rhat","ess_bulk","ess_tail")],
  row.names=FALSE,digits=4)
utils::write.csv(derived_summary[
  match(hazard_names,derived_summary$parameter),],
  file.path(outdir,"AF_slice_mortality_hazard_convergence.csv"),
  row.names=FALSE)
cat("\n=== AF_SLICE THREE-CHAIN HOMOZYGOSITY EFFECTS (PROVISIONAL) ===\n")
keep <- c("acquisition_OR_hom_1SD",
          "progression_OR_hom_1SD",
          "mortality_HR_hom_1SD",
          "SD_frailty_SE","SD_frailty_EI","SD_frailty_D",
          "cor_frailty_SE_EI","cor_frailty_SE_D",
          "cor_frailty_EI_D")
print(derived_summary[match(keep,derived_summary$parameter),
  c("parameter","median","lower95","upper95",
    "rhat","ess_bulk","ess_tail")],row.names=FALSE,digits=4)
cat("\nVERDICT:",verdict,"\n")
cat("Saved diagnostics:",outdir,"\n")
cat("Core trace PDF:",file.path(outdir,
       "core_parameter_traces.pdf"),"\n")
cat("Do not interpret biological associations until convergence,",
    " death provenance, effort assumptions, simulation recovery",
    " and predictive checks have been assessed.\n")
