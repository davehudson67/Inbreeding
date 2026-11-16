# ============================================================================
# Stage 80: DIRECT-DIAGNOSTIC, GENETIC SEI JOINT-LIKELIHOOD SMOKE TEST
#
# S -> E and E -> I probabilities include separate homozygosity effects.
# Latent quarterly states are integrated out by a 3-state HMM forward pass.
# Crucially this FITS RAW BINARY DIAGNOSTICS, not inferred onset times.
#
# INPUT: Stage76 conservative 6-test candidate; only TRAP/NET/DART GUN
# confirmed live encounters, excluding the 317 original-only occasions.
# Keep repeated within-quarter captures and stack their test likelihoods.
# Invalid/contradictory results were NA in Stage76 and stay uninformative.
#
# FULL-ELIGIBLE COMPUTATIONAL ROBUSTNESS CHECK: constant quarterly hazards, fixed Se/Sp
# and exposed-state multipliers from original fitted SEI; no seasonal,
# calendar, age, sex, social-group, survival, trapping or observation model.
# Regularized MAP is NOT a fully Bayesian fit nor a final genetic effect.
#
# Reads NO DB, does NOT edit original data, outputs remain gitignored.
# N_STAGE80 (default 99999) uses ALL eligible genetic individuals; optional subset uses seed 79.
# ============================================================================
options(stringsAsFactors=FALSE)
fail <- function(...)stop(paste0(...),call.=FALSE)
base <- "outputs/sei_transitions"
files <- c(
  diagnostics=file.path(base,"stage76_PRIVATE_CANDIDATE_NOT_VALIDATED.csv"),
  crosswalk=file.path(base,"stage76_private_id_crosswalk.csv"),
  genetics="data/model_inputs/integrated_gm_inputs_group_proxy_2026.rds",
  originalfit="data/infection_posteriors/sei_2026/fits/sei_inferred_accuracy_.rds")
if(!all(file.exists(files)))
  fail("Missing Stage76/GM/original SEI inputs: ",
       paste(names(files)[!file.exists(files)],collapse=", "))
normal <- function(z)toupper(trimws(as.character(z)))
cand <- read.csv(files["diagnostics"],na.strings=c("","NA"),
                 check.names=FALSE)
ids <- read.csv(files["crosswalk"],colClasses=c(tattoo="character"))
x <- readRDS(files["genetics"])
f <- readRDS(files["originalfit"])
cols <- c("time","id","captured",paste0("test",1:6))
if(!all(cols%in%names(cand))||
   !all(c("id","tattoo")%in%names(ids))||
   !all(c("id","hom","sex","birth_q")%in%names(x))||
   !all(c("test","Se","Sp")%in%names(f$sesp)))
  fail("Required input columns changed; do not guess mappings")
stopifnot(nrow(ids)==3214L,length(x$id)==1603L,
          !anyDuplicated(ids$id),!anyDuplicated(ids$tattoo),
          !anyDuplicated(normal(x$id)),all(cand$captured==1L),
          all(cand$time>=1L&cand$time<=200L),
          !anyNA(cand[,c("time","id","captured")]),
          all(cand$id %in% ids$id))
for(nm in paste0("test",1:6))
  if(!all(is.na(cand[[nm]])|cand[[nm]]%in%c(0,1)))
    fail("Nonbinary diagnostic result in ",nm)
if(anyNA(x$hom)||any(!is.finite(as.numeric(x$hom)))||
   any(!x$sex %in% c("MALE","FEMALE"))||
   anyNA(x$birth_q)||sd(as.numeric(x$hom))<=0)
  fail("Genetic input not valid")
if(!identical(c("ELISA_BELT","ELISA_INDIRECT","Culture",
                "DPP","IGRA","StatPak"),
              as.character(f$sesp$test_name)))
  fail("Original six-test ordering changed; do not silently swap assays")
Se <- as.numeric(f$sesp$Se)
Sp <- as.numeric(f$sesp$Sp)
if(anyNA(Se)||anyNA(Sp)||any(Se<=0|Se>=1)||
   any(Sp<=0|Sp>=1))fail("Original Se/Sp invalid")
# Original SEI emission: P(positive | E, assay k) = m_k * Se_k.
# The source has 4-chain posterior DRAWS (Stage75 diagnostics include
# e_mult.4.). In this fit f$parameters does NOT expose the multiplier
# rows, so calculate each m_k from the original draw columns instead.
# Never use diagnostic Rhat/ESS as a multiplier, or guess a default.
# Julia's e_mult[k] can become e_mult.k. after read.csv(check.names=TRUE).
draws <- f$draws
if(!is.data.frame(draws)||nrow(draws)<10L)
  fail("Original SEI posterior draw table unavailable or too small")
names_raw <- names(draws)
where <- grep("^e_mult(\\[|\\.|$)",names_raw)
extract_index <- function(label) {
  # Exact accepted forms: e_mult[1], e_mult.1., e_mult.1
  if(grepl("^e_mult\\[[1-6]\\]$",label))
    return(as.integer(sub("^e_mult\\[([1-6])\\]$","\\1",label)))
  if(grepl("^e_mult\\.[1-6]\\.?$",label))
    return(as.integer(sub("^e_mult\\.([1-6])\\.?$","\\1",label)))
  NA_integer_
}
indices <- vapply(names_raw[where],extract_index,integer(1),
                  USE.NAMES=FALSE)
if(length(where)!=6L||anyNA(indices)||
   !identical(sort(indices),1:6)){
  cat("Source posterior draw columns matching e_mult: ",
      paste(names_raw[where],collapse=", "),"\n",sep="")
  cat("Last 16 posterior draw column names: ",
      paste(tail(names_raw,16L),collapse=", "),"\n",sep="")
  if(is.data.frame(f$parameters)&&"parameter"%in%names(f$parameters))
    cat("Parameter-table e_mult labels: ",
        paste(grep("mult",as.character(f$parameters$parameter),
                   ignore.case=TRUE,value=TRUE),collapse=", "),"\n",sep="")
  fail("Cannot map all 6 exposed assay multipliers from posterior draws.")
}
mapped <- where[match(1:6,indices)]
if(!all(vapply(draws[mapped],is.numeric,logical(1))))
  fail("Exposed multipliers are not numeric posterior columns")
m <- vapply(draws[mapped],function(v) {
  if(anyNA(v)||any(!is.finite(v))||any(v<0|v>1))
    fail("SEI exposed multiplier posterior values outside [0,1]")
  mean(v)
},numeric(1),USE.NAMES=FALSE)
cat("Source exposed-state multipliers recovered from ",
    nrow(draws)," joint posterior draw rows:\n",sep="")
print(data.frame(assay=seq_len(6L),
                 source_column=names_raw[mapped],
                 posterior_mean=signif(m,5L)),row.names=FALSE)
if(anyNA(m)||any(m<=0|m>=1))
  fail("Original exposed-state m multipliers are nonfinite/boundary")
positive <- rbind(1-Sp, m*Se, Se)
rownames(positive)<-c("S","E","I")
colnames(positive)<-paste0("test",1:6)
if(any(positive<=0|positive>=1))
  fail("Nonfinite emission likelihood inputs")
cat("\n=== STAGE 80: DIRECT-DIAGNOSTIC SEI GENETIC SMOKE TEST ===\n")
cat("Stage76 confirmed live captures:",nrow(cand),"\n")
cat("Original-only 317 uncertain capture quarters: EXCLUDED.\n")
cat("Within-quarter repeated captures: RETAINED/STACKED.\n")
cat("Source: inferred-accuracy SEI emission parameters, fixed for pilot.\n")

ids$tattoo <- normal(ids$tattoo)
gm_index <- match(ids$tattoo,normal(x$id))
gen_ids <- ids$id[!is.na(gm_index)]
c_by_id <- split(seq_len(nrow(cand)),as.character(cand$id))
available <- gen_ids[as.character(gen_ids)%in%names(c_by_id)]
# Minimum two distinct quarters and at least one usable assay.
valid <- vapply(available,function(i) {
  ii <- c_by_id[[as.character(i)]]
  gi <- gm_index[match(i,ids$id)]
  length(unique(cand$time[ii]))>=2L &&
    any(!is.na(as.matrix(cand[ii,paste0("test",1:6)]))) &&
    min(cand$time[ii])>=x$birth_q[gi]
},logical(1))
eligible <- available[valid]
n_requested <- suppressWarnings(as.integer(Sys.getenv("N_STAGE80","99999")))
if(is.na(n_requested)||n_requested<40L)
  fail("N_STAGE80 must be an integer >= 40")
if(length(eligible)<40L)
  fail("Too few genotyped badgers with valid repeated diagnostic histories")
set.seed(79L)
chosen <- sort(sample(eligible,min(n_requested,length(eligible)),
                      replace=FALSE))
gis <- gm_index[match(chosen,ids$id)]
hz <- (as.numeric(x$hom)-mean(as.numeric(x$hom)))/
      sd(as.numeric(x$hom))
stopifnot(length(chosen)==length(gis),!anyNA(gis))
cat("Genotyped badgers with >=2 diagnostic quarters:",length(eligible),"\n")
cat("Selected:",length(chosen),"of",length(eligible),
    "eligible genotyped badgers; seed 79 for reproducibility.\n")

# Precompute per-quarter, per-state LOG emissions.
# One positive/negative test contributes a Bernoulli emission; repeated
# captures with different sample results are stacked as original HMM mode.
# Assumes conditional independence of tests/captures given latent state.
log_positive <- log(positive)
log_negative <- log1p(-positive)
obs <- vector("list",length(chosen))
total_captures <- 0L
total_assays <- 0L
for(k in seq_along(chosen)){
  idx <- c_by_id[[as.character(chosen[k])]]
  d <- cand[idx,,drop=FALSE]
  ord <- order(d$time)
  d <- d[ord,,drop=FALSE]
  qs <- seq.int(min(d$time),max(d$time))
  em <- matrix(0,length(qs),3L)
  for(j in seq_len(nrow(d))) {
    q <- match(d$time[j],qs)
    for(assay in 1:6) {
      v <- d[j,paste0("test",assay)]
      if(is.na(v))next
      em[q,] <- em[q,] + if(v==1L)
        log_positive[,assay] else log_negative[,assay]
      total_assays <- total_assays+1L
    }
  }
  obs[[k]] <- list(em=em,nq=nrow(em),
                   h=hz[gis[k]],capture_count=nrow(d))
  total_captures <- total_captures+nrow(d)
}
cat("Pilot capture records:",total_captures,
    "; observed binary tests:",total_assays,"\n")
cat("Median modelled quarters spanning first-to-last capture:",median(vapply(obs,function(z)z$nq,
     integer(1))),"\n")
if(total_assays<100L)
  fail("Too few observed test outcomes for even a smoke test")

# Probability and recursion unit tests, independent of candidate dataset.
# Empty emissions must integrate to likelihood 1 (logL=0), independent
# of starting state and transition probabilities.
hmm_forward <- function(em,lam,tau,pi,rho) {
  if(!is.matrix(em)||ncol(em)!=3L||nrow(em)<1L)
    stop("Invalid emissions")
  a <- c(1-pi, pi*(1-rho), pi*rho)
  logL <- 0
  for(t in seq_len(nrow(em))){
    if(t>1L) {
      a <- c(a[1L]*(1-lam),
             a[1L]*lam+a[2L]*(1-tau),
             a[2L]*tau+a[3L])
    }
    offset <- max(em[t,])
    w <- exp(em[t,]-offset)
    a <- a*w
    scale <- sum(a)
    if(!is.finite(scale)||scale<=0)stop("Impossible latent likelihood")
    logL <- logL+log(scale)+offset
    a <- a/scale
  }
  logL
}
stopifnot(abs(hmm_forward(matrix(0,9L,3L),.12,.2,.1,.4))<1e-9,
          abs(hmm_forward(matrix(0,1L,3L),.12,.2,.1,.4))<1e-9)
# A one-quarter positive test has exact initial-mixture likelihood,
# regardless of the transition hazards (no transition before entry).
pi0 <- .2; rho0 <- .4
p_exact <- (1-pi0)*positive[1L,3L]+
  pi0*(1-rho0)*positive[2L,3L]+
  pi0*rho0*positive[3L,3L]
em1 <- matrix(log_positive[,3L],1L,3L)
stopifnot(abs(hmm_forward(em1,.12,.2,pi0,rho0)-log(p_exact))<1e-10)
cat("Forward-likelihood unit tests: PASSED.\n")

# Simplified HMM parameters:
# (1) S->E baseline logit, (2) E->I baseline logit,
# (3) initial infected prevalence logit, (4) fraction of initial infected
# infectious logit, (5) acquisition hom_z slope, (6) progression hom_z slope.
# Beta(5:6) are per original-cohort SD of multilocus homozygosity.
# This deliberately excludes age/sex/calendar/site and is not publishable.
start <- c(logit_SE=qlogis(.035),logit_EI=qlogis(.05),
           logit_entry_pi=qlogis(.10),
           logit_entry_rho=qlogis(.5),
           beta_hom_SE=0,beta_hom_EI=0)
lower <- c(-9,-9,-9,-7,-4,-4)
upper <- c(0,0,0,7,4,4)
prior_mean <- c(qlogis(.035),qlogis(.05),
                qlogis(.10),qlogis(.5),0,0)
prior_sd <- c(2,2,2,2.5,1.5,1.5)
raw_loglik <- function(theta) {
  pi <- plogis(theta[3L]);rho <- plogis(theta[4L])
  total <- 0
  for(d in obs) {
    lam <- plogis(theta[1L]+theta[5L]*d$h)
    tau <- plogis(theta[2L]+theta[6L]*d$h)
    total <- total+hmm_forward(d$em,lam,tau,pi,rho)
  }
  total
}
neglogpost <- function(theta) {
  -raw_loglik(theta)+.5*sum(((theta-prior_mean)/prior_sd)^2)
}
# For a given nested model, omitted slopes are exactly zero and the
# nuisance parameters are RE-FITTED (not held at joint-model estimates).
# Both penalised objective and unpenalised diagnostic likelihood are saved.
fit_nested <- function(label,active,initial,maxit=250L){
  stopifnot(length(initial)==6L,all(active%in%1:6))
  theta_base <- start
  theta_base[active] <- initial[active]
  objective <- function(par){
    theta_full <- theta_base
    theta_full[active] <- par
    neglogpost(theta_full)
  }
  t0 <- proc.time()[["elapsed"]]
  opt <- stats::optim(theta_base[active],objective,method="L-BFGS-B",
           lower=lower[active],upper=upper[active],
           control=list(maxit=maxit,trace=0))
  theta <- theta_base
  theta[active] <- opt$par
  elapsed <- proc.time()[["elapsed"]]-t0
  ll <- raw_loglik(theta)
  if(!is.finite(opt$value)||!is.finite(ll))
    fail(label,": invalid fitted likelihood")
  cat(sprintf("%-22s convergence=%d, logLik=%8.3f, penaltyObjective=%8.3f, elapsed=%.1f s\n",
              label,opt$convergence,ll,opt$value,elapsed))
  flush.console()
  list(label=label,theta=theta,logLik=ll,
       objective=opt$value,convergence=opt$convergence,
       message=opt$message,elapsed=elapsed,
       active=active)
}

cat("\n=== STAGE 80: ALL-ELIGIBLE NESTED LIKELIHOOD FITS ===\n")
cat("Genotyped individuals fitted:",length(chosen),"of",length(eligible),"eligible.\n")
cat("Initial logLik:",round(raw_loglik(start),3),"\n")
cat("Fitting null, acquisition only, progression only, and full model.\n")
cat("All nested comparisons use the same animals, tests and priors.\n")
t0 <- proc.time()[["elapsed"]]
null <- fit_nested("neither genetic effect",1:4,start)
a0 <- null$theta
a0[5L] <- 0; a0[6L] <- 0
acquisition <- fit_nested("acquisition only",1:5,a0)
progression <- fit_nested("progression only",c(1:4,6),a0)
joint_from_null <- fit_nested("both: null start",1:6,a0)
joint_from_initial <- fit_nested("both: original start",1:6,start)
joint <- if(joint_from_null$objective<=joint_from_initial$objective)
   joint_from_null else joint_from_initial

fits <- list(null=null,acquisition_only=acquisition,
             progression_only=progression,both=joint,
             joint_start2=if(identical(joint$label,"both: null start"))
                joint_from_initial else joint_from_null)
results <- do.call(rbind,lapply(names(fits)[1:4],function(key){
  fit <- fits[[key]]
  data.frame(model=key,parameters=length(fit$active),
    n_badgers=length(chosen),n_captures=total_captures,
    observed_tests=total_assays,
    logLik=fit$logLik,neg_log_posterior=fit$objective,
    convergence=fit$convergence,
    beta_hom_SE=fit$theta[5L],
    beta_hom_EI=fit$theta[6L],
    OR_hom_SE=exp(fit$theta[5L]),
    OR_hom_EI=exp(fit$theta[6L]),
    baseline_SE_q=plogis(fit$theta[1L]),
    baseline_EI_q=plogis(fit$theta[2L]),
    entry_pi=plogis(fit$theta[3L]),
    entry_rho=plogis(fit$theta[4L]))
}))
cat("\n--- COMPARATIVE FITS ---\n")
print(results,row.names=FALSE,digits=5)
cat("\nDescriptive diagnostic log-likelihood improvements vs null:\n")
for(key in c("acquisition_only","progression_only","both")){
  fcur <- fits[[key]]
  cat(key,": 2*deltaLogLik =",
      signif(2*(fcur$logLik-null$logLik),5),"\n")
}
cat("Difference between two full-model optimised objectives:",
    signif(abs(joint_from_initial$objective-
               joint_from_null$objective),5),"\n")
cat("NOTE: these are penalized MAP solutions with fixed Se/Sp.",
    "Likelihood differences are DESCRIPTIVE, NOT calibrated",
    "chi-square hypothesis tests or posterior probabilities.\n")

summary_path <- file.path(base,"stage80_all_eligible_model_comparison.csv")
fit_path <- file.path(base,"stage80_all_eligible_joint_sei.rds")
write.csv(results,summary_path,row.names=FALSE)
saveRDS(list(
  fits=fits,summary=results,
  n_badgers=length(chosen),n_eligible=length(eligible),
  n_captures=total_captures,n_tests=total_assays,
  sample_seed=79L,ids=chosen,
  assay_parameters=list(Se=Se,Sp=Sp,exposed_multiplier=m),
  source_md5=tools::md5sum(files),
  run_elapsed_seconds=proc.time()[["elapsed"]]-t0,
  assumptions=c(
    "Confirmed TRAP/NET/DART GUN captures from Stage76 only",
    "317 original-only uncertain encounter occasions excluded",
    "genotyped with >=2 quarters and at least one readable diagnostic",
    "source exposed sensitivity multiplier means from 4000 posterior draws",
    "fixed source-fitted Se/Sp and exposed multipliers",
    "mixed M.bovis+/negative cultures coded positive at capture",
    "ambiguous DPP and other test results treated NA",
    "repeat within-quarter records stacked as independent emissions",
    "quartely irreversible S->E->I transition model; no direct S->I",
    "constant transition baseline hazards, common initial prevalence",
    "all latent states analytically marginalized with forward HMM",
    "MAP with fixed weak Gaussian priors, nested nuisance refits",
    "both homozygosity slopes are log-odds for 1 SD original genetic sample",
    "omits sex/age/year/group/observation process and test parameter uncertainty",
    "genetic estimates uncalibrated, not confirmatory or causal"
  )),fit_path)
cat("\nSaved:",summary_path,"\n")
cat("Saved:",fit_path,"\n")
if(any(vapply(fits[1:4],function(z)z$convergence!=0L,logical(1))))
  fail("One or more nested fits did not converge; inspect output")
if(joint$objective > null$objective+1e-5)
  fail("Joint fit worse penalised objective than nested null")
cat("STAGE 80 COMPLETED. NEXT: simulation recovery and",
    "time/age/sex diagnostic-bias sensitivity before biology.\n")