# ============================================================================
# Stage 81: DIRECT-DIAGNOSTIC, GENETIC SEI JOINT-LIKELIHOOD SMOKE TEST
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
# SIMULATION-RECOVERY CHECK: constant quarterly hazards, fixed Se/Sp
# and exposed-state multipliers from original fitted SEI; no seasonal,
# calendar, age, sex, social-group, survival, trapping or observation model.
# Regularized MAP is NOT a fully Bayesian fit nor a final genetic effect.
#
# Reads NO DB, does NOT edit original data, outputs remain gitignored.
# N_STAGE81 (default 99999) uses ALL eligible genetic individuals; optional subset uses seed 79.
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
cat("\n=== STAGE 81: DIRECT-DIAGNOSTIC SEI GENETIC SMOKE TEST ===\n")
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
n_requested <- suppressWarnings(as.integer(Sys.getenv("N_STAGE81","250")))
if(is.na(n_requested)||n_requested<40L)
  fail("N_STAGE81 must be an integer >= 40")
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


# ===========================================================================
# The simulation uses the EXACT capture-quarter and six-assay observation
# opportunity pattern in the selected Stage76 candidate records. It preserves:
#   - individual's hom_z and observed study window,
#   - number of actual readable assays of each type in each quarter,
#   - within-quarter repeated captures and test missingness.
# It SIMULATES NEW OUTCOMES (0/1) under specified S/E/I paths.
#
# The generating process intentionally MATCHES the simplified fitted model.
# That makes this a computational/identifiability check conditional on the
# fitted likelihood assumptions, NOT verification of biological truth.
# ===========================================================================
source_stage80 <- file.path(base,"stage80_all_eligible_joint_sei.rds")
if(!file.exists(source_stage80))
  fail("Need completed Stage80 fit: ",source_stage80)
ref <- readRDS(source_stage80)
if(!is.list(ref)||is.null(ref$fits$both$theta)||
   !is.numeric(ref$fits$both$theta)||
   length(ref$fits$both$theta)!=6L ||
   any(!is.finite(ref$fits$both$theta)) ||
   !identical(as.integer(ref$n_eligible),as.integer(length(eligible))))
  fail("Stage80 fit is missing or its genetic eligibility cohort differs")
truth_alt <- as.numeric(ref$fits$both$theta)
names(truth_alt) <- c("logit_SE","logit_EI","logit_entry_pi",
                      "logit_entry_rho","beta_hom_SE","beta_hom_EI")
truth_null <- truth_alt
truth_null[5:6] <- 0
if(any(truth_alt[1:4] < -9) || any(truth_alt[1:2] > 0))
  fail("Unexpected Stage80 baseline parameter values")

# Count the observed assays by quarter for the chosen genotyped badgers.
# Stage76 test values themselves are NOT reused for simulated outcomes.
schedules <- vector("list",length(chosen))
for(k in seq_along(chosen)){
  idx <- c_by_id[[as.character(chosen[k])]]
  d <- cand[idx,,drop=FALSE]
  q <- seq.int(min(d$time),max(d$time))
  counts <- matrix(0L,nrow=length(q),ncol=6L)
  for(j in seq_len(nrow(d))) {
    ii <- match(d$time[j],q)
    v <- unlist(d[j,paste0("test",1:6)],use.names=FALSE)
    counts[ii,] <- counts[ii,]+as.integer(!is.na(v))
  }
  if(!identical(as.integer(rowSums(counts)>0),
                as.integer(rowSums(counts)>0)))
    fail("Internal test mask consistency")
  if(sum(counts)<1L || nrow(counts)!=obs[[k]]$nq)
    fail("Per-animal assay schedule does not match Stage81 preparation")
  schedules[[k]] <- counts
}
if(sum(vapply(schedules,sum,numeric(1)))!=total_assays)
  fail("Simulated assay opportunities differ from real observation counts")
cat("\n=== STAGE 81: SIMULATED GENETIC EFFECT RECOVERY ===\n")
cat("Animals:",length(chosen),
    "observed assay opportunities held fixed:",total_assays,"\n")
cat("Original Stage80 beta_hom_SE:",truth_alt[5L],
    "beta_hom_EI:",truth_alt[6L],"\n")
cat("Initial state parameters and Se/Sp also fixed from Stage80.\n")

simulate_one <- function(truth,seed) {
  set.seed(as.integer(seed))
  pi <- plogis(truth[3L])
  rho <- plogis(truth[4L])
  sim_obs <- vector("list",length(obs))
  latent_counts <- integer(3L)
  for(k in seq_along(obs)) {
    d <- obs[[k]]
    mask <- schedules[[k]]
    nq <- nrow(mask)
    em <- matrix(0,nrow=nq,ncol=3L)
    lam <- plogis(truth[1L]+truth[5L]*d$h)
    tau <- plogis(truth[2L]+truth[6L]*d$h)
    state <- sample.int(3L,1L,
                        prob=c(1-pi,pi*(1-rho),pi*rho))
    for(q in seq_len(nq)) {
      if(q>1L) {
        if(state==1L && runif(1L)<lam)
          state <- 2L
        else if(state==2L && runif(1L)<tau)
          state <- 3L
      }
      latent_counts[state] <- latent_counts[state]+1L
      for(a in 1:6) {
        n <- mask[q,a]
        if(n==0L)next
        npos <- rbinom(1L,size=n,prob=positive[state,a])
        em[q,] <- em[q,]+npos*log_positive[,a]+
          (n-npos)*log_negative[,a]
      }
    }
    sim_obs[[k]] <- list(em=em,h=d$h,nq=nq)
  }
  list(obs=sim_obs,latent_counts=latent_counts)
}

start <- c(logit_SE=qlogis(.035),logit_EI=qlogis(.05),
           logit_entry_pi=qlogis(.10),
           logit_entry_rho=qlogis(.5),
           beta_hom_SE=0,beta_hom_EI=0)
lower <- c(-9,-9,-9,-7,-4,-4)
upper <- c(0,0,0,7,4,4)
prior_mean <- c(qlogis(.035),qlogis(.05),
                qlogis(.10),qlogis(.5),0,0)
prior_sd <- c(2,2,2,2.5,1.5,1.5)

# Rerun the Stage80 penalised MAP using only the newly simulated outcomes.
# Starting all optimisations from the SAME non-genetic initial values
# ensures recovery isn't handed the known true genetic coefficients.
fit_one <- function(sim_obs,maxit=200L) {
  likelihood <- function(theta) {
    pi <- plogis(theta[3L]); rho <- plogis(theta[4L])
    ll <- 0
    for(d in sim_obs) {
      lam <- plogis(theta[1L]+theta[5L]*d$h)
      tau <- plogis(theta[2L]+theta[6L]*d$h)
      ll <- ll+hmm_forward(d$em,lam,tau,pi,rho)
    }
    ll
  }
  objective <- function(theta)
    -likelihood(theta)+.5*sum(((theta-prior_mean)/prior_sd)^2)
  fit <- stats::optim(start,objective,method="L-BFGS-B",
           lower=lower,upper=upper,
           control=list(maxit=maxit,trace=0))
  list(theta=as.numeric(fit$par),
       convergence=fit$convergence,
       objective=fit$value,
       logLik=likelihood(fit$par),
       improvement=objective(start)-fit$value,
       message=fit$message)
}

nsim <- suppressWarnings(as.integer(Sys.getenv("N_SIM","4")))
if(is.na(nsim)||nsim<2L||nsim>100L)
  fail("N_SIM must be between 2 and 100")
out_rows <- vector("list",2L*nsim)
t_begin <- proc.time()[["elapsed"]]
j <- 0L
for(scenario in c("null","stage80_effects")) {
  truth <- if(scenario=="null") truth_null else truth_alt
  for(rep in seq_len(nsim)) {
    j <- j+1L
    seed <- 81000L+if(scenario=="null") rep else 1000L+rep
    sim <- simulate_one(truth,seed)
    fit <- fit_one(sim$obs)
    out_rows[[j]] <- data.frame(
      scenario=scenario,rep=rep,seed=seed,
      n_badgers=length(chosen),n_tests=total_assays,
      true_beta_SE=truth[5L],fit_beta_SE=fit$theta[5L],
      true_beta_EI=truth[6L],fit_beta_EI=fit$theta[6L],
      delta_beta_SE=fit$theta[5L]-truth[5L],
      delta_beta_EI=fit$theta[6L]-truth[6L],
      true_logit_SE=truth[1L],fit_logit_SE=fit$theta[1L],
      true_logit_EI=truth[2L],fit_logit_EI=fit$theta[2L],
      true_entry_pi=plogis(truth[3L]),
      fit_entry_pi=plogis(fit$theta[3L]),
      true_entry_rho=plogis(truth[4L]),
      fit_entry_rho=plogis(fit$theta[4L]),
      convergence=fit$convergence,
      logLik=fit$logLik,
      objective_improvement=fit$improvement,
      n_quarters_S=sim$latent_counts[1L],
      n_quarters_E=sim$latent_counts[2L],
      n_quarters_I=sim$latent_counts[3L])
    cat(sprintf(
      "%s sim %d/%d: convergence=%d, beta(SE) true=%+.3f fitted=%+.3f; beta(EI) true=%+.3f fitted=%+.3f; elapsed=%.1fs\n",
      scenario,rep,nsim,fit$convergence,truth[5L],fit$theta[5L],
      truth[6L],fit$theta[6L],proc.time()[["elapsed"]]-t_begin))
    flush.console()
  }
}
results <- do.call(rbind,out_rows)
if(nrow(results)!=2L*nsim)fail("Simulation count mismatch")
print(results[,c("scenario","rep","true_beta_SE","fit_beta_SE",
                "true_beta_EI","fit_beta_EI","convergence")],
      row.names=FALSE,digits=3)
if(any(results$convergence!=0L))
  cat("WARNING: Some optimisations did not converge.\n")
for(scenario in c("null","stage80_effects")) {
  z <- results[results$scenario==scenario,,drop=FALSE]
  cat("\n",scenario,": fitted genetic-effect ranges:\n",sep="")
  cat("SE beta:",paste(signif(range(z$fit_beta_SE),3),collapse=" to "),
      "; EI beta:",paste(signif(range(z$fit_beta_EI),3),
                         collapse=" to "),"\n")
  cat("Median beta bias (SE,EI):",
      signif(median(z$delta_beta_SE),3),
      signif(median(z$delta_beta_EI),3),"\n")
}
dst <- file.path(base,"stage81_simulation_recovery.csv")
write.csv(results,dst,row.names=FALSE)
saveRDS(list(results=results,truth=list(null=truth_null,
                                        stage80_effects=truth_alt),
  n_individuals=length(chosen),n_tests=total_assays,
  nsim_per_scenario=nsim,model="simple 3-state HMM, fixed assay accuracy",
  source_model_md5=unname(tools::md5sum(source_stage80)),
  source_candidate_md5=unname(tools::md5sum(files["diagnostics"])),
  caveats=c(
    "simulates from the SAME simplified likelihood being fitted",
    "not a test of observation bias or true biological assumptions",
    "only selected genotyped badgers with repeated observable captures",
    "simulated assays have fixed original emission probabilities",
    "preserves assay opportunities, not outcome-dependent test ordering",
    "same weak Gaussian regularisation in generative fit",
    "four replicates per scenario by default insufficient for calibration",
    "no uncertainty intervals or significance testing"
  )),file.path(base,"stage81_simulation_recovery.rds"))
cat("\nSaved:",dst,"\n")
cat("These few simulations are an ENGINEERING AND IDENTIFIABILITY",
    "check only. They cannot establish bias, coverage or reliability",
    "at scientific publication standards.\n")
if(any(results$convergence!=0L))
  fail("Some simulated refits did not converge; see saved diagnostics")
cat("STAGE 81 SIMULATION RECOVERY PILOT COMPLETE.\n")