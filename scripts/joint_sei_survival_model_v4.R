# ============================================================================
# V4 (new files on joint-v3b-benton): copy of joint_sei_survival_model.R that fits
# R/joint_sei_survival_hierarchy_nimble_v4.R (v3 + euthanised/killed
# badgers censored at removal; kernel R/joint_sei_survival_forward_nimble_removal.R). Only hierarchical_compile / full / benchmark
# modes with the default samplers; v1 files, inputs and chains are read-only.
# ============================================================================
# ONE TARGET MODEL: JOINT INFECTION + SURVIVAL
#
# Source contract: docs/one_joint_infection_survival_model_target.md
# Initial MODE=validate prepares the actual full-genetic-cohort interface
# and checks a FOUR-STATE Markov survival/infection transition kernel.
#
# THIS FILE IS THE SINGLE FUTURE JOINT MODEL CODEBASE.
# Do not interpret MODE=validate as a fitted biological model.
#
# FIT CONTRACT — implemented in full hierarchical model:
# - Current-quarter age + sex + standardized homozygosity in ALL three
#   biological hazards: S->E, E->I and S/E/I->death.
# - Persistent, mean-zero, correlated individual random effects (frailties)
#   for S->E, E->I and death; hierarchical covariance and shrinkage priors.
# - Separate capture heterogeneity (individual random intercept with
#   shrinkage if identifiable) to avoid conflating missed captures and death.
# - A four-state forward HMM conditional on the individual's frailties,
#   embedded in hierarchical inference sampling/integrating frailties.
# - Age-dependent Gompertz-Makeham death baseline; age and regularized
#   historical trends in infection transitions. Fixed per-assay Se/Sp.
# - Uncertainty for genetic effects AND random-effect variances/correlations.
#
# Full MCMC is enabled with MODE=full; validation and engineering gates
# remain separate, using this SAME biological model, not variants.
# Mortality inferences are still PROVISIONAL until deliberate management
# deaths, capture effort and parameter recovery are investigated.
#
# No DB writes; no database reads; all exports private and gitignored.
# ============================================================================
options(stringsAsFactors=FALSE)
mode <- Sys.getenv("MODE","validate")
sampler_scheme <- Sys.getenv("SAMPLER_SCHEME","default")
if(!sampler_scheme%in%c("default","blocked","af_slice",
                       "af_slice_calendar"))
  stop("SAMPLER_SCHEME must be default, blocked, af_slice or af_slice_calendar")
if(!mode%in%c("hierarchical_compile","hierarchical_benchmark","full"))
  stop("v4 runs only hierarchical_compile, hierarchical_benchmark or full")
if(sampler_scheme!="default")
  stop("v4 uses the default NIMBLE samplers only")
if(sampler_scheme!="default" && !mode%in%
   c("hierarchical_benchmark","full"))
  stop("Special MCMC samplers apply only to benchmark/full modes")
if(!mode %in% c("validate","likelihood_check","compile",
                  "hierarchical_compile","hierarchical_benchmark",
                  "full"))
  stop("MODE must be validate, likelihood_check, compile,",
       " hierarchical_compile, hierarchical_benchmark or full.")
fail <- function(...)stop(paste(...,collapse=""),call.=FALSE)
normal <- function(z) {
  x <- toupper(trimws(as.character(z)))
  x[is.na(z)|!nzchar(x)] <- NA_character_
  x
}
files <- c(
  survival="data/model_inputs/integrated_gm_inputs_group_proxy_2026.rds",
  diagnostic="outputs/sei_transitions/stage76_PRIVATE_CANDIDATE_NOT_VALIDATED.csv",
  crosswalk="outputs/sei_transitions/stage76_private_id_crosswalk.csv",
  sei_fit="data/infection_posteriors/sei_2026/fits/sei_inferred_accuracy_.rds")
if(!all(file.exists(files)))
  fail("Missing original inputs: ",
       paste(names(files)[!file.exists(files)],collapse=", "))
outdir <- "outputs/joint_sei_survival"
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

# Four-state quarterly transition. Biology:
#  - survive current quarter THEN S->E or E->I, at most one step;
#  - dead is absorbing;
#  - no direct S->I, E->S, I->E.
#   phi_S, phi_E, phi_I: quarterly survival probabilities,
#   lambda: S->E conditional on surviving,
#   tau: E->I conditional on surviving.
# Death probability from each live state is 1-phi(state).
make_transition <- function(phi,lambda,tau) {
  if(length(phi)!=3L || any(!is.finite(phi))||
     any(phi<0|phi>1) ||
     length(lambda)!=1L||length(tau)!=1L||
     !is.finite(lambda)||!is.finite(tau)||
     lambda<0||lambda>1||tau<0||tau>1)
    fail("Invalid four-state quarterly transition probabilities")
  Tmat <- matrix(0,4L,4L,
                 dimnames=list(c("S","E","I","D"),
                               c("S","E","I","D")))
  Tmat[1,1] <- phi[1]*(1-lambda)
  Tmat[1,2] <- phi[1]*lambda
  Tmat[1,4] <- 1-phi[1]
  Tmat[2,2] <- phi[2]*(1-tau)
  Tmat[2,3] <- phi[2]*tau
  Tmat[2,4] <- 1-phi[2]
  Tmat[3,3] <- phi[3]
  Tmat[3,4] <- 1-phi[3]
  Tmat[4,4] <- 1
  if(any(Tmat < -1e-12) || any(Tmat > 1+1e-12)||
     any(abs(rowSums(Tmat)-1)>1e-12))
    fail("Four-state transition rows are not stochastic")
  Tmat
}
# Dynamic biological hazard kernel used by the ONE joint model.
# This is the implementation for age/sex/genotype/individual frailty, but
# not an estimate of those effects. The correlated individual frailty
# PRIOR and time-varying historical baseline will be fitted together with
# all four-state transitions, after capture/death reconciliation.
exprel <- function(z){
  if(abs(z)<1e-6)
    1+z/2+z*z/6+z*z*z/24
  else
    expm1(z)/z
}
gm_quarter_hazard <- function(age_years,log_a,log_b,log_c){
  if(!is.finite(age_years)||age_years<0L)
    fail("Mortality evaluated at negative/nonfinite age")
  a <- exp(log_a)
  b <- exp(log_b)
  c <- exp(log_c)
  H <- .25*(a*exp(b*age_years)*exprel(b*.25)+c)
  if(!is.finite(H)||H<0)
    fail("Nonfinite Gompertz-Makeham integrated hazard")
  H
}
# theta convention: each transition coefficients named
#   intercept, age4, male, hom
# death coefficients: male,hom,E,I
# frailty: named SE,EI,D persistent within an individual.
# calendar: named SE,EI quarter-specific baseline log-odds offsets.
make_covariate_transition <- function(
    age_years,male,hom_z,frailty,calendar,theta) {
  if(!is.list(theta)||!all(c("SE","EI","D","gm")%in%names(theta))||
     !all(c("intercept","age4","male","hom")%in%names(theta$SE))||
     !all(c("intercept","age4","male","hom")%in%names(theta$EI))||
     !all(c("male","hom","E","I")%in%names(theta$D))||
     !all(c("log_a","log_b","log_c")%in%names(theta$gm))||
     !all(c("SE","EI","D")%in%names(frailty))||
     !all(c("SE","EI")%in%names(calendar))||
     !all(is.finite(c(age_years,male,hom_z,frailty,calendar)))||
     !male%in%c(0L,1L)||age_years<0)
    fail("Invalid covariate transition inputs")
  age4 <- age_years/4
  eta_SE <- theta$SE["intercept"]+
       theta$SE["age4"]*age4+theta$SE["male"]*male+
       theta$SE["hom"]*hom_z+frailty["SE"]+calendar["SE"]
  eta_EI <- theta$EI["intercept"]+
       theta$EI["age4"]*age4+theta$EI["male"]*male+
       theta$EI["hom"]*hom_z+frailty["EI"]+calendar["EI"]
  h <- gm_quarter_hazard(
       age_years,theta$gm["log_a"],theta$gm["log_b"],
       theta$gm["log_c"])
  death_shift <- theta$D["male"]*male+
       theta$D["hom"]*hom_z+frailty["D"]
  modifiers <- c(0,theta$D["E"],theta$D["I"])
  phi <- exp(-h*exp(death_shift+modifiers))
  make_transition(phi,unname(plogis(eta_SE)),unname(plogis(eta_EI)))
}
# Fixed demonstration parameters; ONLY for unit tests, not fitted values.
demo <- list(
  SE=c(intercept=-4,age4=-.2,male=.3,hom=.35),
  EI=c(intercept=-2,age4=.2,male=.1,hom=.45),
  D=c(male=.1,hom=.1,E=.2,I=.6),
  gm=c(log_a=log(.005),log_b=log(.15),log_c=log(.003)))
zero_re <- c(SE=0,EI=0,D=0)
zero_cal <- c(SE=0,EI=0)
t0 <- make_covariate_transition(2,0,0,zero_re,zero_cal,demo)
t1 <- make_covariate_transition(2,0,1,zero_re,zero_cal,demo)
t2 <- make_covariate_transition(
  2,0,0,c(SE=.7,EI=0,D=0),zero_cal,demo)
t3 <- make_covariate_transition(
  2,0,0,c(SE=0,EI=0,D=.7),zero_cal,demo)
stopifnot(
  all(abs(rowSums(t0)-1)<1e-12),
  t1["S","E"]>t0["S","E"],
  t2["S","E"]>t0["S","E"],
  t3["S","D"]>t0["S","D"],
  t1["E","I"]>t0["E","I"],
  t0["I","D"]>t0["S","D"])
cat("Age/sex/genotype/frailty dynamic transition kernel tests: PASS\n")

# Forward pass conditional on a known alive first capture;
# first occasion's detection probability cancels in conditional
# likelihood, but its diagnostic likelihood remains.
# obs[t] is 1=live captured, 0=monitored not captured,
# -1=not monitored. No exact recovery dates in this unit kernel.
# logem[t,1:3] is sum of diagnostic log likelihood at any captures.
# 4th state is death; impossible on an observed live capture.
forward_basic <- function(Tmat,obs,p,logem,entry) {
  K <- length(obs)
  if(K<1||length(p)!=K||nrow(logem)!=K||
     ncol(logem)!=3L||length(entry)!=3L||
     any(!obs%in%c(-1L,0L,1L))||
     abs(sum(entry)-1)>1e-10||
     any(entry<0|entry>1)||
     obs[1L]!=1L ||
     any(!is.finite(p)|p<0|p>1) ||
     !all(is.finite(logem)))
    fail("Invalid conditional-entry HMM inputs")
  a <- c(entry,0)
  LL <- 0
  for(t in seq_len(K)){
    if(t>1L)
      a <- drop(a %*% Tmat)
    if(t==1L) {
      # First confirmed capture is conditioned upon.
      em <- c(exp(logem[t,]),0)
    } else if(obs[t]==1L) {
      em <- c(rep(p[t],3L)*exp(logem[t,]),0)
    } else if(obs[t]==0L) {
      em <- c(rep(1-p[t],3L),1)
    } else {
      em <- rep(1,4L)
    }
    a <- a*em
    cscale <- sum(a)
    if(cscale<=0||!is.finite(cscale))
      fail("Zero/nonfinite HMM likelihood during forward recursion")
    LL <- LL+log(cscale)
    a <- a/cscale
  }
  LL
}

# --------------------------------------------------------------------------
# Reference joint disease + death-recovery likelihood.
#
# Within each quarter q:
#   1. Latent state at START of q is alive S/E/I or already dead D.
#   2. Observe live capture / no capture and any diagnostic assays.
#   3. Transition across q to the start of q+1 using state-specific
#      survival, S->E or E->I.
#   4. A recovered death reported IN quarter q is emitted ONLY on
#      the live->D transition, with recovery probability rho.
#
# Thus a live capture and an observed death in the SAME quarter is
# ALLOWED and has positive likelihood: captured before dying.
# The source builder held any recorded death with LATER live capture
# dates; same-date/within-quarter exact sequencing remains audited.
#
# With no reported death, live->D transitions are downweighted
# by (1-rho); already-dead D->D has no extra recovery factor.
# An unreported death can still occur, making later captures impossible.
#
# This is a pure-R REFERENCE IMPLEMENTATION for testing the eventual
# compiled hierarchical joint likelihood. It is NOT a fitted model.
# --------------------------------------------------------------------------
forward_joint_reference <- function(transitions,y,dead,p,logem,entry,rho) {
  K <- length(y)
  if(K<1L||!is.array(transitions)||
     !identical(dim(transitions),c(4L,4L,as.integer(K)))||
     length(dead)!=K||length(p)!=K||
     !is.matrix(logem)||!identical(dim(logem),c(K,3L))||
     length(entry)!=3L||
     any(!y%in%c(-1L,0L,1L))||
     any(!dead%in%c(0L,1L))||
     any(!is.finite(p)|p<0|p>1)||
     any(!is.finite(logem))||
     any(!is.finite(entry)|entry<0)||
     abs(sum(entry)-1)>1e-10||
     length(rho)!=1L||!is.finite(rho)||rho<0||rho>1||
     y[1L]!=1L||sum(dead)>1L)
    fail("Invalid source contract for joint reference HMM")
  if(sum(dead)==1L){
    dq <- which(dead==1L)
    if(dq<K && any(y[seq.int(dq+1L,K)]==1L))
      # This test sequence is IMPOSSIBLE under the disease/death model,
      # not a malformed input to the likelihood. Return log(0) rather
      # than throwing, so the impossibility unit test can succeed.
      return(-Inf)
  }
  a <- c(entry,0)
  logL <- 0
  for(t in seq_len(K)) {
    if(t==1L||y[t]==1L) {
      # First capture conditioned on, not multiplied by p[1].
      off <- max(logem[t,])
      live_em <- exp(logem[t,]-off)
      if(t!=1L)live_em <- live_em*p[t]
      obs_em <- c(live_em,0)
    } else {
      off <- 0
      if(y[t]==0L)
        obs_em <- c(rep(1-p[t],3L),1)
      else
        obs_em <- rep(1,4L)
    }
    a <- a*obs_em
    scale <- sum(a)
    if(!is.finite(scale)||scale<=0)
      return(-Inf)
    logL <- logL+log(scale)+off
    a <- a/scale
    M <- transitions[,,t]
    if(!all(is.finite(M))||any(M<0)||
       any(abs(rowSums(M)-1)>1e-8))
      fail("Invalid transition probability row in joint reference")
    if(dead[t]==1L){
      # Crucially this may happen in a captured-alive quarter.
      b <- c(0,0,0,rho*sum(a[1:3]*M[1:3,4L]))
    } else {
      b <- drop(a %*% M)
      # Subtract recovered live->death events that were not recorded.
      b[4L] <- b[4L]-
          rho*sum(a[1:3]*M[1:3,4L])
    }
    b[abs(b)<1e-14] <- 0
    if(any(b<0))
      fail("Negative transition probability after death report")
    scale <- sum(b)
    if(!is.finite(scale)||scale<=0)
      return(-Inf)
    logL <- logL+log(scale)
    a <- b/scale
  }
  logL
}
# Explicit independent synthetic tests for the 79 observed
# same-quarter live-capture and death cases.
make_reference_sequence <- function(M,K)
  array(rep(unname(M),K),dim=c(4L,4L,as.integer(K)))
unitM <- make_transition(c(.9,.9,.9),.3,.2)
unit1 <- make_reference_sequence(unitM,1L)
unit2 <- make_reference_sequence(unitM,2L)
unit3 <- make_reference_sequence(unitM,3L)
# First quarter: conditioned-on capture, then death observed.
# P(death in quarter and recovered) = 0.1*0.8.
test_same_q <- forward_joint_reference(
  unit1,1L,1L,.5,matrix(0,1,3),c(1,0,0),.8)
stopifnot(is.finite(test_same_q),
          abs(exp(test_same_q)-.08)<1e-12)
# Q1 survives, Q2 captured alive (p=.5), dies and is recovered.
test_later_q <- forward_joint_reference(
  unit2,c(1L,1L),c(0L,1L),c(.5,.5),
  matrix(0,2,3),c(1,0,0),.8)
stopifnot(is.finite(test_later_q),
          abs(exp(test_later_q)-(.9*.5*.1*.8))<1e-12)
# With no death recovery possible (rho=0), and future unmonitored
# intervals, everything after the conditioned first observation
# marginalizes to probability one.
test_unmon <- forward_joint_reference(
  unit3,c(1L,-1L,-1L),rep(0L,3L),
  rep(.5,3L),matrix(0,3,3),c(1,0,0),0)
stopifnot(abs(test_unmon)<1e-12)
# Death in q1 makes a q2 live capture impossible.
test_impossible <- forward_joint_reference(
  unit2,c(1L,1L),c(1L,0L),rep(.5,2L),
  matrix(0,2,3),c(1,0,0),.8)
stopifnot(identical(test_impossible,-Inf))
cat("Joint observation/death recovery reference HMM unit tests: PASS\n")

# Unit tests deliberately independent of research data.
tm <- make_transition(c(.9,.9,.9),.3,.2)
stopifnot(
  identical(dim(tm),c(4L,4L)),
  all(abs(rowSums(tm)-1)<1e-12),
  tm["S","I"]==0,
  tm["I","E"]==0,
  tm["D","D"]==1,
  abs(tm["S","D"]-.1)<1e-12)
# If captured alive in the following quarter with perfect detection,
# and all live states have survival probability .9, total conditional
# likelihood is .9 irrespective of infection transition probabilities.
synthetic <- forward_basic(tm,c(1L,1L),c(1,1),
                           matrix(0,2,3),c(1,0,0))
stopifnot(abs(synthetic-log(.9))<1e-12)
# If the second quarter is unmonitored, S/E/I/D is marginalized;
# the probability of the observed data is exactly 1.
stopifnot(abs(forward_basic(tm,c(1L,-1L),c(1,1),
                          matrix(0,2,3),c(1,0,0)))<1e-12)
# If the second quarter is monitored and NO capture with p=.5,
# P(no-capture)=P(dead)+P(alive)*.5 = .1+.9*.5 = .55.
stopifnot(abs(forward_basic(tm,c(1L,0L),c(1,.5),
                          matrix(0,2,3),c(1,0,0))-
              log(.55))<1e-12)
cat("\n=== ONE JOINT SEI + SURVIVAL MODEL: INPUT VALIDATION ===\n")
cat("Four-state transition and observation unit tests: PASS\n")
cat("State order: S,E,I,D (death absorbing)\n")
cat("Emission assumption: condition on first live capture.\n")

# --------------------------------------------------------------------------
# Engineering check against the PREVIOUSLY VALIDATED PRIVATE data.
# We deliberately do not refit or recompute the source audit in this mode.
# In particular, the 79 live-captured/dead-in-same-quarter histories
# must have finite positive likelihood under the observation/event order.
#
# Demonstration inputs below ARE NOT estimates or scientific results.
# The full hierarchical frailty prior and MCMC are not implemented yet.
# --------------------------------------------------------------------------
if(mode%in%c("likelihood_check","compile","hierarchical_compile",
             "hierarchical_benchmark","full")){
  validated <- file.path(outdir,"joint_PRIVATE_UNRECONCILED_input.rds")
  if(!file.exists(validated))
    fail("Run MODE=validate once to construct the private input")
  inp <- readRDS(validated)
  if(!isTRUE(inp$audit_only)||
     as.integer(inp$N)!=1603L||
     as.integer(inp$T)!=200L||
     !identical(unname(as.character(inp$source_md5)),
                unname(as.character(tools::md5sum(files))))||
     as.integer(inp$clock_discrepancies)!=0L||
     as.integer(inp$hard_chronology_conflicts)!=0L||
     as.integer(inp$same_quarter_capture_death)!=79L)
    fail("Cached private data/audit differs from the validated 1603-badger",
         " 79 same-quarter-death source; revalidate before inference")
  if(!identical(dim(inp$log_diagnostic_emission),
                c(as.integer(inp$N),as.integer(inp$T),3L)))
    fail("Private joint assay tensor has unexpected dimensions")

  death <- as.integer(inp$death_q)
  same <- which(death>0L & vapply(seq_len(inp$N),function(i){
    death[i]>0L && inp$y[i,death[i]]==1L
  },logical(1)))
  other <- which(death>0L & !seq_len(inp$N)%in%same)
  censored <- which(death==0L)
  chosen <- unique(c(head(same,6L),head(other,6L),
                     head(censored,8L)))
  if(length(same)!=79L||length(chosen)<15L)
    fail("Unexpected fatal/live source event representation")

  cat("\n=== JOINT LIKELIHOOD: REAL-HISTORY ENGINEERING CHECK ===\n")
  cat("Reference HMM on",length(chosen),
      "selected badgers (same-quarter death, other death, censored).\n")
  cat("Original assay log likelihoods are used without altering values.\n")
  cat("Illustrative prior/parameter values ONLY; NO FIT.\n")
  demo_entry <- c(.92,.076,.004)
  p_capture <- .23
  rho_recovery <- .4
  rows <- vector("list",length(chosen))
  for(j in seq_along(chosen)) {
    i <- chosen[j]
    qs <- seq.int(as.integer(inp$first_q[i]),as.integer(inp$T))
    K <- length(qs)
    trans <- array(0,c(4L,4L,as.integer(K)))
    sex_male <- as.integer(inp$sex[i]=="MALE")
    for(t in seq_len(K)) {
      age <- (qs[t]-inp$birth_q[i])*.25
      trans[,,t] <- make_covariate_transition(
        age,sex_male,inp$hom_z[i],
        frailty=zero_re,calendar=zero_cal,
        theta=demo)
    }
    em <- matrix(inp$log_diagnostic_emission[i,qs,],
                 nrow=K,ncol=3L)
    ll <- forward_joint_reference(
      trans,as.integer(inp$y[i,qs]),
      as.integer(inp$dead[i,qs]),
      rep(p_capture,K),em,
      demo_entry,rho_recovery)
    if(!is.finite(ll))
      fail("Likelihood not finite for selected individual ",i,
           " (private tattoo omitted)")
    rows[[j]] <- data.frame(
      category=if(i%in%same)"LIVE_AND_DEAD_SAME_QUARTER"
         else if(i%in%other)"RECOVERED_DEATH"
         else "NO_RECORDED_DEATH",
      n_quarters=K,
      n_captures=sum(inp$y[i,qs]==1L),
      n_tests=sum(inp$n_diagnostic_tests[i,qs]),
      ll=ll)
  }
  checks <- do.call(rbind,rows)
  print(stats::aggregate(
    cbind(n_quarters,n_captures,n_tests)~category,
    checks,FUN=sum),row.names=FALSE)
  cat("All",nrow(checks),"illustrative real-history likelihoods: FINITE\n")
  cat("This tests event/observation machinery only.",
      " It is NOT a genotype or mortality estimate.\n")
  cat("JOINT REFERENCE LIKELIHOOD ENGINEERING CHECK: PASS\n")
  if(mode=="compile"){
    if(!requireNamespace("nimble",quietly=TRUE))
      fail("The single joint model requires the installed NIMBLE",
           " package. Ask administrator to install package nimble.")
    suppressPackageStartupMessages(library(nimble))
    source("R/joint_sei_survival_forward_nimble.R",local=TRUE)
    cat("Compiling single-model 4-state SEI + survival NIMBLE kernel...\n")
    compiled <- tryCatch(
      compileNimble(joint_sei_survival_ll),
      error=function(e) {
        try(nimble::printErrors(),silent=TRUE)
        fail("NIMBLE compiled joint kernel failed: ",
             conditionMessage(e))
      })
    betaSE <- unname(demo$SE)
    betaEI <- unname(demo$EI)
    betaD <- unname(demo$D)
    gm <- unname(demo$gm)
    frailty <- as.numeric(c(zero_re,0))
    alpha <- qlogis(p_capture)
    calendarSE <- rep(0,inp$T)
    calendarEI <- rep(0,inp$T)
    per <- rep(0,10L)
    sea <- rep(0,4L)
    period_index <- as.numeric((seq_len(inp$T)-1L)%/%20L+1L)
    season_index <- as.numeric((seq_len(inp$T)-1L)%%4L+1L)
    # This is a simple nimbleFunction without setup(). Both compiled
    # and interpreted objects are directly callable R functions,
    # not reference-class objects exposing $run.
    # See NIMBLE manual Chapter 12 (simple nimbleFunctions).
    if(!is.function(compiled)||!is.function(joint_sei_survival_ll))
      fail("NIMBLE simple function interface is unexpectedly non-callable")
    delta <- numeric(length(chosen))
    for(j in seq_along(chosen)){
      i <- chosen[j]
      qs <- seq.int(inp$first_q[i],inp$T)
      em <- matrix(inp$log_diagnostic_emission[i,,],
                   nrow=inp$T,ncol=3L)
      args <- list(
        y=as.numeric(inp$y[i,]),
        first=as.integer(inp$first_q[i]),
        death=as.integer(ifelse(is.na(inp$death_q[i]),
                               0L,inp$death_q[i])),
        birth=as.integer(inp$birth_q[i]),
        em=em,hom=as.numeric(inp$hom_z[i]),
        male=as.numeric(inp$sex[i]=="MALE"),
        entry=demo_entry,
        betaSE=betaSE,betaEI=betaEI,betaD=betaD,gm=gm,
        frailty=frailty,
        calSE=calendarSE,calEI=calendarEI,
        alphaP=alpha,periodEff=per,seasonEff=sea,
        periodIndex=period_index,seasonIndex=season_index,
        rho=rho_recovery)
      lc <- tryCatch(do.call(compiled,args),
        error=function(e)
          fail("Compiled function invocation failed: ",
               conditionMessage(e)))
      # Exact same inputs in the interpreted nimbleFunction.
      lr <- do.call(joint_sei_survival_ll,args)
      # Compare pure R reference likelihood (all 4 Markov states).
      tr <- array(0,c(4L,4L,length(qs)))
      for(t in seq_along(qs))
        tr[,,t] <- make_covariate_transition(
          (qs[t]-inp$birth_q[i])*.25,
          as.integer(inp$sex[i]=="MALE"),
          inp$hom_z[i],zero_re,zero_cal,demo)
      lreference <- forward_joint_reference(
        tr,as.integer(inp$y[i,qs]),
        as.integer(inp$dead[i,qs]),
        rep(p_capture,length(qs)),
        matrix(inp$log_diagnostic_emission[i,qs,],
               nrow=length(qs),ncol=3L),
        demo_entry,rho_recovery)
      if(!all(is.finite(c(lc,lr,lreference))))
        fail("Nonfinite joint kernel on tested badger #",j)
      delta[j] <- max(abs(c(lc-lr,lc-lreference,lr-lreference)))
      if(delta[j]>1e-7)
        fail("Compiled/R/reference likelihood mismatch on badger #",
             j," (max difference ",signif(delta[j],4),")")
    }
    cat("Compiled NIMBLE vs interpreted NIMBLE vs independent",
        " reference:",length(chosen),"real histories PASS\n")
    cat("Maximum absolute log-likelihood difference:",
        signif(max(delta),5),"\n")
    cat("Conditional on first capture; 79 same-quarter",
        " capture/death records retained.\n")
    cat("NOT YET A BAYESIAN FIT: correlated frailty priors,",
        " smoothing and posterior sampling are still to implement.\n")
    cat("SINGLE JOINT MODEL COMPILED LIKELIHOOD GATE: PASS\n")
  }
  if(mode%in%c("hierarchical_compile","hierarchical_benchmark",
               "full")){
    if(!requireNamespace("nimble",quietly=TRUE))
      fail("NIMBLE is required for the joint model hierarchy")
    suppressPackageStartupMessages(library(nimble))
    source("R/joint_sei_survival_forward_nimble_removal.R",local=TRUE)
    source("R/joint_sei_survival_hierarchy_nimble_v4.R",local=TRUE)

    # The SAME intended scientific model will be used at full N=1603.
    # Only the number of individuals used for this engineering check
    # changes. The scientific regression and prior structure do not.
    full_benchmark <- mode=="hierarchical_benchmark" &&
      identical(Sys.getenv("BENCHMARK_FULL","0"),"1")
    sel <- if(mode=="full"||full_benchmark)
      seq_len(inp$N) else chosen
    n <- length(sel)
    TT <- as.integer(inp$T)
    if(mode=="full" && n!=1603L)
      fail("Full run must include all 1603 genotyped badgers")
    pid <- as.numeric((seq_len(TT)-1L)%/%20L+1L)
    sid <- as.numeric((seq_len(TT)-1L)%%4L+1L)
    yi <- as.integer((seq_len(TT)-1L)%/%20L+1L)
    source("R/joint_removals.R",local=TRUE)
    removed <- joint_removed(inp)
    log_em <- inp$log_diagnostic_emission[sel,,,drop=FALSE]
    dim(log_em) <- c(n,TT,3L)
    if(!identical(dim(log_em),c(n,TT,3L)))
      fail("Hierarchical input assay dimensions incorrect")
    const <- list(
      N=as.integer(n),T=TT,
      y=matrix(as.integer(inp$y[sel,,drop=FALSE]),nrow=n),
      first=as.integer(inp$first_q[sel]),
      death=as.integer(ifelse(is.na(inp$death_q[sel]),0L,
                               inp$death_q[sel])*ifelse(removed[sel],-1L,1L)),
      birth=as.integer(inp$birth_q[sel]),
      em=log_em,
      hom=as.numeric(inp$hom_z[sel]),
      male=as.numeric(inp$sex[sel]=="MALE"),
      yearIndex=as.integer(yi),
      periodIndex=pid,seasonIndex=sid)
    model_data <- list(zeros=rep(0L,n))
    if(any(const$first<1L|const$first>TT)||
       any(const$first<const$birth)||
       any(!is.finite(const$hom))||
       any(!const$male%in%c(0L,1L)))
      fail("Invalid index, age, sex or genotype in hierarchical data")

    # Initialization values only; NOT inferred population parameters.
    init <- list(
      h4=.25,f4=.8,bGM=.18,
      betaSE=c(qlogis(.035),0,0,0),
      betaEI=c(qlogis(.05),0,0,0),
      betaD=rep(0,4),
      entryBeta=c(qlogis(.1),0,0,0,0),
      entryRho=.5,
      sigmaYearSE=.25,sigmaYearEI=.25,
      zYearSE=c(NA_real_,rep(0,9L)),
      zYearEI=c(NA_real_,rep(0,9L)),
      sigmaSE=.3,sigmaEI=.3,load21=0,
      sigmaCapture=.4,
      alphaP=qlogis(.23),rho=.4,
      periodRaw=rep(0,10L),
      seasonEff=c(NA_real_,rep(0,3L)),
      z=cbind(0,0,NA_real_,rep(0,n)))
    # Chain-specific initial values, using the SAME model and priors.
    # Different starting positions and RNG states help reveal poor
    # mixing, unlike identical all-zero latent frailty initializations.
    if(mode=="full"){
      full_chain <- suppressWarnings(as.integer(
        Sys.getenv("CHAIN_ID","1")))
      if(is.na(full_chain)||full_chain<1L||full_chain>99L)
        fail("CHAIN_ID must be 1..99")
      init_seed <- 20261117L+full_chain*13007L
      set.seed(init_seed)
      init$h4 <- exp(rnorm(1,log(.25),.10))
      init$f4 <- plogis(rnorm(1,qlogis(.8),.20))
      init$bGM <- exp(rnorm(1,log(.18),.12))
      init$betaSE[2:4] <- rnorm(3,0,.12)
      init$betaEI[2:4] <- rnorm(3,0,.12)
      init$betaD <- rnorm(4,0,.12)
      init$entryBeta[2:5] <- rnorm(4,0,.12)
      init$entryRho <- plogis(rnorm(1,0,.2))
      init$z <- matrix(rnorm(n*4L,0,.15),nrow=n,ncol=4L)
      init$z[,3] <- NA_real_
      init$zYearSE <- c(NA,rnorm(9,0,.3))
      init$zYearEI <- c(NA,rnorm(9,0,.3))
      init$periodRaw <- rnorm(10,0,.08)
      init$seasonEff[2:4] <- rnorm(3,0,.08)
      init$alphaP <- qlogis(.23)+rnorm(1,0,.15)
      init$rho <- plogis(rnorm(1,qlogis(.4),.12))
    }
    cat("\n=== SINGLE BAYESIAN JOINT INFECTION + SURVIVAL MODEL ===\n")
    cat("Mode:",mode,
        "; subset of genotyped badgers:",n,"/",inp$N,"\n")
    cat("Sex, dynamic age, homozygosity, S/E/I infection state:",
        " all active.\n")
    cat("v4: correlated frailties S->E, E->I; capture frailty;",
        " NO mortality frailty.\n")
    cat("v4: v3 + removals censored (death=-q).\n")
    cat("Same-quarter live capture/death: retained in likelihood.\n")
    cat("Sampler configuration:",sampler_scheme,
        "(same biological model and priors).\n")
    cat("Fixed six-assay accuracy; first live capture conditioned upon.\n")
    cat("Defining Bayesian hierarchy and evaluating initial log density...\n")
    model <- tryCatch(
      nimbleModel(
        joint_sei_survival_model_code,
        constants=const,data=model_data,inits=init,
        check=TRUE,calculate=FALSE),
      error=function(e)
        fail("NIMBLE joint hierarchy setup failed: ",
             conditionMessage(e)))
    lpR <- model$calculate()
    if(!is.finite(lpR))
      fail("Nonfinite initial Bayesian model log probability: ",
           lpR,"; inspect model initialization")
    cat("Interpreted joint posterior log density:",signif(lpR,9),"\n")
    cat("Compiling hierarchical model for ",n,
        " genotyped badgers...\n",sep="")
    cm <- tryCatch(
      compileNimble(model),
      error=function(e){
        try(printErrors(),silent=TRUE)
        fail("NIMBLE hierarchical model compilation failed: ",
             conditionMessage(e))})
    lpC <- cm$calculate()
    if(!is.finite(lpC)||abs(lpR-lpC)>1e-5)
      fail("Compiled hierarchy log probability disagrees with R: ",
           signif(lpR,8)," vs ",signif(lpC,8))
    cat("Compiled joint posterior log density:",
        signif(lpC,9),"\n")
    cat("Difference compiled vs interpreted:",
        signif(abs(lpR-lpC),5),"\n")
    # Require zero diagnostic/entry data mismatch and identical
    # per-badger conditional HMM likelihoods.
    llR <- as.numeric(model$ll)
    llC <- as.numeric(cm$ll)
    if(length(llR)!=n||length(llC)!=n||
       any(!is.finite(llR))||any(!is.finite(llC))||
       max(abs(llR-llC))>1e-7)
      fail("Compiled per-badger HMM likelihoods disagree with R")
    cat("Per-badger forward likelihoods matched for all",
        n,"individuals.\n")
    # Full runs are isolated by chain and immutable; never reuse
    # or overwrite another chain's file or log.
    chain <- suppressWarnings(as.integer(Sys.getenv("CHAIN_ID","1")))
    if(is.na(chain)||chain<1L||chain>99L)
      fail("CHAIN_ID must be an integer in 1..99")
    if(mode=="full") {
      run_dir <- file.path(outdir,Sys.getenv("RUN_NAME","full_joint_v4_removals"),
                           sprintf("chain_%02d",chain))
      dir.create(run_dir,recursive=TRUE,showWarnings=FALSE)
      save_path <- file.path(run_dir,"posterior_samples.rds")
      if(file.exists(save_path))
        fail("Full posterior output already exists: ",save_path,
             " (will never overwrite)")
      if(file.exists(file.path(run_dir,"chain_started.rds")))
        fail("This chain was previously launched: ",run_dir,
             "; inspect it before restarting")
      saveRDS(list(
        chain=chain,n=n,
        sampler_scheme=sampler_scheme,
        input_md5=tools::md5sum(validated),
        source_md5=inp$source_md5,
        model_md5=tools::md5sum(c(
          "scripts/joint_sei_survival_model_v4.R",
          "R/joint_sei_survival_hierarchy_nimble_v4.R",
          "R/joint_sei_survival_forward_nimble_removal.R",
        "R/joint_removals.R")),
        started=Sys.time(),status="started",
        caveat=paste(
          "Exploratory conditional mortality model:",
          "EUTHANISED/KILLED deaths censored at removal (v4);",
          "other recovered deaths (RTA, other/unknown) kept as deaths")
        ),file.path(run_dir,"chain_started.rds"))
      cat("Full chain directory:",run_dir,"\n")
      cat("WARNING: known deaths may include managed removals;",
          " mortality inference is PROVISIONAL.\n")
    } else {
      save_path <- file.path(outdir,
          paste0("joint_hierarchical_v4_",mode,
                 if(sampler_scheme=="blocked")"_blocked" else
                   if(sampler_scheme=="af_slice")"_af_slice" else
                     if(sampler_scheme=="af_slice_calendar")"_af_slice_calendar"
                     else "",
                 if(full_benchmark)"_full_cohort" else "",
                 "_gate.rds"))
      if(file.exists(save_path) &&
         Sys.getenv("OVERWRITE","0")!="1")
        fail("Previous hierarchical gate exists: ",save_path,
             "; use OVERWRITE=1 only after reviewing old output")
    }
    result <- list(
      mode=mode,n=n,total_cohort=inp$N,
      log_posterior_R=lpR,log_posterior_compiled=lpC,
      max_individual_ll_delta=max(abs(llR-llC)),
      input_source_md5=inp$source_md5,
      input_md5=tools::md5sum(validated),
      model_code_md5=tools::md5sum(c(
        "scripts/joint_sei_survival_model_v4.R",
        "R/joint_sei_survival_hierarchy_nimble_v4.R",
        "R/joint_sei_survival_forward_nimble_removal.R",
        "R/joint_removals.R")),
      note="Engineering gate only; pilot subset; not a fitted posterior")

    if(mode%in%c("hierarchical_benchmark","full")){
      cat("Constructing single-model MCMC samplers...\n")
      monitor <- c("betaSE","betaEI","betaD",
         "entryBeta","entryRho","h4","f4","bGM",
         "sigmaSE","sigmaEI","load21","sigmaCapture",
         "sigmaYearSE","sigmaYearEI",
         "yearSE","yearEI","periodEff","seasonEff",
         "alphaP","rho")
      # z (S->E, E->I, capture) kept every 50th iteration for the
      # post-last-capture I-onset check (forward-filter backward-sample).
      monitor2 <- sprintf("z[1:%d, %d]",n,c(1L,2L,4L))
      thin2 <- 50L
      sampler <- configureMCMC(model,monitors=monitor,thin=1,
                               monitors2=monitor2,thin2=thin2)
      block_metadata <- NULL
      if(sampler_scheme=="blocked"){
        source("R/joint_sei_sampler_blocking.R",local=TRUE)
        block_metadata <- configure_joint_sei_blocked(
          sampler,tools::md5sum(validated))
        cat("Blocked proposal tuning based on original 3 chains;",
            " no changes to model likelihood/priors.\n")
      } else if(sampler_scheme%in%c("af_slice","af_slice_calendar")){
        source("R/joint_sei_sampler_af_slice.R",local=TRUE)
        mortality_metadata <- configure_joint_sei_AF_slice(
          sampler,tools::md5sum(validated))
        block_metadata <- mortality_metadata
        if(sampler_scheme=="af_slice_calendar"){
          source("R/joint_sei_sampler_calendar_af_slice.R",local=TRUE)
          calendar_metadata <- configure_joint_sei_calendar_AF_slice(
            sampler,tools::md5sum(validated))
          block_metadata <- list(
            scheme=sampler_scheme,
            mortality=mortality_metadata,
            calendar=calendar_metadata)
        }
        cat("AF_slice proposal directions adapt within burn-in;",
            " unchanged scientific likelihood and priors.\n")
      }
      mcmc <- buildMCMC(sampler)
      cmcmc <- tryCatch(
        compileNimble(mcmc,project=model),
        error=function(e){
          try(printErrors(),silent=TRUE)
          fail("NIMBLE MCMC compilation failed: ",
               conditionMessage(e))})
      # Distinct RNG streams by chain, with SAME intended model/prior.
      seed <- as.integer(20261117L + chain*13007L)
      set.seed(seed)
      if(mode=="hierarchical_benchmark"){
        sweeps <- suppressWarnings(as.integer(
          Sys.getenv("BENCHMARK_SWEEPS","2")))
        if(is.na(sweeps)||sweeps<2L||sweeps>100L)
          fail("BENCHMARK_SWEEPS must be 2..100")
        cat("MCMC timing: cohort N =",n,
            "; sampler =",sampler_scheme,
            "; test sweeps =",sweeps,"\n")
        bench <- system.time(
          draws <- runMCMC(cmcmc,niter=sweeps,nburnin=0L,
              thin=1L,nchains=1L,progressBar=FALSE,
              samplesAsCodaMCMC=FALSE))
        result$sampler_scheme <- sampler_scheme
        result$block_metadata <- block_metadata
        result$benchmark_sweeps <- sweeps
        result$seconds_for_benchmark_sweeps <- unname(bench[["elapsed"]])
        result$seconds_per_sweep <- unname(bench[["elapsed"]])/sweeps
        result$parameter_columns <- colnames(draws$samples)
        cat("Trial MCMC sweeps:",sweeps,
            "; total sampling seconds:",
            signif(result$seconds_for_benchmark_sweeps,6),
            "; seconds per sweep:",
            signif(result$seconds_per_sweep,6),"\n")
        cat("Timing gate only; no biological inference or convergence.\n")
      } else {
        niter <- suppressWarnings(as.integer(
          Sys.getenv("NITER","8000")))
        nburn <- suppressWarnings(as.integer(
          Sys.getenv("NBURN","3000")))
        thin <- suppressWarnings(as.integer(
          Sys.getenv("THIN","5")))
        if(anyNA(c(niter,nburn,thin))||
           niter<2000L||nburn<500L||nburn>=niter||
           thin<1L||((niter-nburn)%/%thin)<100L)
          fail("Invalid full MCMC length/burnin/thinning")
        cat("\n=== FULL COHORT SINGLE JOINT MODEL MCMC ===\n")
        cat("N =",n,"; chain =",chain,"; seed =",seed,
            "; iterations =",niter,"; burnin =",nburn,
            "; thin =",thin,"\n")
        cat("Monitor:",paste(monitor,collapse=", "),"\n")
        cat("Progress follows NIMBLE iteration reporting.\n")
        flush.console()
        t0 <- Sys.time()
        if((niter-nburn)%%thin2!=0L||thin2%%thin!=0L)
          fail("NITER-NBURN must be a multiple of 50 and of THIN")
        run_out <- runMCMC(
          cmcmc,niter=niter,nburnin=nburn,thin=thin,thin2=thin2,
          nchains=1,progressBar=TRUE,
          samplesAsCodaMCMC=FALSE)
        posterior <- run_out$samples
        posterior_z <- run_out$samples2
        wall <- as.numeric(difftime(Sys.time(),t0,units="secs"))
        if(!is.matrix(posterior) ||nrow(posterior)<100L||
           any(!is.finite(posterior)))
          fail("Posterior samples absent or nonfinite")
        result <- list(
          chain=chain,seed=seed,
          sampler_scheme=sampler_scheme,
          block_metadata=block_metadata,
          N=n,iterations=niter,burnin=nburn,thin=thin,
          saved_draws=nrow(posterior),
          elapsed_sampling_seconds=wall,
          model_parameters=colnames(posterior),
          samples=posterior,
          samples_z=posterior_z,thin_z=thin2,
          source_md5=inp$source_md5,
          input_md5=tools::md5sum(validated),
          model_code_md5=tools::md5sum(c(
            "scripts/joint_sei_survival_model_v4.R",
            "R/joint_sei_survival_hierarchy_nimble_v4.R",
            "R/joint_sei_survival_forward_nimble_removal.R",
        "R/joint_removals.R")),
          status="UNVALIDATED_PROVISIONAL",
          mortality_caveat=paste(
            "Natural vs intentional deaths unresolved;",
            "capture effort is a proxy, fixed assay accuracy.")
        )
        # Protect completed result: no overwrites.
        if(file.exists(save_path))
          fail("Posterior file unexpectedly exists; refusing overwrite")
        saveRDS(result,save_path,compress="xz")
        done_path <- file.path(run_dir,"chain_finished.txt")
        writeLines(c(
          paste("chain",chain),paste("N",n),
          paste("niter",niter),paste("burnin",nburn),
          paste("thin",thin),paste("draws",nrow(posterior)),
          paste("model","v4_removals"),
          paste("sampling_seconds",signif(wall,6)),
          paste("sampler_scheme",sampler_scheme),
          paste("status","UNVALIDATED_PROVISIONAL")),
          done_path)
        cat("Saved complete chain posterior:",save_path,"\n")
        cat("Saved chain completion marker:",done_path,"\n")
        cat("FULL JOINT MODEL MCMC CHAIN COMPLETED (PROVISIONAL).\n")
      }
    }
    if(mode!="full"){
      saveRDS(result,save_path)
      cat("Saved:",save_path,"\n")
      cat("SINGLE JOINT HIERARCHICAL MODEL COMPILATION GATE: PASS\n")
    }
  }
  quit(save="no",status=0L)
}

x <- readRDS(files["survival"])
required_x <- c("id","N","T","y","first_q","birth_q","death_q",
                "hom","sex")
if(!all(required_x%in%names(x)))
  fail("Survival input missing: ",
       paste(setdiff(required_x,names(x)),collapse=", "))
N <- as.integer(x$N)
T <- as.integer(x$T)
if(N!=1603L||T!=200L||
   !is.matrix(x$y)||!identical(dim(x$y),c(N,T)))
  fail("Survival input expected 1,603 individuals x 200 quarters")
cohort <- normal(x$id)
if(anyNA(cohort)||anyDuplicated(cohort))
  fail("Missing/duplicate survival tattoo IDs")
if(anyNA(x$y)||!all(x$y%in%c(-1L,0L,1L)))
  fail("Survival observation history contains unrecognized codes")
if(length(x$birth_q)!=N||length(x$first_q)!=N||
   length(x$death_q)!=N||length(x$hom)!=N||
   length(x$sex)!=N)
  fail("Survival demographic and history fields not N-length")
birth <- as.integer(x$birth_q)
first <- as.integer(x$first_q)
death <- as.integer(ifelse(is.na(x$death_q),0L,x$death_q))
if(anyNA(birth)||anyNA(first)||anyNA(death)||
   any(first<1L|first>T)||
   any(first<birth)||
   any(death<0L|death>T)||
   any(!is.finite(as.numeric(x$hom)))||
   sd(as.numeric(x$hom))<=0||
   any(!normal(x$sex)%in%c("MALE","FEMALE")))
  fail("Invalid age, first capture, death or genetic input")
cat("Survival individuals:",N,"quarters:",T,"\n")
cat("Known recorded death quarters:",sum(death>0L),"\n")
cat("Survival observations: captured=",sum(x$y==1L),
    "monitored-not-captured=",sum(x$y==0L),
    "not monitored=",sum(x$y==-1L),"\n")

ids <- utils::read.csv(files["crosswalk"],
          colClasses=c(tattoo="character"))
cand <- utils::read.csv(files["diagnostic"],
                       na.strings=c("","NA"),
                       check.names=FALSE)
required_id <- c("id","tattoo")
required_can <- c("id","time","captured",paste0("test",1:6))
if(!all(required_id%in%names(ids))||
   !all(required_can%in%names(cand)))
  fail("Diagnostic candidate or source ID crosswalk schema changed")
ids$tattoo <- normal(ids$tattoo)
if(nrow(ids)!=3214L||anyDuplicated(ids$id)||
   anyDuplicated(ids$tattoo)||anyNA(ids[,required_id])||
   anyNA(cand[,c("time","id","captured")])||
   any(cand$captured!=1L)||
   any(cand$time<1L|cand$time>T)||
   any(!cand$id%in%ids$id))
  fail("SEI diagnostic candidate/crosswalk integrity failed")
gi <- match(cohort,ids$tattoo)
if(anyNA(gi)||anyDuplicated(gi))
  fail("Genotyped cohort tattoos not all in SEI source crosswalk")
cand$genetic_idx <- match(cand$id,ids$id[gi])
cand <- cand[!is.na(cand$genetic_idx),,drop=FALSE]
for(nm in paste0("test",1:6))
  if(!all(is.na(cand[[nm]])|cand[[nm]]%in%c(0,1)))
    fail("Invalid diagnostic code in ",nm)
cat("Exact genotyped tattoo joins:",sum(!is.na(gi)),"/",N,"\n")
cat("Matched confirmed-live diagnostic capture records:",nrow(cand),"\n")
cat("Observed six-assay outcomes:",
    sum(!is.na(as.matrix(cand[,paste0("test",1:6),drop=FALSE]))),
    "\n")
n_captures <- tabulate(cand$genetic_idx,nbins=N)
quarters <- vector("list",N)
for(k in seq_len(N))
  quarters[[k]] <- sort(unique(cand$time[cand$genetic_idx==k]))
n_quarters <- lengths(quarters)
cat("Genotyped individuals with 0 / 1 / >=2 confirmed-live",
    "diagnostic quarters:",
    sum(n_quarters==0L),sum(n_quarters==1L),sum(n_quarters>=2L),
    "\n")

# Event-clock reconciliation: original y==1 and Stage76 method-confirmed
# diagnostic capture can differ (e.g., death records and earlier cleaning).
# Report differences, NEVER silently replace or merge the two clocks.
y_good <- which(x$y==1L,arr.ind=TRUE)
ykey <- paste(y_good[,1],y_good[,2],sep=":")
ckey <- unique(paste(cand$genetic_idx,cand$time,sep=":"))
source_only <- setdiff(ykey,ckey)
diagnostic_only <- setdiff(ckey,ykey)
cat("Survival live-encounter quarter keys:",length(unique(ykey)),"\n")
cat("Stage76 confirmed live capture quarter keys:",length(ckey),"\n")
cat("Survival-live quarters WITHOUT Stage76 capture:",length(source_only),"\n")
cat("Stage76 confirmed-live quarters WITHOUT survival y=1:",
    length(diagnostic_only),"\n")

record <- data.frame(
  tattoo=cohort,
  genotyped_index=seq_len(N),
  candidate_capture_records=n_captures,
  candidate_capture_quarters=n_quarters,
  survival_live_quarters=rowSums(x$y==1L),
  survival_monitored_no_capture=rowSums(x$y==0L),
  survival_not_monitored=rowSums(x$y==-1L),
  first_q=first,birth_q=birth,death_q=death,
  death_no_after_first=death==0L | death>=first,
  first_capture_marked_live=x$y[cbind(seq_len(N),first)]==1L,
  live_after_recorded_death=vapply(seq_len(N),function(k)
      death[k]>0L&&death[k]<T&&
      any(x$y[k,seq.int(death[k]+1L,T)]==1L),logical(1)),
  diag_before_birth=vapply(seq_len(N),function(k)
      length(quarters[[k]])>0L &&
      any(quarters[[k]]<birth[k]),logical(1)),
  diag_after_recorded_death=vapply(seq_len(N),function(k)
      death[k]>0L&&length(quarters[[k]])>0L &&
      any(quarters[[k]]>death[k]),logical(1)),
  diag_same_quarter_death=vapply(seq_len(N),function(k)
      death[k]>0L&&length(quarters[[k]])>0L &&
      death[k]%in%quarters[[k]],logical(1)))
private <- file.path(outdir,"joint_PRIVATE_cohort_input_audit.csv")
utils::write.csv(record,private,row.names=FALSE)
cat("First capture not marked y=1:",
    sum(!record$first_capture_marked_live),"\n")
cat("Death before first capture:",
    sum(!record$death_no_after_first),"\n")
cat("Survival y=1 later than recorded death quarter:",
    sum(record$live_after_recorded_death),"\n")
cat("Stage76 diagnostic captures before derived birth quarter:",
    sum(record$diag_before_birth),"\n")
cat("Stage76 diagnostic captures after recorded death quarter:",
    sum(record$diag_after_recorded_death),"\n")
cat("Stage76 diagnostic captures SAME quarter as death:",
    sum(record$diag_same_quarter_death),"\n")

# We need paired source-origin encounter keys to decide which encounters
# are appropriate for mortality/detection and how same-quarter death
# observations should enter the likelihood. Export PRIVATE difference
# lists for local reconciliation without disclosing tattoos in logs.
key_table <- function(key,reason) {
  if(!length(key))return(data.frame(
    genetic_idx=integer(),time=integer(),reason=character()))
  parts <- strsplit(key,":",fixed=TRUE)
  data.frame(genetic_idx=as.integer(vapply(parts,function(v) v[[1L]],character(1))),
             time=as.integer(vapply(parts,function(v) v[[2L]],character(1))),
             reason=rep(reason,length(key)))
}
dif <- rbind(
  key_table(source_only,"SURVIVAL_Y_LIVE_NOT_IN_STAGE76"),
  key_table(diagnostic_only,"STAGE76_LIVE_NOT_IN_SURVIVAL_Y"))
if(nrow(dif)>0L){
  dif$tattoo <- cohort[dif$genetic_idx]
  utils::write.csv(dif,
     file.path(outdir,"joint_PRIVATE_discordant_live_quarters.csv"),
     row.names=FALSE)
}

audit_summary <- data.frame(
  total_genotyped=N,
  exact_genotyped_sei_tattoos=sum(!is.na(gi)),
  confirmed_live_diagnostic_records=nrow(cand),
  genotyped_no_diagnostic_quarter=sum(n_quarters==0L),
  genotyped_one_diagnostic_quarter=sum(n_quarters==1L),
  genotyped_two_plus_diagnostic_quarters=sum(n_quarters>=2L),
  survival_live_quarters_not_stage76=length(source_only),
  stage76_live_quarters_not_survival_y=length(diagnostic_only),
  same_quarter_death_and_diagnostic_capture=
      sum(record$diag_same_quarter_death),
  death_before_first_capture=sum(!record$death_no_after_first),
  live_after_recorded_death=sum(record$live_after_recorded_death),
  first_capture_y_not_live=sum(!record$first_capture_marked_live))
utils::write.csv(audit_summary,
  file.path(outdir,"joint_input_readiness_summary.csv"),row.names=FALSE)
cat("Saved aggregated summary: ",
    file.path(outdir,"joint_input_readiness_summary.csv"),"\n",sep="")
cat("Saved private provenance reconciliation: ",private,"\n",sep="")
# --------------------------------------------------------------------------
# Compile a private, source-fingerprinted audit input for this SAME model.
# Preserve all 1,603 genotyped individuals, including zero-test animals.
# This is not yet labelled model-ready because the y vs candidate clocks
# may differ. All observed assays retain their exact q and source ID.
# --------------------------------------------------------------------------
f <- readRDS(files["sei_fit"])
if(!is.list(f) || !is.data.frame(f$sesp)||
   !all(c("test_name","Se","Sp")%in%names(f$sesp)) ||
   !identical(as.character(f$sesp$test_name),
       c("ELISA_BELT","ELISA_INDIRECT","Culture",
         "DPP","IGRA","StatPak")))
  fail("Original diagnostic assay parameter table/order changed")
Se <- as.numeric(f$sesp$Se); Sp <- as.numeric(f$sesp$Sp)
if(anyNA(Se)||anyNA(Sp)||any(Se<=0|Se>=1)||
   any(Sp<=0|Sp>=1))
  fail("Original diagnostic sensitivity or specificity invalid")
draws <- f$draws
if(!is.data.frame(draws)||nrow(draws)<10L)
  fail("Source original SEI exposed-state draws unavailable")
nms <- names(draws)
mapped <- vapply(1:6,function(k){
  names_allowed <- c(paste0("e_mult.",k,"."),
                     paste0("e_mult.",k),
                     paste0("e_mult[",k,"]"))
  which(nms %in% names_allowed)[1L]
},integer(1))
if(anyNA(mapped)||anyDuplicated(mapped))
  fail("Cannot identify six source posterior exposed multipliers")
m <- vapply(draws[mapped],function(z){
  if(!is.numeric(z)||anyNA(z)||any(z<0|z>1))
    fail("Invalid original exposed multiplier posterior draws")
  mean(z)
},numeric(1))
pos <- rbind(1-Sp, m*Se, Se)
rownames(pos) <- c("S","E","I")
colnames(pos) <- paste0("test",1:6)
if(any(pos<=0|pos>=1))
  fail("Invalid probabilities in diagnostic emission matrix")
diag_em <- array(0,c(N,T,3L),
                 dimnames=list(cohort,NULL,c("S","E","I")))
diag_test_count <- matrix(0L,N,T)
for(j in seq_len(nrow(cand))){
  k <- cand$genetic_idx[j]
  q <- cand$time[j]
  for(a in 1:6) {
    z <- cand[[paste0("test",a)]][j]
    if(is.na(z))next
    if(z==1L)
      diag_em[k,q,] <- diag_em[k,q,]+log(pos[,a])
    else
      diag_em[k,q,] <- diag_em[k,q,]+log1p(-pos[,a])
    diag_test_count[k,q] <- diag_test_count[k,q]+1L
  }
}
# No test in a quarter has log likelihood 0 under ALL disease states.
untested <- diag_test_count==0L
if(any(diag_em[,,1L][untested]!=0) ||
   any(diag_em[,,2L][untested]!=0) ||
   any(diag_em[,,3L][untested]!=0))
  fail("Untested quarter inadvertently became a diagnostic observation")
if(sum(diag_test_count)!=sum(!is.na(as.matrix(
  cand[,paste0("test",1:6),drop=FALSE]))))
  fail("Diagnostic log emission aggregation dropped or duplicated tests")
cat("Raw diagnostic log emission tensor:",N,"x",T,"x 3 disease states\n")
cat("Distinct genotyped test-observed badger-quarters:",
    sum(diag_test_count>0L),"\n")
cat("Original assay Se/Sp/multipliers retained as FIXED primary model inputs.\n")

# Historical birth/first capture/death are quarter-indexed.
# Entry-state and quarterly death hazards will be evaluated using age(q).
# Keep raw hom, and separately store the ORIGINAL 1603-cohort scale.
hom <- as.numeric(x$hom)
hom_center <- mean(hom); hom_scale <- sd(hom)
hom_z <- (hom-hom_center)/hom_scale
if(!all(is.finite(hom_z))||abs(mean(hom_z))>1e-10)
  fail("Cannot standardize genotype in original full cohort")
if(!is.matrix(x$dead)||!identical(dim(x$dead),c(N,T))||
   anyNA(x$dead)||any(!x$dead%in%c(0L,1L)))
  fail("Recorded-death matrix absent or has invalid shape/codes")
for(k in seq_len(N)) {
  dq <- death[k]
  if(sum(x$dead[k,])!=as.integer(dq>0L))
    fail("Death matrix/quarter discrepancy in survival source")
  if(dq>0L&&x$dead[k,dq]!=1L)
    fail("Recorded-death quarter disagrees with death matrix")
}
# Known death dates and within-quarter capture chronology still require
# detailed source reconciliation if any same-quarter events occur.
conflicts <- sum(record$live_after_recorded_death|
                 !record$death_no_after_first|
                 record$diag_after_recorded_death|
                 record$diag_before_birth|
                 !record$first_capture_marked_live)
diff_quarters <- length(source_only)+length(diagnostic_only)
cat("Hard event chronology problems:",conflicts,"\n")
cat("Quarter-clock discrepancies:",diff_quarters,"\n")
cat("Important: y=0 uses a monitoring PROXY, not true individual effort.\n")

private_input <- file.path(
  outdir,"joint_PRIVATE_UNRECONCILED_input.rds")
if(file.exists(private_input) && Sys.getenv("OVERWRITE","0")!="1")
  fail("Private joint input already exists; audit the previous results",
       " before intentionally replacing with OVERWRITE=1")
saveRDS(list(
  audit_only=TRUE,
  ready_for_joint_fit=FALSE,
  source_md5=tools::md5sum(files),
  id=cohort,N=N,T=T,
  y=x$y,dead=x$dead,
  first_q=first,birth_q=birth,death_q=death,
  hom=hom,hom_center=hom_center,hom_scale=hom_scale,
  hom_z=hom_z,sex=normal(x$sex),
  log_diagnostic_emission=diag_em,
  n_diagnostic_tests=diag_test_count,
  source_candidate=cand[,c("genetic_idx","time","id",
              "captured",paste0("test",1:6)),drop=FALSE],
  assay=list(Se=Se,Sp=Sp,exposed_multiplier=m,
             positive_probability=pos),
  source_live_only_quarters=source_only,
  diagnostics_only_quarters=diagnostic_only,
  hard_chronology_conflicts=conflicts,
  clock_discrepancies=diff_quarters,
  same_quarter_capture_death=sum(
      record$diag_same_quarter_death),
  monitoring_proxy=x$calendar,
  assumptions=c(
    "individual genotype effects estimated per full-genotyped-cohort SD",
    "sex, time-varying age and biological badger frailties required",
    "monitoring proxy is not a measured per-individual trapping effort",
    "diagnostic positivity conditional on S/E/I with fixed assay accuracy",
    "confirmed live diagnostic encounters only (317 uncertain excluded)",
    "source survival and diagnostic encounter clocks not silently merged",
    "death recoveries are distinct from observed live captures",
    "four-state joint model not fitted; discrepancies require reconciliation"
  )
),private_input)
cat("Saved private audit-only joint input:",private_input,"\n")
cat("NO JOINT MODEL FIT YET. NO SOURCE DATA MODIFIED.\n")
cat("INPUT VALIDATION COMPLETE. Inspect mismatches before enabling fit.\n")
