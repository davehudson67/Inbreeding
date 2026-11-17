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
# Next development will add data-specific observation/death processes,
# compiled HMM marginalisation, and posterior/uncertainty estimation IN
# THIS CODEBASE rather than making further standalone SEI variants.
#
# No DB writes; no database reads; all exports private and gitignored.
# ============================================================================
options(stringsAsFactors=FALSE)
mode <- Sys.getenv("MODE","validate")
if(!mode %in% c("validate"))
  stop("Only MODE=validate is implemented. Joint fitting not yet validated.")
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
  data.frame(genetic_idx=as.integer(vapply(parts,\`[\`,character(1),1L)),
             time=as.integer(vapply(parts,\`[\`,character(1),2L)),
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
cat("NO JOINT MODEL FIT YET. NO SOURCE DATA MODIFIED.\n")
cat("INPUT VALIDATION COMPLETE; reconciliation required before final fitting.\n")
