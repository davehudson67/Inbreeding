# ============================================================================
# Stage 83: QUADRATIC AGE AND CALENDAR SEI SENSITIVITY
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
# QUADRATIC AGE/CALENDAR ROBUSTNESS TEST: time-varying hazards
# and exposed-state multipliers from original fitted SEI; no seasonal,
# calendar, age, sex, social-group, survival, trapping or observation model.
# Regularized MAP is NOT a fully Bayesian fit nor a final genetic effect.
#
# Reads NO DB, does NOT edit original data, outputs remain gitignored.
# N_STAGE83 (default 99999) uses ALL eligible genetic individuals; optional subset uses seed 79.
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
cat("\n=== STAGE 83: DIRECT-DIAGNOSTIC SEI GENETIC SMOKE TEST ===\n")
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
n_requested <- suppressWarnings(as.integer(Sys.getenv("N_STAGE83","99999")))
if(is.na(n_requested)||n_requested<40L)
  fail("N_STAGE83 must be an integer >= 40")
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
cat("Selected confirmed live capture records:",total_captures,
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

# ============================================================================
# STAGE 83: adjust both transition hazards for age, sex, calendar time;
# explicitly model homozygosity in INITIAL prevalence too.
#
# Two nested models on identical individuals/capture test histories:
# (A) age+sex+calendar on S->E, E->I and entry prevalence, PLUS
#     entry homozygosity effect, but NO genetic transition effects;
# (B) same nuisance model PLUS separate hom effects on S->E and E->I.
#
# We let starting prevalence vary with homozygosity in BOTH models,
# since infections predating first capture cannot be assigned to the
# subsequent S->E hazard. Otherwise entry-genetic association could
# masquerade as an ongoing infection hazard effect.
#
# Stage83 will add quadratic calendar/age terms to the linear model;
# the model still lacks seasonal/RW1 temporal baselines and source-
# specific test observation/missingness. This is NOT a final model.
# ============================================================================
ref_path <- file.path(base,"stage80_all_eligible_joint_sei.rds")
if(!file.exists(ref_path))fail("Run Stage80 first: ",ref_path)
ref <- readRDS(ref_path)
if(!is.list(ref)||is.null(ref$fits$both$theta)||
   length(ref$fits$both$theta)!=6L ||
   !identical(unname(as.character(ref$source_md5)),
              unname(as.character(tools::md5sum(files)))))
  fail("Stage80 baseline missing, malformed or from different sources")
if(!isTRUE(all.equal(unname(ref$assay_parameters$exposed_multiplier),
                     unname(m),tolerance=1e-10)))
  fail("Stage80 diagnostic sensitivity multiplier inputs changed")

# Covariates are quarter-specific. logit hazards are predicted for EACH
# interval q>first quarter; no retroactive infection attribution.
# Ages use original birth_q, so do not substitute inferred first capture.
for(k in seq_along(obs)){
  gi <- gis[k]
  idx <- c_by_id[[as.character(chosen[k])]]
  q <- seq.int(min(cand$time[idx]),max(cand$time[idx]))
  if(length(q)!=obs[[k]]$nq)
    fail("Stage82 q sequence and Stage76 emissions disagree")
  age_years <- (q-as.integer(x$birth_q[gi]))*.25
  if(any(!is.finite(age_years))||any(age_years<0))
    fail("Nonfinite or pre-birth age in selected cohort")
  age_scaled <- (age_years-4)/4
  cal_scaled <- (1976+(q-1L)%/%4L-2000)/10
  if(!x$sex[gi]%in%c("MALE","FEMALE"))
    fail("Unknown sex in genotyped individual")
  obs[[k]]$age <- age_scaled
  obs[[k]]$calendar <- cal_scaled
  obs[[k]]$male <- as.integer(x$sex[gi]=="MALE")
  obs[[k]]$entry_age <- age_scaled[1L]
  obs[[k]]$entry_calendar <- cal_scaled[1L]
}
cat("\n=== STAGE 83: TIME/SEX ADJUSTED JOINT GENETIC SEI ===\n")
cat("Badgers:",length(chosen),"test results:",total_assays,"\n")
cat("Time-varying transition effects: age (per 4 years),",
    "calendar (per 10 years), sex (male vs female).\n")
cat("Initial prevalence adjusts for age at first capture, calendar",
    "year, sex and homozygosity.\n")

# ============================================================================
# STAGE 83: Curved age/calendar time sensitivity.
# Add a squared age term and squared calendar time term separately to
# S->E hazard, E->I hazard and initial infection prevalence.
# Everything else (including sex, entry hom and fixed Se/Sp) remains
# the same as Stage82. No categorical year effects or random walk yet.
#
# Keep Stage82 as a nested model: quadratic coefficient=0 exactly
# recovers Stage82; thus large genetic shifts indicate sensitivity
# to the choice of age/calendar functional form.
#
# This is an exploratory, penalized MAP robustness fit and not a
# hypothesis test with calibrated uncertainty.
# ============================================================================
adj_path <- sprintf(
  "outputs/sei_transitions/stage82_age_sex_calendar_adjustment_%04d.rds",
  length(chosen))
if(!file.exists(adj_path))
  fail("Run matching cohort-size Stage82 first; missing ",adj_path)
adj <- readRDS(adj_path)
if(!is.list(adj) || is.null(adj$fits$null$theta)||
   is.null(adj$fits$genetics$theta)||
   length(adj$fits$null$theta)!=16L||
   length(adj$fits$genetics$theta)!=16L||
   !identical(as.integer(adj$selected),as.integer(length(chosen)))||
   !identical(unname(as.character(adj$source_md5)),
              unname(as.character(tools::md5sum(files)))))
  fail("Stage82 source inputs/cohort differ or fit incomplete")
if(!identical(unname(as.character(adj$source_stage80_md5)),
              unname(as.character(tools::md5sum(ref_path)))))
  fail("Stage82 is not based on the currently loaded Stage80 fit")
for(k in seq_along(obs)){
  obs[[k]]$age2 <- obs[[k]]$age^2
  obs[[k]]$calendar2 <- obs[[k]]$calendar^2
  obs[[k]]$entry_age2 <- obs[[k]]$entry_age^2
  obs[[k]]$entry_calendar2 <- obs[[k]]$entry_calendar^2
}
cat("\n=== STAGE 83: QUADRATIC AGE/CALENDAR JOINT SEI SENSITIVITY ===\n")
cat("Individuals:",length(chosen),
    "; observed diagnostic outcomes:",total_assays,"\n")
cat("Baseline: Stage82 fitted linear age, sex, calendar and entry hom.\n")
cat("Extra terms: age^2 and calendar^2 for acquisition, progression, entry.\n")

# First 16 parameters exactly match Stage82 layout.
names_theta <- c(
 "logit_SE","logit_EI","logit_entry_pi","logit_entry_rho",
 "beta_hom_SE","beta_hom_EI",
 "beta_age_SE","beta_age_EI","beta_male_SE","beta_male_EI",
 "beta_calendar_SE","beta_calendar_EI",
 "beta_entry_age","beta_entry_male","beta_entry_calendar",
 "beta_entry_hom",
 "beta_age2_SE","beta_age2_EI","beta_calendar2_SE",
 "beta_calendar2_EI","beta_entry_age2","beta_entry_calendar2")
npar <- length(names_theta)
if(!identical(as.character(adj$parameter_names),names_theta[1:16]))
  fail("Stage82 model parameter ordering changed; stop")
start <- setNames(c(as.numeric(adj$fits$genetics$theta),
                    rep(0,6)),names_theta)
lower <- c(-9,-9,-9,-7,rep(-4,12),rep(-2,6))
upper <- c(0,0,0,7,rep(4,12),rep(2,6))
prior_mean <- c(qlogis(.035),qlogis(.05),qlogis(.1),
                qlogis(.5),rep(0,18))
prior_sd <- c(2,2,2,2.5,1.5,1.5,
              1,1,1,1,.8,.8,
              1,1,1,1,
              .65,.65,.55,.55,.65,.55)
stopifnot(length(start)==npar,length(lower)==npar,
          length(upper)==npar,length(prior_mean)==npar,
          length(prior_sd)==npar,all(start>=lower),
          all(start<=upper))
hmm_vary <- function(em,lam,tau,pi,rho) {
  if(!is.matrix(em)||ncol(em)!=3L||nrow(em)!=length(lam)||
     length(tau)!=length(lam))
    fail("Stage83 emission/transition sizes differ")
  a <- c(1-pi,pi*(1-rho),pi*rho)
  LL <- 0
  for(t in seq_len(nrow(em))) {
    if(t>1L) {
      l <- lam[t];u <- tau[t]
      a <- c(a[1L]*(1-l),
             a[1L]*l+a[2L]*(1-u),
             a[2L]*u+a[3L])
    }
    offset <- max(em[t,])
    a <- a*exp(em[t,]-offset)
    zz <- sum(a)
    if(!is.finite(zz)||zz<=0)
      fail("Stage83 nonfinite HMM emission probability")
    LL <- LL+log(zz)+offset
    a <- a/zz
  }
  LL
}
# Reassure ourselves the time-varying recursion reduces to Stage80 when
# all hazards are constant and there are no covariates.
stopifnot(abs(hmm_vary(matrix(0,7,3),rep(.1,7),rep(.2,7),
                       .07,.02))<1e-10)
set.seed(82L)
em_test <- matrix(runif(21,-3,0),7,3)
stopifnot(abs(hmm_vary(em_test,rep(.1,7),rep(.2,7),.07,.02)-
              hmm_forward(em_test,.1,.2,.07,.02))<1e-10)

# Stage82 likelihood is reproduced if square terms are zero.
logLik_theta <- function(theta) {
  LL <- 0
  for(d in obs) {
    pi <- plogis(
      theta[3]+theta[13]*d$entry_age+
      theta[14]*d$male+theta[15]*d$entry_calendar+
      theta[16]*d$h+
      theta[21]*d$entry_age2+theta[22]*d$entry_calendar2)
    rho <- plogis(theta[4])
    lam <- plogis(
      theta[1]+theta[5]*d$h+
      theta[7]*d$age+theta[9]*d$male+
      theta[11]*d$calendar+
      theta[17]*d$age2+theta[19]*d$calendar2)
    tau <- plogis(
      theta[2]+theta[6]*d$h+
      theta[8]*d$age+theta[10]*d$male+
      theta[12]*d$calendar+
      theta[18]*d$age2+theta[20]*d$calendar2)
    LL <- LL+hmm_vary(d$em,lam,tau,pi,rho)
  }
  LL
}
nll_penalty <- function(theta)
  -logLik_theta(theta)+.5*sum(((theta-prior_mean)/prior_sd)^2)

# Exact Stage82 identity test: the full data log likelihood should
# match to numerical precision with the 6 quadratic terms fixed zero.
for(nm in c("null","genetics")) {
  a <- adj$fits[[nm]]
  theta <- c(as.numeric(a$theta),rep(0,6))
  ll_new <- logLik_theta(theta)
  cat("Stage82 likelihood identity",nm,": stored=",
      round(a$logLik,6)," recalculated=",round(ll_new,6),"\n")
  if(!isTRUE(all.equal(ll_new,as.numeric(a$logLik),
                      tolerance=1e-7)))
    fail("Stage82 HMM likelihood is NOT reproduced; stop before fitting")
}
cat("Stage82 nested-likelihood consistency checks: PASSED\n")

fit_model <- function(label,include_genetics,initial,maxit=250L){
  active <- if(include_genetics) seq_len(npar) else
    setdiff(seq_len(npar),5:6)
  anchor <- start
  anchor[active] <- initial[active]
  if(!include_genetics)anchor[5:6]<-0
  obj <- function(par){
    th <- anchor
    th[active] <- par
    nll_penalty(th)
  }
  t0 <- proc.time()[["elapsed"]]
  opt <- stats::optim(anchor[active],obj,method="L-BFGS-B",
    lower=lower[active],upper=upper[active],
    control=list(maxit=maxit,trace=0))
  theta <- anchor
  theta[active] <- opt$par
  elapsed <- proc.time()[["elapsed"]]-t0
  LL <- logLik_theta(theta)
  if(!is.finite(LL)||!is.finite(opt$value))
    fail(label,": nonfinite likelihood or fit objective")
  cat(sprintf("%-22s convergence=%d logLik=%.3f objective=%.3f elapsed=%.1fs\n",
      label,opt$convergence,LL,opt$value,elapsed))
  flush.console()
  list(label=label,theta=theta,logLik=LL,
       objective=opt$value,convergence=opt$convergence,
       elapsed=elapsed,message=opt$message,active=active)
}
# Both fits have the SAME age2/calendar2 nuisance terms; only
# transition homozygosity coefficients differ between the models.
null_initial <- c(as.numeric(adj$fits$null$theta),rep(0,6))
fit_null <- fit_model("quadratic_null",FALSE,null_initial)
gen_initial <- c(as.numeric(adj$fits$genetics$theta),rep(0,6))
fit_genetics <- fit_model("quadratic_genetics",TRUE,gen_initial)
fit_genetics_start2 <- fit_model("quadratic_alt_start",TRUE,
  fit_null$theta)
best <- if(fit_genetics$objective<=fit_genetics_start2$objective)
  fit_genetics else fit_genetics_start2

make_row <- function(key,fit){
  data.frame(model=key,n_badgers=length(chosen),
    n_tests=total_assays,convergence=fit$convergence,
    logLik=fit$logLik,penalised_objective=fit$objective,
    beta_hom_SE=fit$theta[5],beta_hom_EI=fit$theta[6],
    OR_hom_SE=exp(fit$theta[5]),OR_hom_EI=exp(fit$theta[6]),
    beta_entry_hom=fit$theta[16],OR_entry_hom=exp(fit$theta[16]),
    beta_age2_SE=fit$theta[17],
    beta_age2_EI=fit$theta[18],
    beta_calendar2_SE=fit$theta[19],
    beta_calendar2_EI=fit$theta[20],
    beta_entry_age2=fit$theta[21],
    beta_entry_calendar2=fit$theta[22])
}
results <- rbind(make_row("quadratic_null",fit_null),
                 make_row("quadratic_genetics",best))
cat("\n--- STAGE83 SENSITIVITY COMPARISON ---\n")
print(results,row.names=FALSE,digits=4)
cat("Stage82 adjusted SE and EI betas:",
    paste(signif(adj$fits$genetics$theta[5:6],4),collapse=", "),"\n")
cat("Stage83 quadratic SE and EI betas:",
    paste(signif(best$theta[5:6],4),collapse=", "),"\n")
cat("Adjusted model improvement (2*deltaLogLik):",
    round(2*(best$logLik-fit_null$logLik),4),"\n")
cat("Stage83 vs Stage82 genetic model 2*deltaLogLik:",
    round(2*(best$logLik-adj$fits$genetics$logLik),4),
    "(6 additional quadratic coefficients).\n")
cat("Agreement of alternate Stage83 starts:",
    signif(abs(fit_genetics$objective-
               fit_genetics_start2$objective),4),"\n")
cat("These are descriptive likelihood differences using regularised",
    " MAP estimates, NOT calibrated likelihood-ratio tests.\n")
cat("No diagnostic accuracy, capture probability, seasonality,",
    " survival or spatial social-group dependence is refitted.\n")

stem <- sprintf("stage83_quadratic_age_calendar_%04d",length(chosen))
dst <- file.path(base,paste0(stem,".csv"))
fit_path <- file.path(base,paste0(stem,".rds"))
if((file.exists(dst)||file.exists(fit_path)) &&
   Sys.getenv("OVERWRITE","0")!="1")
  fail("Existing Stage83 fit would be overwritten; keep it or explicitly",
       " use OVERWRITE=1.")
write.csv(results,dst,row.names=FALSE)
saveRDS(list(
    comparison=results,
    fits=list(null=fit_null,genetics=best,
      genetics_start1=fit_genetics,
      genetics_start2=fit_genetics_start2),
    parameter_names=names_theta,
    source_stage82_md5=unname(tools::md5sum(adj_path)),
    source_md5=tools::md5sum(files),
    n_eligible=length(eligible),n_selected=length(chosen),
    n_captures=total_captures,n_tests=total_assays,
    caveats=c(
      "same Stage76 confirmed live captures and test code",
      "317 original-only uncertain encounter occasions excluded",
      "same 1224 genetic eligibility and prior input",
      "six quadratic age/calendar terms augment the Stage82 linear model",
      "both genetic models adjust entry prevalence for homozygosity",
      "fixed diagnostic sensitivities/specificities and exposed multipliers",
      "simple age/calendar polynomial not a seasonal/annual random-walk",
      "observed capture and assay missingness are conditioned upon",
      "regularized MAP only, no standard errors/credible intervals",
      "no joint survival process or informative capture selection"
    )),fit_path)
cat("Saved:",dst,"\n")
cat("Saved:",fit_path,"\n")
if(any(vapply(list(fit_null,fit_genetics,fit_genetics_start2),
             function(z)z$convergence!=0L,logical(1))))
  fail("One or more Stage83 optimisations did not converge.")
cat("STAGE 83 QUADRATIC AGE/CALENDAR SENSITIVITY COMPLETE.\n")
