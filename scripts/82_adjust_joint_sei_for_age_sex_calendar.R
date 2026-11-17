# ============================================================================
# Stage 82: DIRECT-DIAGNOSTIC, GENETIC SEI JOINT-LIKELIHOOD SMOKE TEST
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
# AGE/SEX/CALENDAR ADJUSTED ROBUSTNESS PILOT: time-varying transition hazards
# and exposed-state multipliers from original fitted SEI; no seasonal,
# calendar, age, sex, social-group, survival, trapping or observation model.
# Regularized MAP is NOT a fully Bayesian fit nor a final genetic effect.
#
# Reads NO DB, does NOT edit original data, outputs remain gitignored.
# N_STAGE82 (default 99999) uses ALL eligible genetic individuals; optional subset uses seed 79.
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
cat("\n=== STAGE 82: DIRECT-DIAGNOSTIC SEI GENETIC SMOKE TEST ===\n")
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
n_requested <- suppressWarnings(as.integer(Sys.getenv("N_STAGE82","99999")))
if(is.na(n_requested)||n_requested<40L)
  fail("N_STAGE82 must be an integer >= 40")
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

# ============================================================================
# STAGE 82: adjust both transition hazards for age, sex, calendar time;
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
# The calendar trend is deliberately first-order linear (per decade);
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
cat("\n=== STAGE 82: TIME/SEX ADJUSTED JOINT GENETIC SEI ===\n")
cat("Badgers:",length(chosen),"test results:",total_assays,"\n")
cat("Time-varying transition effects: age (per 4 years),",
    "calendar (per 10 years), sex (male vs female).\n")
cat("Initial prevalence adjusts for age at first capture, calendar",
    "year, sex and homozygosity.\n")

# Parameter layout:
# 1 SE intercept, 2 EI intercept, 3 entry prevalence intercept,
# 4 P(initial I|already infected) intercept,
# 5 hom SE, 6 hom EI,
# 7 age SE, 8 age EI, 9 male SE, 10 male EI,
# 11 calendar SE, 12 calendar EI,
# 13 entry age, 14 entry male, 15 entry calendar, 16 entry hom.
names_theta <- c("logit_SE","logit_EI","logit_entry_pi",
    "logit_entry_rho","beta_hom_SE","beta_hom_EI",
    "beta_age_SE","beta_age_EI","beta_male_SE","beta_male_EI",
    "beta_calendar_SE","beta_calendar_EI",
    "beta_entry_age","beta_entry_male",
    "beta_entry_calendar","beta_entry_hom")
npar <- length(names_theta)
start <- setNames(rep(0,npar),names_theta)
start[1:4] <- as.numeric(ref$fits$both$theta[1:4])
start[5:6] <- 0
lower <- c(-9,-9,-9,-7,rep(-4,12))
upper <- c(0,0,0,7,rep(4,12))
# Same priors on shared parameters as Stage80.
prior_mean <- c(qlogis(.035),qlogis(.05),qlogis(.1),
                qlogis(.5),0,0,rep(0,10))
prior_sd <- c(2,2,2,2.5,1.5,1.5,
              1,1,1,1,.8,.8,
              1,1,1,1)
if(length(lower)!=npar||length(upper)!=npar||
   length(prior_sd)!=npar||length(prior_mean)!=npar)
  fail("Stage82 parameter/prior dimension error")
if(any(start<lower)|any(start>upper))
  fail("Stage80 initial estimates outside Stage82 bounds")

# For variable transition hazards, every interval has its own lambda/tau.
# Emissions at q are observed conditional on state at q, as Stage80.
hmm_vary <- function(em,lam,tau,pi,rho) {
  if(!is.matrix(em)||ncol(em)!=3L||nrow(em)!=length(lam)||
     length(tau)!=length(lam))
    fail("Stage82 emission/transition sizes differ")
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
      fail("Stage82 nonfinite HMM emission probability")
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

logLik_theta <- function(theta) {
  LL <- 0
  for(d in obs) {
    pi <- plogis(theta[3]+theta[13]*d$entry_age+
       theta[14]*d$male+theta[15]*d$entry_calendar+
       theta[16]*d$h)
    rho <- plogis(theta[4])
    lam <- plogis(theta[1]+theta[5]*d$h+
       theta[7]*d$age+theta[9]*d$male+
       theta[11]*d$calendar)
    tau <- plogis(theta[2]+theta[6]*d$h+
       theta[8]*d$age+theta[10]*d$male+
       theta[12]*d$calendar)
    LL <- LL+hmm_vary(d$em,lam,tau,pi,rho)
  }
  LL
}
nll_penalty <- function(theta) {
  -logLik_theta(theta)+.5*sum(((theta-prior_mean)/prior_sd)^2)
}
fit_model <- function(label,fix_genetic,initial,maxit=160L) {
  active <- if(fix_genetic) setdiff(seq_len(npar),5:6) else seq_len(npar)
  anchor <- start
  anchor[active] <- initial[active]
  if(fix_genetic)anchor[5:6]<-0
  objective <- function(v){
    theta <- anchor
    theta[active] <- v
    nll_penalty(theta)
  }
  t0 <- proc.time()[["elapsed"]]
  opt <- stats::optim(anchor[active],objective,method="L-BFGS-B",
             lower=lower[active],upper=upper[active],
             control=list(maxit=maxit,trace=0))
  theta <- anchor
  theta[active] <- opt$par
  LL <- logLik_theta(theta)
  elapsed <- proc.time()[["elapsed"]]-t0
  if(!is.finite(opt$value)||!is.finite(LL))
    fail("Nonfinite Stage82 model objective/log likelihood: ",label)
  cat(sprintf("%-22s convergence=%d, logLik=%.3f, penalisedObjective=%.3f, elapsed=%.1fs\n",
              label,opt$convergence,LL,opt$value,elapsed))
  flush.console()
  list(label=label,theta=theta,convergence=opt$convergence,
       logLik=LL,objective=opt$value,message=opt$message,
       elapsed=elapsed,active=active)
}
cat("\nFitting sex/age/time/entry-adjusted NO-transition-genetics model.\n")
fit_null <- fit_model("adjusted_null",TRUE,start)
cat("Fitting adjusted model WITH both homozygosity transitions.\n")
initial <- fit_null$theta
initial[5:6] <- as.numeric(ref$fits$both$theta[5:6])
fit_genetics <- fit_model("adjusted_genetics",FALSE,initial)
cat("Refitting from zero genetic transition effects to check start.\n")
fit_genetics_start2 <- fit_model("adjusted_start2",FALSE,
                                  fit_null$theta)
best <- if(fit_genetics$objective<=fit_genetics_start2$objective)
  fit_genetics else fit_genetics_start2

summary <- do.call(rbind,lapply(list(null=fit_null,
                                    genetics=best),function(z){
  data.frame(model=z$label,convergence=z$convergence,
             n_badgers=length(chosen),n_tests=total_assays,
             logLik=z$logLik,penalised_objective=z$objective,
             beta_hom_SE=z$theta[5],
             beta_hom_EI=z$theta[6],
             OR_hom_SE=exp(z$theta[5]),
             OR_hom_EI=exp(z$theta[6]),
             beta_entry_hom=z$theta[16],
             OR_entry_hom=exp(z$theta[16]),
             beta_age_SE=z$theta[7],
             beta_age_EI=z$theta[8],
             beta_male_SE=z$theta[9],
             beta_male_EI=z$theta[10],
             beta_calendar_SE=z$theta[11],
             beta_calendar_EI=z$theta[12])
}))
cat("\n--- STAGE82 FIT COMPARISON ---\n")
print(summary,row.names=FALSE,digits=4)
cat("\nOriginal Stage80 genetic coefficients: ",
    paste(round(ref$fits$both$theta[5:6],4),collapse=", "),"\n")
cat("Adjusted transition genetic coefficients: ",
    paste(round(best$theta[5:6],4),collapse=", "),"\n")
cat("Descriptive 2*logLik improvement:",
    round(2*(best$logLik-fit_null$logLik),4),"\n")
cat("Agreement of alternate starting values (penalised objective):",
    signif(abs(fit_genetics$objective-fit_genetics_start2$objective),4),
    "\n")
cat("Model is conditional on observed captures and tests;",
    "Se/Sp fixed. It is not a survival model.\n")
dst <- file.path(base,"stage82_age_sex_calendar_adjustment.csv")
fit_path <- file.path(base,"stage82_age_sex_calendar_adjustment.rds")
write.csv(summary,dst,row.names=FALSE)
saveRDS(list(summary=summary,
   fits=list(null=fit_null,genetics=best,genetics_start1=fit_genetics,
             genetics_start2=fit_genetics_start2),
   parameter_names=names_theta,source_md5=tools::md5sum(files),
   source_stage80_md5=unname(tools::md5sum(ref_path)),
   eligible=length(eligible),selected=length(chosen),
   captures=total_captures,diagnostic_tests=total_assays,
   assumptions=c(
      "confirmed live captures only; 317 original-only quarters excluded",
      "genotyped selected by >=2 diagnostic quarters and readable assays",
      "underlying age via prior GM birth_q, sex via original genetic dataset",
      "age continuous linear per 4 years; calendar linear per decade",
      "each transition has age, male, calendar coefficients",
      "initial infection prevalence has age, sex, calendar, hom effects",
      "initial infectious share common, no direct S->I",
      "raw six-test binary diagnostic likelihood with multiple captures stacked",
      "source Se/Sp and exposed multiplier posterior means treated as exact",
      "MAP weak Gaussian regularisation; optimise nested models separately",
      "SEI disease states marginalized over forward HMM",
      "not yet seasonal, year random-walk, capture/missingness/survival joint",
      "no calibrated interval estimates or confirmatory inference"
   )),fit_path)
cat("\nSaved:",dst,"\n")
cat("Saved:",fit_path,"\n")
if(fit_null$convergence!=0L||fit_genetics$convergence!=0L||
   fit_genetics_start2$convergence!=0L)
  fail("Stage82 optimisation did not converge; inspect fitted outputs")
if(best$objective > fit_null$objective+1e-5)
  fail("Stage82 genetic model worse penalised objective than nested null")
cat("STAGE 82 ADJUSTED DIRECT-DIAGNOSTIC SEI COMPLETE.\n")