# ============================================================================
# Stage 79: DIRECT-DIAGNOSTIC, GENETIC SEI JOINT-LIKELIHOOD SMOKE TEST
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
# PILOT/PROOF OF COMPUTATION ONLY: constant quarterly hazards, fixed Se/Sp
# and exposed-state multipliers from original fitted SEI; no seasonal,
# calendar, age, sex, social-group, survival, trapping or observation model.
# Regularized MAP is NOT a fully Bayesian fit nor a final genetic effect.
#
# Reads NO DB, does NOT edit original data, outputs remain gitignored.
# N_PILOT (default 250) randomly selects genetic individuals using seed 79.
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
cat("\n=== STAGE 79: DIRECT-DIAGNOSTIC SEI GENETIC SMOKE TEST ===\n")
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
n_requested <- suppressWarnings(as.integer(Sys.getenv("N_PILOT","250")))
if(is.na(n_requested)||n_requested<40L)
  fail("N_PILOT must be an integer >= 40")
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
cat("Pilot randomly selected:",length(chosen),
    "using fixed seed 79 (no outcome-based enrichment).\n")

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
cat("Median observed quarters:",median(vapply(obs,function(z)z$nq,
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

# theta = logit(S->E baseline), logit(E->I baseline),
#         logit(prevalence at first observed capture),
#         logit(infectious fraction among prevalent infected),
#         genetic log-odds effect on S->E, on E->I.
# Penalisation is deliberately weak but avoids separated/unidentified
# estimates in a tiny computational pilot.
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
  pi <- plogis(theta[3L])
  rho <- plogis(theta[4L])
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
cat("Initial observed test log likelihood:",
    round(raw_loglik(start),3),"\n")
init_time <- proc.time()[["elapsed"]]
fit <- stats::optim(start,neglogpost,method="L-BFGS-B",
                    lower=lower,upper=upper,
                    control=list(maxit=200,trace=0))
elapsed <- proc.time()[["elapsed"]]-init_time
cat("\n--- STAGE 79 PENALIZED JOINT DIAGNOSTIC FIT ---\n")
cat("Optim convergence code:",fit$convergence,
    "(0 means optimizer converged)\n")
cat("Message:",fit$message,"\n")
cat("Elapsed seconds:",round(elapsed,1),"\n")
cat("Final observed diagnostic log likelihood:",
    round(raw_loglik(fit$par),3),"\n")
cat("Starting log posterior objective:",round(neglogpost(start),3),
    "optimized:",round(fit$value,3),"\n")
result <- data.frame(
  parameter=names(start),initial=as.numeric(start),
  fit=as.numeric(fit$par),lower=lower,upper=upper)
result$exp_fit <- exp(result$fit)
print(result[,c("parameter","initial","fit","exp_fit")],
      row.names=FALSE,digits=4)
cat("Estimated baseline quarterly S->E:",
    signif(plogis(fit$par[1L]),3),"\n")
cat("Estimated baseline quarterly E->I:",
    signif(plogis(fit$par[2L]),3),"\n")
cat("Exploratory homozygosity acquisition OR:",
    signif(exp(fit$par[5L]),3),"\n")
cat("Exploratory homozygosity progression OR:",
    signif(exp(fit$par[6L]),3),"\n")
cat("These are preliminary regularized MAP odds ratios,",
    "not credible intervals, and NOT valid publishable effects.\n")
dest <- file.path(base,"stage79_joint_sei_genetic_smoke.rds")
saveRDS(list(
  params=fit$par,convergence=fit$convergence,
  convergence_message=fit$message,
  optim_value=fit$value,elapsed=elapsed,
  n_genotyped=length(chosen),n_captures=total_captures,
  n_tests=total_assays,all_eligible=length(eligible),
  ids=chosen,pilot_seed=79L,
  fitted_emissions=list(Se=Se,Sp=Sp,exposed_multiplier=m),
  source_md5=tools::md5sum(files),
  assumptions=c(
    "Stage76: excludes 317 original-only nonconfirmed capture quarters",
    "Stage76: 301 DPP within-capture conflicts encoded missing",
    "Stage76: any M.bovis-positive culture plate yields culture positive",
    "within-quarter repeated captures and observed assays stacked",
    "no diagnostic information treated as a negative test",
    "genotyped pilot random among >=2-quarter/assay eligibility",
    "fixed original SEI inferred accuracy and exposed multipliers",
    "S/E/I latent states marginalized via HMM forward algorithm",
    "constant quarterly S->E and E->I baseline hazards",
    "common initial prevalence and initial infectious fraction",
    "genetic effects enter BOTH quarterly transition hazards directly",
    "regularized MAP with Gaussian parameter penalties",
    "no age, sex, calendar, territory, survival, trapping effort model",
    "not an exact reconstruction of the original input matrix",
    "NOT confirmatory or causal; no interval estimates"
  )),dest)
write.csv(result,file.path(base,"stage79_smoke_parameter_summary.csv"),
          row.names=FALSE)
cat("Saved:",dest,"\n")
if(fit$convergence!=0L)
  stop("STAGE 79: optimizer did not converge. See saved diagnostics.")
if(!is.finite(fit$value)||fit$value>neglogpost(start)+1e-6)
  stop("STAGE 79: optimization did not improve objective")
cat("STAGE 79 COMPUTATIONAL PILOT PASSED. NO ORIGINAL DATA CHANGED.\n")