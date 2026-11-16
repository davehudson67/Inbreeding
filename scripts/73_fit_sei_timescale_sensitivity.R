# ============================================================================
# Stage 73: EXPLORATORY SEI transition model sensitivity.
# Changes relative to Stage 69:
#
#   S->E: actual calendar quarter is a TIME-VARYING covariate, not
#         first-capture/entry decade. Age remains Cox baseline time.
#   E->I: TIME SINCE EXPOSURE is the Cox time scale; actual calendar
#         quarter is time-varying; age at exposure is an adjustment.
# IMPORTANT: with age as baseline, calendar quarter at an age risk set
# also contrasts BIRTH COHORTS. With exposure-duration baseline,
# calendar quarter at a duration risk set also contrasts EXPOSURE
# COHORTS. These fitted calendar coefficients are therefore NOT
# independently identified estimates of secular calendar trends.
#
# Cohort: 1603 genotyped badgers, shared SEI posterior draw columns.
# Left-truncate at first observed live capture, SEI history start, birth;
# censor at last SEI capture or before death quarter. Individuals
# inferred exposed before birth are excluded from E->I (QC reported).
# Do not extrapolate unobserved post-last-capture disease stage.
#
# MODEL IS MODULAR / EXPLORATORY, NOT A JOINT DIAGNOSTIC LIKELIHOOD FIT.
# It does NOT address changing diagnostic accuracy, biased observation
# windows, social groups or uncertainty in the upstream SEI model.
#
# Pilot: MODE=pilot (25 spread-out joint draws). Full: MODE=full
# (1000). Optional SHARD_ID/N_SHARDS to distribute full draw columns.
# Outputs in a new stage73_* namespace. Safe rerun with OVERWRITE=1.
# ============================================================================
options(stringsAsFactors=FALSE)
if(!requireNamespace("survival",quietly=TRUE))
  stop("R package 'survival' required")
suppressPackageStartupMessages(library(survival))
mode <- Sys.getenv("MODE","pilot")
if(!mode %in% c("pilot","full"))
  stop("MODE must be 'pilot' or 'full'")
shard <- suppressWarnings(as.integer(Sys.getenv("SHARD_ID","1")))
nshards <- suppressWarnings(as.integer(Sys.getenv("N_SHARDS","1")))
if(!is.finite(shard)||!is.finite(nshards)||
   nshards<1L||nshards>1000L||shard<1L||shard>nshards)
  stop("Invalid SHARD_ID/N_SHARDS")
sei_path <- "data/model_inputs/integrated_gm_sei_aligned_2026.rds"
model_path <- "data/model_inputs/integrated_gm_inputs_group_proxy_2026.rds"
stopifnot(file.exists(sei_path),file.exists(model_path))
z <- readRDS(sei_path); x <- readRDS(model_path)
N <- as.integer(z$N)
D <- length(z$draw_ids)
stopifnot(
  identical(N,1603L),
  identical(as.integer(z$draw_ids),1:1000),
  identical(dim(z$exposed_onset_q),c(N,D)),
  identical(dim(z$infectious_onset_q),c(N,D)),
  identical(z$id,toupper(trimws(as.character(x$id)))),
  identical(z$survival_input_md5,unname(tools::md5sum(model_path)))
)
sei_start <- as.integer(z$sei_start_q)
first_live <- as.integer(x$first_q)
birth <- as.integer(x$birth_q)
last <- as.integer(z$sei_last_capture_q)
death <- as.integer(z$death_q)
hom <- as.numeric(x$hom)
male <- as.integer(x$sex=="MALE")
if(anyNA(sei_start)||anyNA(first_live)||anyNA(birth)||
   anyNA(last)||anyNA(death)||anyNA(hom)||anyNA(male)||
   any(sei_start<1L)||any(last>200L)||any(sei_start>last)||
   any(first_live<birth)||any(first_live>last)||
   any(death<0L)||any(!x$sex %in% c("FEMALE","MALE"))||
   !is.finite(sd(hom))||sd(hom)<=0)
  stop("Input chronology, sex or genetic quality gate failed")
hom_z <- as.numeric((hom-mean(hom))/sd(hom))
landmark <- pmax(sei_start,first_live,birth)
last_risk <- pmin(last,ifelse(death>0L,death-1L,last))
if(mode=="pilot"){
  draw_ids <- unique(as.integer(round(seq.int(1,D,length.out=25L))))
} else {
  draw_ids <- seq_len(D)
}
draw_ids <- draw_ids[(draw_ids-1L)%%nshards==(shard-1L)]
if(!length(draw_ids))stop("No posterior draw assigned to this shard")
out <- "outputs/sei_transitions"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
stem <- sprintf("stage73_%s_shard_%02d_of_%02d",mode,shard,nshards)
resfile <- file.path(out,paste0(stem,".rds"))
csvfile <- file.path(out,paste0(stem,"_per_draw.csv"))
if(file.exists(resfile) && Sys.getenv("OVERWRITE","0")!="1")
  stop("Already fitted: ",resfile)
cat("\n=== STAGE 73: CORRECTED CALENDAR / EXPOSURE-TIME ANALYSIS ===\n")
cat("Cohort:",N,"; draws in this job:",length(draw_ids),
    "; mode:",mode,"shard:",shard,"/",nshards,"\n")
cat("Time scales: S->E attained age; E->I time since exposure.\n")
cat("Calendar covariate: actual at-risk quarter, updated quarterly.\n")
cat("Original ZIP integrity verified:",
    isTRUE(z$original_zip_checksum_verified),"\n")
cat("Post-first-capture follow-up starts: before birth prohibited.\n")
cat("SEI start before estimated birth:",sum(sei_start<birth),"\n")
cat("No post-entry risk-time:",sum(landmark>last_risk),"\n")
cat("These are inferred disease histories, not observed transition dates.\n")

# One row per at-risk quarter. The corresponding event, if any, is at
# the END of its quarter, giving strictly positive (start, stop] time.
# Time-varying calendar decade per row = (quarter - 100)/40.
# We do not use stage labels after the final diagnostic capture.
create_risk <- function(draw,transition) {
  E <- as.integer(z$exposed_onset_q[,draw])
  I <- as.integer(z$infectious_onset_q[,draw])
  if(anyNA(E)||anyNA(I)||any(E<0L)||any(I<0L)||
     any(E>last & E>0L)||any(I>last & I>0L)||
     any(I>0L & (E==0L|I<E)))
    stop("Invalid disease history in draw ",draw)
  if(transition=="SE") {
    eligible <- landmark<=last_risk & !(E>0L & E<landmark)
    entry <- landmark
    finish <- ifelse(E>0L & E>=entry & E<=last_risk,
                     E,last_risk)
    eventq <- E
    age_at_exposure <- rep(NA_real_,N)
  } else if(transition=="EI") {
    # Prevalent E before first capture is retained if the infectious
    # transition has NOT already happened at first capture. The
    # exposure onset must still be biologically after birth.
    eligible <- landmark<=last_risk &
      E>0L & E>=birth & E<=last_risk &
      (I==0L|I>=landmark)
    entry <- pmax(E,landmark)
    finish <- ifelse(I>0L & I>=entry & I<=last_risk,
                     I,last_risk)
    eventq <- I
    age_at_exposure <- (E-birth)*.25
  } else stop("Unknown transition")
  idx <- which(eligible & finish>=entry)
  if(!length(idx))stop("Empty transition risk set")
  n_quarters <- as.integer(finish[idx]-entry[idx]+1L)
  if(any(n_quarters<=0L))stop("Invalid length of at-risk interval")
  ii <- rep(idx,times=n_quarters)
  q <- as.integer(entry[ii]+sequence(n_quarters)-1L)
  event <- as.integer(eventq[ii]>0L & q==eventq[ii])
  if(transition=="SE") {
    t0 <- (q-birth[ii])*.25
    t1 <- t0+.25
  } else {
    t0 <- (q-E[ii])*.25
    t1 <- t0+.25
  }
  d <- data.frame(
    badger=ii,event=event,start=t0,stop=t1,
    hom_z=hom_z[ii],male=male[ii],
    calendar_decade=(q-100)/40,
    age_exposure=age_at_exposure[ii])
  if(anyNA(d[,c("start","stop","event","hom_z","male",
               "calendar_decade")])||
     any(d$start<0)||any(d$stop<=d$start)||
     sum(d$event)>length(idx))
    stop("Invalid quarter-level risk history in draw ",draw)
  list(data=d,n_at_risk=length(idx),n_person_quarters=nrow(d),
       events=sum(d$event),
       E_before_birth=sum(E>0L & E<birth),
       I_before_birth=sum(I>0L & I<birth),
       E_before_landmark=sum(E>0L & E<landmark),
       I_before_landmark=sum(I>0L & I<landmark))
}
fit_one <- function(draw,transition) {
  dat <- create_risk(draw,transition)
  d <- dat$data
  if(dat$events<15L||dat$n_at_risk<40L)
    stop("Too few events/persons in ",transition," draw ",draw,
         ": ",dat$events," / ",dat$n_at_risk)
  f <- if(transition=="SE") {
    survival::Surv(start,stop,event) ~ hom_z + male + calendar_decade
  } else {
    if(anyNA(d$age_exposure))
      stop("Unexpected NA in age at exposure")
    survival::Surv(start,stop,event) ~
      hom_z + male + calendar_decade + age_exposure
  }
  warnings <- character()
  fit <- withCallingHandlers(
    survival::coxph(
      f,data=d,ties="efron",singular.ok=FALSE,
      control=survival::coxph.control(iter.max=50,eps=1e-9)
    ),
    warning=function(w){
      warnings <<- c(warnings,conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  co <- stats::coef(fit)
  V <- stats::vcov(fit)
  terms <- if(transition=="SE")
    c("hom_z","male","calendar_decade") else
    c("hom_z","male","calendar_decade","age_exposure")
  if(length(warnings)||!identical(names(co),terms)||
     any(!is.finite(co))||any(!is.finite(V))||
     any(diag(V)<=0))
    stop("Unstable Cox regression in ",transition," draw ",draw,
         "; warnings: ",paste(unique(warnings),collapse=";"))
  out <- data.frame(
    draw=draw,transition=transition,
    n_risk=dat$n_at_risk,person_quarters=dat$n_person_quarters,
    events=dat$events,
    E_before_birth=dat$E_before_birth,
    I_before_birth=dat$I_before_birth,
    E_before_landmark=dat$E_before_landmark,
    I_before_landmark=dat$I_before_landmark)
  # Keep an identical output schema for SE and EI, including
  # explicitly missing age-at-exposure coefficients in SE.
  for(j in c("hom_z","male","calendar_decade","age_exposure")){
    out[[paste0("beta_",j)]] <-
      if(j %in% terms) unname(co[j]) else NA_real_
    out[[paste0("se_",j)]] <-
      if(j %in% terms) sqrt(V[j,j]) else NA_real_
  }
  out
}
rows <- vector("list",2L*length(draw_ids))
tic <- proc.time()[3]
for(k in seq_along(draw_ids)) {
  draw <- draw_ids[[k]]
  rows[[2L*k-1L]] <- fit_one(draw,"SE")
  rows[[2L*k]] <- fit_one(draw,"EI")
  if(k==1L||k%%25L==0L||k==length(draw_ids)){
    cat(sprintf(
      "Done %d/%d draws (draw %d): acquisition %d events; progression %d; elapsed %.1f seconds\n",
      k,length(draw_ids),draw,
      rows[[2L*k-1L]]$events,rows[[2L*k]]$events,
      proc.time()[3]-tic))
    flush.console()
  }
}
res <- do.call(rbind,rows)
res <- res[order(res$draw,res$transition),,drop=FALSE]
if(nrow(res)!=2L*length(draw_ids)||
   anyDuplicated(res[,c("draw","transition")]))
  stop("Output coverage check failed")
pool <- function(d,key) {
  betas <- d[[paste0("beta_",key)]]
  ses <- d[[paste0("se_",key)]]
  if(any(!is.finite(betas))||any(!is.finite(ses))||any(ses<=0))
    stop("Nonfinite pooled variable: ",key)
  m <- length(betas); b <- if(m>1)stats::var(betas) else 0
  w <- mean(ses^2); tot <- w+(1+1/m)*b
  df <- if(b==0) Inf else (m-1)*(1+w/((1+1/m)*b))^2
  crit <- if(is.finite(df))qt(.975,df) else qnorm(.975)
  data.frame(
    transition=d$transition[1],predictor=key,draws=m,
    HR=exp(mean(betas)),low95=exp(mean(betas)-crit*sqrt(tot)),
    high95=exp(mean(betas)+crit*sqrt(tot)),
    between_fraction=(1+1/m)*b/tot)
}
pooled <- do.call(rbind,lapply(c("SE","EI"),function(t){
  dat <- res[res$transition==t,,drop=FALSE]
  parameters <- if(t=="SE")
    c("hom_z","male","calendar_decade") else
    c("hom_z","male","calendar_decade","age_exposure")
  do.call(rbind,lapply(parameters,function(p)pool(dat,p)))
}))
cat("\n=== EXPLORATORY TIME-SCALE SENSITIVITY RESULTS ===\n")
print(pooled,row.names=FALSE,digits=4)
cat("\nEvent and source chronology QC across draws in this run:\n")
for(t in c("SE","EI")) {
  a <- res[res$transition==t,,drop=FALSE]
  cat(t,": median events ",median(a$events),
      ", range ",paste(range(a$events),collapse="–"),
      ", median risk set ",median(a$n_risk),"\n",sep="")
}
cat("E onset before estimated birth, median individuals per draw:",
    median(res$E_before_birth),"\n")
cat("I onset before estimated birth, median individuals per draw:",
    median(res$I_before_birth),"\n")
cat("Explicitly excluded E-before-birth histories from progression.\n")
saveRDS(list(results=res,pooled=pooled,mode=mode,
             shard=shard,nshards=nshards,draw_ids=draw_ids,
             model_input_md5=unname(tools::md5sum(model_path)),
             sei_input_md5=unname(tools::md5sum(sei_path)),
             script_md5=unname(tools::md5sum(
               "scripts/73_fit_sei_timescale_sensitivity.R")),
             assumptions=c(
               "post-first-capture left truncation",
               "last SEI capture and pre-death censoring",
               "S->E age baseline with time-varying real calendar quarter",
               "E->I time since exposed baseline with calendar quarter and age at exposure",
               "exposure before biological birth excluded from progression",
               "SEI posterior kept modular, diagnostic observations not refitted",
               "approximate Rubin pooling; no causal interpretation",
               "social-group and observation effort not adjusted"
             )),resfile)
write.csv(res,csvfile,row.names=FALSE)
cat("Saved:",resfile,"\nSaved:",csvfile,"\n")
cat("Calendar HR interpretation: age-time-scale S->E contrasts birth",
    "cohorts, while exposure-duration E->I contrasts exposure cohorts.",
    "Do NOT interpret either HR as a pure secular calendar-period effect.\n")
cat("STAGE 73 COMPLETE. This does not resolve diagnostic observation bias.\n")
