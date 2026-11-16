# ============================================================================
# Stage 69: FIRST MODULAR SEI TRANSITION REGRESSION (1603 genotyped badgers)
#
# Separate age-time-scale left-truncated Cox models for:
#  S -> E : susceptible to exposed/infected onset;
#  E -> I : exposure to infectious onset.
# Model coefficients: standardized homozygosity (+1 SD), male sex,
# and calendar decade of entry. Flexible age dependence is absorbed by
# the Cox baseline hazard, not forced into a linear age coefficient.
#
# One WHOLE coherent shared posterior draw is used at a time. Models
# fitted per draw are exploratory CUT/MODULAR posterior imputations:
# the SEI posterior was learned from diagnostic histories without
# conditioning on genotype, sex or these survival-regression effects.
# Consequently these pooled coefficients must NOT be described as a
# fully joint posterior, causal effects or validated frequentist CIs.
#
# Observation windows: enter S at SEI start; ignore prevalent E before
# start; enter E at inferred onset; censor at SEI last capture and
# before any death quarter to avoid assuming later disease states.
# Transitions occurring in the death quarter are excluded (unknown order).
# Same-quarter E and I is allowed as one-quarter E -> I duration.
#
# MODE=pilot (25 spread-out draws) or MODE=full (all 1000).
# SHARD_ID/N_SHARDS partition whole draw columns for independent jobs.
# No alteration of source data or previous survival fits.
# Requires 'survival' R package (base recommended in most R installs).
# ============================================================================
options(stringsAsFactors=FALSE)
if(!requireNamespace("survival",quietly=TRUE))
  stop("Please install the 'survival' R package on this server.")
suppressPackageStartupMessages(library(survival))
mode <- Sys.getenv("MODE","pilot")
if(!mode %in% c("pilot","full"))stop("MODE must be pilot or full")
shard <- suppressWarnings(as.integer(Sys.getenv("SHARD_ID","1")))
nshards <- suppressWarnings(as.integer(Sys.getenv("N_SHARDS","1")))
if(!is.finite(shard)||!is.finite(nshards)||
   shard<1L||nshards<1L||shard>nshards||nshards>1000L)
  stop("Invalid SHARD_ID or N_SHARDS")
sourcefile <- "data/model_inputs/integrated_gm_sei_aligned_2026.rds"
original <- "data/model_inputs/integrated_gm_inputs_group_proxy_2026.rds"
if(!file.exists(sourcefile)||!file.exists(original))
  stop("Run Stages 66-67 first: SEI or original survival data missing")
z <- readRDS(sourcefile)
x <- readRDS(original)
stopifnot(
  identical(as.integer(z$N),1603L),
  identical(as.integer(z$T),200L),
  identical(as.integer(z$draw_ids),1:1000),
  identical(dim(z$exposed_onset_q),c(1603L,1000L)),
  identical(dim(z$infectious_onset_q),c(1603L,1000L)),
  identical(z$id,toupper(trimws(as.character(x$id)))),
  identical(z$survival_input_md5,unname(tools::md5sum(original)))
)
if(is.null(x$birth_q)||is.null(x$hom)||is.null(x$sex))
  stop("Original input lacks required birth_q/hom/sex fields")
N <- length(z$id); D <- length(z$draw_ids)
start <- as.integer(z$sei_start_q)
last <- as.integer(z$sei_last_capture_q)
death <- as.integer(z$death_q)
birth <- as.integer(x$birth_q)
male <- as.integer(x$sex=="MALE")
rawhom <- as.numeric(x$hom)
if(length(start)!=N||length(last)!=N||length(death)!=N||
   length(birth)!=N||length(male)!=N||length(rawhom)!=N||
   anyNA(start)||anyNA(last)||anyNA(death)||anyNA(birth)||
   anyNA(male)||anyNA(rawhom)||
   any(start<1L)||any(start>last)||any(last>z$T)||
   any(death<0L)||any(!x$sex %in% c("MALE","FEMALE"))||
   sd(rawhom)<=0)
  stop("Invalid SEI observation windows, birth times, sex or homozygosity")
hom_z <- as.numeric((rawhom-mean(rawhom))/sd(rawhom))
# Quarter q represents the [q,q+1) age interval; event in quarter q
# occurs at its end, avoiding zero-duration same-quarter E->I events.
# Person-time after a recorded death is never included.
follow_end <- pmin(last,ifelse(death>0L,death-1L,last))
if(any(start<birth))
  stop("SEI history begins before birth for at least one individual; audit times")
if(any(!is.finite(hom_z)))stop("Nonfinite homozygosity covariate")
source_md5 <- unname(tools::md5sum(sourcefile))
original_md5 <- unname(tools::md5sum(original))
if(mode=="pilot") {
  draw_ids <- unique(as.integer(round(seq.int(1,D,length.out=25L))))
} else {
  draw_ids <- seq_len(D)
}
draw_ids <- draw_ids[(draw_ids-1L)%%nshards==(shard-1L)]
if(!length(draw_ids))stop("Shard was assigned zero draws")
out <- "outputs/sei_transitions"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
stem <- sprintf("%s_shard_%02d_of_%02d",mode,shard,nshards)
dest <- file.path(out,paste0(stem,".rds"))
tabdest <- file.path(out,paste0(stem,"_per_draw.csv"))
if(file.exists(dest) && Sys.getenv("OVERWRITE","0")!="1")
  stop("Output already exists: ",dest,
       ". Set OVERWRITE=1 only when intentionally replacing it.")

cat("\n=== STAGE 69: SEI ACQUISITION AND PROGRESSION FIT ===\n")
cat("Mode:",mode,"shard",shard,"of",nshards,"\n")
cat("Cohort:",N,"matched genotyped badgers; posterior draws:",
    D,"\n")
cat("Draws in this job:",length(draw_ids),"; first:",
    head(draw_ids,5),"\n")
cat("SEI source mode:",z$source_mode,"\n")
cat("Survival input MD5:",original_md5,"\n")
cat("Death before/within SEI end (censor prior to death quarter):",
    sum(death>0L & death<=last),"\n")
cat("Mean homozygosity:",mean(rawhom),"; SD:",sd(rawhom),
    "; unique values:",length(unique(rawhom)),"\n")
cat("OUTCOMES are INFERRED disease trajectories, not directly",
    "observed events.\n")

build_data <- function(draw, transition) {
  e <- as.integer(z$exposed_onset_q[,draw])
  i <- as.integer(z$infectious_onset_q[,draw])
  if(anyNA(e)||anyNA(i)||any(e<0L)||any(i<0L)||
     any(i>0L & (e==0L | i<e)) )
    stop("Invalid paired onset history in draw ",draw)
  # All event times limited to the person's observed window,
  # with the additional terminal death-quarter censoring applied.
  if(transition=="SE") {
    # Individuals exposed before the start are left-prevalent and
    # cannot be treated as susceptible incident cases.
    eligible <- start<=follow_end & !(e>0L & e<start)
    event <- e>0L & e>=start & e<=follow_end
    entry <- start
    exit <- ifelse(event,e,follow_end)
  } else if(transition=="EI") {
    # Incident E within follow-up only: no prevalent E histories.
    eligible <- e>0L & e>=start & e<=follow_end
    entry <- e
    event <- eligible & i>0L & i>=e & i<=follow_end
    exit <- ifelse(event,i,follow_end)
  } else stop("Unknown transition")
  keep <- which(eligible & exit>=entry)
  if(!length(keep))return(NULL)
  data.frame(
    id=keep,
    age_start=.25*(entry[keep]-birth[keep]),
    age_stop=.25*(exit[keep]+1L-birth[keep]),
    event=as.integer(event[keep]),
    hom_z=hom_z[keep],
    male=male[keep],
    # +1 unit corresponds to 10 calendar years later.
    entry_decade=(entry[keep]-100)/40
  )
}
terms <- c("hom_z","male","entry_decade")
fit_one <- function(draw,transition) {
  d <- build_data(draw,transition)
  if(is.null(d)||nrow(d)<40L||sum(d$event)<15L||
     anyNA(d)||any(d$age_start<0)||
     any(d$age_stop<=d$age_start))
    stop("Invalid or insufficient at-risk data for draw ",
         draw," transition ",transition)
  warnings <- character()
  fit <- withCallingHandlers(
    survival::coxph(
      survival::Surv(age_start,age_stop,event) ~
        hom_z + male + entry_decade,
      data=d,ties="efron",singular.ok=FALSE,
      control=survival::coxph.control(iter.max=35,
                                      eps=1e-9)
    ),
    warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  beta <- stats::coef(fit)
  covariance <- stats::vcov(fit)
  if(length(beta)!=length(terms)||!identical(names(beta),terms)||
     any(!is.finite(beta))||any(!is.finite(covariance))||
     any(diag(covariance)<=0))
    stop("Nonfinite/singular Cox coefficients in draw ",draw,
         " transition ",transition)
  if(length(warnings))
    stop("Cox model warning in draw ",draw," transition ",
         transition,": ",paste(unique(warnings),collapse="; "))
  data.frame(
    draw=draw,transition=transition,
    n_risk=nrow(d),events=sum(d$event),
    person_quarters=sum(d$age_stop-d$age_start)/.25,
    beta_hom_z=unname(beta["hom_z"]),
    se_hom_z=sqrt(covariance["hom_z","hom_z"]),
    beta_male=unname(beta["male"]),
    se_male=sqrt(covariance["male","male"]),
    beta_entry_decade=unname(beta["entry_decade"]),
    se_entry_decade=sqrt(covariance["entry_decade","entry_decade"]),
    row.names=NULL
  )
}

results <- vector("list",length(draw_ids)*2L)
idx <- 0L
start_time <- Sys.time()
for(j in seq_along(draw_ids)) {
  draw <- draw_ids[[j]]
  for(transition in c("SE","EI")) {
    idx <- idx + 1L
    results[[idx]] <- fit_one(draw,transition)
  }
  if(j==1L || j%%25L==0L || j==length(draw_ids)) {
    lasttwo <- results[(idx-1L):idx]
    cat(sprintf("Finished %d / %d shared posterior draws; draw %d; ",
                j,length(draw_ids),draw))
    cat(sprintf("S->E events %d; E->I events %d; elapsed %.1f min\n",
                lasttwo[[1]]$events,lasttwo[[2]]$events,
                as.numeric(difftime(Sys.time(),start_time,units="mins"))))
    flush.console()
  }
}
result <- do.call(rbind,results)
if(nrow(result)!=2L*length(draw_ids)||
   anyDuplicated(result[,c("draw","transition")]))
  stop("Unexpected draw-wise result count or duplicate fits")

# A concise exploratory pooled output; not a fully joint SEI model.
pool <- function(d,variable) {
  b <- d[[paste0("beta_",variable)]]
  se <- d[[paste0("se_",variable)]]
  M <- length(b)
  mu <- mean(b)
  within <- mean(se^2)
  between <- if(M>1L) stats::var(b) else 0
  total <- within+(1+1/M)*between
  # Rubin's large-sample approximation; especially tentative for
  # inferred-stage/posterior inputs and informative follow-up.
  df <- if(between<=0) Inf else
    (M-1L)*(1+within/((1+1/M)*between))^2
  crit <- if(is.finite(df))qt(.975,df) else qnorm(.975)
  data.frame(transition=d$transition[1],parameter=variable,
             n_draws=M,logHR=mu,logHR_se=sqrt(total),
             HR=exp(mu),HR_low95=exp(mu-crit*sqrt(total)),
             HR_high95=exp(mu+crit*sqrt(total)),
             within_var=within,between_var=between,
             df=if(is.finite(df))df else NA_real_)
}
pooled <- do.call(rbind,lapply(c("SE","EI"),function(stage){
  d <- result[result$transition==stage,,drop=FALSE]
  do.call(rbind,lapply(terms,function(name)pool(d,name)))
}))
cat("\n=== EXPLORATORY TRANSITION HAZARD RATIOS ===\n")
print(pooled[,c("transition","parameter","n_draws","HR",
                "HR_low95","HR_high95")],
      row.names=FALSE,digits=4)
cat("\nTransition risk and events (across draws in this shard):\n")
for(t in c("SE","EI")) {
  d <- result[result$transition==t,,drop=FALSE]
  cat(t,": median risk set ",median(d$n_risk),
      "; median event count ",median(d$events),
      "; range ",paste(range(d$events),collapse=".."),"\n",sep="")
}
cat("\nCAUTION: This is an exploratory modular (cut) analysis.",
    "It uses inferred SEI histories from another fitted model,",
    "not raw diagnostic likelihood; pooled intervals are approximate.",
    "SEI follow-up may be informative. Do not infer causation.\n")
saveRDS(list(
  results=result,pooled=pooled,
  model="Age-timescale left-truncated Cox PH, Efron ties",
  transition_states=c("S->E","E->I"),
  mode=mode,shard=shard,nshards=nshards,
  ids=z$id,draw_ids=draw_ids,N=N,
  source_mode=z$source_mode,
  model_input_md5=original_md5,sei_input_md5=source_md5,
  assumptions=c(
    "1000 shared SEI posterior draw columns preserved across individuals",
    "modular cut analysis; SEI diagnostics are not refitted",
    "SEI history begins at individual start quarter and ends at last capture",
    "no assumption of disease state after last SEI capture",
    "exclude incident infection prior to modelled SEI observation start",
    "exclude exposure-to-infectiousness prevalent E at observation start",
    "stop risk at quarter before any known death",
    "quarterly event times treated as occurring at end of event quarter",
    "same-quarter E->I allowed as quarter-long at-risk period",
    "flexible baseline age dependence via Cox age timescale",
    "calendar decade covariate uses entry quarter",
    "hom effect per 1 original-cohort homozygosity standard deviation",
    "approximate Rubin pooling not a joint posterior or causal inference",
    "post-last-capture censoring may be informative"
  )
),dest)
write.csv(result,tabdest,row.names=FALSE)
cat("\nSaved:",dest,"\n")
cat("Saved:",tabdest,"\n")
cat("STAGE 69",toupper(mode),"SHARD",shard,"OF",nshards,"COMPLETE\n")
