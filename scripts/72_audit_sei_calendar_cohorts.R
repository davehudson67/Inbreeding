# ===========================================================================
# Stage 72: Read-only audit of apparent SEI calendar effects.
# Summarizes:
#   - Cohort entry dates and available SEI follow-up;
#   - S->E and E->I inferred events by first-capture cohort;
#   - Crude events per 100 animal-years with person-time denominators;
#   - Calendar time of inferred events (not the same as entry decade).
#
# Uses the SAME post-first-capture risk definitions as Stage 69.
# 1000 shared posterior draws are summarized without treating them as
# 1000 different badgers. Intervals quantify *posterior trajectory*
# variation only; they do NOT reflect complete regression uncertainty.
# No fitting, no modification of data, no Supabase or MCMC.
# ===========================================================================
options(stringsAsFactors=FALSE)
path <- "data/model_inputs/integrated_gm_sei_aligned_2026.rds"
input <- "data/model_inputs/integrated_gm_inputs_group_proxy_2026.rds"
if(!file.exists(path)||!file.exists(input))
  stop("Missing existing aligned SEI inputs; run Stage 66 first")
z <- readRDS(path); x <- readRDS(input)
stopifnot(z$N==1603L,
          identical(as.integer(z$draw_ids),1:1000),
          identical(z$id,toupper(trimws(as.character(x$id)))),
          identical(z$survival_input_md5,unname(tools::md5sum(input))))
N <- as.integer(z$N); D <- length(z$draw_ids)
start <- as.integer(z$sei_start_q)
last <- as.integer(z$sei_last_capture_q)
first_live <- as.integer(x$first_q)
birth <- as.integer(x$birth_q)
death <- as.integer(z$death_q)
E <- z$exposed_onset_q
I <- z$infectious_onset_q
stopifnot(identical(dim(E),c(N,D)),identical(dim(I),c(N,D)),
          length(start)==N,length(last)==N,
          length(first_live)==N,length(birth)==N,length(death)==N,
          all(is.finite(start)),all(is.finite(last)),
          all(is.finite(first_live)),all(is.finite(birth)),
          all(is.finite(death)),all(first_live>=birth),
          all(first_live<=last),all(last<=200L))
landmark <- pmax(start,first_live,birth)
end <- pmin(last,ifelse(death>0L,death-1L,last))
bins <- seq_len(5L)
labels <- sprintf("%d-%d",1976+10L*(bins-1L),
                  1985+10L*(bins-1L))
quarter_bin <- function(q)
  pmin(5L,pmax(1L,as.integer((as.integer(q)-1L)%/%40L+1L)))
cohort <- quarter_bin(first_live)
range_checks <- c(
  first_capture_before_sei_start=sum(first_live<start),
  sei_start_before_first_capture=sum(start<first_live),
  sei_start_before_birth=sum(start<birth),
  no_post_entry_followup=sum(landmark>end))
cat("\n=== STAGE 72: SEI CALENDAR AND RISK-SET AUDIT ===\n")
cat("Full cohort:",N,"joint draws:",D,"\n")
cat("SEI SHA256 original ZIP verified:",isTRUE(z$original_zip_checksum_verified),"\n")
cat("Start/entry reconciliation:\n")
print(range_checks)

capture_q <- rowSums(x$y==1L)
cohort_summary <- do.call(rbind,lapply(bins,function(k){
  ix <- which(cohort==k)
  if(length(ix)==0L) {
    # Empty recruitment decades are valid and must remain visible;
    # return NA for undefined demographic summaries, not a model error.
    return(data.frame(
      capture_cohort=labels[k],
      n_badgers=0L,
      median_age_first_capture=NA_real_,
      median_SEI_followup_quarters=NA_real_,
      mean_homozygosity=NA_real_,
      male_fraction=NA_real_,
      median_live_captures=NA_real_,
      proportion_without_post_entry_followup=NA_real_
    ))
  }
  data.frame(
    capture_cohort=labels[k],
    n_badgers=length(ix),
    median_age_first_capture=median((first_live[ix]-birth[ix])*.25),
    median_SEI_followup_quarters=median(pmax(0,end[ix]-landmark[ix]+1L)),
    mean_homozygosity=mean(x$hom[ix]),
    male_fraction=mean(x$sex[ix]=="MALE"),
    median_live_captures=median(capture_q[ix]),
    proportion_without_post_entry_followup=mean(landmark[ix]>end[ix])
  )
}))
stopifnot(sum(cohort_summary$n_badgers)==N)
if(any(cohort_summary$n_badgers==0L))
  cat("Decades with zero first captures:",
      paste(cohort_summary$capture_cohort[
        cohort_summary$n_badgers==0L],collapse=", "), "\n")
cat("\nFirst-capture cohort composition:\n")
print(cohort_summary,row.names=FALSE,digits=4)

# Per draw and cohort, count events and total person-time at risk.
by_draw <- vector("list",D*5L*2L)
by_event_calendar <- vector("list",D*5L*2L)
idx <- 0L
for(draw in seq_len(D)){
  e <- E[,draw]; inf <- I[,draw]
  if(anyNA(e)||anyNA(inf)||
     any(inf>0L & (e==0L|inf<e))||
     any(e>last & e>0L)||any(inf>last & inf>0L))
    stop("Unusable SEI history in draw ",draw)
  for(trans in c("SE","EI")){
    if(trans=="SE"){
      eligible <- landmark<=end & !(e>0L & e<landmark)
      entry <- landmark
      event <- eligible & e>0L & e>=entry & e<=end
      exit <- ifelse(event,e,end)
      event_q <- e
    } else {
      eligible <- landmark<=end & e>0L & e<=end &
        (inf==0L | inf>=landmark)
      entry <- pmax(e,landmark)
      event <- eligible & inf>0L & inf>=entry & inf<=end
      exit <- ifelse(event,inf,end)
      event_q <- inf
    }
    duration <- ifelse(eligible,exit-entry+1L,0L)
    if(any(duration<0L)|anyNA(duration))
      stop("Negative/NA risk time at draw ",draw," transition ",trans)
    for(k in bins){
      ix <- cohort==k
      events <- sum(event[ix])
      years <- sum(duration[ix])*.25
      idx <- idx+1L
      by_draw[[idx]] <- data.frame(
        draw=draw,transition=trans,
        capture_cohort=labels[k],n_at_risk=sum(eligible[ix]),
        events=events,risk_years=years,
        crude_events_per_100_years=if(years>0)events/years*100 else NA_real_)
      by_event_calendar[[idx]] <- data.frame(
        draw=draw,transition=trans,
        event_calendar_period=labels[k],
        events=sum(event & quarter_bin(event_q)==k))
    }
  }
  if(draw%%250L==0L)cat("Audited",draw,"/",D,"joint trajectories\n")
}
risk <- do.call(rbind,by_draw)
events <- do.call(rbind,by_event_calendar)
out <- "outputs/sei_transitions"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(cohort_summary,file.path(out,"stage72_cohort_composition.csv"),
          row.names=FALSE)
write.csv(risk,file.path(out,"stage72_crude_rates_by_draw.csv"),
          row.names=FALSE)
write.csv(events,file.path(out,"stage72_event_calendar_by_draw.csv"),
          row.names=FALSE)
q3 <- function(v){
  # A zero-risk cohort cannot have a rate: retain NA for that decade.
  # If only some posterior histories have risk time, report intervals
  # over defined rates and separately count how many were defined.
  v <- v[is.finite(v)]
  if(!length(v)) return(rep(NA_real_,3))
  unname(quantile(v,c(.025,.5,.975)))
}
rates <- do.call(rbind,lapply(c("SE","EI"),function(t){
  do.call(rbind,lapply(labels,function(k){
    v <- risk[risk$transition==t & risk$capture_cohort==k,]
    ev <- q3(v$events)
    ps <- q3(v$risk_years)
    rr <- q3(v$crude_events_per_100_years)
    data.frame(transition=t,capture_cohort=k,
               median_events=ev[2],
               events_2.5pct=ev[1],events_97.5pct=ev[3],
               median_risk_years=ps[2],
               median_crude_events_per_100_years=rr[2],
               rate_2.5pct=rr[1],rate_97.5pct=rr[3],
               draws_with_nonzero_risk=sum(v$risk_years>0),
               row.names=NULL)
  }))
}))
events_summary <- do.call(rbind,lapply(c("SE","EI"),function(t){
  do.call(rbind,lapply(labels,function(k){
    v <- events[events$transition==t &
                events$event_calendar_period==k,"events"]
    q <- q3(v)
    data.frame(transition=t,event_calendar_period=k,
               median_events=q[2],lo=q[1],hi=q[3])
  }))
}))
write.csv(rates,file.path(out,"stage72_cohort_crude_rate_summary.csv"),
          row.names=FALSE)
write.csv(events_summary,file.path(out,"stage72_event_calendar_summary.csv"),
          row.names=FALSE)
cat("\nCRUDE events per 100 years, stratified by FIRST CAPTURE decade\n")
print(rates,row.names=FALSE,digits=4)
cat("\nInferred event calendar quarter counts (no person-time denominator)\n")
print(events_summary,row.names=FALSE,digits=4)
cat("\nThese are DESCRIPTIVE rates, not age-/sex-adjusted models.\n")
cat("Calendar-entry covariate in Stage 69 differs between SE and EI:\n")
cat("SE uses FIRST at-risk quarter; EI uses EXPOSURE/first at-risk quarter.\n")
cat("Changing test sensitivity, effort, cohort selection, age effects,",
    "and nonrandom censoring can affect the patterns.\n")
cat("E->I here uses age as Cox time scale; duration since exposure",
    "is a different biological time scale and requires sensitivity work.\n")
cat("No new conclusions about causal time trends or homozygosity.\n")
cat("STAGE 72 COMPLETE. Outputs:",out,"\n")
