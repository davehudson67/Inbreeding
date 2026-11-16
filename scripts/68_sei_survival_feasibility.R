# =====================================================================
# Stage 68: Is there adequate SEI-stage and death evidence to motivate
# an integrated SEI survival model? READ ONLY, no model fitting.
#
# Uses already aligned coherent 1000 shared posterior trajectory draws.
# Stage 67's 85 deaths within the SEI end quarter are NOT assumed to
# be valid stage-at-death observations until we also check SEI start.
#
# Caveats: stages are inferred from diagnostic histories, not direct
# measured disease states; posterior summaries alone do not identify
# mortality contrasts, and deaths after last capture are unclassified.
# =====================================================================
options(stringsAsFactors=FALSE)
file <- "data/model_inputs/integrated_gm_sei_aligned_2026.rds"
input <- "data/model_inputs/integrated_gm_inputs_group_proxy_2026.rds"
out <- "outputs/sei_survival"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(file) || !file.exists(input))
  stop("Run Stage 66 to align SEI trajectories before Stage 68.")
z <- readRDS(file)
x <- readRDS(input)
stopifnot(x$N==z$N, z$N==1603L,
          identical(as.character(z$id),
                    toupper(trimws(as.character(x$id)))),
          identical(z$survival_input_md5,unname(tools::md5sum(input))))
E <- z$exposed_onset_q
I <- z$infectious_onset_q
N <- nrow(E); D <- ncol(E)
stopifnot(identical(dim(E),dim(I)),D==1000L,
          length(z$sei_start_q)==N,length(z$sei_last_capture_q)==N)
start <- as.integer(z$sei_start_q)
last <- as.integer(z$sei_last_capture_q)
death <- as.integer(z$death_q)
stopifnot(all(!is.na(start)),all(!is.na(last)),
          all(!is.na(death)),all(death>=0),
          all(start<=last),all(last<=200),
          all(I==0L | (I>=E & E>0L)))
cat("\n=== STAGE 68: IS SEI-STAGE MORTALITY FEASIBLE? ===\n")
cat("Badgers:",N,"coherent joint SEI draws:",D,"\n")
cat("Dataset source:",z$source_mode,"\n")
cat("Known deaths:",sum(death>0L),"\n")

# Do not conflate "death at/before SEI last" with "death during a
# disease state-supported SEI observation window".
d_in <- which(death>0L & death>=start & death<=last)
d_before <- which(death>0L & death<start)
d_after <- which(death>last)
cat("\nDeaths by relationship to SEI observation window:\n")
cat("  Inside start..last (stage may be inferred):",length(d_in),"\n")
cat("  Before SEI start (state unavailable):",length(d_before),"\n")
cat("  After SEI last capture (state unavailable):",length(d_after),"\n")
stopifnot(length(d_in)+length(d_before)+length(d_after)==sum(death>0L))

# A quarter-specific state is determined FROM EACH draw's paired E/I
# onset quarters; columns stay coherent across individuals.
stage_at <- function(indices,quarters){
  n <- length(indices)
  if(n==0L) return(list(
    S=matrix(FALSE,nrow=0,ncol=D),
    E=matrix(FALSE,nrow=0,ncol=D),
    I=matrix(FALSE,nrow=0,ncol=D)))
  q <- matrix(as.integer(quarters),nrow=n,ncol=D)
  e <- E[indices,,drop=FALSE]
  inf <- I[indices,,drop=FALSE]
  active_I <- inf>0L & inf<=q
  active_E <- e>0L & e<=q & !active_I
  active_S <- !active_I & !active_E
  stopifnot(all(active_S+active_E+active_I==1))
  list(S=active_S,E=active_E,I=active_I)
}
summ <- function(v){
  qs <- quantile(v,c(.025,.5,.975),names=FALSE,na.rm=FALSE)
  c(median=qs[2],lower95=qs[1],upper95=qs[3])
}
count_stage <- function(label,indices,quarters){
  if(!length(indices)) return(NULL)
  ss <- stage_at(indices,quarters)
  vals <- lapply(ss,colSums)
  ans <- do.call(rbind,lapply(names(vals),function(k){
    q <- summ(vals[[k]])
    data.frame(window=label,stage=k,n_animals=length(indices),
               median=q[1],lower95=q[2],upper95=q[3],
               posterior_mean=mean(vals[[k]]))
  }))
  cat("\n",label," (",length(indices)," badgers):\n",sep="")
  print(ans[,c("stage","median","lower95","upper95")],
        row.names=FALSE,digits=4)
  list(aggregate=ans,by_draw=vals,by_animal=ss)
}
inside <- count_stage("DEATH quarter inside SEI window",d_in,death[d_in])
# IMPORTANT: stage at LAST capture for recovered deaths afterwards
# is not stage at death and cannot be used as if it were.
after <- count_stage("LAST SEI capture prior to a later death",
                     d_after,last[d_after])
others <- which(death==0L)
alive <- count_stage("LAST SEI capture (no known death)",
                     others,last[others])
summaries <- Filter(Negate(is.null),list(inside,after,alive))
write.csv(do.call(rbind,lapply(summaries,function(a)a$aggregate)),
          file.path(out,"sei_stage_death_feasibility.csv"),
          row.names=FALSE)

# Save posterior per-animal stage probabilities for deaths that lie
# inside the valid SEI observation window, not just pooled counts.
if(length(d_in)){
  states <- inside$by_animal
  person <- data.frame(
    tattoo=z$id[d_in],
    disease_start_q=start[d_in],
    disease_last_q=last[d_in],
    death_q=death[d_in],
    p_susceptible_at_death=rowMeans(states$S),
    p_exposed_at_death=rowMeans(states$E),
    p_infectious_at_death=rowMeans(states$I))
  write.csv(person,file.path(out,"sei_within_window_death_stage.csv"),
            row.names=FALSE)
  # How many individuals strongly support an infectious death state?
  cat("\nWithin-window deaths with posterior Pr(I at death) >= 0.5:",
      sum(person$p_infectious_at_death>=0.5),"\n")
  cat("Within-window deaths with posterior Pr(I at death) >= 0.9:",
      sum(person$p_infectious_at_death>=0.9),"\n")
  cat("Within-window deaths with posterior Pr(E or I at death) >= 0.5:",
      sum(person$p_exposed_at_death+
          person$p_infectious_at_death>=0.5),"\n")
}

# Disease progression among *observed* progressors, not eventual
# transition risk or time-to-progression for all exposed animals.
progressed <- I>0L & E>0L
delays <- (I[progressed]-E[progressed])*.25
cat("\nSEI progression timing (only E->I within follow-up):\n")
cat("E->I observed transitions per draw:\n")
print(summ(colSums(progressed)),digits=4)
if(length(delays)){
  cat("Time from E to I among observed transitions, years:\n")
  print(summ(delays),digits=4)
  cat("Proportion with same-quarter E and I onset:",
      mean(I[progressed]==E[progressed]),"\n")
}
# Stage 67's split at median included every tie in the high group.
hom <- as.numeric(x$hom)
cat("\nHomozygosity distribution QC:\n")
cat("Unique homozygosity values:",length(unique(hom)),"\n")
cat("Proportion of badgers exactly at the median:",
    mean(hom==median(hom)),"\n")
cat("Median-split group sizes (>= median, < median):",
    sum(hom>=median(hom)),sum(hom<median(hom)),"\n")
cat("Do not interpret imbalanced median-split contrasts as adjusted",
    "effects of homozygosity on disease progression.\n")
cat("\nINTERPRETATION GATE: No SEI survival association has been",
    "estimated. Disease stage after the SEI end quarter is UNKNOWN.\n")
cat("Saved Stage 68 CSVs to",out,"\n")
