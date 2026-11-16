# ====================================================================
# Stage 67: Biological reconnaissance of the SEI joint posterior
# WITHOUT fitting a survival model or double-counting diagnostic tests.
# All 1000 draw columns are coherent across all matched badgers.
#
# Progression by last capture is observationally window-limited.
# Comparisons by homozygosity are DESCRIPTIVE, not causal effects.
# ====================================================================
options(stringsAsFactors=FALSE)
file <- "data/model_inputs/integrated_gm_sei_aligned_2026.rds"
original <- "data/model_inputs/integrated_gm_inputs_group_proxy_2026.rds"
out <- "outputs/sei_survival"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(file))stop("Missing SEI alignment: run Stage 66")
if(!file.exists(original))stop("Missing integrated survival input")
z <- readRDS(file)
x <- readRDS(original)
stopifnot(identical(z$id,toupper(trimws(x$id))),
          identical(as.integer(z$draw_ids),1:1000),
          all(dim(z$exposed_onset_q)==c(x$N,1000)),
          all(dim(z$infectious_onset_q)==c(x$N,1000)),
          identical(z$survival_input_md5,unname(tools::md5sum(original))))
E <- z$exposed_onset_q
I <- z$infectious_onset_q
N <- x$N
D <- ncol(E)
cat("\n=== STAGE 67: BIOLOGICAL INFORMATION IN SEI TRAJECTORIES ===\n")
cat("Matched badgers:",N,"joint disease-history draws:",D,"\n")
cat("All onsets only define progression through SEI last capture.\n")

ever_E <- E>0L
ever_I <- I>0L
n_E <- colSums(ever_E)
n_I <- colSums(ever_I)
n_EI <- colSums(ever_E & ever_I)
prob_E <- rowMeans(ever_E)
prob_I <- rowMeans(ever_I)
years <- (z$sei_last_capture_q-z$sei_start_q+1L)*.25
mean_follow <- mean(years)
cat("Mean documented disease-history window in years:",
    round(mean_follow,2),"\n")
fmt <- function(v) {
  q <- quantile(v,c(.025,.5,.975))
  sprintf("median %0.2f (95%% interval %0.2f-%0.2f)",
          q[2],q[1],q[3])
}
cat("Animals with E onset by last capture, posterior:",fmt(n_E),"\n")
cat("Animals with I onset by last capture, posterior:",fmt(n_I),"\n")
cat("Among observed onsets, share progressing to I by last capture:",
    fmt(n_I/pmax(1,n_E)),"\n")
cat("The last percentage is not an eventual progression probability.\n")

# Report descriptive associations with homozygosity while carrying
# the SHARED 1000 draw columns through all cross-badger summaries.
hom <- (x$hom-mean(x$hom))/sd(x$hom)
hi <- hom >= median(hom)
lo <- !hi
capture_span <- z$sei_last_capture_q-z$sei_start_q+1L
cat("\nHomozygosity strata (descriptive only):\n")
cat("High-hom badgers:",sum(hi),
    "; lower-hom badgers:",sum(lo),"\n")
cat("Median documented quarters, high/lower hom:",
    median(capture_span[hi]),"/",median(capture_span[lo]),"\n")
diff_E <- colMeans(ever_E[hi,,drop=FALSE]) -
  colMeans(ever_E[lo,,drop=FALSE])
diff_I <- colMeans(ever_I[hi,,drop=FALSE]) -
  colMeans(ever_I[lo,,drop=FALSE])
cat("High minus lower hom; E onset by LAST capture:",fmt(diff_E),"\n")
cat("High minus lower hom; I onset by LAST capture:",fmt(diff_I),"\n")
cat("These contrasts conflate infection, progression, different",
    "capture windows, age and differential survival.\n")

# Death-recovery evidence is often after last capture, when SEI
# trajectories provide NO valid stage assignment.
known_death <- z$death_q>0L
post_SEI_death <- known_death & z$death_q>z$sei_last_capture_q
cat("\n=== MORTALITY/DISEASE FOLLOW-UP GAP ===\n")
cat("Known deaths in matched survival cohort:",sum(known_death),"\n")
cat("Deaths after SEI last capture:",sum(post_SEI_death),"\n")
cat("Deaths within recorded SEI state window:",
    sum(known_death & !post_SEI_death),"\n")
cat("Fraction of deaths beyond SEI state window:",
    if(sum(known_death)>0L) {
      round(mean(post_SEI_death[known_death]),3)
    } else {
      NA_real_
    },"\n")
cat("Do NOT treat zero onset as uninfected after last capture.\n")

# Current SEI-state posterior probabilities at the disease end
# are meaningful; disease states at later death are not.
state_at_last <- list(
  susceptible=1-prob_E,
  exposed=rowMeans(ever_E & !ever_I),
  infectious=prob_I
)
stopifnot(all(abs(state_at_last$susceptible+
                  state_at_last$exposed+
                  state_at_last$infectious-1)<1e-12))
per_badger <- data.frame(
  tattoo=z$id,
  sei_last_q=z$sei_last_capture_q,
  followup_years=years,
  known_death=known_death,
  death_after_sei_last=post_SEI_death,
  p_ever_exposed_by_last=prob_E,
  p_ever_infectious_by_last=prob_I,
  p_susceptible_at_last=state_at_last$susceptible,
  p_exposed_at_last=state_at_last$exposed,
  p_infectious_at_last=state_at_last$infectious,
  standardized_homozygosity=hom
)
write.csv(per_badger,file.path(out,"sei_individual_state_summary.csv"),
          row.names=FALSE)
by_draw <- data.frame(draw=seq_len(D),
                      n_exposed_by_last=n_E,
                      n_infectious_by_last=n_I,
                      fraction_infectious_among_exposed=n_I/pmax(1,n_E),
                      high_minus_lower_hom_E=diff_E,
                      high_minus_lower_hom_I=diff_I)
write.csv(by_draw,file.path(out,"sei_joint_draw_summary.csv"),
          row.names=FALSE)
cat("\nSaved individual-level and draw-level SEI summaries in ",out,"\n",sep="")
cat("STAGE 67 COMPLETE: descriptive only, NO SEI mortality model fitted.\n")
