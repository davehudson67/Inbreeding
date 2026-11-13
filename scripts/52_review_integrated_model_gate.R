# ============================================================
# Stage 52: Final integrated model-input adjudication
# No fitting. No database writes. Audits facts that change
# scientific likelihood, rather than exploratory development.
# ============================================================
options(stringsAsFactors = FALSE)
stopifnot(
  file.exists("data/model_inputs/integrated_gm_inputs_2026.rds"),
  file.exists("outputs/model_audit/integrated_accepted_death_sources.csv")
)
x <- readRDS("data/model_inputs/integrated_gm_inputs_2026.rds")
published <- readRDS("data/badger/legacy_published/published_2023_analytical_inputs.rds")
hist <- published$CH
hist$ID <- toupper(trimws(as.character(hist$ID)))
hist <- hist[!duplicated(hist$ID), , drop = FALSE]
orig_ids <- hist$ID
hist_dead <- setNames(!is.na(hist$death), hist$ID)
modern_dead <- setNames(!is.na(x$death_q), x$id)

source <- read.csv(
  "outputs/model_audit/integrated_accepted_death_sources.csv",
  colClasses = c(capture_id = "character"),
  na.strings = c("", "NA")
)
source$capture_date <- as.Date(source$capture_date)
source <- source[source$id %in% x$id, , drop = FALSE]

cat("\n=== HISTORICAL COHORT ATTRITION ===\n")
absent <- setdiff(orig_ids, x$id)
cat("Original unique CH IDs:", length(orig_ids), "\n")
cat("Included in current input:", length(x$id), "\n")
cat("Excluded historical IDs:", length(absent), "\n")
print(data.frame(
  id = absent,
  infection_posterior = absent %in%
    toupper(trimws(readRDS(
      "data/infection_posteriors/badger_infection_trajectories_all_tests_inferred.rds"
    )$tattoo)),
  historical_death = unname(hist_dead[absent]),
  stringsAsFactors = FALSE
), row.names = FALSE)

cat("\n=== DEATH CONCORDANCE BY ID (NOT JUST TOTALS) ===\n")
both <- x$id[hist_dead[x$id] & modern_dead[x$id]]
historical_only <- x$id[hist_dead[x$id] & !modern_dead[x$id]]
updated_only <- x$id[!hist_dead[x$id] & modern_dead[x$id]]
neither <- x$id[!hist_dead[x$id] & !modern_dead[x$id]]
concordance <- data.frame(
  comparison = c("historical_and_updated",
                 "historical_only",
                 "updated_only",
                 "neither"),
  n = c(length(both), length(historical_only),
        length(updated_only), length(neither))
)
print(concordance, row.names = FALSE)
cat("Historical-only IDs:",
    if(length(historical_only)) paste(historical_only, collapse=", ") else "(none)", "\n")
cat("Updated-only IDs:",
    if(length(updated_only)) paste(updated_only, collapse=", ") else "(none)", "\n")

cat("\n=== ACCEPTED TERMINAL DEATH CAUSES ===\n")
print(sort(table(source$cause_class, useNA="ifany"), decreasing=TRUE))
cat("Selected death source rows:", nrow(source), "\n")
cat("Distinct selected death IDs:", length(unique(source$id)), "\n")
if (nrow(source) != sum(modern_dead))
  stop("Selected death evidence count and cohort terminal deaths disagree")

# Audit single-quarter live trap + terminal death. Such an event can
# legitimately contain both an alive encounter and a subsequent death.
death_quarter_detected <- which(!is.na(x$death_q) &
  x$y[cbind(seq_len(x$N), ifelse(is.na(x$death_q), 1L, x$death_q))] == 1L)
cat("Died in same quarter as a live capture:",
    length(death_quarter_detected), "\n")
if(length(death_quarter_detected))
  print(data.frame(id=x$id[death_quarter_detected],
                   death_q=x$death_q[death_quarter_detected]),
        row.names=FALSE)

cat("\n=== HELD/REPEATED TERMINAL RECORDS ===\n")
held <- read.csv(
  "outputs/model_audit/integrated_held_death_records.csv",
  colClasses=c(capture_id="character")
)
repeat_deaths <- read.csv(
  "outputs/model_audit/integrated_repeat_death_candidates.csv",
  colClasses=c(capture_id="character")
)
print(held, row.names=FALSE)
cat("Repeated candidate rows:", nrow(repeat_deaths),
    "across", length(unique(repeat_deaths$id)), "IDs\n")
if(nrow(repeat_deaths))
  print(repeat_deaths[, c("id","capture_id","capture_date","moc","cause_class")],
        row.names=FALSE)

cat("\n=== TIME-VARYING INFECTION COMPATIBILITY ===\n")
dead_i <- which(!is.na(x$death_q))
n_later <- 0
n_positive <- 0
case_list <- list()
for(i in dead_i) {
  times <- x$infection_time[i, ]
  positives <- sum(times > 0L)
  later <- sum(times > x$death_q[i])
  n_later <- n_later + later
  n_positive <- n_positive + positives
  if(later > 0L)
    case_list[[length(case_list)+1L]] <-
      data.frame(
        id=x$id[i], death_q=x$death_q[i],
        posterior_draws_infection_after_death=later,
        positive_infection_draws=positives
      )
}
cat("Infection-positive draws among known-death badgers:", n_positive, "\n")
cat("Positive infection times AFTER recorded death:", n_later, "\n")
cat("Interpretation: times after recorded death cannot contribute",
    "infected person-time.\n")
if(length(case_list)) {
  postmortem <- do.call(rbind, case_list)
  write.csv(postmortem,
     "outputs/model_audit/integrated_postmortem_infection_draws.csv",
     row.names=FALSE)
  print(head(postmortem[order(-postmortem$posterior_draws_infection_after_death),], 12),
        row.names=FALSE)
}

cat("\n=== OBSERVATION-PROCESS CHECK ===\n")
noncaptured <- rowSums(x$y == 0L)
captured <- rowSums(x$y == 1L)
cat("Informative noncapture cells under current global proxy:",
    sum(noncaptured), "\n")
cat("Recapture cells:", sum(captured), "\n")
cat("Median noncaptures per animal:", median(noncaptured), "\n")
cat("Median detected quarters per animal:", median(captured), "\n")
cat("Nondetection/observed-detection ratio:",
    round(sum(noncaptured)/sum(captured), 2), "\n")
cat("NOTE: these noncaptures are only meaningful if animal",
    "was in monitored and available area.\n")

cat("\n=== FINAL REVIEW OUTPUT ===\n")
write.csv(
  concordance,
  "outputs/model_audit/integrated_historical_death_concordance.csv",
  row.names=FALSE
)
write.csv(
  data.frame(
    id=x$id,
    historical_death=unname(hist_dead[x$id]),
    updated_death=unname(modern_dead[x$id]),
    first_q=x$first_q,
    death_q=x$death_q,
    n_live_quarters=captured,
    n_proxy_noncaptures=noncaptured
  ),
  "outputs/model_audit/integrated_individual_model_gate.csv",
  row.names=FALSE
)
cat("Reports saved in outputs/model_audit\n")
cat("No model fitted, no data modified.\n")
