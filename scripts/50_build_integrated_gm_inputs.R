# ==============================================================
# V2026 integrated survival/recapture/infection inputs
# Read-only from Supabase; NEVER modifies original datasets.
# Output is the input to ONE intended integrated GM model.
# This builder is an audit gate, not a fitted development model.
# ==============================================================

options(stringsAsFactors = FALSE)
root <- normalizePath(getwd(), mustWork = TRUE)
stopifnot(file.exists(file.path(root, "R/dist_gompertz_makeham.R")))
dir.create("data/model_inputs", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/model_audit", recursive = TRUE, showWarnings = FALSE)
out <- "outputs/model_audit"
fail <- function(...) stop(paste0(...), call. = FALSE)
normalise <- function(x) toupper(trimws(as.character(x)))
q_index <- function(date) {
  date <- as.Date(date)
  (as.integer(format(date, "%Y")) - 1976L) * 4L +
    (as.integer(format(date, "%m")) - 1L) %/% 3L + 1L
}

stopifnot(
  requireNamespace("DBI", quietly = TRUE),
  requireNamespace("RPostgres", quietly = TRUE)
)

gene_path <- "data/badger/genetics/legacy_workspace_gene.rds"
ch_path <- "data/badger/legacy_published/published_2023_analytical_inputs.rds"
inf_path <- "data/infection_posteriors/badger_infection_trajectories_all_tests_inferred.rds"
if (!all(file.exists(c(gene_path, ch_path, inf_path))))
  fail("Required audited source RDS files are missing")

gene <- readRDS(gene_path)
CH <- readRDS(ch_path)$CH
inf <- readRDS(inf_path)
stopifnot(
  all(c("ID", "mlh") %in% names(gene)),
  all(c("ID", "birth", "last_seen", "death", "censored") %in% names(CH)),
  all(c("tattoo", "draws", "infection_time", "start_year") %in% names(inf)),
  is.matrix(inf$infection_time),
  identical(rownames(inf$infection_time), inf$tattoo),
  identical(colnames(inf$infection_time), as.character(inf$draws)),
  ncol(inf$infection_time) == 500L,
  inf$start_year == 1976L,
  !anyNA(inf$infection_time),
  all(inf$infection_time >= 0)
)
gene$ID <- normalise(gene$ID)
CH$ID <- normalise(CH$ID)
infection_ids <- normalise(inf$tattoo)
if (anyDuplicated(gene$ID) || anyDuplicated(infection_ids))
  fail("Duplicated genetic/infection IDs: manual identity review required")
if (nrow(gene) != 1933L)
  fail("Genetic cohort changed unexpectedly: ", nrow(gene))

# Historical duplicated IDs have identical survival histories;
# deliberately deduplicate only after checking that fact.
ch_vars <- c("birth", "last_seen", "death", "censored")
for (id in unique(CH$ID[duplicated(CH$ID)])) {
  z <- CH[CH$ID == id, ch_vars, drop = FALSE]
  if (nrow(unique(z)) != 1L)
    fail("Conflicting historical survival inputs for ", id)
}
CH <- CH[!duplicated(CH$ID), , drop = FALSE]
if (nrow(CH) != 1606L)
  fail("Historical unique-individual count unexpectedly changed")

if (!nzchar(Sys.getenv("BADGER_DB_USER")) ||
    !nzchar(Sys.getenv("BADGER_DB_PASSWORD")))
  fail("Set BADGER_DB_USER and BADGER_DB_PASSWORD in the R session")

con <- DBI::dbConnect(
  RPostgres::Postgres(),
  host = Sys.getenv("BADGER_DB_HOST", "db.klpvihubzlyzubyyubtw.supabase.co"),
  port = as.integer(Sys.getenv("BADGER_DB_PORT", "5432")),
  dbname = "postgres",
  user = Sys.getenv("BADGER_DB_USER"),
  password = Sys.getenv("BADGER_DB_PASSWORD"),
  sslmode = "require"
)
on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

# CAST bigint IDs TO TEXT IN SQL: coercing bit64 afterwards can corrupt IDs.
sql <- "
SELECT h.capture_id::text AS capture_id,
       upper(trim(h.tattoo)) AS id,
       h.capture_date,
       h.moc, h.is_pm, h.pm_cause,
       h.sex, h.age_first_capture,
       e.encounter_type, e.alive_when_encountered,
       e.encounter_outcome, e.resolved_episode_outcome,
       e.episode_record_role
FROM analysis.capture_history_v1 AS h
LEFT JOIN analysis.encounters_with_episodes_v1 AS e
  ON e.capture_id = h.capture_id
WHERE h.capture_date >= DATE '1976-01-01'
  AND h.capture_date < DATE '2026-01-01'
ORDER BY h.capture_id
"
cat("Querying historical-to-2025 encounters (read only)...\n")
events <- DBI::dbGetQuery(con, sql)
DBI::dbDisconnect(con)
if (anyDuplicated(events$capture_id))
  fail("Join duplicated capture IDs; fix source join before modelling")
if (!nrow(events))
  fail("No capture records returned")

events$id <- normalise(events$id)
events$moc_clean <- normalise(events$moc)
events$moc_clean[is.na(events$moc)] <- ""
events$episode_clean <- normalise(events$resolved_episode_outcome)
events$episode_clean[is.na(events$resolved_episode_outcome)] <- ""
events$state_clean <- normalise(events$alive_when_encountered)
events$state_clean[is.na(events$alive_when_encountered)] <- ""
events$q <- q_index(events$capture_date)
events$is_pm[is.na(events$is_pm)] <- FALSE

# Source-method overrides explicitly distinguish life status from the
# inherited view's UNKNOWN label (notably RTA, SHOT, FOUND EXT).
death_methods <- c(
  "RTA", "ON ROAD", "FOUND DEAD", "FOUND EXT",
  "SHOT", "SHOT (GP)", "REMOVAL", "PUT DOWN", "DEAD SN/TR"
)
death_outcomes <- c(
  "FOUND_DEAD", "KILLED", "EUTHANISED", "DIED",
  "DIED_AFTER_CAPTURE"
)
events$source_dead <- events$moc_clean %in% death_methods |
  events$episode_clean %in% death_outcomes |
  events$is_pm |
  events$state_clean == "DEAD"

# Only live traps/confirmed live handling count for recapture;
# SNARE alone and missing MOC do not get silently recoded alive.
live_methods <- c("TRAP", "NET", "DART GUN")
events$source_live <- (
  events$moc_clean %in% live_methods & !events$is_pm
) | (
  events$state_clean == "ALIVE" & !events$is_pm
)
events$source_live[is.na(events$source_live)] <- FALSE
events$source_dead[is.na(events$source_dead)] <- FALSE

# U005 death is contradicted by a later confirmed-live trap:
# do not destroy the source record or assign an observed death.
events$held_death <- events$id == "U005" &
  events$capture_id == "16968"
events$source_dead[events$held_death] <- FALSE

# Death evidence and an alive observation can occur in the SAME quarter
# (animal captured and later euthanised/shot). Retain both.
events$cause_class <- ifelse(
  events$moc_clean %in% c("RTA", "ON ROAD"), "RTA",
  ifelse(
    events$moc_clean %in% c("SHOT", "SHOT (GP)", "REMOVAL"),
    "KILLED",
    ifelse(
      events$moc_clean %in% c("PUT DOWN") |
        events$episode_clean == "EUTHANISED", "EUTHANISED",
      ifelse(events$source_dead, "DEAD_OTHER_OR_UNKNOWN", NA_character_)
    )
  )
)

# Any death followed by a live record is a held identity contradiction.
# Record rather than deciding silently which original event is false.
live <- events[events$source_live, , drop = FALSE]
dead <- events[events$source_dead, , drop = FALSE]
first_live_dates <- tapply(live$capture_date, live$id, min)
last_live_dates <- tapply(live$capture_date, live$id, max)
dead$later_live <- !is.na(last_live_dates[dead$id]) &
  as.Date(last_live_dates[dead$id], origin = "1970-01-01") >
  as.Date(dead$capture_date)
dead$later_live[is.na(dead$later_live)] <- FALSE
holds <- dead[dead$later_live, , drop = FALSE]
dead <- dead[!dead$later_live, , drop = FALSE]
events$held_death <- events$held_death |
  events$capture_id %in% holds$capture_id
write.csv(
  events[events$held_death, c("id","capture_id","capture_date","moc",
                               "resolved_episode_outcome")],
  file.path(out, "integrated_held_death_records.csv"),
  row.names = FALSE
)

# Multiple candidate deaths for one animal must remain visible in audit.
repeated_dead <- dead[duplicated(dead$id) |
  duplicated(dead$id, fromLast = TRUE), , drop = FALSE]
write.csv(
  repeated_dead[, c("id","capture_id","capture_date","moc","cause_class")],
  file.path(out, "integrated_repeat_death_candidates.csv"),
  row.names = FALSE
)
dead <- dead[order(dead$id, dead$capture_date), , drop = FALSE]
dead_first <- dead[!duplicated(dead$id), , drop = FALSE]


# Record the exact source record selected for each accepted death.
# This is needed to distinguish genuine mortality concordance from
# accidentally matching the historical total number of deaths.
selected_deaths <- dead_first[, c(
  "id", "capture_id", "capture_date", "moc", "is_pm",
  "pm_cause", "cause_class", "resolved_episode_outcome"
), drop = FALSE]
write.csv(
  selected_deaths,
  file.path(out, "integrated_accepted_death_sources.csv"),
  row.names = FALSE, na = ""
)

# Calendar evidence is a proxy for observation effort, NOT independently
# recorded trap effort. Do not infer a noncapture in a quarter with no
# qualifying trapping evidence anywhere in the study.
Tmax <- 200L  # 1976 Q1 -> 2025 Q4; avoid partial 2026 year.
calendar <- data.frame(
  q = seq_len(Tmax),
  year = 1976L + (seq_len(Tmax) - 1L) %/% 4L,
  season = (seq_len(Tmax) - 1L) %% 4L + 1L
)
effort_count <- tabulate(
  events$q[events$source_live & events$q >= 1L &
             events$q <= Tmax],
  nbins = Tmax
)
calendar$n_source_live_records <- effort_count
calendar$observed_trapping_proxy <- effort_count > 0L
write.csv(calendar, file.path(out, "integrated_effort_proxy.csv"),
          row.names = FALSE)

# Main age-verified cohort: original 1606 unique histories, whose
# historical birth-quarter estimates are available.
# All 1933 genetic IDs are retained in the attrition/audit records.
birth <- setNames(as.integer(CH$birth), CH$ID)
g_hom <- setNames(1 - as.numeric(gene$mlh), gene$ID)
unmatched_infection <- setdiff(gene$ID, infection_ids)
write.csv(
  data.frame(id = unmatched_infection),
  file.path(out, "integrated_missing_infection_ids.csv"),
  row.names = FALSE
)
ids <- CH$ID[CH$ID %in% gene$ID &
               CH$ID %in% infection_ids &
               CH$ID %in% names(first_live_dates)]
N <- length(ids)
if (N < 1200L)
  fail("Too many historical cohort IDs lost: N=", N)
inf_row <- match(ids, infection_ids)
first_dates <- as.Date(first_live_dates[ids], origin = "1970-01-01")
first_q <- q_index(first_dates)
birth_q <- unname(birth[ids])
if (anyNA(first_q) || anyNA(birth_q) || any(first_q < birth_q))
  fail("Missing or negative age at entry")

death_q <- q_index(as.Date(
  dead_first$capture_date[match(ids, dead_first$id)]
))
death_q <- as.integer(death_q)
if (any(!is.na(death_q) & death_q < first_q))
  fail("Death before first live encounter: manual resolution required")
if (any(!is.na(death_q) & death_q > Tmax))
  fail("Unexpected death beyond study horizon")

# y=1 confirmed-live quarter; y=0 observed non-capture in globally
# surveyed quarter; y=-1 no capture information / after known death.
# first observation is conditioned upon during model fitting.
y <- matrix(-1L, N, Tmax, dimnames = list(ids, NULL))
d <- matrix(0L, N, Tmax, dimnames = list(ids, NULL))
for (i in seq_len(N)) {
  active <- seq.int(first_q[i], Tmax)
  y[i, active[calendar$observed_trapping_proxy[active]]] <- 0L
  if (!is.na(death_q[i])) {
    d[i, death_q[i]] <- 1L
    if (death_q[i] < Tmax)
      y[i, (death_q[i] + 1L):Tmax] <- -1L
  }
}
live_cohort <- live[live$id %in% ids &
                      live$q >= 1L & live$q <= Tmax, , drop = FALSE]
ii <- match(live_cohort$id, ids)
bad_after_death <- !is.na(death_q[ii]) & live_cohort$q > death_q[ii]
if (any(bad_after_death))
  fail("Live observation after selected death: source reconciliation required")
y[cbind(ii, live_cohort$q)] <- 1L
if (!all(y[cbind(seq_len(N), first_q)] == 1L))
  fail("Entry detections missing")

# Keep the 500 draws AND their original aligned column order.
infection_time <- inf$infection_time[inf_row, , drop = FALSE]
rownames(infection_time) <- ids
if (!identical(colnames(infection_time), as.character(inf$draws)))
  fail("Infection posterior columns lost their original alignment")

# Individual sex data, only if source entries are consistent.
sex_value <- function(id) {
  x <- unique(normalise(events$sex[events$id == id]))
  x <- x[!is.na(x) & x != "" & x != "NA"]
  if (!length(x)) return(NA_character_)
  if (length(x) > 1L) return("CONFLICT")
  x
}
sex <- vapply(ids, sex_value, character(1))
sex_conflict <- ids[sex == "CONFLICT"]
write.csv(data.frame(id = sex_conflict),
          file.path(out, "integrated_sex_conflicts.csv"),
          row.names = FALSE)

# Audits: retained individual-by-individual provenance and attrition.
audit <- data.frame(
  id = gene$ID,
  has_original_birth = gene$ID %in% CH$ID,
  has_infection_draws = gene$ID %in% infection_ids,
  has_confirmed_live = gene$ID %in% names(first_live_dates),
  included_age_verified = gene$ID %in% ids
)
write.csv(audit, file.path(out, "integrated_cohort_attrition.csv"),
          row.names = FALSE)
saveRDS(
  list(
    id = ids, N = N, T = Tmax, hom = unname(g_hom[ids]),
    sex = sex, birth_q = birth_q, first_q = first_q,
    death_q = death_q, y = y, dead = d,
    calendar = calendar, infection_time = infection_time,
    infection_draw_ids = inf$draws,
    source_hashes = tools::md5sum(c(gene_path, ch_path, inf_path)),
    data_provenance = list(
      genetics = gene_path,
      historical_birth = ch_path,
      infection = inf_path,
      encounters = "analysis.capture_history_v1 + analysis.encounters_with_episodes_v1",
      study_end = "2025 Q4",
      effort = "global live-trapping encounter proxy (NOT independent effort)",
      source_method_death_override = death_methods,
      held_death = "U005 capture 16968 and any confirmed live-after-death records"
    )
  ),
  "data/model_inputs/integrated_gm_inputs_2026.rds"
)
cat("\n=== INTEGRATED DATASET BUILD ===\n")
cat("Genotyped individuals:", nrow(gene), "\n")
cat("Age-verified historical IDs:", nrow(CH), "\n")
cat("Included with infection and live observation:", N, "\n")
cat("Matched infection posterior draws:", ncol(infection_time), "\n")
cat("Live observations (individual-quarter):", sum(y == 1L), "\n")
cat("Known deaths retained:", sum(d == 1L), "\n")
cat("Death records held:", sum(events$held_death), "\n")
cat("Repeated death candidate rows:", nrow(repeated_dead), "\n")
cat("Sex conflicts:", length(sex_conflict), "\n")
cat("Quarters with study-wide live-trapping evidence:",
    sum(calendar$observed_trapping_proxy), "/", Tmax, "\n")
cat("Genetic IDs with no infection posterior:",
    paste(unmatched_infection, collapse = ", "), "\n")
cat("CHECK: study-wide observed trapping is only an effort PROXY.\n")
cat("CHECK: age-uncertain additional genetic IDs require separate inclusion rule.\n")
cat("READY FOR MODEL GATE REVIEW; model not yet fitted.\n")
