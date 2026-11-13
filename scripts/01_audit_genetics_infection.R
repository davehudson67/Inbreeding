# ============================================================
# INBREEDING REANALYSIS 2026
# Stage 01: Genetic and infection trajectory audit
#
# Read-only with respect to original research data.
# Writes reports to outputs/data_audit (Git-ignored).
#
# Assumptions:
#   infection index 0 = never infected during model period
#   infection index 1 = 1976 Q1
#   4 quarters per year
#   posterior columns must retain original order and IDs
# ============================================================

options(stringsAsFactors = FALSE)

root <- path.expand("~/Inbreeding")
out <- file.path(root, "outputs/data_audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================\n")
cat("STAGE 01: GENETICS + INFECTION AUDIT\n")
cat("============================================\n")

# ------------------------------------------------------------
# 1. Verify retained files
# ------------------------------------------------------------

genetic_dir <- file.path(root, "data/badger/genetics")

historical_file <- file.path(
  root, "data/badger/legacy_published",
  "published_2023_analytical_inputs.rds"
)

gene_file <- file.path(
  genetic_dir, "legacy_workspace_gene.rds"
)

infection_file <- file.path(
  root, "data/infection_posteriors",
  "badger_infection_trajectories_all_tests_inferred.rds"
)

required <- c(historical_file, gene_file, infection_file)

missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Missing required files:\n",
       paste(missing, collapse = "\n"))
}

# ------------------------------------------------------------
# 2. Load historical genetic summaries
# ------------------------------------------------------------

gene <- readRDS(gene_file)
published <- readRDS(historical_file)

stopifnot(
  is.data.frame(gene),
  all(c("ID", "mlh") %in% names(gene)),
  "mlh" %in% names(published)
)

cat("\nGENETIC SUMMARY\n")
cat("Historical gene rows:", nrow(gene), "\n")
cat("Historical MLH rows:", nrow(published$mlh), "\n")

cat("\nGene columns:\n")
print(names(gene))

cat("\nPublished MLH columns:\n")
print(names(published$mlh))

cat("\nPublished CH columns:\n")
print(names(published$CH))

# ------------------------------------------------------------
# 3. Inspect original genotype file candidates
# ------------------------------------------------------------

genetic_files <- list.files(
  genetic_dir,
  full.names = TRUE,
  recursive = TRUE
)

genetic_files <- genetic_files[
  grepl(
    "gene_data|newgenedata|gene.info|hetinbreed|loci",
    basename(genetic_files),
    ignore.case = TRUE
  )
]

genetic_inventory <- data.frame(
  filename = basename(genetic_files),
  size_bytes = file.info(genetic_files)$size,
  extension = tolower(tools::file_ext(genetic_files))
)

write.csv(
  genetic_inventory,
  file.path(out, "genetic_file_inventory.csv"),
  row.names = FALSE
)

cat("\nGENOTYPE FILE CANDIDATES\n")
print(genetic_inventory, row.names = FALSE)

# ------------------------------------------------------------
# 4. Load full posterior infection trajectories
# ------------------------------------------------------------

infection <- readRDS(infection_file)

required_names <- c(
  "tattoo", "draws", "infection_time", "start_year"
)

stopifnot(
  all(required_names %in% names(infection)),
  is.matrix(infection$infection_time),
  nrow(infection$infection_time) == length(infection$tattoo),
  ncol(infection$infection_time) == length(infection$draws),
  identical(
    rownames(infection$infection_time),
    infection$tattoo
  ),
  identical(
    colnames(infection$infection_time),
    as.character(infection$draws)
  ),
  infection$start_year == 1976L,
  !anyNA(infection$infection_time),
  all(infection$infection_time >= 0),
  !anyDuplicated(infection$tattoo),
  !anyDuplicated(infection$draws)
)

times <- infection$infection_time

cat("\nINFECTION TRAJECTORIES\n")
cat("Badgers:", nrow(times), "\n")
cat("Posterior draws:", ncol(times), "\n")
cat("Minimum index:", min(times), "\n")
cat("Maximum index:", max(times), "\n")

cat(
  "Original posterior columns preserved:",
  identical(colnames(times), as.character(infection$draws)),
  "\n"
)

# ------------------------------------------------------------
# 5. Time-index conversion
# ------------------------------------------------------------

index_to_year <- function(index, start_year = 1976L) {
  ifelse(index == 0L, NA_integer_,
         start_year + (index - 1L) %/% 4L)
}

index_to_quarter <- function(index) {
  ifelse(index == 0L, NA_integer_,
         ((index - 1L) %% 4L) + 1L)
}

indices <- seq_len(max(times))

time_lookup <- data.frame(
  infection_index = indices,
  year = index_to_year(indices),
  quarter = index_to_quarter(indices)
)

write.csv(
  time_lookup,
  file.path(out, "infection_time_lookup.csv"),
  row.names = FALSE
)

cat("\nTIME INDEX CHECK\n")
print(head(time_lookup, 8))
print(tail(time_lookup, 8))

# ------------------------------------------------------------
# 6. Per-individual infection probability
# ------------------------------------------------------------

p_infected <- rowMeans(times > 0L)

cat("\nPOSTERIOR PROBABILITY EVER INFECTED\n")
print(summary(p_infected))

probability_bins <- cut(
  p_infected,
  breaks = c(-0.001, 0, 0.1, 0.5, 0.9, 0.999, 1),
  include.lowest = TRUE
)

print(table(probability_bins))

# ------------------------------------------------------------
# 7. Compare individual identifiers
# ------------------------------------------------------------

normalise_id <- function(x) {
  toupper(trimws(as.character(x)))
}

genetic_ids <- normalise_id(gene$ID)
infection_ids <- normalise_id(infection$tattoo)

valid_genetic <- !is.na(genetic_ids) &
  nzchar(genetic_ids) & genetic_ids != "NA"

valid_infection <- !is.na(infection_ids) &
  nzchar(infection_ids) & infection_ids != "NA"

matched <- intersect(
  unique(genetic_ids[valid_genetic]),
  unique(infection_ids[valid_infection])
)

cat("\nIDENTIFIER MATCHING\n")
cat("Unique valid genetic IDs:",
    length(unique(genetic_ids[valid_genetic])), "\n")

cat("Unique valid infection IDs:",
    length(unique(infection_ids[valid_infection])), "\n")

cat("Matched IDs:", length(matched), "\n")

cat("Duplicate genetic IDs:",
    sum(duplicated(genetic_ids[valid_genetic])), "\n")

# Never modify the original posterior matrix.
# Row selection retains ALL columns and their ordering.

matched_original <- infection$tattoo[
  infection_ids %in% matched
]

matched_times <- times[
  matched_original, ,
  drop = FALSE
]

stopifnot(
  identical(colnames(matched_times), colnames(times))
)

cat("Matched matrix dimensions:",
    paste(dim(matched_times), collapse = " x "), "\n")

# ------------------------------------------------------------
# 8. Write a small, non-identifying audit summary
# ------------------------------------------------------------

audit <- data.frame(
  metric = c(
    "n_genetic_rows",
    "n_published_mlh_rows",
    "n_infection_badgers",
    "n_posterior_draws",
    "n_matched_badgers",
    "n_duplicate_genetic_ids",
    "infection_start_year",
    "max_infection_index",
    "mean_posterior_infection_probability"
  ),
  value = as.character(c(
    nrow(gene),
    nrow(published$mlh),
    nrow(times),
    ncol(times),
    length(matched),
    sum(duplicated(genetic_ids[valid_genetic])),
    infection$start_year,
    max(times),
    mean(p_infected)
  ))
)

write.csv(
  audit,
  file.path(out, "stage01_summary.csv"),
  row.names = FALSE
)

cat("\n============================================\n")
cat("STAGE 01 COMPLETED\n")
cat("Reports:", out, "\n")
cat("Original genetic and infection files unchanged\n")
cat("No MCMC fitting performed\n")
cat("============================================\n")
