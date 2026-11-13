# ============================================================
# INBREEDING REANALYSIS 2026
# Stage 02: Historical genotype and MLH audit
#
# No original data modified.
# No genetic measures recalculated yet.
# ============================================================

options(stringsAsFactors = FALSE)

root <- path.expand("~/Inbreeding")
out <- file.path(root, "outputs/genetic_audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

normalise_id <- function(x) {
  toupper(trimws(as.character(x)))
}

gene <- readRDS(file.path(
  root, "data/badger/genetics",
  "legacy_workspace_gene.rds"
))

published <- readRDS(file.path(
  root, "data/badger/legacy_published",
  "published_2023_analytical_inputs.rds"
))

infection <- readRDS(file.path(
  root, "data/infection_posteriors",
  "badger_infection_trajectories_all_tests_inferred.rds"
))

CH <- published$CH
mlh <- published$mlh

stopifnot(
  "ID" %in% names(gene),
  "ID" %in% names(CH),
  "ID" %in% names(mlh),
  is.matrix(infection$infection_time)
)

# ------------------------------------------------------------
# 1. Individual identifiers
# ------------------------------------------------------------

ids_gene <- normalise_id(gene$ID)
ids_ch <- normalise_id(CH$ID)
ids_mlh <- normalise_id(mlh$ID)
ids_inf <- normalise_id(infection$tattoo)

cat("\n=== INDIVIDUAL IDENTIFIERS ===\n")

counts <- data.frame(
  dataset = c(
    "Genetic summary",
    "Published CH",
    "Published MLH",
    "Infection trajectories"
  ),
  rows = c(
    length(ids_gene),
    length(ids_ch),
    length(ids_mlh),
    length(ids_inf)
  ),
  unique_ids = c(
    length(unique(ids_gene)),
    length(unique(ids_ch)),
    length(unique(ids_mlh)),
    length(unique(ids_inf))
  )
)

print(counts, row.names = FALSE)

unmatched <- setdiff(unique(ids_gene), unique(ids_inf))

cat("\nGenetic IDs without infection trajectories:\n")
print(unmatched)

write.csv(
  data.frame(ID = unmatched),
  file.path(out, "unmatched_genetic_ids.csv"),
  row.names = FALSE
)

cat("\nMatching summary:\n")
cat("Genetic / infection:",
    length(intersect(ids_gene, ids_inf)), "\n")
cat("Genetic / published CH:",
    length(intersect(ids_gene, ids_ch)), "\n")
cat("Genetic / published MLH:",
    length(intersect(ids_gene, ids_mlh)), "\n")

# ------------------------------------------------------------
# 2. Allele columns in published CH
# ------------------------------------------------------------

allele_cols <- grep(
  "\\.[ab]$",
  names(CH),
  value = TRUE
)

loci <- unique(sub("\\.[ab]$", "", allele_cols))

complete_loci <- loci[
  paste0(loci, ".a") %in% names(CH) &
  paste0(loci, ".b") %in% names(CH)
]

cat("\n=== GENOTYPE STRUCTURE ===\n")
cat("Published CH rows:", nrow(CH), "\n")
cat("Allele columns:", length(allele_cols), "\n")
cat("Complete locus pairs:", length(complete_loci), "\n")

cat("\nLoci:\n")
print(complete_loci)

locus_summary <- do.call(
  rbind,
  lapply(complete_loci, function(locus) {

    a <- as.character(CH[[paste0(locus, ".a")]])
    b <- as.character(CH[[paste0(locus, ".b")]])

    a <- trimws(a)
    b <- trimws(b)

    # This counts physical presence only.
    # It does NOT yet interpret zero or other missing codes.
    present <- !is.na(a) & !is.na(b) &
      nzchar(a) & nzchar(b)

    data.frame(
      locus = locus,
      n_rows = length(a),
      n_both_nonblank = sum(present),
      n_equal_nonblank = sum(a[present] == b[present]),
      n_different_nonblank = sum(a[present] != b[present]),
      n_zero_in_either = sum(
        (!is.na(a) & a == "0") |
        (!is.na(b) & b == "0")
      )
    )
  })
)

print(locus_summary, row.names = FALSE)

write.csv(
  locus_summary,
  file.path(out, "locus_structure.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. Most frequent allele codes
# ------------------------------------------------------------

allele_values <- unlist(
  lapply(allele_cols, function(n) {
    trimws(as.character(CH[[n]]))
  }),
  use.names = FALSE
)

cat("\n=== COMMON ALLELE VALUES ===\n")

print(
  head(
    sort(table(allele_values, useNA = "ifany"),
         decreasing = TRUE),
    25
  )
)

# ------------------------------------------------------------
# 4. Compare historical MLH measurements
# ------------------------------------------------------------

cat("\n=== HISTORICAL GENETIC MEASURES ===\n")

cat("\nGenetic summary MLH:\n")
print(summary(gene$mlh))

for (variable in c("IR", "SH", "HL", "Fbar", "hom")) {
  if (variable %in% names(mlh)) {
    cat("\n", variable, ":\n", sep = "")
    print(summary(mlh[[variable]]))
  }
}

# Compare numeric genetic measures between the
# published CH and published MLH tables by ID.

if (!anyDuplicated(ids_ch) && !anyDuplicated(ids_mlh)) {

  idx <- match(ids_ch, ids_mlh)

  comparison <- data.frame(
    ID = ids_ch,
    matched_mlh = !is.na(idx)
  )

  common <- intersect(
    c("IR", "SH", "HL", "Fbar"),
    intersect(names(CH), names(mlh))
  )

  for (variable in common) {

    x <- suppressWarnings(as.numeric(CH[[variable]]))
    y <- suppressWarnings(
      as.numeric(mlh[[variable]][idx])
    )

    ok <- is.finite(x) & is.finite(y)

    cat("\n", variable, "comparison:\n", sep = "")
    cat("Comparable rows:", sum(ok), "\n")

    if (sum(ok) > 1) {
      cat("Maximum absolute difference:",
          max(abs(x[ok] - y[ok])), "\n")
    }

    comparison[[paste0(variable, "_difference")]] <- x - y
  }

  write.csv(
    comparison,
    file.path(out, "published_genetic_measure_comparison.csv"),
    row.names = FALSE
  )

} else {
  cat("\nDuplicate IDs detected; direct comparisons skipped.\n")
}

# ------------------------------------------------------------
# 5. Preserve infection draw alignment
# ------------------------------------------------------------

times <- infection$infection_time

matched_ids <- infection$tattoo[
  normalise_id(infection$tattoo) %in% ids_gene
]

matched_matrix <- times[
  matched_ids, ,
  drop = FALSE
]

stopifnot(
  identical(colnames(matched_matrix), colnames(times)),
  ncol(matched_matrix) == 500L
)

cat("\n=== POSTERIOR ALIGNMENT ===\n")
cat("Matched badgers:", nrow(matched_matrix), "\n")
cat("Preserved posterior columns:", ncol(matched_matrix), "\n")
cat("Column alignment verified: TRUE\n")

cat("\n=== STAGE 02 COMPLETED ===\n")
cat("Reports:", out, "\n")
cat("No input datasets modified\n")
