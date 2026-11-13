# ============================================================
# Integrated GM data and mathematics: pre-fit validation gate.
# Run after scripts/50_build_integrated_gm_inputs.R.
# Does NOT fit a model, change data or access another repository.
# ============================================================
options(stringsAsFactors = FALSE)
infile <- "data/model_inputs/integrated_gm_inputs_2026.rds"
if (!file.exists(infile)) stop("Run stage 50 first.")
x <- readRDS(infile)
fail <- function(...) stop(paste0(...), call. = FALSE)
cat("\n=== FINAL MODEL INPUT CHECKS ===\n")
stopifnot(
  length(x$id) == x$N,
  length(unique(x$id)) == x$N,
  x$T == 200L,
  is.matrix(x$y), identical(dim(x$y), c(x$N, x$T)),
  is.matrix(x$dead), identical(dim(x$dead), c(x$N, x$T)),
  is.matrix(x$infection_time),
  identical(dim(x$infection_time), c(x$N, 500L)),
  identical(colnames(x$infection_time), as.character(x$infection_draw_ids)),
  identical(rownames(x$infection_time), x$id),
  all(!is.na(x$hom) & x$hom >= 0 & x$hom <= 1),
  all(x$birth_q <= x$first_q),
  all(x$first_q >= 1L & x$first_q <= x$T),
  all(x$y %in% c(-1L, 0L, 1L)),
  all(x$dead %in% 0:1),
  all(x$infection_time >= 0)
)
cat("Basic ID, genetics, calendar and 500-draw alignment: PASS\n")
stopifnot(all(x$y[cbind(seq_len(x$N), x$first_q)] == 1L))
cat("First capture observed for every included badger: PASS\n")
stopifnot(all(rowSums(x$dead) <= 1L))
cat("At most one accepted terminal death per badger: PASS\n")
for (i in seq_len(x$N)) {
  if (x$first_q[i] > 1L &&
      any(x$y[i, seq_len(x$first_q[i] - 1L)] != -1L))
    fail("Pre-entry detection history incorrectly coded for ", x$id[i])
  dq <- x$death_q[i]
  if (!is.na(dq)) {
    if (x$dead[i, dq] != 1L ||
        (dq < x$T && any(x$y[i, (dq + 1L):x$T] == 1L)))
      fail("Post-death detection or missing event in ", x$id[i])
  }
}
cat("Entry conditioning and terminal-death chronology: PASS\n")

# An absence is informative only on a quarter with study-wide
# trapping evidence; this proxy is not individual-specific effort.
no_effort <- which(!x$calendar$observed_trapping_proxy)
if (length(no_effort) && any(x$y[, no_effort, drop = FALSE] == 0L))
  fail("False noncaptures on quarters with no trapping evidence")
cat("No false noncaptures in no-evidence quarters: PASS\n")
cat("WARNING: quarter-level global trapping activity is an effort PROXY.\n")
cat("It does not establish individual/group-specific availability.\n")

# Test exact quarterly hazard integration with multiplicative,
# age-dependent genetic and time-varying infection effects.
# GM baseline h(age) = a*exp(b*age)+c.
# multiplier = exp(gamma0 + gammaAge*age), constant in a quarter.
# Integral is analytic with a stable exprel() at zero slope.
exprel <- function(z) {
  ifelse(abs(z) < 1e-6,
         1 + z/2 + z^2/6 + z^3/24,
         expm1(z)/z)
}
gm_interval <- function(age, dt, a, b, c, gamma0, gammaAge) {
  if (a < 0 || b < 0 || c < 0 || age < 0 || dt <= 0)
    stop("Invalid mortality parameter")
  exp(gamma0 + gammaAge*age) * dt *
    (a*exp(b*age)*exprel((b+gammaAge)*dt) +
       c*exprel(gammaAge*dt))
}
test_cases <- data.frame(
  age = c(0, 2, 7, 12, 4),
  a = c(.02, .015, .03, .01, .01),
  b = c(.03, .1, .05, .08, .01),
  c = c(.002, .001, .01, .005, .001),
  gamma0 = c(0, -.5, .4, 1, 0),
  gammaAge = c(0, .06, -.04, .1, -.01)
)
for (j in seq_len(nrow(test_cases))) {
  o <- test_cases[j, ]
  numeric_truth <- integrate(
    function(t) (o$a*exp(o$b*t) + o$c) *
      exp(o$gamma0 + o$gammaAge*t),
    o$age, o$age + .25, rel.tol = 1e-12
  )$value
  analytic <- gm_interval(
    o$age, .25, o$a, o$b, o$c, o$gamma0, o$gammaAge
  )
  if (!is.finite(analytic) ||
      abs(analytic - numeric_truth) > 1e-9)
    fail("Quarterly hazard integration mismatch: test ", j)
}
cat("Piecewise GM hazards including age-dependent multipliers: PASS\n")
cat("Exact integration uses age in YEARS and quarter width = .25 years.\n")

# Audit the infection-onset interpretation without averaging draws.
# draw 0 = not infected by time horizon; positive q = onset quarter.
n_ever_draw <- sum(x$infection_time > 0L)
n_preentry <- sum(
  sweep(x$infection_time, 1L, x$first_q, FUN = "<=") &
    x$infection_time > 0L
)
cat("Positive infection-timing posterior entries:", n_ever_draw, "\n")
cat("Of these, on/before first capture:", n_preentry, "\n")
cat("Global draw columns preserved:", ncol(x$infection_time), "\n")
cat("Never-infected code remains zero.\n")

# Audit observed survival and effort support by calendar year.
year_summary <- aggregate(
  cbind(n_source_live_records, observed_trapping_proxy) ~ year,
  data = x$calendar, FUN = sum
)
names(year_summary) <- c(
  "year","live_records","quarters_with_trapping_evidence"
)
write.csv(
  year_summary,
  "outputs/model_audit/integrated_effort_year_summary.csv",
  row.names = FALSE
)
print(year_summary, row.names = FALSE)
cat("\n=== MODEL-READY COHORT SUMMARY ===\n")
cat("Age-verified, infection-matched badgers:", x$N, "\n")
cat("Confirmed alive badger-quarters:", sum(x$y == 1), "\n")
cat("Recorded deaths:", sum(x$dead == 1), "\n")
cat("Detected quarters with proxy effort:",
    sum(x$calendar$observed_trapping_proxy), "\n")
cat("Sex UNKNOWN:", sum(is.na(x$sex)), " CONFLICT:",
    sum(x$sex == "CONFLICT", na.rm = TRUE), "\n")
cat("\nNot yet validated by this script:\n")
cat("- individual/social-group trapping effort and spatial availability\n")
cat("- full-cohort birth-age uncertainty outside historical CH\n")
cat("- NIMBLE compiled custom likelihood and simulated parameter recovery\n")
cat("- infection timing-draw uncertainty integration and chain mixing\n")
cat("\nThis gate validates model inputs and hazard maths, not inference.\n")
