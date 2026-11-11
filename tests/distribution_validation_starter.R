# Standalone phase-1 distribution validation. Run from any working directory:
# Rscript inbreeding_validation_tests.R /path/to/SimulationStudy
# Requires: R, lamW. NIMBLE test optional via option below.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop('Usage: Rscript inbreeding_validation_tests.R /path/to/SimulationStudy')
base <- normalizePath(args[[1]], mustWork = TRUE)
stopifnot(requireNamespace('lamW', quietly = TRUE))
source(file.path(base, 'Distributions', 'Gompertz', 'Dist_Gompertz.R'))
source(file.path(base, 'Distributions', 'GompertzMakeham', 'Dist_GompzMake.R'))
source(file.path(base, 'Distributions', 'Siler', 'Dist_Siler.R'))

checks <- list(
 list(name='Gompertz', dens=function(x) dGompertz(x,.08,.2), cdf=function(x) pGompertz(x,.08,.2), quant=function(p) qGompertz(p,.08,.2)),
 list(name='Gompertz-Makeham', dens=function(x) dGompzMake(x,.1,.3,.03), cdf=function(x) pGompzMake(x,.1,.3,.03), quant=function(p) qGompzMake(p,.1,.3,.03)),
 list(name='Siler', dens=function(x) dSiler(x,.3,.09,.8,.35,.06), cdf=function(x) pSiler(x,.3,.09,.8,.35,.06), quant=function(p) qSiler(p,.3,.09,.8,.35,.06))
)
failed <- FALSE
for (spec in checks) {
  cat('\nTESTING', spec$name, '\n')
  result <- tryCatch({
    area <- integrate(spec$dens, lower=0, upper=150, rel.tol=1e-7)$value
    stopifnot(is.finite(area), abs(area-1) < 1e-5)
    for (p in c(.01, .25, .5, .9, .95, .99)) {
      q <- spec$quant(p)
      value <- spec$cdf(q)
      stopifnot(is.finite(q), q>=0, is.finite(value), abs(value-p)<2e-4)
    }
    cat('PASS: density integrates to 1; quantile-CDF round trip\n')
    TRUE
  }, error=function(e) {cat('FAIL:',conditionMessage(e),'\n');FALSE})
  if (!result) failed <- TRUE
}

# Validate the known NIMBLE Siler quantile bug separately; this is only a static regression marker.
pfile <- file.path(base,'Distributions','Siler','Dist_SilerNim.R')
code <- paste(readLines(pfile, warn=FALSE),collapse='\n')
if (grepl('while\\s*\\(qU\\s*<\\s*p\\)',code)) {
  cat('KNOWN BUG: Dist_SilerNim.R uses while(qU < p) instead of while(pU < p)\n')
  failed <- TRUE
}
if (failed) quit(status=1L)
cat('\nAll checks passed. Note: this does NOT establish correct likelihood, MCMC recovery or WAIC.\n')
