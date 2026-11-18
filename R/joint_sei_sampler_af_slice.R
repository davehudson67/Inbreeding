# ============================================================================
# NIMBLE automated factor slice sampler (AF_slice) for the SAME joint
# infection-survival Bayesian scientific model and prior.
#
# Replace default RW updates for SIX highly correlated mortality parameters
# with ONE joint AF_slice sampler: h4/f4/bGM/sigmaD/load31/load32.
# It learns its factor directions from the chain. No biological model changes.
#
# Initial per-dimension slice widths come from within-chain standard
# deviations in the original three full 1603-badger MCMC runs. These draws
# tune proposals only; they are not data for the biological likelihood.
#
# Unlike the earlier separate RW_block groups this single six-node block can
# follow dependencies between GM mortality parameters and frailty loadings.
# Keep adaptation confined to burn-in by choosing the max iteration
# from NBURN (default 3000) when compiling the sampler.
# ============================================================================
configure_joint_sei_AF_slice <- function(
    conf, input_md5,
    root="outputs/joint_sei_survival/full_joint_20261117") {
  fail_af <- function(...)stop(paste0(...),call.=FALSE)
  paths <- file.path(root,sprintf("chain_%02d",1:3),
                     "posterior_samples.rds")
  if(!all(file.exists(paths)))
    fail_af("AF_slice tuning needs the original 3 chain files: ",
            paste(paths[!file.exists(paths)],collapse=", "))
  chains <- lapply(paths,readRDS)
  full_md5 <- function(z)
    paste(unname(as.character(z)),collapse=":")
  source_keys <- vapply(chains,function(z)
    full_md5(z$source_md5),character(1))
  input_keys <- vapply(chains,function(z)
    full_md5(z$input_md5),character(1))
  science_keys <- vapply(chains,function(z){
    h <- z$model_code_md5
    need <- c("R/joint_sei_survival_hierarchy_nimble.R",
              "R/joint_sei_survival_forward_nimble.R")
    if(!all(need%in%names(h)))
      fail_af("Old chain has no exact scientific model hashes")
    full_md5(h[need])
  },character(1))
  current_input <- full_md5(input_md5)
  if(length(unique(source_keys))!=1L ||
     length(unique(input_keys))!=1L ||
     length(unique(science_keys))!=1L ||
     !identical(input_keys[[1L]],current_input))
    fail_af("AF_slice source data/model provenance differs from original chains")
  targets <- c("h4","f4","bGM","sigmaD","load31","load32")
  oldM <- lapply(chains,function(z)as.matrix(z$samples))
  if(!all(vapply(oldM,function(z)
      all(targets%in%colnames(z)),logical(1))))
    fail_af("Original chains do not monitor every AF_slice target")
  # Robust initial proposal widths; AF_slice will adapt rotations/widths.
  scales <- vapply(oldM,function(x)
    vapply(targets,function(nm)stats::sd(x[,nm]),numeric(1)),
    numeric(length(targets)))
  if(any(!is.finite(scales)) || any(scales<=0))
    fail_af("Invalid posterior scales for AF_slice initial slice widths")
  widths <- pmax(0.015,pmin(0.5,.7*rowMeans(scales)))
  names(widths) <- targets

  burn <- suppressWarnings(as.integer(Sys.getenv("NBURN","3000")))
  if(is.na(burn) || burn<200L)
    fail_af("AF_slice needs NBURN >= 200 for its factor adaptation")
  # A pilot can have a shorter burn-in. Adapt factor directions only
  # before the saved draw interval, never throughout the posterior sample.
  max_factor_iterations <- as.integer(min(15000L,burn))
  controls <- list(
    sliceWidths=unname(as.numeric(widths)),
    sliceAdaptFactorMaxIter=max_factor_iterations,
    sliceAdaptFactorInterval=200L,
    sliceAdaptWidthMaxIter=512L,
    sliceAdaptWidthTolerance=0.10,
    sliceMaxSteps=100L,
    maxContractions=100L,
    maxContractionsWarning=TRUE
  )
  conf$removeSamplers(targets,print=FALSE)
  conf$addSampler(target=targets,type="AF_slice",
                  control=controls,print=TRUE)

  cat("Added single NIMBLE AF_slice sampler across SIX mortality nodes:\n")
  cat("  ",paste(targets,collapse=" + "),"\n")
  cat("Initial widths:",
      paste(sprintf("%s=%.4g",targets,widths),collapse="; "),"\n")
  cat("Factor direction adaptation ends by iteration",
      max_factor_iterations,"(burn-in =",burn,").\n")
  cat("These are proposal controls ONLY; model and priors unchanged.\n")

  list(scheme="af_slice",targets=targets,
       widths=widths,controls=controls,
       source_rds=paths,
       source_rds_md5=tools::md5sum(paths),
       original_scientific_model_md5=science_keys[[1L]])
}
