# ============================================================================
# NEW MCMC SCHEDULING OPTION: TWO AF_slice calendar blocks, preserving the
# IDENTICAL joint infection-survival likelihood, parameterization and priors.
#
# Mortality AF_slice from R/joint_sei_sampler_af_slice.R is still used.
# Each period block updates the infection intercept jointly with the nine
# stochastic 5-year period levels (anchored period 1 remains deterministic).
#
# S->E: betaSE[1] + yearSE[2:10]
# E->I: betaEI[1] + yearEI[2:10]
#
# Old posterior samples initialize sampler widths ONLY, never likelihood data.
# The new design is an EXPLORATORY sampler, not a second biological model.
# Requires a separate compilation/timing benchmark before any new long chain.
# ============================================================================
configure_joint_sei_calendar_AF_slice <- function(conf,input_md5,
   root="outputs/joint_sei_survival/full_joint_20261117"){
  fail_cal <- function(...)stop(paste0(...),call.=FALSE)
  ids <- 6:8
  files <- file.path(root,sprintf("chain_%02d",ids),
                     "posterior_samples.rds")
  if(!all(file.exists(files)))
    fail_cal("Missing AF_slice convergence source chains 6-8")
  sources <- lapply(files,readRDS)
  canonical_md5 <- function(z)
    paste(unname(as.character(z)),collapse=":")
  model_files <- c("R/joint_sei_survival_hierarchy_nimble.R",
                   "R/joint_sei_survival_forward_nimble.R")
  model_hashes_now <- tools::md5sum(model_files)
  ref_input <- canonical_md5(input_md5)
  for(j in seq_along(sources)){
    z <- sources[[j]]
    if(as.integer(z$chain)!=ids[j] || as.integer(z$N)!=1603L ||
       !identical(z$sampler_scheme,"af_slice") ||
       !identical(z$status,"UNVALIDATED_PROVISIONAL") ||
       canonical_md5(z$input_md5)!=ref_input)
      fail_cal("AF reference chain ",ids[j]," input/status mismatch")
    h <- z$model_code_md5
    if(!all(model_files%in%names(h)) ||
       !identical(canonical_md5(h[model_files]),
                  canonical_md5(model_hashes_now)))
      fail_cal("AF reference chain ",ids[j],
               " scientific model code differs from current version")
  }
  mats <- lapply(sources,function(z)as.matrix(z$samples))
  targets <- list(
    acquisition=c("betaSE[1]",paste0("yearSE[",2:10,"]")),
    progression=c("betaEI[1]",paste0("yearEI[",2:10,"]")))
  if(!all(vapply(mats,function(x)
    all(unlist(targets,use.names=FALSE)%in%colnames(x)),logical(1))))
    fail_cal("AF source chains missing calendar/intercept coefficients")
  burn <- suppressWarnings(as.integer(Sys.getenv("NBURN","3000")))
  if(is.na(burn) || burn<1000L)
    fail_cal("NBURN must be >= 1000 for calendar AF_slice adaptation")
  stop_at <- as.integer(min(15000L,burn-600L))
  metadata <- list(
    scheme="af_slice_calendar",
    source_ids=ids,
    source_rds=files,
    source_rds_md5=tools::md5sum(files),
    science_md5=model_hashes_now,
    groups=targets,
    proposals=list())
  for(k in names(targets)){
    target <- targets[[k]]
    scales <- vapply(mats,function(x)
      vapply(target,function(name)sd(x[,name]),numeric(1)),
      numeric(length(target)))
    if(any(!is.finite(scales))||any(scales<=0))
      fail_cal("Invalid calendar posterior scale for ",k)
    # Robust widths for high-dimensional, correlated factor slice.
    # Width adaptation occurs in burn-in, with a bounded initial range.
    widths <- pmax(.03,pmin(.80,.7*rowMeans(scales)))
    names(widths) <- target
    ctl <- list(
      sliceWidths=unname(as.numeric(widths)),
      sliceAdaptFactorMaxIter=stop_at,
      sliceAdaptFactorInterval=200L,
      sliceAdaptWidthMaxIter=512L,
      sliceAdaptWidthTolerance=.10,
      sliceMaxSteps=100L,
      maxContractions=100L,
      maxContractionsWarning=TRUE)
    conf$removeSamplers(target,print=FALSE)
    conf$addSampler(target=target,type="AF_slice",
                    control=ctl,print=TRUE)
    cat("Calendar AF_slice (",k,") target nodes: ",
        paste(target,collapse=", "),"\n",sep="")
    cat("Calendar initial widths: ",
        paste(sprintf("%s=%.4g",target,widths),collapse="; "),"\n",sep="")
    cat("Calendar factor adaptations end by iteration ",stop_at,
        " (NBURN=",burn,").\n",sep="")
    metadata$proposals[[k]] <- list(widths=widths,controls=ctl)
  }
  cat("Calendar AF_slice blocks change only MCMC proposals;",
      " biological model code and priors are IDENTICAL.\n")
  metadata
}
