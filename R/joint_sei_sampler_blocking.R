# MCMC algorithm tuning ONLY for the single joint SEI-survival model.
# Does not alter the model likelihood, biological coefficients or priors.
#
# Initial block proposal covariance is estimated from within-chain
# covariance of THREE previously completed, nonconverged chains.
# Covariances are regularized so the proposal is positive definite.
# Old posterior draws are used to tune proposals, NOT as new data.
#
# Separate blocks: (h4,f4,bGM,sigmaD) and (load31,load32).
# This directly targets the observed GM-mortality-frailty ridge without
# duplicating samplers on any target variable.
configure_joint_sei_blocked <- function(conf, input_md5,
   root="outputs/joint_sei_survival/full_joint_20261117") {
  fail <- function(...)stop(paste0(...),call.=FALSE)
  paths <- file.path(root,sprintf("chain_%02d",1:3),
                     "posterior_samples.rds")
  if(!all(file.exists(paths)))
    fail("Blocked sampler tuning requires all three original chains")
  ch <- lapply(paths,readRDS)
  sourcekeys <- vapply(ch,function(z)
    paste(unname(as.character(z$source_md5)),collapse=":"),character(1))
  inputkeys <- vapply(ch,function(z)
    paste(unname(as.character(z$input_md5)),collapse=":"),character(1))
  modelkeys <- vapply(ch,function(z)
    paste(unname(as.character(z$model_code_md5)),collapse=":"),character(1))
  if(length(unique(sourcekeys))!=1L||
     length(unique(inputkeys))!=1L||
     length(unique(modelkeys))!=1L||
     !identical(inputkeys[[1L]],paste(
       unname(as.character(input_md5)),collapse=":")))
    fail("Original 3 chains disagree on data/model or current input")
  mats <- lapply(ch,function(z)as.matrix(z$samples))
  groups <- list(
    mortality_baseline=c("h4","f4","bGM","sigmaD"),
    mortality_frailty_loadings=c("load31","load32"))
  if(!all(vapply(mats,function(m)
    all(unlist(groups,use.names=FALSE)%in%colnames(m)),logical(1))))
    fail("MCMC covariances missing required parameters")
  output <- list(
    scheme="blocked",
    source_rds=paths,
    source_rds_md5=tools::md5sum(paths),
    blocks=groups,
    proposals=list())
  for(k in names(groups)){
    target <- groups[[k]]
    # Within-chain covariances avoid artificially inflated covariance
    # from between-chain disagreement. A modest diagonal shrink
    # regularizes near-linear trade-offs.
    covs <- lapply(mats,function(m)
      stats::cov(m[,target,drop=FALSE]))
    cov0 <- Reduce("+",covs)/length(covs)
    d <- diag(cov0)
    if(any(!is.finite(cov0))||any(!is.finite(d))||any(d<=0))
      fail("Nonfinite/degenerate empirical covariance for block ",k)
    prop <- .90*cov0+.10*diag(d,nrow(cov0))
    diag(prop) <- diag(prop)+1e-6
    if(min(eigen(prop,symmetric=TRUE,
                 only.values=TRUE)$values)<=0)
      fail("Non-positive sampler proposal covariance")
    conf$removeSamplers(target,print=FALSE)
    conf$addSampler(target=target,type="RW_block",
      control=list(
        propCov=unname(prop),
        scale=.6,
        adaptive=TRUE,
        adaptInterval=100L))
    output$proposals[[k]] <- prop
    cat("Added adaptive RW_block:",paste(target,collapse="+"),
        "(initial proposal scale 0.6; adapt interval 100).\n")
    print(round(prop,5))
  }
  output
}
