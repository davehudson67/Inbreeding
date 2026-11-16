# ==============================================================
# Stage 63: Recapture/detection inference from ANCHORED GM chains.
#
# Tests whether slowly mixing alphaP and calendar-period effects
# nevertheless yield stable predicted recapture probabilities.
#
# All quantities are conditional on the unvalidated SG monitoring
# proxy; not a validation of observation effort or emigration.
# No model refitting or source-data changes.
# ==============================================================
options(stringsAsFactors=FALSE)
stopifnot(requireNamespace("coda",quietly=TRUE))
out <- "outputs/integrated_gm/reparam"
ff <- file.path(out,sprintf("full_chain_%02d.rds",1:3))
stopifnot(all(file.exists(ff)))
x <- lapply(ff,readRDS)
data_hash <- vapply(x,function(z)as.character(z$input_md5),character(1))
code_hash <- vapply(x,function(z)paste(z$code_md5,collapse=":"),character(1))
if(length(unique(data_hash))!=1L || length(unique(code_hash))!=1L)
  stop("Chains have different inputs or code")
M <- lapply(x,function(z)as.matrix(z$samples))
stopifnot(all(vapply(M,function(z)identical(colnames(z),colnames(M[[1]])),
                     logical(1))))
need <- c("alphaP","sigmaU",sprintf("periodEff[%d]",1:10),
          sprintf("seasonEff[%d]",1:4))
missing <- setdiff(need,colnames(M[[1]]))
if(length(missing))stop("Missing monitors: ",paste(missing,collapse=", "))
stopifnot(all(vapply(M,function(z)all(is.finite(z[,need])),logical(1))))
cat("\n=== ANCHORED GM RECAPTURE DIAGNOSTICS ===\n")
cat("Individuals:",x[[1]]$n,"; retained draws:",
    paste(vapply(M,nrow,integer(1)),collapse=", "),"\n")
cat("Input MD5:",data_hash[1],"\n")

# A quadrature-like, deterministic equally weighted normal-quantile
# approximation for the marginal recapture probability in a
# randomly sampled badger: E_u[logit^-1(alpha + ui + period + season)].
# This integration is not a posterior estimate of any individual ui.
normal_quantiles <- qnorm(((1:41)-.5)/41)
marginal_p <- function(eta,sigma){
  p <- numeric(length(eta))
  for(zk in normal_quantiles)
    p <- p + plogis(eta + sigma*zk)
  p/length(normal_quantiles)
}
predict <- function(z){
  alpha <- z[,"alphaP"]; sig <- z[,"sigmaU"]
  out <- matrix(nrow=nrow(z),ncol=0)
  for(j in 1:10){
    pj <- z[,sprintf("periodEff[%d]",j)]
    for(k in 1:4){
      sk <- z[,sprintf("seasonEff[%d]",k)]
      eta <- alpha+pj+sk
      a <- marginal_p(eta,sig)
      out <- cbind(out,a)
      colnames(out)[ncol(out)] <- sprintf("p_period%02d_season%d",j,k)
    }
  }
  # The difference between marginal detection in the first
  # and last periods with season 1, a useful identifiable contrast.
  out <- cbind(out,
    p_period10_minus_period01=
      out[,"p_period10_season1"]-out[,"p_period01_season1"],
    p_period05_minus_period01=
      out[,"p_period05_season1"]-out[,"p_period01_season1"])
  out
}
P <- lapply(M,predict)
summarise <- function(a) {
  ml <- coda::mcmc.list(lapply(a,coda::mcmc))
  psrf <- coda::gelman.diag(ml,autoburnin=FALSE,
                            multivariate=FALSE)$psrf[,1]
  ess <- coda::effectiveSize(ml)
  all <- do.call(rbind,a)
  qs <- apply(all,2,quantile,c(.025,.5,.975))
  nn <- colnames(all)
  data.frame(quantity=nn,mean=colMeans(all),
             median=unname(qs[2,]),
             lower95=unname(qs[1,]),
             upper95=unname(qs[3,]),
             rhat=unname(psrf[nn]),
             ESS=unname(ess[nn]),row.names=NULL)
}
post <- summarise(P)
write.csv(post,file.path(out,"recapture_predicted_probabilities.csv"),
          row.names=FALSE)
cat("\nPredicted recapture probability, season 1 by five-year period:\n")
keep <- sprintf("p_period%02d_season1",1:10)
print(post[match(keep,post$quantity),
           c("quantity","median","lower95","upper95","rhat","ESS")],
      row.names=FALSE,digits=4)
cat("\nRecapture probability contrasts:\n")
print(tail(post,2),row.names=FALSE,digits=4)

# Distinguish weakly-mixing intercept/period coefficients from
# the marginal predictions they jointly produce.
raws <- lapply(M,function(z){
  z[,c("alphaP","sigmaU",sprintf("periodEff[%d]",2:10),
       sprintf("seasonEff[%d]",2:4)),drop=FALSE]
})
rawsummary <- summarise(raws)
write.csv(rawsummary,file.path(out,"recapture_coefficients_convergence.csv"),
          row.names=FALSE)
cat("\nRecapture coefficient diagnostics:\n")
print(rawsummary[,c("quantity","rhat","ESS")],
      row.names=FALSE,digits=4)
cat("\nMinimum ESS of 40 derived period-season probabilities:",
    min(post$ESS[1:40]),"\n")
cat("Maximum Rhat of 40 derived period-season probabilities:",
    max(post$rhat[1:40]),"\n")
if(any(post$rhat[1:40]>1.05) || any(post$ESS[1:40]<400)){
  cat("RECAPTURE GATE: INSUFFICIENT MIXING in some predicted",
      "recapture probabilities.\n")
} else {
  cat("RECAPTURE GATE: period-season predicted probabilities",
      "have adequate initial mixing diagnostics.\n")
}
cat("Regardless of MCMC mixing, these predictions rely on",
    "social-group-specific capture-activity proxy, not measured effort.\n")
