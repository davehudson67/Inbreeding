# ============================================================================
# READ-ONLY diagnosis of five-year infection-pressure sampling in the
# SAME full 1,603-badger joint SEI-survival model.
#
# Compares complete original RW chains (1:3) with complete AF_slice
# chains (6:8) WITHOUT pooling samplers or refitting the model.
#
# Computes S->E and E->I quarterly transition probabilities exactly as
# the forward likelihood uses them, conditional on:
# female (male=0), age 4 (age/4=1), cohort-mean hom (hom=0),
# zero individual frailty; these are reference risks, NOT population
# averages and NOT observed incident rates. Death happens separately.
# ============================================================================
options(stringsAsFactors=FALSE)
fail <- function(...)stop(paste0(...),call.=FALSE)
if(!requireNamespace("posterior",quietly=TRUE))
  fail("Install the posterior package for rank-normalized diagnostics")
root <- "outputs/joint_sei_survival/full_joint_20261117"
groups <- list(original_RW=1:3,AF_slice=6:8)
read_group <- function(ids,scheme){
  paths <- file.path(root,sprintf("chain_%02d",ids),
                     "posterior_samples.rds")
  if(any(!file.exists(paths)))
    fail("Missing chain files: ",paste(paths[!file.exists(paths)],collapse=", "))
  chains <- lapply(paths,readRDS)
  for(i in seq_along(chains)){
    z <- chains[[i]]
    s <- if(is.null(z$sampler_scheme)) "default" else z$sampler_scheme
    if(as.integer(z$N)!=1603L||as.integer(z$chain)!=ids[i]||
       !identical(s,scheme)||
       !identical(z$status,"UNVALIDATED_PROVISIONAL"))
      fail("Unexpected model/sampler provenance in chain ",ids[i])
  }
  hashes <- vapply(chains,function(z)
    paste(as.character(unname(z$input_md5)),collapse=":"),
    character(1))
  sci <- vapply(chains,function(z){
    v <- z$model_code_md5[c(
      "R/joint_sei_survival_hierarchy_nimble.R",
      "R/joint_sei_survival_forward_nimble.R")]
    if(anyNA(v))fail("Missing scientific model hash")
    paste(unname(as.character(v)),collapse=":")
  },character(1))
  if(length(unique(hashes))>1L || length(unique(sci))>1L)
    fail("Different input or scientific model within sampler group")
  ms <- lapply(chains,function(z)as.matrix(z$samples))
  if(length(unique(vapply(ms,nrow,integer(1))))!=1L||
     !all(vapply(ms,function(x)
       identical(colnames(x),colnames(ms[[1]])),logical(1))))
    fail("Different draws/monitors within sampler group")
  list(id=ids,chains=chains,samples=ms,
       input_hash=hashes[1],science_hash=sci[1])
}
G <- list(original_RW=read_group(groups$original_RW,"default"),
          AF_slice=read_group(groups$AF_slice,"af_slice"))
if(!identical(G[[1]]$input_hash,G[[2]]$input_hash)||
   !identical(G[[1]]$science_hash,G[[2]]$science_hash))
  fail("Scientific model or input differs between original and AF chains")
vars <- c("betaSE[1]","betaSE[2]","betaEI[1]","betaEI[2]",
          "sigmaYearSE","sigmaYearEI",
          paste0("yearSE[",2:10,"]"),
          paste0("yearEI[",2:10,"]"))
if(!all(vapply(G,function(g)all(vars%in%colnames(g$samples[[1]])),
               logical(1))))
  fail("Not all calendar/intercept parameters are monitored")

period_labels <- sprintf("%d-%d",seq(1976,2021,5),seq(1980,2025,5))
derive <- function(x){
  se <- ei <- matrix(NA_real_,nrow=nrow(x),ncol=10)
  for(j in 1:10){
    ys <- if(j==1) 0 else x[,paste0("yearSE[",j,"]")]
    yi <- if(j==1) 0 else x[,paste0("yearEI[",j,"]")]
    se[,j] <- plogis(x[,"betaSE[1]"] + x[,"betaSE[2]"] + ys)
    ei[,j] <- plogis(x[,"betaEI[1]"] + x[,"betaEI[2]"] + yi)
  }
  colnames(se) <- paste0("p_SE_",period_labels)
  colnames(ei) <- paste0("p_EI_",period_labels)
  out <- cbind(
    betaSE_intercept=x[,"betaSE[1]"],
    betaEI_intercept=x[,"betaEI[1]"],
    sigmaYearSE=x[,"sigmaYearSE"],
    sigmaYearEI=x[,"sigmaYearEI"],
    se,ei)
  out
}
diagnostics <- function(group,label){
  d <- lapply(group$samples,derive)
  nm <- colnames(d[[1]])
  a <- array(NA_real_,
    dim=c(nrow(d[[1]]),3,length(nm)),
    dimnames=list(NULL,paste0("chain",group$id),nm))
  for(j in 1:3)a[,j,] <- d[[j]]
  ps <- as.data.frame(posterior::summarise_draws(
    posterior::as_draws_array(a)))
  ii <- match(nm,as.character(ps$variable))
  if(anyNA(ii))fail("Missing posterior diagnostic labels")
  all <- do.call(rbind,d)
  out <- data.frame(
    sampler=label,
    quantity=nm,
    median=unname(apply(all,2,median)),
    lower95=unname(apply(all,2,quantile,.025)),
    upper95=unname(apply(all,2,quantile,.975)),
    rhat=ps$rhat[ii],
    ess_bulk=ps$ess_bulk[ii],
    ess_tail=ps$ess_tail[ii])
  by_chain <- do.call(rbind,lapply(seq_along(d),function(j){
    x<-d[[j]];h<-nrow(x)%/%2L
    data.frame(sampler=label,chain=group$id[j],
      quantity=colnames(x),
      median=unname(apply(x,2,median)),
      first_half=unname(apply(x[seq_len(h),,drop=FALSE],2,median)),
      second_half=unname(apply(
        x[seq.int(h+1L,nrow(x)),,drop=FALSE],2,median)))
  }))
  # Correlation of intercept with calendar effects. Strong negative
  # correlations can explain a shift-ridge with the anchored first year.
  corrs <- do.call(rbind,lapply(seq_along(group$samples),function(j){
    x <- group$samples[[j]]
    do.call(rbind,lapply(c("SE","EI"),function(type){
      b <- paste0("beta",type,"[1]")
      y <- paste0("year",type,"[",2:10,"]")
      data.frame(sampler=label,chain=group$id[j],
        process=type,period=period_labels[2:10],
        intercept_year_corr=vapply(y,function(yy)
           cor(x[,b],x[,yy]),numeric(1)))
    }))
  }))
  list(stats=out,by_chain=by_chain,corr=corrs)
}
res <- Map(diagnostics,G,names(G))
stats <- do.call(rbind,lapply(res,function(z)z$stats))
by_chain <- do.call(rbind,lapply(res,function(z)z$by_chain))
corrs <- do.call(rbind,lapply(res,function(z)z$corr))

cat("\n=== JOINT SEI MODEL: FIVE-YEAR CALENDAR CONVERGENCE ===\n")
cat("Independent RW 1-3 vs AF_slice 6-8, same exact data and model: PASS\n")
cat("Reference risk = female, age 4, mean homozygosity, frailty zero;\n")
cat("quarterly conditional S->E / E->I transition probability.\n\n")
for(process in c("SE","EI")){
 cat("\n=== ",process," REFERENCE QUARTERLY TRANSITION RISK ===\n",sep="")
 x <- stats[grepl(paste0("^p_",process,"_"),stats$quantity),]
 print(x[,c("sampler","quantity","median","rhat",
            "ess_bulk","ess_tail")],row.names=FALSE,digits=4)
}
cat("\n=== CALENDAR INTERCEPT AND RW SCALE DIAGNOSTICS ===\n")
print(stats[stats$quantity%in%c(
  "betaSE_intercept","betaEI_intercept",
  "sigmaYearSE","sigmaYearEI"),],
  row.names=FALSE,digits=4)
for(label in names(G)){
  x<-stats[stats$sampler==label,]
  for(process in c("SE","EI")){
    z <- x[grepl(paste0("^p_",process,"_"),x$quantity),]
    cat(label,process,
      "reference-risk max Rhat:",signif(max(z$rhat),4),
      " min bulk ESS:",signif(min(z$ess_bulk),4),"\n")
  }
}
cat("\n=== MEDIAN INTERCEPT-CALENDAR POSTERIOR CORRELATIONS ===\n")
agg<-aggregate(intercept_year_corr~sampler+process,
  data=corrs,FUN=function(x)median(x))
print(agg,row.names=FALSE,digits=4)

out <- file.path(root,
   paste0("calendar_risk_review_",format(Sys.time(),"%Y%m%dT%H%M%S")))
if(dir.exists(out))fail("Output folder collision")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(stats,file.path(out,"calendar_probability_convergence.csv"),
          row.names=FALSE)
write.csv(by_chain,file.path(out,"calendar_probabilities_by_chain_half.csv"),
          row.names=FALSE)
write.csv(corrs,file.path(out,"intercept_period_correlations.csv"),
          row.names=FALSE)
writeLines(c(
 "Original RW vs AF_slice calendar progression/acquisition comparison.",
 "Scientific model identical; original five-year effects are RW with first period anchored to zero.",
 "Reference risks are conditional age-4 female, mean homozygosity, frailty zero, not population-average effects.",
 "Original chain 1-3 and AF slice chain 6-8 each analysed independently.",
 "Do not infer biological/causal effects while model convergence, prior sensitivity, and death-provenance issues remain."
),file.path(out,"limitations.txt"))
cat("\nSaved read-only calendar audit:",out,"\n")
