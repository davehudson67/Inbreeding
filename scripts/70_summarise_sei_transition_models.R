# ==========================================================================
# Stage 70: combine ALL Stage 69 full shared-draw SEI transition regressions.
# Read only. Aggregate exploratory Rubin estimates across the 1000
# distinct posterior draw columns, checking coverage and code/input hashes.
# ==========================================================================
options(stringsAsFactors=FALSE)
out <- "outputs/sei_transitions"
nshards <- as.integer(Sys.getenv("N_SHARDS","3"))
if(!is.finite(nshards)||nshards<1L||nshards>1000L)
  stop("Bad N_SHARDS")
paths <- file.path(out,sprintf("full_shard_%02d_of_%02d.rds",
                              seq_len(nshards),nshards))
if(!all(file.exists(paths))) {
  cat("Missing shard files:\n")
  print(paths[!file.exists(paths)])
  stop("All full Stage 69 shards must finish first.")
}
Z <- lapply(paths,readRDS)
if(any(vapply(Z,function(a)a$mode!="full",logical(1))) ||
   !identical(vapply(Z,function(a)as.integer(a$shard),integer(1)),
              seq_len(nshards)) ||
   any(vapply(Z,function(a)as.integer(a$nshards)!=nshards,logical(1))))
  stop("Shard IDs, modes, or total shard count do not match")
first <- Z[[1]]
for(j in seq_along(Z)) {
  other <- Z[[j]]
  if(!identical(other$model,first$model)||
     !identical(other$N,first$N)||
     !identical(other$model_input_md5,first$model_input_md5)||
     !identical(other$sei_input_md5,first$sei_input_md5)||
     !identical(other$fit_source_md5,first$fit_source_md5)||
     !identical(other$ids,first$ids)||
     !identical(other$assumptions,first$assumptions))
    stop("Different cohort, model definition, assumption set or input MD5 at shard ",j)
}
d <- do.call(rbind,lapply(Z,function(z)z$results))
all_draws <- sort(unique(d$draw))
if(!identical(as.integer(all_draws),1:1000))
  stop("Full-run draw coverage is incomplete or unexpected")
if(anyDuplicated(d[,c("draw","transition")]) ||
   nrow(d)!=2000L ||
   !all(table(d$transition)==1000L) ||
   !all(d$transition %in% c("SE","EI")))
  stop("Missing/duplicated transition-draw combinations")
d <- d[order(d$transition,d$draw),,drop=FALSE]
params <- c("hom_z","male","entry_decade")
pool <- function(a,par) {
  beta <- a[[paste0("beta_",par)]]
  se <- a[[paste0("se_",par)]]
  if(length(beta)!=1000L||any(!is.finite(beta))||
     any(!is.finite(se))||any(se<=0))
    stop("Invalid draw coefficients or SE for ",par)
  m <- length(beta)
  u <- mean(se^2)
  b <- var(beta)
  total <- u+(1+1/m)*b
  df <- if(b<=0) Inf else (m-1L)*(1+u/((1+1/m)*b))^2
  critical <- if(is.finite(df))qt(.975,df) else qnorm(.975)
  v <- sqrt(total)
  data.frame(
    transition=a$transition[1],parameter=par,
    logHR=mean(beta),totalSE=v,
    HR=exp(mean(beta)),HR_low95=exp(mean(beta)-critical*v),
    HR_high95=exp(mean(beta)+critical*v),
    mean_conditional_variance=u,
    between_draw_variance=b,
    fraction_between=(1+1/m)*b/total,
    pooled_df=if(is.finite(df))df else NA_real_,
    mean_events=mean(a$events),
    median_events=median(a$events),
    mean_risk_set=mean(a$n_risk),
    row.names=NULL)
}
tab <- do.call(rbind,lapply(c("SE","EI"),function(t){
  a <- d[d$transition==t,,drop=FALSE]
  do.call(rbind,lapply(params,function(p)pool(a,p)))
}))
write.csv(tab,file.path(out,"full_pooled_transition_effects.csv"),
          row.names=FALSE)
write.csv(d,file.path(out,"full_1000_draws_transition_coefficients.csv"),
          row.names=FALSE)
cat("\n=== STAGE 70: ALL 1000 SHARED SEI DRAWS ===\n")
cat("Input MD5:",first$model_input_md5,"\n")
cat("SEI aligned MD5:",first$sei_input_md5,"\n")
cat("Shards:",nshards,"; 1000 unique joint draws; 2,000 models\n")
cat("Age-time-scale Cox regressions; continuous hom (+1 SD),",
    "male, entry calendar decade.\n\n")
print(tab[,c("transition","parameter","HR","HR_low95",
             "HR_high95","fraction_between","median_events")],
      row.names=FALSE,digits=4)
cat("\nEvent counts and risk-set ranges across imputed draws:\n")
for(t in c("SE","EI")) {
  a <- d[d$transition==t,,drop=FALSE]
  cat(t,": Events median ",median(a$events),
      " range ",paste(range(a$events),collapse="-"),
      "; risk set median ",median(a$n_risk),
      "; risk set range ",paste(range(a$n_risk),collapse="-"),"\n",sep="")
}
# Empirical imputation stability: effects estimated using first 250/500
# posterior columns vs all 1000, not an MCMC diagnostic for SEI fit.
for(t in c("SE","EI")){
  a <- d[d$transition==t,,drop=FALSE]
  cat("\n",t,"homozygosity log-HR mean by number of SEI draws:\n",sep="")
  for(k in c(100L,250L,500L,1000L))
    cat("  ",k," draws: ",
        round(mean(a$beta_hom_z[a$draw<=k]),5),"\n",sep="")
}
cat("\nCRITICAL INTERPRETATION:\n")
cat("This is a modular exploratory analysis of an EXISTING SEI",
    "posterior, not a fully joint disease-history regression.\n")
cat("No diagnostic likelihood refit, no adjustment for social-group",
    "clustering or informative capture/follow-up.\n")
cat("Calendar is entry-decade covariate, age a flexible baseline.",
    "Censoring occurs before known death, after last capture.\n")
cat("Reported intervals are APPROXIMATE Rubin pooled intervals;",
    "do not treat them as definitive causal or calibrated uncertainty.\n")
saveRDS(list(results=d,pooled=tab,source_shards=paths,
             model=first$model,assumptions=first$assumptions,
             cohort_N=first$N,draws=1000L,
             sei_input_md5=first$sei_input_md5,
             model_input_md5=first$model_input_md5),
        file.path(out,"full_1000_draws_combined.rds"))
cat("\nSTAGE 70 COMPLETE: outputs in ",out,"\n",sep="")
