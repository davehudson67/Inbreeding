# =============================================================================
# Read-only investigation of NONCONVERGED single joint SEI-survival posterior.
# Contrasts mortality-baseline/frailty chains; audits accepted death provenance.
# Never changes RDS chain files or model assumptions.
# =============================================================================
options(stringsAsFactors=FALSE)
fail <- function(...)stop(paste0(...),call.=FALSE)
root <- "outputs/joint_sei_survival/full_joint_20261117"
p <- file.path(root,sprintf("chain_%02d",1:3),"posterior_samples.rds")
if(!all(file.exists(p)))fail("Missing MCMC chains: ",paste(p[!file.exists(p)],collapse=", "))
chains <- lapply(p,readRDS)
M <- lapply(chains,function(z)as.matrix(z$samples))
needed <- c("h4","f4","bGM","sigmaD","load31","load32",
            "sigmaSE","sigmaEI","betaD[2]","betaSE[4]","betaEI[4]",
            "rho","alphaP")
if(any(!needed%in%colnames(M[[1]])))
  fail("Missing required posterior monitors: ",
       paste(setdiff(needed,colnames(M[[1]])),collapse=", "))
if(!all(vapply(M,function(m)identical(colnames(m),colnames(M[[1]])),
               logical(1))))
  fail("Incompatible posterior monitor order")
cat("\n=== SINGLE JOINT MODEL: NONCONVERGENCE INVESTIGATION ===\n")
cat("Three chains, ",nrow(M[[1]]),
    " retained draws per chain, n=1603 animals.\n",sep="")
print(group_data <- data.frame(
 chain=seq_along(chains),
 seed=vapply(chains,function(z)as.integer(z$seed),integer(1)),
 niter=vapply(chains,function(z)as.integer(z$iterations),integer(1)),
 elapsed_hours=round(vapply(chains,function(z)
    as.numeric(z$elapsed_sampling_seconds)/3600,numeric(1)),2)),
 row.names=FALSE)

key <- c("h4","f4","bGM","sigmaD","load31","load32",
         "sigmaSE","sigmaEI","rho","alphaP",
         "betaSE[4]","betaEI[4]","betaD[2]")
tab <- do.call(rbind,lapply(seq_along(M),function(i){
  x <- M[[i]]
  half <- floor(nrow(x)/2)
  data.frame(
   chain=i,parameter=key,
   median=vapply(key,function(k)median(x[,k]),numeric(1)),
   first_half=vapply(key,function(k)median(x[seq_len(half),k]),numeric(1)),
   second_half=vapply(key,function(k)
       median(x[seq.int(half+1L,nrow(x)),k]),numeric(1)),
   q025=vapply(key,function(k)as.numeric(quantile(x[,k],.025)),
                numeric(1)),
   q975=vapply(key,function(k)as.numeric(quantile(x[,k],.975)),
                numeric(1))
  )
}))
cat("\n=== BASELINE/FRAILTY CHAIN MEDIANS AND DRIFT ===\n")
print(tab[tab$parameter%in%c("h4","f4","bGM","sigmaD",
                            "load31","load32"),],
      row.names=FALSE,digits=4)

pairs <- list(
 c("h4","f4"),
 c("h4","bGM"),
 c("h4","sigmaD"),
 c("f4","sigmaD"),
 c("sigmaD","load31"),
 c("sigmaD","load32"),
 c("load31","load32"),
 c("betaSE[4]","sigmaSE"),
 c("betaEI[4]","sigmaEI"))
corr <- do.call(rbind,lapply(seq_along(M),function(i)
  do.call(rbind,lapply(pairs,function(pair) {
    data.frame(chain=i,pair=paste(pair,collapse=" ~ "),
      correlation=cor(M[[i]][,pair[1]],M[[i]][,pair[2]]))
  }))))
cat("\n=== WITHIN-CHAIN POSTERIOR CORRELATIONS ===\n")
print(corr,row.names=FALSE,digits=4)

# Frailty scales have uniform priors with maxima 1.5 (biological)
# and 2.0 (capture). Near-boundary posterior mass is a warning for
# possible prior sensitivity, not direct proof of misspecification.
boundary <- do.call(rbind,lapply(seq_along(M),function(i)
  do.call(rbind,lapply(
    list(c("sigmaSE",1.5),c("sigmaEI",1.5),c("sigmaD",1.5),
         c("sigmaCapture",2.0)),
    function(k){
      x <- M[[i]][,k[1]]
      cap <- as.numeric(k[2])
      data.frame(chain=i,parameter=k[1],
        upper_prior=cap,
        frac_above_90pct_limit=mean(x>cap*.90),
        frac_above_95pct_limit=mean(x>cap*.95))
    }))))
cat("\n=== RANDOM-EFFECT SCALE PRIOR BOUNDARY AUDIT ===\n")
print(boundary,row.names=FALSE,digits=4)

# Accepted source death records were saved locally during the original
# audited read-only database extraction. Work with aggregate cause counts,
# not identity/capture IDs.
infile <- "outputs/joint_sei_survival/joint_PRIVATE_UNRECONCILED_input.rds"
deathfile <- "outputs/model_audit/integrated_accepted_death_sources.csv"
if(file.exists(infile)&&file.exists(deathfile)){
  x <- readRDS(infile)
  d <- utils::read.csv(deathfile,stringsAsFactors=FALSE)
  req <- c("id","capture_date","moc","cause_class")
  if(!all(req%in%names(d)))
    fail("Accepted death provenance CSV missing: ",
         paste(setdiff(req,names(d)),collapse=", "))
  norm <- function(z)toupper(trimws(as.character(z)))
  d$id <- norm(d$id)
  id <- norm(x$id)
  dead_idx <- which(!is.na(x$death_q) & x$death_q>0)
  if(anyDuplicated(d$id))fail("Duplicate accepted death source IDs")
  matched <- d[match(id[dead_idx],d$id),,drop=FALSE]
  ok <- !is.na(matched$id)
  cat("\n=== ACCEPTED DEATH SOURCE PROVENANCE ===\n")
  cat("Known death quarters in full joint cohort:",length(dead_idx),"\n")
  cat("Accepted death source records found:",sum(ok),"\n")
  if(any(!ok))
    cat("WARNING: Death provenance absent for",sum(!ok),
        "badgers; do not infer natural mortality.\n")
  matched <- matched[ok,,drop=FALSE]
  dead_idx <- dead_idx[ok]
  if(nrow(matched)){
    dates <- as.Date(matched$capture_date)
    q <- (as.integer(format(dates,"%Y"))-1976L)*4L+
          (as.integer(format(dates,"%m"))-1L)%/%3L+1L
    mismatch <- sum(is.na(q)|q!=x$death_q[dead_idx])
    cat("Death quarter/date provenance mismatches:",mismatch,"\n")
    cause <- norm(matched$cause_class)
    cause[is.na(cause)|cause==""] <- "UNCLASSIFIED"
    counts <- as.data.frame(table(cause),
                            stringsAsFactors=FALSE)
    names(counts) <- c("cause_class","n")
    counts$share <- round(counts$n/sum(counts$n),4)
    print(counts,row.names=FALSE)
    managed <- cause%in%c("KILLED","EUTHANISED")
    cat("Explicitly coded KILLED/EUTHANISED:",
        sum(managed),"/",nrow(matched),"\n")
    cat("NOTE: RTA and other recovered deaths are not synonymous",
        " with natural biological mortality.\n")
    cat("NOTE: No death classes have been excluded or recoded.\n")
  }
} else {
  cat("\nDeath-provenance source unavailable at expected path:\n",
      deathfile,"\n",
      "No death causes inferred; original chain data unchanged.\n",sep="")
}

out <- file.path(root,
  paste0("postfit_mixing_death_audit_",format(Sys.time(),"%Y%m%dT%H%M%S")))
if(dir.exists(out))fail("Audit path exists: ",out)
dir.create(out,recursive=TRUE)
write.csv(tab,file.path(out,"parameter_chain_drift.csv"),row.names=FALSE)
write.csv(corr,file.path(out,"posterior_parameter_correlations.csv"),
          row.names=FALSE)
write.csv(boundary,file.path(out,"frailty_prior_boundary.csv"),
          row.names=FALSE)
if(exists("counts",inherits=FALSE))
  write.csv(counts,file.path(out,"accepted_death_causes_aggregate.csv"),
            row.names=FALSE)
cat("\nSaved non-destructive diagnostic audit:",out,"\n")
cat("Interpretation remains provisional; adjust sampling within",
    " existing scientific model, not new biological variants.\n")
