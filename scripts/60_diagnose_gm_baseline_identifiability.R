# =============================================================
# Stage 60: Diagnose baseline GM identifiability across the
# completed three-chain real-data fits.
#
# Computes age-specific hazards and interval survival per chain,
# checks R-hat/ESS of DERIVED quantities, and compares traces.
# No model refit, database access, or data modifications.
# =============================================================
options(stringsAsFactors=FALSE)
stopifnot(requireNamespace("coda",quietly=TRUE))
out <- "outputs/integrated_gm"
files <- file.path(out,sprintf("full_chain_%02d.rds",1:3))
stopifnot(all(file.exists(files)))
z <- lapply(files,readRDS)
if(!all(vapply(z,function(x)x$n==1603L,logical(1))))
  stop("Unexpected cohort sizes")
hashes <- vapply(z,function(x)as.character(x$input_md5),character(1))
if(length(unique(hashes))!=1L)stop("Data differ between chains")
codehash <- vapply(z,function(x)paste(x$code_md5,collapse=";"),character(1))
if(length(unique(codehash))!=1L)stop("Code differs between chains")
m <- lapply(z,function(x)as.matrix(x$samples))
nms <- colnames(m[[1]])
if(!all(vapply(m,function(x)identical(colnames(x),nms),logical(1))))
  stop("Inconsistent monitor columns")
required <- c("a","b","c1","bh","bha","bi","bia","bhi","bs",
              "rho","alphaP","sigmaU")
if(!all(required %in% nms))
  stop("Missing required monitors: ",paste(setdiff(required,nms),collapse=", "))
if(!all(vapply(m,function(x)all(is.finite(x[,required])),logical(1))))
  stop("Nonfinite core parameter samples")
if(!length(unique(vapply(m,nrow,integer(1))))==1L)
  stop("Chains have different posterior lengths")

# Cumulative hazard and survival between exact ages (years).
# expm1 form ensures numerical precision near b = 0.
integral <- function(A,B,C,age0,age1){
  delta <- age1-age0
  term <- ifelse(abs(B)<1e-7,
                 A*exp(B*age0)*delta*(1+B*delta/2),
                 A*exp(B*age0)*expm1(B*delta)/B)
  term+C*delta
}
ages <- c(0,2,4,6,8,10,12,14)
derived_for <- function(x){
  A <- x[,"a"];B <- x[,"b"];C <- x[,"c1"]
  ans <- cbind(a=A,b=B,c1=C,
               baseline_hazard_at_age_2=A*exp(2*B)+C,
               baseline_hazard_at_age_4=A*exp(4*B)+C,
               baseline_hazard_at_age_6=A*exp(6*B)+C,
               baseline_hazard_at_age_8=A*exp(8*B)+C,
               baseline_hazard_at_age_10=A*exp(10*B)+C,
               baseline_hazard_at_age_12=A*exp(12*B)+C,
               baseline_S_age_2_to_6=exp(-integral(A,B,C,2,6)),
               baseline_S_age_2_to_8=exp(-integral(A,B,C,2,8)),
               relative_hazard_age_8_over_2=
                 (A*exp(8*B)+C)/(A*exp(2*B)+C),
               infection_HR_age_4=exp(x[,"bi"]),
               infection_HR_age_8=exp(x[,"bi"]+4*x[,"bia"]),
               homozygosity_plus_1SD_HR_age_4=exp(x[,"bh"]),
               homozygosity_plus_1SD_HR_age_8=
                 exp(x[,"bh"]+4*x[,"bha"]),
               male_HR=exp(x[,"bs"]))
  ans
}
d <- lapply(m,derived_for)
nm <- colnames(d[[1]])
ml <- coda::mcmc.list(lapply(d,coda::mcmc))
gr <- coda::gelman.diag(ml,autoburnin=FALSE,multivariate=FALSE)$psrf[,1]
ess <- coda::effectiveSize(ml)
pool <- do.call(rbind,d)
qs <- apply(pool,2,quantile,probs=c(.025,.5,.975))
summary <- data.frame(
  quantity=nm,
  mean=colMeans(pool),
  median=unname(qs[2,]),
  lower95=unname(qs[1,]),
  upper95=unname(qs[3,]),
  rhat=as.numeric(gr[nm]),
  ESS=as.numeric(ess[nm]),
  stringsAsFactors=FALSE
)
bychain <- do.call(rbind,lapply(seq_along(d),function(j) {
  x <- d[[j]]
  data.frame(chain=j,quantity=colnames(x),
             median=apply(x,2,median),
             mean=colMeans(x),
             low=apply(x,2,quantile,.025),
             high=apply(x,2,quantile,.975),
             row.names=NULL)
}))
cat("\n=== BASELINE GM IDENTIFIABILITY REVIEW ===\n")
cat("Chains:",length(d),"; draws per chain:",nrow(d[[1]]),"\n")
cat("Input MD5:",unique(hashes),"\n")
cat("\nPosterior pooled summaries of baseline and derived survival:\n")
print(summary,row.names=FALSE,digits=4)
cat("\nChain-specific reference hazards at age 4, age 8, and survival 2 to 8:\n")
important <- c("a","b","c1","baseline_hazard_at_age_4",
               "baseline_hazard_at_age_8","baseline_S_age_2_to_8")
print(bychain[bychain$quantity %in% important,
              c("chain","quantity","median","low","high")],
      row.names=FALSE,digits=4)

# Report early-half/late-half drift to detect lack of stationarity.
drift <- do.call(rbind,lapply(seq_along(d),function(j){
  x <- d[[j]]
  n <- nrow(x)
  k <- floor(n/2)
  if(k<40) return(NULL)
  measure <- c("a","b","c1","baseline_hazard_at_age_4",
               "baseline_hazard_at_age_8","baseline_S_age_2_to_8")
  data.frame(chain=j,quantity=measure,
             mean_first=colMeans(x[seq_len(k),measure,drop=FALSE]),
             mean_second=colMeans(x[(k+1):n,measure,drop=FALSE]),
             stringsAsFactors=FALSE)
}))
drift$change <- drift$mean_second-drift$mean_first
cat("\nFirst-half vs second-half posterior means:\n")
print(drift,row.names=FALSE,digits=4)
write.csv(summary,file.path(out,"derived_age_hazards_convergence.csv"),
          row.names=FALSE)
write.csv(bychain,file.path(out,"derived_age_hazards_by_chain.csv"),
          row.names=FALSE)
write.csv(drift,file.path(out,"derived_age_hazards_drift.csv"),
          row.names=FALSE)

# Does the available data include sufficient older living badgers?
input_path <- z[[1]]$input
if(!is.null(input_path) && file.exists(input_path)){
  inp <- readRDS(input_path)
  if(all(c("birth_q","y","id") %in% names(inp))) {
    indices <- which(inp$y==1L,arr.ind=TRUE)
    age <- (indices[,2]-inp$birth_q[indices[,1]])*.25
    sample_age <- data.frame(
      threshold_years=c(2,4,6,8,10,12,14),
      observed_badger_quarters=vapply(c(2,4,6,8,10,12,14),
        function(a)sum(age>=a),integer(1)),
      unique_detected_badgers=vapply(c(2,4,6,8,10,12,14),
        function(a)length(unique(indices[age>=a,1])),integer(1))
    )
    cat("\nAge coverage among confirmed-live encounters:\n")
    print(sample_age,row.names=FALSE)
    write.csv(sample_age,
              file.path(out,"empirical_live_age_coverage.csv"),
              row.names=FALSE)
    cat("Median live-observed age:",median(age),
        "; maximum:",max(age),"\n")
  }
}

# One concise figure: chain-specific median baseline hazard vs age.
png(file.path(out,"baseline_mortality_by_chain.png"),
    width=1150,height=700,res=135)
cols <- c("#1F77B4","#D62728","#2CA02C")
ylim <- range(unlist(lapply(d,function(q) {
  vapply(ages,function(age)
    median(q[,"a"]*exp(age*q[,"b"])+q[,"c1"]),numeric(1))
})))
plot(ages,ages*0,type="n",ylim=ylim,xlab="Age (years)",
     ylab="Baseline annual mortality hazard",
     main="GM baseline mortality: chain-specific posterior medians")
for(j in seq_along(d)){
  v <- vapply(ages,function(age)
    median(d[[j]][,"a"]*exp(age*d[[j]][,"b"])+d[[j]][,"c1"]),
    numeric(1))
  lines(ages,v,col=cols[j],lwd=2.5,type="b",pch=15+j)
}
legend("topleft",legend=paste("Chain",1:3),
       col=cols,lwd=2.5,pch=16:18,bty="n")
dev.off()

png(file.path(out,"baseline_parameter_traces.png"),
    width=1200,height=900,res=130)
op <- par(mfrow=c(3,1),mar=c(3.6,4.5,2.2,1))
for(p in c("a","b","c1")) {
  yy <- unlist(lapply(m,function(v)v[,p]))
  plot(seq_len(nrow(m[[1]])),m[[1]][,p],type="l",
       ylim=range(yy),col=cols[1],xlab="Retained sample",
       ylab=p,main=p)
  for(j in 2:3)lines(seq_len(nrow(m[[j]])),m[[j]][,p],col=cols[j])
}
par(op)
dev.off()
cat("\nSaved diagnostics CSVs and two PNG plots to",out,"\n")
cat("INTERPRETATION GATE:\n")
if(any(summary$rhat[summary$quantity %in% c(
  "baseline_hazard_at_age_4","baseline_hazard_at_age_8",
  "baseline_S_age_2_to_8")]>1.05)) {
  cat("Derived age-specific mortality/survival ALSO NOT CONVERGED.\n")
  cat("Do not interpret baseline survival or senescence yet.\n")
} else {
  cat("Derived mortality/survival diagnostics look more stable than",
      "raw a,b,c1, but parameters are NOT globally converged.\n")
  cat("Do not treat these pooled predictions as final inference.\n")
}
