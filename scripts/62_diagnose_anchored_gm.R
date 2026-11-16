# ==============================================================
# Stage 62: Anchored GM three-chain posterior and identifiability
# evaluation. Read-only; does NOT rerun or change model.
# Requires successfully completed Stage 61 full_chain_01..03.rds.
# ==============================================================
options(stringsAsFactors = FALSE)
stopifnot(requireNamespace("coda", quietly = TRUE))
out <- "outputs/integrated_gm/reparam"
files <- file.path(out, sprintf("full_chain_%02d.rds",1:3))
if(!all(file.exists(files))) {
  stop("Missing anchored chain output(s): ",
       paste(basename(files[!file.exists(files)]), collapse=", "))
}

chains <- lapply(files,readRDS)
cat("\n=== STAGE 62: ANCHORED GM THREE-CHAIN DIAGNOSTICS ===\n")
n <- vapply(chains,function(x)as.integer(x$n),integer(1))
datahash <- vapply(chains,function(x)as.character(x$input_md5),character(1))
codehash <- vapply(chains,function(x)paste(x$code_md5,collapse=";"),character(1))
if(length(unique(n)) != 1L || length(unique(datahash)) != 1L ||
   length(unique(codehash)) != 1L)
  stop("Three anchored chains do not share exactly the same input and model source")
cat("Badgers:",n[1],"\n")
cat("Input MD5:",datahash[1],"\n")
exitpaths <- file.path(out,sprintf("full_chain_%02d.exit",1:3))
if(all(file.exists(exitpaths))) {
  exitcodes <- vapply(exitpaths,function(p)trimws(readLines(p,warn=FALSE)[1]),character(1))
  cat("Process exit codes:",paste(exitcodes,collapse=", "),"\n")
  if(any(exitcodes != "0")) stop("One or more anchored chains exited with errors")
} else {
  cat("NOTE: Not all process exit files exist; verifying saved samples instead.\n")
}
M <- lapply(chains,function(x)as.matrix(x$samples))
if(!all(vapply(M,function(m)identical(colnames(m),colnames(M[[1]])),logical(1))))
  stop("Monitors differ between chains")
cat("Retained draws per chain:",paste(vapply(M,nrow,integer(1)),collapse=", "),"\n")
if(any(vapply(M,nrow,integer(1)) < 100L))
  stop("Insufficient retained posterior draws")
if(any(vapply(M,function(m)any(!is.finite(m)),logical(1))))
  stop("Nonfinite monitored posterior samples")
required <- c("h4","f4","b","a","c1","h8",
              "bh","bha","bi","bia","bhi","bs","rho","sigmaU","alphaP")
if(!all(required %in% colnames(M)))
  stop("Missing expected monitored parameters: ",
       paste(setdiff(required,colnames(M)),collapse=", "))

# Remove deterministic reference-category columns before Rhat.
vary <- vapply(seq_along(colnames(M[[1]])),function(j) {
  v <- unlist(lapply(M,function(x)x[,j]))
  sd(v) > 1e-12
},logical(1))
cat("Constant monitors excluded:",paste(colnames(M[[1]])[!vary],collapse=", "),"\n")
D <- lapply(M,function(x)x[,vary,drop=FALSE])
summarise_samples <- function(X) {
  ml <- coda::mcmc.list(lapply(X,coda::mcmc))
  rhat <- coda::gelman.diag(ml,autoburnin=FALSE,multivariate=FALSE)$psrf[,1]
  ess <- coda::effectiveSize(ml)
  all <- do.call(rbind,X)
  q <- apply(all,2,quantile,probs=c(.025,.5,.975),na.rm=FALSE)
  nm <- colnames(all)
  data.frame(parameter=nm,mean=colMeans(all),median=unname(q[2,]),
             lower95=unname(q[1,]),upper95=unname(q[3,]),
             rhat=as.numeric(rhat[nm]),effective_n=as.numeric(ess[nm]),
             row.names=NULL)
}
parameter_summary <- summarise_samples(D)
write.csv(parameter_summary,file.path(out,"anchored_posterior_summary.csv"),
          row.names=FALSE)
cat("\nCore parameter estimates and convergence:\n")
print(parameter_summary[match(required,parameter_summary$parameter),
      c("parameter","median","lower95","upper95","rhat","effective_n")],
      row.names=FALSE,digits=4)

# Derived annual baseline hazard h(t)=a exp(bt)+c1. Integration
# uses exprel to be stable for b approximately zero.
exprel <- function(z) {
  ans <- numeric(length(z))
  small <- abs(z)<1e-6
  ans[small] <- 1 + z[small]/2 + z[small]^2/6 + z[small]^3/24
  ans[!small] <- expm1(z[!small])/z[!small]
  ans
}
derived <- lapply(M,function(m) {
  A <- m[,"a"];B <- m[,"b"];C <- m[,"c1"]
  hazard <- function(age)A*exp(B*age)+C
  survival <- function(x,y)
    exp(-((y-x)*A*exp(B*x)*exprel((y-x)*B)+C*(y-x)))
  cbind(
    hazard_age_2=hazard(2),hazard_age_4=hazard(4),
    hazard_age_6=hazard(6),hazard_age_8=hazard(8),
    hazard_age_10=hazard(10),hazard_age_12=hazard(12),
    survival_age_2_to_6=survival(2,6),
    survival_age_2_to_8=survival(2,8),
    hazard_ratio_age_8_to_2=hazard(8)/hazard(2),
    infection_HR_age_4=exp(m[,"bi"]),
    infection_HR_age_8=exp(m[,"bi"]+4*m[,"bia"]),
    homozygosity_plus_1SD_HR_age_4=exp(m[,"bh"]),
    homozygosity_plus_1SD_HR_age_8=exp(m[,"bh"]+4*m[,"bha"]),
    male_HR=exp(m[,"bs"])
  )
})
derived_summary <- summarise_samples(derived)
write.csv(derived_summary,file.path(out,"anchored_derived_survival_convergence.csv"),
          row.names=FALSE)
cat("\nDerived mortality, survival and biological effects:\n")
print(derived_summary,row.names=FALSE,digits=4)

bychain <- do.call(rbind,lapply(seq_along(derived),function(j) {
  m <- derived[[j]]
  data.frame(chain=j,quantity=colnames(m),
             median=apply(m,2,median),
             lower95=apply(m,2,quantile,probs=.025),
             upper95=apply(m,2,quantile,probs=.975),
             row.names=NULL)
}))
write.csv(bychain,file.path(out,"anchored_derived_by_chain.csv"),
          row.names=FALSE)
cat("\nChain-specific median hazards and survival:\n")
keep <- c("hazard_age_4","hazard_age_8","survival_age_2_to_8",
          "hazard_age_12")
print(bychain[bychain$quantity %in% keep,],
      row.names=FALSE,digits=4)

# Within-chain drift of the baseline shape.
drift <- do.call(rbind,lapply(seq_along(M),function(j) {
  z <- cbind(M[[j]][,c("h4","f4","b","a","c1"),drop=FALSE],
             derived[[j]][,c("hazard_age_4","hazard_age_8",
                            "survival_age_2_to_8"),drop=FALSE])
  n <- nrow(z); half <- floor(n/2)
  data.frame(chain=j,parameter=colnames(z),
             first_half_mean=colMeans(z[1:half,,drop=FALSE]),
             second_half_mean=colMeans(z[(half+1):n,,drop=FALSE]),
             row.names=NULL)
}))
drift$change <- drift$second_half_mean-drift$first_half_mean
write.csv(drift,file.path(out,"anchored_chain_drift.csv"),
          row.names=FALSE)
cat("\nBaseline first-half/second-half drift:\n")
print(drift,row.names=FALSE,digits=4)

colours <- c("#1F77B4","#E45756","#2A9D8F")
png(file.path(out,"anchored_mortality_curves_by_chain.png"),
    width=1100,height=700,res=130)
ages <- seq(0,13,by=.25)
haz <- lapply(M,function(m)
  vapply(ages,function(age)median(m[,"a"]*exp(age*m[,"b"])+m[,"c1"]),
         numeric(1)))
plot(ages,haz[[1]],type="l",lwd=3,col=colours[1],
     ylim=range(unlist(haz)),xlab="Age (years)",
     ylab="Annual baseline mortality hazard",
     main="Anchored GM: mortality by posterior chain")
for(j in 2:3)lines(ages,haz[[j]],col=colours[j],lwd=3)
legend("topleft",legend=paste("Chain",1:3),
       col=colours,lwd=3,bty="n")
dev.off()

png(file.path(out,"anchored_baseline_parameter_traces.png"),
    width=1200,height=1000,res=120)
op <- par(mfrow=c(3,1),mar=c(3.7,4.2,2.1,1))
for(parname in c("h4","f4","b")) {
  allvalues <- unlist(lapply(M,function(m)m[,parname]))
  plot(seq_len(nrow(M[[1]])),M[[1]][,parname],type="l",
       col=colours[1],ylim=range(allvalues),
       xlab="Retained posterior draw",ylab=parname,main=parname)
  for(j in 2:3)lines(seq_len(nrow(M[[j]])),M[[j]][,parname],
                      col=colours[j])
}
par(op)
dev.off()

important <- c("h4","f4","b","a","c1",
               "hazard_age_4","hazard_age_8","survival_age_2_to_8")
rhat <- c(setNames(parameter_summary$rhat,parameter_summary$parameter),
          setNames(derived_summary$rhat,derived_summary$parameter))
ESS <- c(setNames(parameter_summary$effective_n,parameter_summary$parameter),
         setNames(derived_summary$effective_n,derived_summary$parameter))
cat("\n=== CONVERGENCE GATE ===\n")
cat("Maximum Rhat across all monitored parameters:",
    round(max(parameter_summary$rhat),4),"\n")
cat("Number of parameters with Rhat > 1.05:",
    sum(parameter_summary$rhat > 1.05),"\n")
cat("Most problematic monitored parameters:\n")
print(head(parameter_summary[order(-parameter_summary$rhat),
     c("parameter","rhat","effective_n")],12),row.names=FALSE,digits=4)
cat("\nPriority biological and baseline quantities:\n")
print(data.frame(parameter=important,
                 rhat=round(rhat[important],4),
                 ESS=round(ESS[important],1)),
      row.names=FALSE)
if(any(!is.finite(rhat[important])) ||
   any(rhat[important]>1.05) || any(ESS[important]<400)) {
  cat("\nNOT CONVERGED: do not interpret baseline ageing",
      "or pooled survival estimates as final.\n")
} else {
  cat("\nPRIORITY DIAGNOSTICS ACCEPTABLE (<1.05 Rhat, >=400 ESS).",
      "\nInspect all parameters and trace plots before interpretation.\n")
}
cat("NOTE: new anchored baseline prior differs from the old fit;",
    "this is a model sensitivity comparison, not just a sampler fix.\n")
cat("Saved Stage 62 CSVs and PNGs in",out,"\n")
