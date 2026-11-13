# ===============================================================
# Summarise completed integrated-GM chains (posterior diagnostics)
# Checks exact data/code provenance before any pooling.
# ===============================================================
options(stringsAsFactors=FALSE)
stopifnot(requireNamespace("coda",quietly=TRUE))
out <- "outputs/integrated_gm"
files <- file.path(out,sprintf("full_chain_%02d.rds",1:3))
if(!all(file.exists(files)))
  stop("Require all three completed chain RDS files:\n",
       paste(files[!file.exists(files)],collapse="\n"))
chains <- lapply(files,readRDS)
n <- vapply(chains, function(z) z$n, integer(1))
hashes <- vapply(chains,function(z) as.character(z$input_md5),character(1))
code_hashes <- vapply(chains,function(z)
  paste(z$code_md5,collapse=":"),character(1))
if(length(unique(n))!=1L || length(unique(hashes))!=1L ||
   length(unique(code_hashes))!=1L)
  stop("Chains refer to different data/code; do not combine")
cat("\n=== INTEGRATED GM: THREE-CHAIN ANALYSIS ===\n")
cat("Badgers:",unique(n),"\n")
cat("Input MD5:",unique(hashes),"\n")
m <- lapply(chains,function(z) as.matrix(z$samples))
nameslist <- lapply(m,colnames)
if(!all(vapply(nameslist,identical,logical(1),nameslist[[1L]])))
  stop("Monitor names/order differ across chain outputs")
if(any(vapply(m,nrow,integer(1)) < 100L))
  warning("Fewer than 100 retained draws in at least one chain")
# Fixed reference period and reference season are zero; drop all
# constant columns before applying Gelman-Rubin.
allcols <- colnames(m[[1]])
vary <- vapply(seq_along(allcols),function(j) {
  z <- unlist(lapply(m,function(a) a[,j]))
  all(is.finite(z)) && sd(z)>1e-10
},logical(1))
cat("Excluded constant or invalid monitors:",
    paste(allcols[!vary],collapse=", "),"\n")
if(sum(vary)<3L)stop("No variable monitored parameters")
m <- lapply(m,function(a) a[,vary,drop=FALSE])
pnames <- colnames(m[[1]])
stopifnot(all(vapply(m,function(a) all(is.finite(a)),logical(1))))
ml <- coda::mcmc.list(lapply(m,coda::mcmc))
gr <- tryCatch(coda::gelman.diag(
  ml,autoburnin=FALSE,multivariate=FALSE)$psrf[,1],
  error=function(e) {
    warning("Gelman diagnostics failed: ",conditionMessage(e))
    setNames(rep(NA_real_,length(pnames)),pnames)
  }
)
ess <- coda::effectiveSize(ml)
combined <- do.call(rbind,m)
q <- apply(combined,2,quantile,probs=c(.025,.5,.975),na.rm=TRUE)
summary <- data.frame(
  parameter=pnames,
  mean=colMeans(combined),
  median=q[2,],
  lower95=q[1,],
  upper95=q[3,],
  rhat=as.numeric(gr[pnames]),
  effective_n=as.numeric(ess[pnames]),
  row.names=NULL
)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(summary,file.path(out,"full_posterior_summary.csv"),
          row.names=FALSE)
cat("\nKey biological parameters:\n")
key <- summary$parameter %in% c("a","b","c1","bh","bs","bi",
                                "bhi","bha","bia",
                                "alphaP","sigmaU","rho")
print(summary[key,],row.names=FALSE,digits=4)
cat("\nMaximum Rhat:",max(summary$rhat,na.rm=TRUE),"\n")
cat("Minimum bulk-like effective size (coda ESS):",
    min(summary$effective_n,na.rm=TRUE),"\n")
bad <- summary[is.finite(summary$rhat) & summary$rhat>1.05,]
cat("Parameters with Rhat > 1.05:",nrow(bad),"\n")
if(nrow(bad))print(bad[,c("parameter","rhat","effective_n")],
                    row.names=FALSE)
cat("\nImportant caveats:\n")
cat("- Survival estimates conditional on social-group monitoring proxy\n")
cat("- Genetic hom was standardized across included animals\n")
cat("- Infection uncertainty integrated per animal; shared-draw",
    "dependence not represented\n")
cat("- Recovery probability assumes common rho across years/causes\n")
cat("- All-cause apparent survival; unmodelled emigration/availability\n")
cat("Full posterior summaries saved; convergence must be reviewed",
    "before biological interpretation.\n")
