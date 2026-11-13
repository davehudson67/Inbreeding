# ===================================================================
# Stage 56: FINAL GENETIC-COVARIATE ORIENTATION GATE
# An independent check against the historical published hom covariate.
# Stage 55 simulations cannot detect an inverted real genetic variable.
# No database or source file modifications.
# ===================================================================
options(stringsAsFactors=FALSE)
gene <- readRDS("data/badger/genetics/legacy_workspace_gene.rds")
pub <- readRDS("data/badger/legacy_published/published_2023_analytical_inputs.rds")
x <- readRDS("data/model_inputs/integrated_gm_inputs_group_proxy_2026.rds")
norm <- function(z) toupper(trimws(as.character(z)))
stopifnot(all(c("ID","mlh") %in% names(gene)),
          is.numeric(gene$mlh),
          !anyNA(gene$mlh))
gene$ID <- norm(gene$ID)
if(anyDuplicated(gene$ID))
  stop("Genetic IDs duplicated: cannot verify orientation")
names(x$hom) <- x$id
if(any(abs(unname(x$hom)-
         (1-gene$mlh[match(x$id,gene$ID)]))>1e-12))
  stop("Current fitted hom values disagree with current builder")
cat("\n=== GENETIC COVARIATE ORIENTATION ===\n")
cat("Current builder uses: hom = 1 - gene$mlh\n")
cat("Original legacy modelling code uses published CH$hom\n")
cat("Gene mlh summary:\n")
print(summary(gene$mlh))
refs <- list()
for(section in c("CH","mlh")) {
  v <- pub[[section]]
  if(!is.data.frame(v) || !"ID" %in% names(v)) next
  for(col in c("hom","hom.x","hom.y")) {
    if(col %in% names(v) && is.numeric(v[[col]])) {
      z <- data.frame(ID=norm(v$ID),value=v[[col]])
      z <- z[!is.na(z$ID) & !is.na(z$value),]
      if(!nrow(z)) next
      # Drop duplicated identities only if source values agree.
      inconsistent <- vapply(split(z$value,z$ID),function(u)
        diff(range(u))>1e-8, logical(1))
      if(any(inconsistent)) {
        cat(section,"$",col,": source has",
            sum(inconsistent),"IDs with conflicting values;",
            "using median as reference only\n")
        z <- aggregate(value~ID,z,median)
      } else {
        z <- z[!duplicated(z$ID),]
      }
      refs[[paste(section,col,sep="$")]] <- z
    }
  }
}
if(!length(refs))
  stop("No historical hom covariate available; identity must be verified")
tab <- do.call(rbind,lapply(names(refs),function(nm) {
  z <- refs[[nm]]
  k <- match(z$ID,gene$ID)
  ok <- !is.na(k)
  if(sum(ok) < 20) return(NULL)
  a <- gene$mlh[k[ok]]
  h <- z$value[ok]
  data.frame(
    historical_source=nm,n_matched=sum(ok),
    mean_abs_diff_if_mlh=mean(abs(a-h)),
    mean_abs_diff_if_1_minus_mlh=mean(abs(1-a-h)),
    correlation_mlh=if(sd(a)>0 && sd(h)>0) cor(a,h) else NA_real_,
    correlation_1_minus_mlh=if(sd(a)>0 && sd(h)>0) cor(1-a,h) else NA_real_,
    stringsAsFactors=FALSE
  )
}))
if(is.null(tab)||!nrow(tab))
  stop("Insufficient historical overlap to validate genetic orientation")
print(tab,row.names=FALSE,digits=6)
dir.create("outputs/model_audit",recursive=TRUE,showWarnings=FALSE)
write.csv(tab,
  "outputs/model_audit/integrated_genetic_orientation.csv",
  row.names=FALSE)
row <- which.max(tab$n_matched)
if(tab$n_matched[row] < 500L)
  stop("GENETIC ORIENTATION UNCERTAIN: fewer than 500 matched records")
err_same <- tab$mean_abs_diff_if_mlh[row]
err_flip <- tab$mean_abs_diff_if_1_minus_mlh[row]
if(err_same + 1e-4 < err_flip) {
  stop(
    "CRITICAL: gene$mlh matches historical hom better than 1-gene$mlh. ",
    "Current builder may invert the inbreeding variable. ",
    "DO NOT RUN REAL MCMC until corrected."
  )
}
if(err_flip + 1e-4 < err_same) {
  cat("ORIENTATION GATE PASSED: current 1 - mlh matches the",
      "historical hom direction better.\n")
} else {
  stop("GENETIC ORIENTATION AMBIGUOUS: manual genotype-scale check needed.")
}
cat("\nThis gate verifies direction, not absolute equality; inspect",
    "the differences for historical transformations/rounding.\n")
