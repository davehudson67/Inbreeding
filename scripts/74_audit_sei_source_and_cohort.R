# ===========================================================================
# Stage 74: SAFE, READ-ONLY inventory of the source SEI diagnostic data
# and recruitment/monitoring time coverage of the 1603 genetic cohort.
#
# Purpose: identify which original SEI bundle tables can actually support
# a future joint diagnostic-observation + S/E/I transition model.
#
# NO model fitting, NO Supabase, NO writes to inputs, NO huge CSV loading.
# Only CSV headers and file sizes are sampled; raw per-badger records
# are not printed. Structured output is saved in outputs/sei_transitions.
# ===========================================================================
options(stringsAsFactors=FALSE)
root <- "data/infection_posteriors/sei_2026"
input <- "data/model_inputs/integrated_gm_inputs_group_proxy_2026.rds"
sei <- "data/model_inputs/integrated_gm_sei_aligned_2026.rds"
out <- "outputs/sei_transitions"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(!dir.exists(root)||!file.exists(input)||!file.exists(sei))
  stop("Missing existing local SEI source or Stage 66 model inputs.")
x <- readRDS(input)
z <- readRDS(sei)
stopifnot(length(x$id)==1603L,z$N==1603L,
          identical(z$id,toupper(trimws(as.character(x$id)))),
          identical(z$survival_input_md5,unname(tools::md5sum(input))))
cat("\n=== STAGE 74: ORIGINAL SEI SOURCE / OBSERVATION AUDIT ===\n")
cat("Matched age-verified genetic cohort:",length(x$id),"\n")
cat("Original SEI bundle: extracted locally; original ZIP verified:",
    isTRUE(z$original_zip_checksum_verified),"\n")
cat("This audit inventories source data, it DOES NOT infer new effects.\n")
cat("\n--- ORIGINAL BUNDLE README (up to 120 lines) ---\n")
readme <- file.path(root,"README_SEI_ONLY.md")
if(file.exists(readme)){
  lines <- readLines(readme,warn=FALSE)
  cat(paste(head(lines,120L),collapse="\n"),"\n")
  writeLines(lines,file.path(out,"stage74_original_sei_readme_copy.txt"))
} else {
  cat("README_SEI_ONLY.md not found.\n")
}
cat("\n--- CSV FILE INVENTORY (headers, sizes, no data rows) ---\n")
paths <- sort(list.files(root,pattern="[.]csv$",recursive=TRUE,
                         full.names=TRUE))
if(!length(paths))stop("No CSV tables found in local SEI source")
headers <- lapply(paths,function(path){
  # Reading nrows=0 only parses the header (no giant draw tables).
  h <- tryCatch(
    names(utils::read.csv(path,nrows=0,check.names=FALSE)),
    error=function(e)paste0("HEADER_ERROR: ",conditionMessage(e))
  )
  data.frame(
    file=sub(paste0("^",root,"/"),"",path),
    bytes=as.numeric(file.info(path)$size),
    columns=paste(h,collapse=" | "),
    row.names=NULL)
})
inventory <- do.call(rbind,headers)
write.csv(inventory,file.path(out,"stage74_sei_csv_inventory.csv"),
          row.names=FALSE)
for(k in seq_len(nrow(inventory))){
  cat(sprintf("%s [%s bytes]: %s\n",
              inventory$file[k],
              format(inventory$bytes[k],scientific=FALSE,trim=TRUE),
              inventory$columns[k]))
}
cat("\n--- SOURCE BUNDLE MANIFEST (up to 80 lines) ---\n")
manifest <- file.path(root,"bundle_manifest.json")
if(file.exists(manifest)){
  lines <- readLines(manifest,warn=FALSE)
  cat(paste(head(lines,80L),collapse="\n"),"\n")
} else {
  cat("bundle_manifest.json not found.\n")
}

# Monitor independent cohort selection relative to the 2025 study end.
# Assumes calendar indexing from original integrated input T=200:
# 1976 Q1 = q1, 2016 Q1 = q161.
first <- as.integer(x$first_q)
birth <- as.integer(x$birth_q)
last <- as.integer(z$sei_last_capture_q)
death <- as.integer(z$death_q)
stopifnot(x$T==200L,all(is.finite(first)),
          all(first>=1L & first<=200L),
          all(first>=birth),all(last>=first),
          length(death)==length(first))
latest_live <- vapply(seq_along(first),function(i){
  q <- which(x$y[i,]==1L)
  if(length(q))as.integer(max(q))else NA_integer_
},integer(1))
if(anyNA(latest_live))stop("An included badger has no live encounter")
periods <- data.frame(
  period=c("1976-1985","1986-1995","1996-2005",
           "2006-2015","2016-2025"),
  first=seq.int(1L,161L,by=40L),
  last=seq.int(40L,200L,by=40L))
obs <- do.call(rbind,lapply(seq_len(nrow(periods)),function(k){
  lo <- periods$first[k];hi <- periods$last[k]
  data.frame(
    decade=periods$period[k],
    n_first_live=sum(first>=lo & first<=hi),
    n_last_live=sum(latest_live>=lo & latest_live<=hi),
    n_any_live_capture=sum(rowSums(x$y[,lo:hi,drop=FALSE]==1L)>0L),
    total_live_capture_quarters=sum(x$y[,lo:hi,drop=FALSE]==1L),
    n_sei_last_capture=sum(last>=lo & last<=hi),
    n_death=sum(death>0L & death>=lo & death<=hi)
  )
}))
cat("\n--- COHORT AGE-VERIFICATION / DATA-COVERAGE SUMMARY ---\n")
print(obs,row.names=FALSE)
write.csv(obs,file.path(out,"stage74_cohort_time_coverage.csv"),
          row.names=FALSE)
cat("\nCHECK: 2016-2025 may still include repeat captures despite",
    "zero first-capture entrants in this selected cohort.\n")
cat("Interpret model-based infection changes only after comparing",
    "diagnostic methods, test data availability, capture opportunities,",
    "cohort selection and SEI prior assumptions by historical period.\n")
cat("STAGE 74 COMPLETE; all source files unchanged.\n")
