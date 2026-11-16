# ======================================================================
# Stage 66: Audit and align coherent S -> E -> I posterior trajectories
# with the existing 1603-badger integrated GM survival input.
#
# Input: checksum-verified private Supabase Storage SEI bundle fetched by
#        scripts/65_fetch_sei_from_supabase.py
# Output: safe, reproducible, model-ready draw-indexed joint SEI histories
#         + comprehensive overlap and post-last-capture audit.
#
# NO MCMC fits. NO replacement of existing infection input or Supabase
# database. Crucially, unobserved disease state after last capture is NA.
# ======================================================================
options(stringsAsFactors=FALSE)
if(!requireNamespace("data.table",quietly=TRUE))
  stop("This stage needs data.table; use install.packages('data.table') once.")
suppressPackageStartupMessages(library(data.table))
root <- "data/infection_posteriors/sei_2026"
archive <- file.path(root,"SEI_only_research_bundle.zip")
source_manifest <- file.path(root,"sei_storage_source_verified.json")
input <- Sys.getenv("MODEL_INPUT",
 "data/model_inputs/integrated_gm_inputs_group_proxy_2026.rds")
out <- "outputs/sei_survival"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(source_manifest))
  stop("Missing SEI source manifest. First run scripts/65_fetch_sei_from_supabase.py")
if(!file.exists(input))
  stop("Missing integrated survival inputs: ",input)
x <- readRDS(input)
stopifnot(x$N==length(x$id),x$T==200L,
          all(x$y %in% c(-1L,0L,1L)))
expected_sha256 <- "65a1f1a7c84e2f9bfbbe149a3d92078cebd185e72281d7e5452cf0eadaaf9f98"
# Accept either a checksum-verified original ZIP or an explicitly
# fingerprinted extracted source. Do NOT attribute the ZIP's SHA256
# to extracted files when the original ZIP is unavailable.
zip_verified <- file.exists(archive)
sha256file <- function(filename) {
  if(!nzchar(Sys.which("sha256sum")))
    stop("sha256sum is required to audit SEI input files")
  val <- system2("sha256sum",shQuote(filename),
                 stdout=TRUE,stderr=FALSE)
  if(length(val)!=1L)
    stop("Failed to calculate SHA256 for ",filename)
  strsplit(val," ",fixed=TRUE)[[1]][1]
}
if(zip_verified) {
  if(file.info(archive)$size!=23662164)
    stop("SEI ZIP has unexpected size: halt")
  if(!identical(sha256file(archive),expected_sha256))
    stop("SEI ZIP checksum changed after verification: halt")
}

norm <- function(z) toupper(trimws(as.character(z)))
ids <- norm(x$id)
if(anyNA(ids) || anyDuplicated(ids))stop("Bad or duplicate survival IDs")
cat("\n=== STAGE 66: SUPABASE SEI -> SURVIVAL AUDIT ===\n")
cat("Survival IDs:",length(ids),"quarterly occasions:",x$T,"\n")
cat("Supabase private source: research-model-outputs/SEI_only_research_bundle.zip\n")
cat("Original ZIP SHA256 verified:",zip_verified,"\n")
if(zip_verified) {
  cat("Original ZIP SHA256:",expected_sha256,"\n")
} else {
  cat("Original ZIP unavailable: checking extracted CSV fingerprints.\n")
}
cat("Input MD5:",unname(tools::md5sum(input)),"\n")

# In verified-ZIP mode, always freshly extract these two files from
# the verified archive (never trust potentially stale extracted CSVs).
# Otherwise locate the ONE uploaded directory containing both files.
required <- c("badgers.csv","trajectories_by_draw.csv")
if(zip_verified) {
  exdir <- file.path(root,"extracted")
  dir.create(exdir,recursive=TRUE,showWarnings=FALSE)
  cat("Extracting verified source tables into:",exdir,"\n")
  utils::unzip(archive,files=required,exdir=exdir,overwrite=TRUE)
} else {
  matches <- list.files(root,pattern="^badgers[.]csv$",
                        recursive=TRUE,full.names=TRUE)
  sources <- unique(dirname(matches[
    file.exists(file.path(dirname(matches),"trajectories_by_draw.csv"))
  ]))
  if(length(sources)!=1L)
    stop("Need exactly one folder containing badgers.csv AND ",
         "trajectories_by_draw.csv. Found ",length(sources),
         ". Resolve duplicate extracted copies.")
  exdir <- sources[[1L]]
  cat("Using already-extracted SEI tables from:",exdir,"\n")
  manifest_fingerprints <- file.path(root,"sei_extracted_files_sha256.tsv")
  if(!file.exists(manifest_fingerprints))
    stop("Missing extracted input SHA256 manifest: rerun Stage 65")
  fingerprint <- read.delim(manifest_fingerprints,
                            stringsAsFactors=FALSE)
  if(!all(c("filename","sha256") %in% names(fingerprint)) ||
     anyDuplicated(fingerprint$filename) ||
     !all(required %in% fingerprint$filename))
    stop("Invalid extracted source fingerprint manifest")
  for(nm in required) {
    actual <- sha256file(file.path(exdir,nm))
    expected <- fingerprint$sha256[match(nm,fingerprint$filename)]
    if(!identical(actual,expected))
      stop("Extracted file changed since Stage 65: ",nm)
  }
  cat("Extracted CSV SHA256 fingerprints: PASS\n")
}
if(!all(file.exists(file.path(exdir,required))))
  stop("Missing essential SEI CSV table after source resolution")
source_file_sha256 <- setNames(vapply(required,function(nm)
  sha256file(file.path(exdir,nm)),character(1)),required)
badgers <- data.table::fread(file.path(exdir,"badgers.csv"),
                            colClasses=list(character="tattoo"))
traj <- data.table::fread(file.path(exdir,"trajectories_by_draw.csv"),
                         colClasses=list(character=c("model","tattoo"),
                           integer=c("draw","infection_time","infectious_time")))
needb <- c("tattoo","start_time","first_capture_time","last_capture_time")
needt <- c("model","tattoo","draw","infection_time","infectious_time")
stopifnot(all(needb %in% names(badgers)),
          all(needt %in% names(traj)))
badgers[,tattoo:=norm(tattoo)]
traj[,tattoo:=norm(tattoo)]
if(nrow(badgers)!=3214L || anyDuplicated(badgers$tattoo))
  stop("Badgers identity inventory disagrees with SEI research bundle")
if(nrow(traj)!=3214000L ||
   data.table::uniqueN(traj$tattoo)!=3214L ||
   data.table::uniqueN(traj$draw)!=1000L)
  stop("The SEI trajectory table is not the expected 3214 x 1000 draws")
if(!identical(sort(unique(traj$draw)),1:1000))
  stop("Expected 1000 draw IDs labelled 1..1000")
if(anyNA(traj[, ..needt]))
  stop("Missing values in SEI onset/draw index columns")
if(!all(traj$model=="SEI (inferred accuracy)"))
  stop("Unexpected model labels: SEI-only bundle required")
if(any(traj$infection_time<0L | traj$infection_time>200L |
       traj$infectious_time<0L | traj$infectious_time>200L))
  stop("SEI posterior transition quarters outside 0..200")
if(any(traj$infectious_time>0L &
       (traj$infection_time==0L |
        traj$infectious_time<traj$infection_time)))
  stop("SEI pathology: infectiousness before infection/exposure")
if(anyNA(badgers[,..needb]))
  stop("Missing disease-history capture-window data")
if(any(badgers$start_time>badgers$last_capture_time))
  stop("SEI start time later than SEI last capture")

badgers[,present_in_survival:=tattoo %chin% ids]
coverage <- data.table(tattoo=ids,
                       sei_index=match(ids,badgers$tattoo))
coverage[,in_sei:=!is.na(sei_index)]
data.table::fwrite(coverage[,.(tattoo,in_sei)],
  file.path(out,"sei_survival_id_overlap.csv"))
cat("SEI IDs:",nrow(badgers),"\n")
cat("Survival IDs matched:",sum(coverage$in_sei),
    "of",nrow(coverage),"\n")
if(!all(coverage$in_sei)) {
  cat("Unmatched IDs:\n")
  print(coverage[!in_sei,.(tattoo)],nrows=50L)
  stop("Unmatched survival identities; exported ID report. Resolve before model fitting.")
}

N <- length(ids)
sei_start <- as.integer(badgers$start_time[coverage$sei_index])
sei_first <- as.integer(badgers$first_capture_time[coverage$sei_index])
sei_last <- as.integer(badgers$last_capture_time[coverage$sei_index])
stopifnot(all(is.finite(sei_start)),all(is.finite(sei_last)))
# Restrict to model cohort before reshaping 1.6 million observations.
selection <- traj[tattoo %chin% ids]
if(nrow(selection)!=N*1000L ||
   data.table::uniqueN(selection,by=c("tattoo","draw"))!=N*1000L)
  stop("Duplicate or missing per-individual/per-draw SEI histories")
selection[, individual_index:=match(tattoo,ids)]
setorder(selection,draw,individual_index)
if(!identical(as.integer(selection$individual_index),
              rep(seq_len(N),times=1000L)) ||
   !identical(as.integer(selection$draw),
              rep(seq_len(1000L),each=N)))
  stop("Posterior draw alignment failed: do not transpose/mix individual draws")

E <- matrix(as.integer(selection$infection_time),nrow=N,ncol=1000L,
            dimnames=list(ids,as.character(1:1000)))
I <- matrix(as.integer(selection$infectious_time),nrow=N,ncol=1000L,
            dimnames=list(ids,as.character(1:1000)))
rm(traj,selection)
gc(verbose=FALSE)
if(any(E>matrix(sei_last,nrow=N,ncol=1000L)) ||
   any(I>matrix(sei_last,nrow=N,ncol=1000L)))
  stop("Disease transition after last SEI capture; inspect source times")
if(any(I>0L & (E==0L | I<E)))
  stop("Draw-specific SEI paths violate S->E->I ordering")
if(any(E>0L & E<matrix(sei_start,nrow=N,ncol=1000L))) {
  warning("Some inferred onset times precede the modelled SEI start window")
}

# Compare independent data reconciliation windows.
live_last <- vapply(seq_len(N),function(i) {
  q <- which(x$y[i,]==1L)
  if(!length(q)) return(NA_integer_)
  max(q)
},integer(1))
death <- as.integer(ifelse(is.na(x$death_q),0L,x$death_q))
supported_end <- pmin(sei_last,x$T)
post_last_monitored <- vapply(seq_len(N),function(i){
  end <- supported_end[i]
  if(end>=x$T)return(0L)
  sum(x$y[i,seq.int(end+1L,x$T)]==0L)
},integer(1))
audit <- data.table(
  tattoo=ids,
  model_first_q=as.integer(x$first_q),
  model_last_live_q=live_last,
  sei_start_q=sei_start,
  sei_first_capture_q=sei_first,
  sei_last_capture_q=sei_last,
  death_q=death,
  model_last_live_after_sei_last=live_last>sei_last,
  death_after_sei_last=death>sei_last,
  n_monitored_nondetections_post_sei_last=post_last_monitored,
  fraction_draws_exposed=as.numeric(rowMeans(E>0L)),
  fraction_draws_infectious=as.numeric(rowMeans(I>0L))
)
data.table::fwrite(audit,file.path(out,"sei_vs_survival_followup.csv"))
cat("\n=== LIVE CAPTURE AND DEATH FOLLOW-UP ===\n")
cat("Survival live encounters after SEI last capture:",
    sum(audit$model_last_live_after_sei_last,na.rm=TRUE),"\n")
cat("Known deaths after SEI last capture:",
    sum(audit$death_after_sei_last),"\n")
cat("Animals with monitored noncaptures after SEI last capture:",
    sum(post_last_monitored>0),"\n")
cat("Total monitored noncaptures AFTER SEI last capture:",
    sum(post_last_monitored),"\n")
cat("Median SEI last-capture quarter:",median(sei_last),"\n")
cat("Fraction of 1.603m aligned draws with E onset recorded:",
    mean(E>0L),"\n")
cat("Fraction with I onset recorded:",mean(I>0L),"\n")

# Preserve the 1000 SHARED columns across badgers. Using independent
# draw resampling per badger would destroy joint dependence.
# Zero is only "not transitioned by last capture", not "never".
result <- list(
  id=ids, N=N, draw_ids=seq_len(1000L), T=x$T,
  exposed_onset_q=E, infectious_onset_q=I,
  sei_start_q=sei_start,sei_last_capture_q=sei_last,
  model_first_q=as.integer(x$first_q),
  model_last_live_q=live_last, death_q=death,
  source_mode=if(zip_verified) "verified_original_zip" else "extracted_csvs",
  original_zip_checksum_verified=zip_verified,
  source_archive=if(zip_verified) archive else NA_character_,
  source_sha256=if(zip_verified) expected_sha256 else NA_character_,
  source_extracted_directory=exdir,
  source_csv_sha256=source_file_sha256,
  survival_input=input,
  survival_input_md5=unname(tools::md5sum(input)),
  warnings=c(
    "infection_time=0 means not infected through SEI last capture, not never afterwards",
    "infectious_time=0 means not infectious through SEI last capture",
    "State is UNKNOWN before SEI start and after SEI last capture",
    "1000 columns are coherent joint posterior draws across individuals",
    "SEI disease posterior does not itself estimate disease progress after last capture",
    "First capture is conditioned upon in the current survival model",
    "No individual survival model has been fitted to these SEI histories",
    if(!zip_verified)
      "Original Supabase ZIP checksum could not be verified; essential CSV files fingerprinted and schema/chronology audited"
    else
      "Original Supabase ZIP checksum verified"
  )
)
dir.create("data/model_inputs",recursive=TRUE,showWarnings=FALSE)
dest <- "data/model_inputs/integrated_gm_sei_aligned_2026.rds"
saveRDS(result,dest,compress="xz")
cat("\nAligned SEI research input saved:",dest,"\n")
cat("Dimensions of E and I onset matrices:",
    paste(dim(E),collapse=" x "),"\n")
cat("\n=== BIOLOGICAL INTERPRETATION GATE ===\n")
cat("SEI transition times ending at last capture CANNOT be extrapolated",
    "as if 0 meant permanently uninfected.\n")
cat("A survival+SEI fit must model/extrapolate post-last disease",
    "progression, or explicitly censor/condition the survival data.\n")
cat("Do not treat averaged quarterly S/E/I probabilities as full",
    "joint disease-history uncertainty.\n")
cat("STAGE 66 ALIGNMENT AND AUDIT COMPLETED; no SEI survival fit yet.\n")
