# ===========================================================================
# Stage 75: Assess original SEI mode, diagnostic-test evidence and whether
# raw capture-by-test observations are available inside the fitted R object.
# STRICTLY READ ONLY, no fitting, no Supabase, no modifying source files.
# Avoid outputting thousands of individual IDs or posterior columns.
# ===========================================================================
options(stringsAsFactors=FALSE)
root <- "data/infection_posteriors/sei_2026"
out <- "outputs/sei_transitions"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
paths <- setNames(file.path(root,c(
  "diagnostics.csv",
  "diagnostics_by_parameter.csv",
  "init_optima.csv",
  "sensitivity_specificity.csv",
  "badgers.csv",
  "fits/sei_inferred_accuracy_.rds"
)),c("diagnostics","parameters","starts","tests","badgers","fit"))
if(!all(file.exists(paths)))
  stop("Missing local SEI original source files: ",
       paste(names(paths)[!file.exists(paths)],collapse=", "))
read_table <- function(path)
  utils::read.csv(path,check.names=FALSE,stringsAsFactors=FALSE)
cat("\n=== STAGE 75: ORIGINAL SEI MODEL GATES ===\n")
cat("Analysed source: extracted local SEI inferred-accuracy bundle\n")
cat("No new biological model is fitted.\n")

# Original fit summary, not posterior draws.
diagnostics <- read_table(paths["diagnostics"])
cat("\n--- FIT DIAGNOSTICS ---\n")
if(nrow(diagnostics)>0L) {
  print(diagnostics,row.names=FALSE,digits=5)
} else {
  cat("No diagnostic rows found\n")
}

parameters <- read_table(paths["parameters"])
cat("\n--- WORST-MIXING PARAMETER DIAGNOSTICS ---\n")
if(all(c("parameter","rhat","ess_bulk")%in%names(parameters))) {
  parameters$rhat <- suppressWarnings(as.numeric(parameters$rhat))
  parameters$ess_bulk <- suppressWarnings(as.numeric(parameters$ess_bulk))
  cat("Number of parameters:",nrow(parameters),"\n")
  cat("Finite Rhat values:",sum(is.finite(parameters$rhat)),"\n")
  cat("Rhat > 1.01:",sum(parameters$rhat>1.01,na.rm=TRUE),"\n")
  cat("Rhat > 1.05:",sum(parameters$rhat>1.05,na.rm=TRUE),"\n")
  cat("Bulk ESS below 100:",sum(parameters$ess_bulk<100,na.rm=TRUE),"\n")
  cat("Bulk ESS below 400:",sum(parameters$ess_bulk<400,na.rm=TRUE),"\n")
  ord <- order(-ifelse(is.finite(parameters$rhat),parameters$rhat,-Inf),
               ifelse(is.finite(parameters$ess_bulk),
                      parameters$ess_bulk,Inf))
  print(utils::head(parameters[ord,c(
    intersect(c("model","parameter","rhat","ess_bulk"),
              names(parameters))),drop=FALSE],12L),
    row.names=FALSE,digits=4)
} else {
  cat("Unrecognized parameter diagnostics schema:\n")
  print(names(parameters))
}

starts <- read_table(paths["starts"])
cat("\n--- OPTIMIZED STARTS / ALTERNATIVE MODE ---\n")
if(nrow(starts)>0) print(starts,row.names=FALSE,digits=6)
cat("IMPORTANT: A smaller mode is described in the source README.",
    "Optimization results alone do not establish posterior mode masses",
    "or mixing between modes.\n")
write.csv(diagnostics,file.path(out,"stage75_fit_diagnostics.csv"),row.names=FALSE)
write.csv(starts,file.path(out,"stage75_init_optima.csv"),row.names=FALSE)

test <- read_table(paths["tests"])
cat("\n--- ORIGINAL TEST SENSITIVITY / SPECIFICITY ---\n")
if(nrow(test)) print(test,row.names=FALSE,digits=5)
write.csv(test,file.path(out,"stage75_test_characteristics.csv"),row.names=FALSE)
cat("Test estimates may be model-inferred; consult 'estimated' column.\n")
cat("These rows describe test parameters, NOT observed",
    "capture-by-test histories.\n")

# Full SEI cohort recruitment and follow-up (3214, not just the
# age-verified 1603 genetic sample).
b <- read_table(paths["badgers"])
req <- c("tattoo","first_capture_time","last_capture_time",
         "start_time","n_captures")
if(!all(req %in% names(b))) stop("Unexpected badgers.csv schema")
stopifnot(nrow(b)==3214L)
b$first_capture_time <- as.integer(b$first_capture_time)
b$last_capture_time <- as.integer(b$last_capture_time)
b$start_time <- as.integer(b$start_time)
b$n_captures <- as.integer(b$n_captures)
stopifnot(!anyNA(b[,req[-1L]]),
          all(b$first_capture_time>=1L),
          all(b$first_capture_time<=b$last_capture_time),
          all(b$last_capture_time<=200L))
decades <- data.frame(
  decade=c("1976-1985","1986-1995","1996-2005",
           "2006-2015","2016-2025"),
  lo=seq.int(1L,161L,by=40L),
  hi=seq.int(40L,200L,by=40L))
full_cohort <- do.call(rbind,lapply(seq_len(nrow(decades)),function(k){
  lo <- decades$lo[k];hi <- decades$hi[k]
  ix <- which(b$first_capture_time>=lo&
              b$first_capture_time<=hi)
  data.frame(
    decade=decades$decade[k],
    first_captures=length(ix),
    last_capture_records=sum(b$last_capture_time>=lo&
                             b$last_capture_time<=hi),
    total_capture_occasions=sum(b$n_captures[ix]),
    median_captures_per_badger=if(length(ix))
      median(b$n_captures[ix])else NA_real_,
    median_followup_quarters=if(length(ix))
      median(b$last_capture_time[ix]-b$first_capture_time[ix]+1L)
      else NA_real_
  )
}))
cat("\n--- ALL 3214 SEI BADGERS: RECRUITMENT DECADES ---\n")
print(full_cohort,row.names=FALSE)
write.csv(full_cohort,file.path(out,"stage75_full_sei_cohort_decades.csv"),
          row.names=FALSE)
cat("\nSEI start earlier than first capture:",
    sum(b$start_time<b$first_capture_time),"\n")
cat("SEI start earlier than recorded birth (when present): ")
if("birth_time"%in%names(b)) {
  valid_birth <- !is.na(b$birth_time)&nzchar(as.character(b$birth_time))
  birth_time <- suppressWarnings(as.integer(b$birth_time))
  cat(sum(b$start_time<birth_time,na.rm=TRUE),
      "of",sum(!is.na(birth_time)),"with birth_time\n")
} else {
  cat("Not recorded\n")
}

# Inspect the fitted object's STRUCTURE without printing an array of
# posterior draws, sensitive tattoo-level records or enormous fits.
inspect_structure <- function(v,max_names=35L){
  typ <- paste(class(v),collapse=", ")
  d <- dim(v)
  n <- if(is.null(d)) length(v) else paste(d,collapse="x")
  nm <- if(isS4(v)) methods::slotNames(v) else
    if(is.environment(v))ls(v,all.names=TRUE) else names(v)
  if(is.null(nm))nm <- character()
  list(
    class=typ,typeof=typeof(v),shape=as.character(n),
    n_named_fields=length(nm),
    first_field_names=paste(head(nm,max_names),collapse=", "),
    possible_observation_fields=paste(
      grep("data|test|obs|capture|history|infect|diagnos|likeli|prior|param",
           nm,ignore.case=TRUE,value=TRUE),
      collapse=", ")
  )
}
cat("\n--- FITTED SEI R OBJECT: STRUCTURAL INSPECTION ---\n")
cat("Fitted RDS bytes:",
    format(as.numeric(file.info(paths["fit"])$size),
           scientific=FALSE,trim=TRUE),"\n")
fit <- tryCatch(readRDS(paths["fit"]),error=function(e)e)
if(inherits(fit,"error")) {
  cat("Cannot read fitted object:",conditionMessage(fit),"\n")
  cat("Source raw diagnostic availability remains UNKNOWN.\n")
} else {
  top <- inspect_structure(fit)
  for(key in names(top))
    cat(key,":",top[[key]],"\n")
  if(is.list(fit) && !is.object(fit)) {
    found <- grep("data|test|obs|capture|history|infect|diagnos|model|fit",
                  names(fit),ignore.case=TRUE,value=TRUE)
    for(key in head(found,12L)) {
      if(!is.null(fit[[key]])){
        cat("\n  FIT TOP-LEVEL FIELD:",key,"\n")
        info <- inspect_structure(fit[[key]],max_names=18L)
        for(nm in names(info))
          cat("   ",nm,":",info[[nm]],"\n")
      }
    }
  }
  cat("Only structural metadata was inspected; this does not confirm",
      "that raw diagnostic results can be recovered from the object.\n")
}
cat("\nSOURCE GATE:\n")
cat("Bundle CSVs include inferred test characteristics and state",
    "probabilities, not an identified raw capture-by-test table.\n")
cat("We need original diagnostic outcomes, test availability by era,",
    "SEI likelihood and priors, and evidence on the second mode",
    "before a confirmatory joint S/E/I + genetic analysis.\n")
cat("STAGE 75 COMPLETE; no inputs modified.\n")
