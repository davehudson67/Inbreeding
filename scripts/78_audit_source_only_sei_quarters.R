# ============================================================================
# STAGE 78: Identify why 317 original SEI badger-quarter occasions
# are absent from the strict method-confirmed live-capture reconstruction.
#
# Reads source p_infected_per_capture.csv and private Stage76 capture
# crosswalk. Queries ONLY capture ID, tattoo, date, MOC, is_pm from
# analysis.capture_history_v1. No writes to Supabase; no new model fits.
#
# Distinguishes:
#   (a) original occasion present in DB with another capture method,
#   (b) original occasion present only as PM/dead record,
#   (c) no exact-quarter encounter in Supabase,
#   (d) unexpected match to a strict-live record (integrity failure).
# Also reports neighbouring-quarter encounters but never reassigns dates.
#
# Saves individual IDs to gitignored outputs/sei_transitions/ ONLY.
# ============================================================================
options(stringsAsFactors=FALSE)
if(!requireNamespace("DBI",quietly=TRUE)||
   !requireNamespace("RPostgres",quietly=TRUE))
  stop("R packages DBI and RPostgres required")
origdir <- "data/infection_posteriors/sei_2026"
out <- "outputs/sei_transitions"
source_file <- file.path(origdir,"p_infected_per_capture.csv")
candidate_file <- file.path(out,"stage76_private_capture_crosswalk.csv")
fit_file <- file.path(origdir,"fits/sei_inferred_accuracy_.rds")
if(!all(file.exists(c(source_file,candidate_file,fit_file))))
  stop("Run Stages 74 and 76 first; source or crosswalk missing")
qindex <- function(date) {
  d <- as.Date(date)
  (as.integer(format(d,"%Y"))-1976L)*4L+
    (as.integer(format(d,"%m"))-1L)%/%3L+1L
}
normal <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z[is.na(x)|!nzchar(z)] <- NA_character_
  z
}
old <- utils::read.csv(source_file,check.names=FALSE,
                       colClasses=c(tattoo="character"))
candidates <- utils::read.csv(candidate_file,
                              colClasses=c(capture_id="character",
                                           tattoo="character"))
if(!all(c("tattoo","time") %in% names(old))||
   !all(c("tattoo","time","capture_id") %in% names(candidates)))
  stop("Original source or crosswalk schema changed")
old$tattoo <- normal(old$tattoo)
candidates$tattoo <- normal(candidates$tattoo)
old$time <- as.integer(old$time)
candidates$time <- as.integer(candidates$time)
if(nrow(old)!=14944L ||
   anyNA(old[,c("tattoo","time")])||
   anyNA(candidates[,c("tattoo","time")])||
   any(old$time<1L|old$time>200L)||
   any(candidates$time<1L|candidates$time>200L))
  stop("Source/candidate count or index gate failed")
k <- function(a,b) paste(a,as.integer(b),sep="::")
old_keys <- k(old$tattoo,old$time)
candidate_keys <- k(candidates$tattoo,candidates$time)
if(anyDuplicated(old_keys))
  stop("Original SEI source should contain one row per badger-quarter")
if(!all(candidate_keys %in% old_keys))
  stop("Candidate includes quarters missing from original source; Stage77 assumption changed")
missing <- old[!old_keys%in%candidate_keys,
               c("tattoo","time"),drop=FALSE]
if(!nrow(missing)) {
  cat("Nothing to reconcile: all original quarters in candidate.\n")
  quit(save="no",status=0)
}
cat("\n=== STAGE 78: ORIGINAL SEI-ONLY OBSERVATION AUDIT ===\n")
cat("Original-only badger-quarter occasions:",nrow(missing),"\n")
cat("Strict-live candidate records:",nrow(candidates),"\n")

# Report original SEI fit settings relevant to time/capture construction.
f <- tryCatch(readRDS(fit_file),error=function(e)e)
cat("\n--- ORIGINAL SEI MODEL SETTINGS ---\n")
if(inherits(f,"error")) {
  cat("Could not read fit settings:",conditionMessage(f),"\n")
} else {
  settings <- f$settings
  keys <- c("model","seasons","start_year","repeat_captures",
            "year_process","year_weight","accuracy_prior",
            "method","draws","warmup","traj_draws")
  for(key in keys) {
    val <- settings[[key]]
    if(!is.null(val)&&length(val)<=5L)
      cat(key,":",paste(val,collapse=" | "),"\n")
  }
}
rm(f)

if(!nzchar(Sys.getenv("BADGER_DB_USER"))||
   !nzchar(Sys.getenv("BADGER_DB_PASSWORD")))
  stop("BADGER_DB_USER and BADGER_DB_PASSWORD must be defined")
con <- DBI::dbConnect(
  RPostgres::Postgres(),
  host=Sys.getenv("BADGER_DB_HOST",
                  "db.klpvihubzlyzubyyubtw.supabase.co"),
  port=as.integer(Sys.getenv("BADGER_DB_PORT","5432")),
  dbname="postgres",
  user=Sys.getenv("BADGER_DB_USER"),
  password=Sys.getenv("BADGER_DB_PASSWORD"),sslmode="require")
sql <- "
SELECT capture_id::text AS capture_id,
       upper(btrim(tattoo)) AS tattoo,
       capture_date, moc, coalesce(is_pm,false) AS is_pm
FROM analysis.capture_history_v1
WHERE capture_date >= DATE '1976-01-01'
  AND capture_date < DATE '2026-01-01'
  AND tattoo IS NOT NULL
"
allcap <- DBI::dbGetQuery(con,sql)
DBI::dbDisconnect(con)
allcap$tattoo <- normal(allcap$tattoo)
allcap$time <- qindex(allcap$capture_date)
allcap$moc <- normal(allcap$moc)
if(anyNA(allcap[,c("capture_id","tattoo","time","is_pm")])||
   anyDuplicated(allcap$capture_id))
  stop("Supabase capture provenance/identity gate failed")
allcap$method_confirmed_live <- !allcap$is_pm &
  !is.na(allcap$moc) &
  allcap$moc %in% c("TRAP","NET","DART GUN")
allcap$key <- k(allcap$tattoo,allcap$time)
exact_indices <- split(seq_len(nrow(allcap)),allcap$key)
individual_indices <- split(seq_len(nrow(allcap)),allcap$tattoo)
period_labels <- c("1976-1985","1986-1995","1996-2005",
                   "2006-2015","2016-2025")
audit <- vector("list",nrow(missing))
for(i in seq_len(nrow(missing))) {
  tattoo <- missing$tattoo[i]
  q <- missing$time[i]
  idx <- exact_indices[[k(tattoo,q)]]
  if(is.null(idx))idx <- integer()
  dat <- allcap[idx,,drop=FALSE]
  catg <- if(nrow(dat)==0L) {
    "NO_EXACT_QUARTER_DATABASE_ENCOUNTER"
  } else if(any(dat$method_confirmed_live)) {
    "UNEXPECTED_STRICT_LIVE_MATCH"
  } else if(any(!dat$is_pm)) {
    "NON_PM_ENCOUNTER_OTHER_OR_UNKNOWN_METHOD"
  } else {
    "PM_OR_DEAD_ONLY_EXACT_QUARTER"
  }
  all_i <- individual_indices[[tattoo]]
  if(is.null(all_i))all_i <- integer()
  nearest <- if(length(all_i)) {
    min(abs(allcap$time[all_i]-q))
  } else {
    NA_integer_
  }
  recent <- if(length(all_i)) {
    sum(abs(allcap$time[all_i]-q)==1L)
  } else 0L
  audit[[i]] <- data.frame(
    tattoo=tattoo, time=q,
    year=1976L+(q-1L)%/%4L,
    quarter=(q-1L)%%4L+1L,
    era=period_labels[pmin(5L,(q-1L)%/%40L+1L)],
    category=catg,
    exact_records=nrow(dat),
    exact_pm_records=sum(dat$is_pm),
    exact_non_pm_records=sum(!dat$is_pm),
    exact_moc=paste(sort(unique(c(ifelse(is.na(dat$moc),
                     "<NULL>",dat$moc)))),collapse=";"),
    exact_capture_ids=paste(dat$capture_id,collapse=";"),
    nearest_database_encounter_quarters=nearest,
    neighbouring_quarter_records=recent,
    any_database_encounter_for_tattoo=length(all_i)>0L
  )
}
audit <- do.call(rbind,audit)
if(nrow(audit)!=nrow(missing))
  stop("Audit lost unmatched observations")
write.csv(audit,file.path(out,
  "stage78_PRIVATE_original_only_occasion_reconciliation.csv"),
  row.names=FALSE)
aggregate <- as.data.frame(table(
  era=factor(audit$era,levels=period_labels),
  category=audit$category),stringsAsFactors=FALSE)
names(aggregate)[3] <- "n"
aggregate <- aggregate[aggregate$n>0,,drop=FALSE]
cat("\n--- 317-QUARTER CLASSIFICATION BY ERA ---\n")
print(aggregate,row.names=FALSE)
write.csv(aggregate,file.path(out,
  "stage78_original_only_by_era_category.csv"),row.names=FALSE)

cat("\n--- OVERALL SOURCE-ONLY CLASSIFICATION ---\n")
print(as.data.frame(table(audit$category)),row.names=FALSE)
cat("\nQuarters without matching encounter but with +/-1-quarter encounter:",
    sum(audit$exact_records==0L &
        audit$neighbouring_quarter_records>0L),"\n")
cat("Original-only occasions whose tattoo is absent from ALL database captures:",
    sum(!audit$any_database_encounter_for_tattoo),"\n")
cat("Original-only occasions with no same or adjacent-quarter capture:",
    sum(audit$exact_records==0L &
        audit$neighbouring_quarter_records==0L),"\n")

# Whole-cohort checks to distinguish 1,284 additional same-quarter captures.
candidate_tab <- table(candidate_keys)
cat("\n--- REPEAT-QUARTER RECORDS IN CANDIDATE ---\n")
cat("Candidate rows:",length(candidate_keys),
    "distinct animal-quarters:",length(candidate_tab),"\n")
cat("Extra repeat records:",sum(as.integer(candidate_tab)-1L),"\n")
if(any(audit$category=="UNEXPECTED_STRICT_LIVE_MATCH"))
  stop("Unexpected original-only row matches strict-live method. Review join.")
cat("\nNo dates, capture methods, outcomes or diagnostic results",
    "were modified. This is NOT a model input approval.\n")
cat("STAGE 78 COMPLETE\n")