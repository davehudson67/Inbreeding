# ============================================================================
# STAGE 77: Original SEI per-capture vs candidate capture-quarter reconciliation
#
# No database calls. No fitting. Reads outputs from Stage 76 plus
# original source p_infected_per_capture.csv and badgers.csv.
#
# Important: p_infected_per_capture.csv is NOT the raw diagnostic matrix.
# Its tattoo/time pairs allow a much stronger audit of WHICH capture
# occasions the original fit used, but not how the six assays were coded.
#
# Privacy: all tattoo-level comparisons stay in outputs/sei_transitions/
# (gitignored). Terminal output is aggregated only.
# ============================================================================
options(stringsAsFactors=FALSE)
dir <- "outputs/sei_transitions"
orig <- "data/infection_posteriors/sei_2026"
files <- c(
  original_events=file.path(orig,"p_infected_per_capture.csv"),
  original_badgers=file.path(orig,"badgers.csv"),
  candidate=file.path(dir,"stage76_PRIVATE_CANDIDATE_NOT_VALIDATED.csv"),
  crosswalk=file.path(dir,"stage76_private_capture_crosswalk.csv"),
  id_crosswalk=file.path(dir,"stage76_private_id_crosswalk.csv"),
  ambiguity=file.path(dir,"stage76_private_assay_ambiguities.csv"))
if(!all(file.exists(files)))
  stop("Missing original bundle or Stage 76 files: ",
       paste(names(files)[!file.exists(files)],collapse=", "))
fail <- function(...)stop(paste0(...),call.=FALSE)
normal <- function(x) {
  toupper(trimws(as.character(x)))
}
read <- function(path)utils::read.csv(path,check.names=FALSE,
                                      stringsAsFactors=FALSE)
old <- read(files["original_events"])
b <- read(files["original_badgers"])
candidate <- read(files["candidate"])
cw <- read(files["crosswalk"])
ids <- read(files["id_crosswalk"])
amb <- read(files["ambiguity"])
req <- function(df,columns,label){
  if(!all(columns%in%names(df)))
    fail(label,": missing required columns: ",
         paste(setdiff(columns,names(df)),collapse=", "))
}
req(old,c("tattoo","time","year","quarter"),"original events")
req(b,c("tattoo","n_captures","first_capture_time",
        "last_capture_time"),"original badgers")
req(candidate,c("time","id","captured",paste0("test",1:6)),
    "candidate diagnostic data")
req(cw,c("capture_id","tattoo","time","numeric_id"),
    "Stage 76 capture crosswalk")
req(ids,c("id","tattoo"),"Stage 76 ID crosswalk")
req(amb,c("assay","source_rows","unreadable","conflict"),
    "Stage 76 ambiguity report")
if(nrow(b)!=3214L||nrow(old)!=14944L||
   nrow(ids)!=3214L||nrow(candidate)!=nrow(cw))
  fail("Cohort sizes or Stage76 output lengths changed unexpectedly")
old$tattoo <- normal(old$tattoo)
b$tattoo <- normal(b$tattoo)
cw$tattoo <- normal(cw$tattoo)
ids$tattoo <- normal(ids$tattoo)
old$time <- suppressWarnings(as.integer(old$time))
old$year <- suppressWarnings(as.integer(old$year))
old$quarter <- suppressWarnings(as.integer(old$quarter))
cw$time <- suppressWarnings(as.integer(cw$time))
candidate$time <- suppressWarnings(as.integer(candidate$time))
cw$numeric_id <- suppressWarnings(as.integer(cw$numeric_id))
candidate$id <- suppressWarnings(as.integer(candidate$id))
ids$id <- suppressWarnings(as.integer(ids$id))
b$n_captures <- suppressWarnings(as.integer(b$n_captures))
if(anyNA(old[,c("tattoo","time","year","quarter")])||
   anyNA(cw[,c("tattoo","time","numeric_id")])||
   anyNA(candidate[,c("time","id","captured")])||
   anyDuplicated(ids$tattoo)||anyDuplicated(ids$id)||
   anyDuplicated(cw$capture_id)||
   any(candidate$captured!=1L)||
   any(old$time<1L|old$time>200L)||
   any(candidate$time<1L|candidate$time>200L)||
   !identical(as.integer(cw$numeric_id),as.integer(candidate$id))||
   !identical(as.integer(cw$time),as.integer(candidate$time))||
   any(ids$tattoo[match(candidate$id,ids$id)]!=cw$tattoo))
  fail("Original/candidate integrity or crosswalk validation failed")
if(any(old$year != 1976L+(old$time-1L)%/%4L)||
   any(old$quarter!=(old$time-1L)%%4L+1L))
  fail("Original per-capture time/year/quarter are not aligned: investigate")
if(!all(old$tattoo%in%b$tattoo)||
   !all(cw$tattoo%in%b$tattoo))
  fail("Unexpected tattoos not in original badgers.csv")
if(sum(b$n_captures)!=nrow(old))
  fail("Original source badgers$n_captures does not sum to 14,944")
cat("\n=== STAGE 77: ORIGINAL VS RECONSTRUCTED CAPTURE OCCASIONS ===\n")
cat("Original SEI per-capture rows:",nrow(old),"\n")
cat("Stage76 method-confirmed candidate rows:",nrow(candidate),"\n")
cat("Original cohort:",nrow(b),"badgers\n")

# Aggregate event records by tattoo and QUARTER, retaining repeats separately.
count_group <- function(tattoo,time,what){
  a <- stats::aggregate(rep.int(1L,length(time)),
       by=list(tattoo=tattoo,time=time),FUN=sum)
  names(a)[3] <- what
  a
}
s <- count_group(old$tattoo,old$time,"source_count")
c <- count_group(cw$tattoo,cw$time,"candidate_count")
g <- merge(s,c,by=c("tattoo","time"),all=TRUE,sort=TRUE)
g$source_count[is.na(g$source_count)] <- 0L
g$candidate_count[is.na(g$candidate_count)] <- 0L
g$delta <- g$candidate_count - g$source_count
g$period <- c("1976-1985","1986-1995","1996-2005",
              "2006-2015","2016-2025")[
              pmin(5L,(g$time-1L)%/%40L+1L)]
if(sum(g$source_count)!=nrow(old)||
   sum(g$candidate_count)!=nrow(candidate)||
   anyDuplicated(g[,c("tattoo","time")]))
  fail("Quarter-based reconciliation lost capture records")
summary <- data.frame(
  source_rows=nrow(old),
  candidate_rows=nrow(candidate),
  source_unique_animal_quarters=nrow(s),
  candidate_unique_animal_quarters=nrow(c),
  source_repeated_quarter_rows=sum(pmax(0L,s$source_count-1L)),
  candidate_repeated_quarter_rows=sum(pmax(0L,c$candidate_count-1L)),
  source_only_quarters=sum(g$source_count>0L & g$candidate_count==0L),
  candidate_only_quarters=sum(g$candidate_count>0L & g$source_count==0L),
  shared_quarters_equal_counts=sum(g$source_count==g$candidate_count),
  shared_quarters_unequal_counts=sum(g$source_count>0L &
                                     g$candidate_count>0L &
                                     g$delta!=0L),
  net_candidate_excess=sum(g$delta),
  total_source_records_missing_from_candidate=sum(pmax(0L,-g$delta)),
  total_candidate_records_exceeding_original=sum(pmax(0L,g$delta))
)
print(summary,row.names=FALSE)
write.csv(summary,file.path(dir,"stage77_reconciliation_totals.csv"),
          row.names=FALSE)
write.csv(g[g$delta!=0L,,drop=FALSE],
          file.path(dir,"stage77_PRIVATE_quarter_discrepancies.csv"),
          row.names=FALSE)
cat("Distinct badger-quarters: original=",nrow(s),
    ", candidate=",nrow(c),"\n",sep="")
cat("Badger-quarters with matching count:",
    sum(g$delta==0L),"/",nrow(g),"\n")
cat("Original-only badger-quarters:",
    sum(g$source_count>0L & g$candidate_count==0L),"\n")
cat("Candidate-only badger-quarters:",
    sum(g$candidate_count>0L & g$source_count==0L),"\n")

# Identify individual animals whose captures are missing or differ.
all_tattoo <- sort(unique(b$tattoo))
tab_original <- tabulate(match(old$tattoo,all_tattoo),
                         nbins=length(all_tattoo))
tab_candidate <- tabulate(match(cw$tattoo,all_tattoo),
                          nbins=length(all_tattoo))
cohort <- data.frame(
  tattoo=all_tattoo,
  original_events=tab_original,
  candidate_events=tab_candidate,
  difference=tab_candidate-tab_original,
  original_inventory_n=b$n_captures[match(all_tattoo,b$tattoo)])
stopifnot(all(cohort$original_events==cohort$original_inventory_n))
write.csv(cohort,file.path(dir,
  "stage77_PRIVATE_individual_count_reconciliation.csv"),
  row.names=FALSE)
cat("Original badgers absent from candidate:",
    sum(cohort$candidate_events==0L),"\n")
cat("Badgers with exactly equal number of capture events:",
    sum(cohort$difference==0L),"\n")
cat("Badgers with fewer reconstructed events:",
    sum(cohort$difference<0L),"\n")
cat("Badgers with more reconstructed events:",
    sum(cohort$difference>0L),"\n")
cat("Capture-count delta frequencies (candidate minus original):\n")
dist <- as.data.frame(table(cohort$difference))
names(dist) <- c("delta","n_badgers")
print(dist,row.names=FALSE)
write.csv(dist,file.path(dir,"stage77_count_difference_distribution.csv"),
          row.names=FALSE)

# Era totals separately for original and candidate capture occasions.
periods <- c("1976-1985","1986-1995","1996-2005",
             "2006-2015","2016-2025")
era <- do.call(rbind,lapply(periods,function(k){
  z <- g[g$period==k,,drop=FALSE]
  data.frame(
    period=k,
    original_records=sum(z$source_count),
    candidate_records=sum(z$candidate_count),
    net_excess=sum(z$delta),
    source_only_badger_quarters=sum(z$source_count>0 & z$candidate_count==0),
    candidate_only_badger_quarters=sum(z$candidate_count>0 & z$source_count==0),
    counts_differ_badger_quarters=sum(z$delta!=0)
  )
}))
cat("\n--- CAPTURE DIFFERENCES BY ERA ---\n")
print(era,row.names=FALSE)
write.csv(era,file.path(dir,"stage77_era_capture_comparison.csv"),
          row.names=FALSE)

# Existing assay ambiguity report from Stage76: print aggregate only.
if(nrow(amb)>0L) {
  ambiguity <- do.call(rbind,lapply(1:6,function(k){
    d <- amb[amb$assay==k,,drop=FALSE]
    data.frame(assay=paste0("test",k),
       flagged_capture_assays=nrow(d),
       conflicting_positive_negative=sum(d$conflict),
       contains_unreadable_result=sum(d$unreadable>0),
       source_rows_in_flagged_assays=sum(d$source_rows))
  }))
  cat("\n--- ASSAY AMBIGUITIES (from Stage 76) ---\n")
  print(ambiguity,row.names=FALSE)
  write.csv(ambiguity,file.path(dir,"stage77_assay_ambiguity_summary.csv"),
            row.names=FALSE)
}

cat("\nINTERPRETATION GATE:\n")
cat("p_infected_per_capture.csv establishes original capture-time",
    "membership, NOT six-test observation outcomes.\n")
if(all(g$delta==0L)) {
  cat("ALL ORIGINAL CAPTURE-TIME COUNTS AGREE; still need",
      "raw diagnostic-input verification before fitting.\n")
} else {
  cat("CAPTURE-TIME RECONCILIATION FAILED; prioritize the",
      "private quarter-level disagreements.\n")
}
cat("Do not automatically drop duplicate quarters or force",
    "ambiguous assays to 0/1. No source was modified.\n")
cat("STAGE 77 COMPLETE\n")
