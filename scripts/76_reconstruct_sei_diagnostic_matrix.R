# ============================================================================
# Stage 76: RECONSTRUCT / AUDIT a candidate BadgerInfeR SEI test input matrix.
#
# Source:
#   local original SEI cohort inventory (3214 tattoos / b$n_captures)
#   Supabase analysis.capture_history_v1
#   Supabase analysis.diagnostic_results
# Read-only SELECT. No database changes or infection model fits.
#
# This is a CANDIDATE reconstruction, NOT necessarily the exact original
# test matrix. A missing/ambiguous test -> NA, never an invented negative.
# Mixed culture plates are treated as positive if any M.BOVIS is present,
# and reported explicitly for source-protocol validation.
#
# Original HMM format: time,id,captured,test1..test6
#   1 BELT ELISA; 2 INDIRECT ELISA; 3 Culture; 4 DPP;
#   5 IGRA (IFN-gamma); 6 StatPak
#
# Outputs: outputs/sei_transitions/stage76_*
# WARNING: local outputs include research diagnostic histories; keep
# private and NEVER git add these CSVs.
# ============================================================================
options(stringsAsFactors=FALSE)
needed <- c("DBI","RPostgres")
if(!all(vapply(needed,requireNamespace,logical(1),quietly=TRUE)))
  stop("Need R packages: DBI and RPostgres")
src <- "data/infection_posteriors/sei_2026"
src_file <- file.path(src,"badgers.csv")
stopifnot(file.exists(src_file))
# Source is an extracted research bundle; verify that the inventory has
# not changed since Stage 65 before linking a different database source.
fingerprint_path <- file.path(src,"sei_extracted_files_sha256.tsv")
if(!file.exists(fingerprint_path)||
   !nzchar(Sys.which("sha256sum")))
  stop("Stage 65 SHA256 manifest or sha256sum is unavailable")
fingerprints <- utils::read.delim(fingerprint_path,stringsAsFactors=FALSE)
if(!all(c("filename","sha256")%in%names(fingerprints))||
   anyDuplicated(fingerprints$filename)||
   !"badgers.csv"%in%fingerprints$filename)
  stop("Original SEI fingerprint manifest is invalid")
expected <- fingerprints$sha256[
  match("badgers.csv",fingerprints$filename)]
digest_line <- system2("sha256sum",shQuote(src_file),
                       stdout=TRUE,stderr=FALSE)
if(length(digest_line)!=1L||
   !identical(strsplit(digest_line," ",fixed=TRUE)[[1]][1],
              expected))
  stop("SEI badgers.csv SHA256 mismatch: do not combine datasets")
out <- "outputs/sei_transitions"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
normal <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z[is.na(x)|z==""] <- NA_character_
  z
}
q_index <- function(x){
  d <- as.Date(x)
  (as.integer(format(d,"%Y"))-1976L)*4L+
    (as.integer(format(d,"%m"))-1L)%/%3L+1L
}
fail <- function(...)stop(paste0(...),call.=FALSE)

orig <- utils::read.csv(src_file,colClasses=c(tattoo="character"),
                        check.names=FALSE)
required <- c("tattoo","first_capture_time","last_capture_time",
              "start_time","n_captures")
if(!all(required %in% names(orig)))fail("Missing cohort inventory fields")
orig$tattoo <- normal(orig$tattoo)
if(nrow(orig)!=3214L||anyNA(orig$tattoo)||
   anyDuplicated(orig$tattoo))
  fail("Original SEI 3214-tattoo inventory mismatch")
for(k in setdiff(required,"tattoo"))
  orig[[k]] <- as.integer(orig[[k]])
if(anyNA(orig[,required])||
   any(orig$first_capture_time<1L)||
   any(orig$first_capture_time>orig$last_capture_time)||
   any(orig$last_capture_time>200L)||
   any(orig$n_captures<1L))fail("Invalid original cohort metadata")
orig$id <- seq_len(nrow(orig))
write.csv(orig[,c("id","tattoo","first_capture_time",
                  "last_capture_time","start_time","n_captures")],
          file.path(out,"stage76_private_id_crosswalk.csv"),
          row.names=FALSE)
cat("\n=== STAGE 76: READ-ONLY DIAGNOSTIC INPUT RECONCILIATION ===\n")
cat("Original cohort:",nrow(orig),
    "animals; declared source captures:",sum(orig$n_captures),"\n")
cat("HMM test columns: BELT, INDIRECT, CULTURE, DPP, IGRA, STATPAK\n")

if(!nzchar(Sys.getenv("BADGER_DB_USER"))||
   !nzchar(Sys.getenv("BADGER_DB_PASSWORD")))
  fail("Set BADGER_DB_USER and BADGER_DB_PASSWORD (same as Stage 50).")
con <- DBI::dbConnect(
  RPostgres::Postgres(),
  host=Sys.getenv("BADGER_DB_HOST",
                  "db.klpvihubzlyzubyyubtw.supabase.co"),
  port=as.integer(Sys.getenv("BADGER_DB_PORT","5432")),
  dbname="postgres",user=Sys.getenv("BADGER_DB_USER"),
  password=Sys.getenv("BADGER_DB_PASSWORD"),sslmode="require")
on.exit(try(DBI::dbDisconnect(con),silent=TRUE),add=TRUE)
# Plain SELECTs only. Do not pull location, comments, genetic data or PMs.
sql_captures <- "
SELECT capture_id::text AS capture_id,
       upper(btrim(tattoo)) AS tattoo,
       capture_date,moc,coalesce(is_pm,false) AS is_pm
FROM analysis.capture_history_v1
WHERE capture_date >= DATE '1976-01-01'
  AND capture_date < DATE '2026-01-01'
  AND tattoo IS NOT NULL
"
sql_tests <- "
SELECT source_table,source_record_id::text AS source_record_id,
       capture_id::text AS capture_id,upper(btrim(tattoo)) AS tattoo,
       capture_date,is_pm,test_type,result,result_code,
       visual_line_1_result,visual_line_2_result
FROM analysis.diagnostic_results
WHERE capture_date >= DATE '1976-01-01'
  AND capture_date < DATE '2026-01-01'
"
cat("Reading capture and six-assay diagnostic views...\n")
captures <- DBI::dbGetQuery(con,sql_captures)
tests <- DBI::dbGetQuery(con,sql_tests)
DBI::dbDisconnect(con)
if(anyDuplicated(captures$capture_id)||
   anyNA(captures$capture_id))fail("Duplicate/null capture IDs")
captures$tattoo <- normal(captures$tattoo)
tests$tattoo <- normal(tests$tattoo)
captures$q <- q_index(captures$capture_date)
tests$q <- q_index(tests$capture_date)
if(anyNA(captures$q)||anyNA(tests$q))fail("Unparseable capture/test dates")
if(anyDuplicated(paste(tests$source_table,tests$source_record_id,
                       sep="|")))
  fail("Duplicate source-table/source-record ID in diagnostic view")

# Prevent phantom tests from joining to the wrong animal or date.
ix <- match(tests$capture_id,captures$capture_id)
bad_link <- is.na(ix) |
  (tests$tattoo != captures$tattoo[ix]) |
  (tests$q != captures$q[ix])
bad_link[is.na(bad_link)] <- TRUE
cat("Diagnostics with missing/inconsistent capture linkage:",
    sum(bad_link),"of",nrow(tests),"\n")
# Some source PM records may be validly absent from the live capture set;
# audit and exclude but do not silently recode their outcomes.
if(any(bad_link)){
  dangling <- tests[bad_link,c("source_table","source_record_id",
                               "capture_id","tattoo","capture_date",
                               "test_type"),drop=FALSE]
  write.csv(dangling,file.path(out,
    "stage76_private_unmatched_diagnostic_sources.csv"),row.names=FALSE)
}
tests <- tests[!bad_link,,drop=FALSE]

method <- normal(captures$moc)
live <- !captures$is_pm & method %in% c("TRAP","NET","DART GUN")
live[is.na(live)] <- FALSE
captures <- captures[live,,drop=FALSE]
captures$id <- match(captures$tattoo,orig$tattoo)
cat("Database method-confirmed live captures:",nrow(captures),
    "; among original SEI tattoos:",sum(!is.na(captures$id)),"\n")
outside <- captures[is.na(captures$id),,drop=FALSE]
if(nrow(outside))
  cat("Live captures excluded for no original SEI ID:",nrow(outside),"\n")
captures <- captures[!is.na(captures$id),,drop=FALSE]
captures$id <- as.integer(captures$id)
captures$orig_first <- orig$first_capture_time[captures$id]
captures$orig_last <- orig$last_capture_time[captures$id]
in_window <- captures$q>=captures$orig_first &
             captures$q<=captures$orig_last
cat("Matched live captures outside original SEI first/last quarter:",
    sum(!in_window),"\n")
write.csv(captures[!in_window,
  c("capture_id","tattoo","q","orig_first","orig_last"),drop=FALSE],
  file.path(out,"stage76_private_outside_sei_windows.csv"),
  row.names=FALSE)
captures <- captures[in_window,,drop=FALSE]
captures <- captures[order(captures$id,captures$q,
                           captures$capture_id),,drop=FALSE]
stopifnot(!anyDuplicated(captures$capture_id))
in_selected <- tests$capture_id %in% captures$capture_id
tests <- tests[in_selected,,drop=FALSE]
cat("Diagnostic source rows attached to matched captures:",nrow(tests),"\n")

# Map observation labels, preserving unknown codes as NA. Do not collapse
# 'I', 'BIT', 'ITT', 'NO SERA', pending, contaminated, 'X' flags to 0.
assay <- function(source,type){
  src <- normal(source); typ <- normal(type)
  if(src=="ELISA_TESTS" & typ=="BELT") return(1L)
  if(src=="ELISA_TESTS" & typ=="INDIRECT") return(2L)
  if(src=="CULTURE") return(3L)
  if(src=="DPP_TESTS" & typ %in% c("DOUBLE","SINGLE"))return(4L)
  if(src=="IFN_GAMMA_ELISA")return(5L)
  if(src=="STAT_PAK_TESTS")return(6L)
  if(src=="BROCK_TESTS")return(0L)
  NA_integer_
}
test_index <- mapply(assay,tests$source_table,tests$test_type,
                     USE.NAMES=FALSE)
if(anyNA(test_index)){
  write.csv(unique(tests[is.na(test_index),
    c("source_table","test_type"),drop=FALSE]),
    file.path(out,"stage76_unrecognised_assays.csv"),row.names=FALSE)
  fail("Unknown diagnostic assay labels: inspect stage76_unrecognised_assays.csv")
}
cat("Brock rows not among the 6 fitted assays:",
    sum(test_index==0L),"\n")
tests <- tests[test_index!=0L,,drop=FALSE]
test_index <- test_index[test_index!=0L]
normval <- function(x){
  v <- normal(x)
  v
}
map_test <- function(j,r,rc,v1,v2,type){
  r <- normval(r);rc <- normval(rc)
  v1 <- normval(v1);v2 <- normval(v2)
  type <- normval(type)
  if(j %in% c(1L,2L)){
    if(!is.na(r) && r %in% c("P","POSITIVE")) return(1L)
    if(!is.na(r) && r %in% c("N","NEGATIVE")) return(0L)
    return(NA_integer_)
  }
  if(j==3L){
    if(!is.na(r)&&r=="M.BOVIS")return(1L)
    if(!is.na(r)&&r=="NEGATIVE")return(0L)
    return(NA_integer_)
  }
  if(j==4L){
    valid <- c("POSITIVE","NEGATIVE")
    if(is.na(v1)||!v1%in%valid)return(NA_integer_)
    if(type=="DOUBLE"){
      if(is.na(v2)||!v2%in%valid)return(NA_integer_)
      return(as.integer(v1=="POSITIVE"||v2=="POSITIVE"))
    }
    if(type=="SINGLE")return(as.integer(v1=="POSITIVE"))
    return(NA_integer_)
  }
  if(j==5L){
    if(!is.na(r)&&r=="POSITIVE")return(1L)
    if(!is.na(r)&&r=="NEGATIVE")return(0L)
    return(NA_integer_)
  }
  if(j==6L){
    if(!is.na(rc)&&rc=="P")return(1L)
    if(!is.na(rc)&&rc=="N")return(0L)
    return(NA_integer_)
  }
  NA_integer_
}
# Small internal regression tests to guard against accidental relabelling.
stopifnot(
  identical(map_test(1L,"P",NA,NA,NA,"BELT"),1L),
  identical(map_test(2L,"I",NA,NA,NA,"INDIRECT"),NA_integer_),
  identical(map_test(3L,"FAILED",NA,NA,NA,"Culture"),NA_integer_),
  identical(map_test(3L,"M.BOVIS",NA,NA,NA,"Culture"),1L),
  identical(map_test(4L,NA,NA,"Positive","Negative","Double"),1L),
  identical(map_test(4L,NA,NA,"NegativeX","NegativeX","Double"),NA_integer_),
  identical(map_test(5L,"Negative",NA,NA,NA,"IGRA"),0L),
  identical(map_test(6L,"No Clotted Sample","C",NA,NA,"StatPak"),
            NA_integer_)
)
tests$assay <- test_index
tests$binary <- mapply(map_test,tests$assay,tests$result,
                       tests$result_code,tests$visual_line_1_result,
                       tests$visual_line_2_result,tests$test_type,
                       USE.NAMES=FALSE)
# Includes repeat assays: Culture may have positive/negative plates
# at one capture. A capture-level positive if ANY M.BOVIS; otherwise
# only all valid NEGATIVE -> negative. Unknown culture results are
# excluded unless a positive was also observed.
# For non-culture: conflicting positive/negative same capture -> NA.
split_key <- paste(tests$capture_id,tests$assay,sep="|")
groups <- split(seq_len(nrow(tests)),split_key)
collapsed <- lapply(groups,function(ii){
  v <- tests$binary[ii]
  assay_id <- tests$assay[ii[1]]
  valid <- v[!is.na(v)]
  has_pos <- any(valid==1L)
  has_neg <- any(valid==0L)
  has_unreadable <- anyNA(v)
  conflict <- has_pos && has_neg
  value <- if(assay_id==3L) {
    if(has_pos) 1L
    else if(has_neg&&!has_unreadable) 0L
    else NA_integer_
  } else {
    if(conflict || !length(valid)) NA_integer_
    else if(length(unique(valid))==1L) valid[1L]
    else NA_integer_
  }
  data.frame(capture_id=tests$capture_id[ii[1]],
             assay=assay_id,value=value,
             source_rows=length(ii),unreadable=sum(is.na(v)),
             conflict=conflict,
             contains_unreadable=has_unreadable)
})
collapsed <- do.call(rbind,collapsed)
if(is.null(collapsed))fail("No assay rows after filtering")
cat("Non-culture assay contradictions (assigned NA):",
    sum(collapsed$conflict & collapsed$assay!=3L),"\n")
cat("Mixed positive/negative culture plates (positive-rule candidates):",
    sum(collapsed$conflict & collapsed$assay==3L),"\n")
write.csv(collapsed[collapsed$conflict |
                     collapsed$contains_unreadable,
                     c("capture_id","assay","source_rows",
                       "unreadable","conflict","value")],
  file.path(out,"stage76_private_assay_ambiguities.csv"),row.names=FALSE)

# Numeric IDs correspond to original badgers.csv order, not legacy
# tattoo numeric substrings. Preserve private crosswalk above.
candidate <- data.frame(time=captures$q,id=captures$id,captured=1L,
                        test1=NA_integer_,test2=NA_integer_,
                        test3=NA_integer_,test4=NA_integer_,
                        test5=NA_integer_,test6=NA_integer_)
pos <- match(collapsed$capture_id,captures$capture_id)
if(anyNA(pos))fail("Collapsed tests not found in capture inventory")
for(k in 1:6){
  sel <- which(collapsed$assay==k)
  candidate[pos[sel],paste0("test",k)] <- collapsed$value[sel]
}
if(anyDuplicated(paste(captures$capture_id,candidate$id))||
   any(candidate$time<1|candidate$time>200)||
   any(candidate$captured!=1L))
  fail("Candidate test input validation failed")
name <- file.path(out,"stage76_PRIVATE_CANDIDATE_NOT_VALIDATED.csv")
write.csv(candidate,name,row.names=FALSE,na="",quote=FALSE)
write.csv(data.frame(capture_id=captures$capture_id,
                     tattoo=captures$tattoo,time=captures$q,
                     numeric_id=candidate$id),
          file.path(out,"stage76_private_capture_crosswalk.csv"),
          row.names=FALSE)

# Reconcile source per-individual capture counts and observation windows.
count <- tabulate(candidate$id,nbins=nrow(orig))
first <- rep(NA_integer_,nrow(orig))
final <- rep(NA_integer_,nrow(orig))
for(i in which(count>0)){
  qq <- candidate$time[candidate$id==i]
  first[i] <- min(qq)
  final[i] <- max(qq)
}
rec <- data.frame(
  id=orig$id,tattoo=orig$tattoo,
  source_n_captures=orig$n_captures,
  candidate_n_captures=count,
  source_first=orig$first_capture_time,candidate_first=first,
  source_last=orig$last_capture_time,candidate_last=final,
  match_count=count==orig$n_captures,
  match_first=!is.na(first)&first==orig$first_capture_time,
  match_last=!is.na(final)&final==orig$last_capture_time)
rec$pass <- rec$match_count&rec$match_first&rec$match_last
write.csv(rec,file.path(out,
  "stage76_private_original_capture_reconciliation.csv"),row.names=FALSE)
cat("\n--- RECONCILIATION ---\n")
cat("Original 3214-badger declared captures:",sum(orig$n_captures),"\n")
cat("Candidate captures:",nrow(candidate),"\n")
cat("Matched original animals with >=1 capture:",sum(count>0),"\n")
cat("Individuals matching ALL count/first/last:",
    sum(rec$pass),"/",nrow(orig),"\n")
cat("Count disagreements:",sum(!rec$match_count),
    "first-quarter disagreements:",sum(!rec$match_first),
    "last-quarter disagreements:",sum(!rec$match_last),"\n")
cat("Assay coverage among candidate captures:\n")
for(k in 1:6){
  val <- candidate[[paste0("test",k)]]
  cat(sprintf("  test%d: tested=%d, positive=%d, negative=%d, NA=%d\n",
              k,sum(!is.na(val)),sum(val==1L,na.rm=TRUE),
              sum(val==0L,na.rm=TRUE),sum(is.na(val))))
}
# By-era test coverage (descriptive; no individual raw records printed).
era <- pmin(5L,pmax(1L,(candidate$time-1L)%/%40L+1L))
labels <- c("1976-85","1986-95","1996-05","2006-15","2016-25")
coverage <- do.call(rbind,lapply(1:5,function(i){
  ix <- which(era==i)
  data.frame(era=labels[i],captures=length(ix),
    t1=sum(!is.na(candidate$test1[ix])),
    t2=sum(!is.na(candidate$test2[ix])),
    t3=sum(!is.na(candidate$test3[ix])),
    t4=sum(!is.na(candidate$test4[ix])),
    t5=sum(!is.na(candidate$test5[ix])),
    t6=sum(!is.na(candidate$test6[ix])))
}))
write.csv(coverage,file.path(out,"stage76_test_coverage_by_era.csv"),
          row.names=FALSE)
print(coverage,row.names=FALSE)
cat("\nFINAL GATE: ")
if(all(rec$pass) &&
   !any(collapsed$conflict & collapsed$assay!=3L)){
  cat("CAPTURE COUNTS MATCH, but original binary test matrix still",
      "NOT independently verified. Do not claim exact reproduction.\n")
} else {
  cat("FAILED exact-input reconciliation: diagnostic candidate ONLY.\n")
}
cat("Retain source test matrix and analysis protocol for manual",
    "record-level reconciliation before a joint model refit.\n")
cat("Saved candidate:",name,"\n")
cat("NO DATABASE RECORDS WERE MODIFIED. NO SEI MODEL WAS FITTED.\n")
cat("STAGE 76 COMPLETE\n")
