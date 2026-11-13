# ============================================================
# Stage 53: Replace study-wide noncapture proxy with a
# social-group-specific MONITORING proxy.
# This is NOT measured trap effort: distinct trapping records
# are a lower bound on real deployment/effort.
# No source records or original integrated input overwritten.
# ============================================================
options(stringsAsFactors = FALSE)
stopifnot(requireNamespace("DBI", quietly=TRUE),
          requireNamespace("RPostgres", quietly=TRUE))
x <- readRDS("data/model_inputs/integrated_gm_inputs_2026.rds")
if (x$T != 200L) stop("Expected 1976-2025 inclusive")
normalise <- function(s) {
  s <- toupper(trimws(as.character(s)))
  s[is.na(s) | s=="NA" | s==""] <- NA_character_
  s
}
qidx <- function(d) {
  d <- as.Date(d)
  (as.integer(format(d,"%Y"))-1976L)*4L +
    (as.integer(format(d,"%m"))-1L)%/%3L+1L
}
con <- DBI::dbConnect(
  RPostgres::Postgres(),
  host=Sys.getenv("BADGER_DB_HOST","db.klpvihubzlyzubyyubtw.supabase.co"),
  port=as.integer(Sys.getenv("BADGER_DB_PORT","5432")),
  dbname="postgres",
  user=Sys.getenv("BADGER_DB_USER"),
  password=Sys.getenv("BADGER_DB_PASSWORD"),
  sslmode="require"
)
on.exit(try(DBI::dbDisconnect(con),silent=TRUE),add=TRUE)
e <- DBI::dbGetQuery(con,"
SELECT capture_id::text AS capture_id,
       upper(trim(tattoo)) AS id, capture_date, moc, is_pm,
       normalized_social_group_name AS sg
FROM analysis.capture_history_v1
WHERE capture_date >= DATE '1976-01-01'
  AND capture_date < DATE '2026-01-01'
ORDER BY capture_id
")
DBI::dbDisconnect(con)
if (anyDuplicated(e$capture_id)) stop("Duplicate source capture_id")
e$id <- normalise(e$id)
e$sg <- normalise(e$sg)
e$q <- qidx(e$capture_date)
e$moc <- normalise(e$moc)
is_live <- e$moc %in% c("TRAP","NET","DART GUN") &
  !is.na(e$moc) & (is.na(e$is_pm)|!e$is_pm)
is_live[is.na(is_live)] <- FALSE
e <- e[is_live & e$q >= 1L & e$q <= x$T, , drop=FALSE]

# Source capture in a group establishes some study activity there.
# This is NOT true effort, and cannot prove zero trapping elsewhere.
activity <- unique(e[
  !is.na(e$sg), c("sg","q"), drop=FALSE
])
activity_key <- paste(activity$sg, activity$q, sep="@@")

# Only overwrite the 0s currently inferred from global proxy.
# Keep confirmed detections and explicitly missing occasions unchanged.
y_global <- x$y
y_group <- y_global
y_group[y_group == 0L] <- -1L
location_known <- logical(x$N)
last_sg <- rep(NA_character_, x$N)
group_switched <- logical(x$N)

# An animal's last observed SG is carried forward ONLY as a
# provisional group reference for sensitivity comparisons.
# Movement/temporary emigration still require a model component.
for(i in seq_len(x$N)) {
  obs <- e[e$id == x$id[i] & !is.na(e$sg),c("q","sg")]
  obs <- obs[order(obs$q), ,drop=FALSE]
  if (!nrow(obs)) next
  location_known[i] <- TRUE
  group_switched[i] <- length(unique(obs$sg)) > 1L
  last_sg[i] <- obs$sg[nrow(obs)]
  qobs <- unique(obs$q)
  for(q in seq.int(x$first_q[i],x$T)) {
    k <- which(obs$q <= q)
    if (!length(k)) next
    sg <- obs$sg[max(k)]
    active <- paste(sg,q,sep="@@") %in% activity_key
    if(active && y_global[i,q] == 0L) y_group[i,q] <- 0L
  }
}
if (!identical(which(y_global==1L),which(y_group==1L)))
  stop("Confirmed-live observation altered")
if (any(y_group[y_global == -1L] != -1L))
  stop("Missing/no-observation cells became noncaptures")
if (any(y_group == 0L & y_global != 0L))
  stop("Invented noncaptures")

x$y <- y_group
x$data_provenance$effort <- paste(
  "SG quarter activity from any live trapping in same group",
  "with last observed group carried forward: PROXY ONLY"
)
x$data_provenance$original_global_proxy_file <-
  "data/model_inputs/integrated_gm_inputs_2026.rds"
x$data_provenance$group_switch_warning <-
  "Movement/dispersal makes last observed SG potentially incorrect"
saveRDS(x,"data/model_inputs/integrated_gm_inputs_group_proxy_2026.rds")
individual <- data.frame(
  id=x$id,
  n_global_noncaptures=rowSums(y_global == 0L),
  n_sg_proxy_noncaptures=rowSums(y_group == 0L),
  had_sg_location=location_known,
  group_switched_observed=group_switched,
  last_observed_sg=last_sg
)
write.csv(individual,
  "outputs/model_audit/group_proxy_detection_comparison.csv",
  row.names=FALSE,na=""
)
cat("\n=== GROUP-SPECIFIC MONITORING PROXY ===\n")
cat("Individuals:",x$N,"\n")
cat("With recorded social-group name:",sum(location_known),"\n")
cat("With observed multiple groups:",sum(group_switched),"\n")
cat("Global proxy noncaptures:",sum(y_global==0L),"\n")
cat("SG-specific proxy noncaptures:",sum(y_group==0L),"\n")
cat("Noncaptures removed as unverified:",
    sum(y_global==0L)-sum(y_group==0L),"\n")
cat("Confirmed-live badger-quarters unchanged:",sum(y_group==1L),"\n")
cat("Deaths unchanged:",sum(x$dead==1L),"\n")
cat("\nThis proxy remains insufficient for estimating true individual effort.\n")
cat("Keep BOTH versions for model sensitivity to monitoring assumptions.\n")
