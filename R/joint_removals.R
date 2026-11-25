# Removal flags for the joint model (v4+): badgers whose accepted death
# record is EUTHANISED or KILLED are censored at that quarter, not counted
# as deaths. Source: outputs/model_audit/integrated_accepted_death_sources.csv
# (written by scripts/50_build_integrated_gm_inputs.R), matched on id.
joint_removed <- function(inp,
    file="outputs/model_audit/integrated_accepted_death_sources.csv") {
  if(!file.exists(file)) stop("Missing death-cause file: ",file,call.=FALSE)
  d <- utils::read.csv(file,stringsAsFactors=FALSE)
  nm <- function(z) toupper(trimws(as.character(z)))
  if(anyDuplicated(nm(d$id))) stop("Duplicate ids in ",file,call.=FALSE)
  cause <- nm(d$cause_class)[match(nm(inp$id),nm(d$id))]
  dead <- !is.na(inp$death_q) & inp$death_q>0
  out <- dead & cause%in%c("EUTHANISED","KILLED")
  cat("Recorded deaths:",sum(dead),"; with cause:",sum(dead&!is.na(cause)),
      "; censored as removals (EUTHANISED/KILLED):",sum(out),"\n")
  print(table(cause[dead],useNA="ifany"))
  out
}
