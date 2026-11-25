# ============================================================================
# Read-only review of the three joint v3d chains (homozygosity x infection
# on senescence rate). Adds: b multiplier per SD homozygosity, and the total mortality
# hazard ratio for +1 SD homozygosity at ages 2, 4, 6, 8, 10 (S state, female).
#  1. provenance: exit 0, identical hashes/settings, distinct seeds
#  2. rank-normalised split R-hat, bulk and tail ESS (posterior package)
#  3. MLH and sex effects as ORs (S->E, E->I) and HRs (mortality);
#     E and I vs S mortality HRs; frailty SDs; derived period risks
#  4. posterior check: share of inferred I onsets after each badger's last
#     live capture (forward-filter backward-sample per draw, using z draws)
# Writes outputs/joint_sei_survival/full_joint_v3d_infsen/review_<UTC>/.
# Usage: Rscript scripts/summarise_joint_v3d_chains.R  [NDRAW_Z=100]
# ============================================================================
options(stringsAsFactors=FALSE)
fail <- function(...) stop(paste0(...),call.=FALSE)
root <- file.path("outputs/joint_sei_survival",
  Sys.getenv("RUN_NAME","full_joint_v3d_infsen"))
ids <- 1:3
dirs <- file.path(root,sprintf("chain_%02d",ids))
for(d in dirs){
  st <- file.path(d,"process_exit_status.txt")
  if(!file.exists(file.path(d,"posterior_samples.rds"))||!file.exists(st)||
     !any(readLines(st)=="EXIT_STATUS=0"))
    fail("Chain not finished cleanly: ",d)
}
ch <- lapply(file.path(dirs,"posterior_samples.rds"),readRDS)
same <- function(k) length(unique(vapply(ch,function(z)
  paste(unname(as.character(z[[k]])),collapse=":"),"")))==1L
for(k in c("source_md5","input_md5","model_code_md5","N","iterations",
           "burnin","thin","thin_z","sampler_scheme"))
  if(!same(k)) fail("Chains disagree on ",k,"; do not pool")
if(anyDuplicated(vapply(ch,function(z)z$seed,1)))fail("Duplicate seeds")
cat("Provenance: 3 chains, identical hashes/settings, distinct seeds: PASS\n")
cat("Draws per chain:",paste(vapply(ch,function(z)nrow(z$samples),1L)),
    "; sampler:",ch[[1]]$sampler_scheme,"\n")

# ---- derived quantities, per draw -----------------------------------------
per <- sprintf("%d-%02d",seq(1976,2021,5),seq(80,125,5)%%100)
derive <- function(x){
  q <- function(n) x[,n]
  rSE <- sapply(1:10,function(j)plogis(q("betaSE[1]")+q("betaSE[2]")+
                                         q(sprintf("yearSE[%d]",j))))
  rEI <- sapply(1:10,function(j)plogis(q("betaEI[1]")+q("betaEI[2]")+
                                         q(sprintf("yearEI[%d]",j))))
  colnames(rSE) <- paste0("qrisk_SE_",per)
  colnames(rEI) <- paste0("qrisk_EI_",per)
  sEI <- sqrt(q("load21")^2+q("sigmaEI")^2)
  cbind(OR_acquisition_MLH=exp(q("betaSE[4]")),
        OR_progression_MLH=exp(q("betaEI[4]")),
        HR_mortality_MLH_constant=exp(q("betaD[2]")),
        b_multiplier_MLH=exp(q("gammaB")),
        b_multiplier_MLH_infected=exp(q("gammaB")+q("deltaB")),
        b_MLH_x_infection=exp(q("deltaB")),
        sapply(c(2,4,6,8,10),function(x){
          b0 <- q("bGM"); b1 <- b0*exp(q("gammaB")); f <- q("f4"); h4 <- q("h4")
          hz <- function(b) h4*(1-f)*exp(b*(x-4))*expm1(b/4)/(b/4)+h4*f
          exp(q("betaD[2]"))*hz(b1)/hz(b0)}) |>
          (\(m) {colnames(m) <- paste0("HR_mortality_MLH_age",c(2,4,6,8,10)); m})(),
        sapply(c(2,4,6,8,10),function(x){
          b0 <- q("bGM"); b1 <- b0*exp(q("gammaB")+q("deltaB")); f <- q("f4"); h4 <- q("h4")
          hz <- function(b) h4*(1-f)*exp(b*(x-4))*expm1(b/4)/(b/4)+h4*f
          exp(q("betaD[2]"))*hz(b1)/hz(b0)}) |>
          (\(m) {colnames(m) <- paste0("HR_mortality_MLH_infected_age",c(2,4,6,8,10)); m})(),
        OR_acquisition_male=exp(q("betaSE[3]")),
        OR_progression_male=exp(q("betaEI[3]")),
        HR_mortality_male=exp(q("betaD[1]")),
        HR_mortality_E_vs_S=exp(q("betaD[3]")),
        HR_mortality_I_vs_S=exp(q("betaD[4]")),
        SD_frailty_SE=q("sigmaSE"),SD_frailty_EI=sEI,
        cor_frailty_SE_EI=q("load21")/sEI,
        qrisk_SE_longrun=plogis(q("betaSE[1]")+q("betaSE[2]")),
        qrisk_EI_longrun=plogis(q("betaEI[1]")+q("betaEI[2]")),
        rSE,rEI)
}
M <- lapply(ch,function(z)z$samples)
keep <- apply(do.call(rbind,M),2,sd)>1e-12
M <- lapply(M,function(m)m[,keep,drop=FALSE])
D <- lapply(ch,function(z)derive(z$samples))
summ <- function(L){
  arr <- array(NA_real_,c(nrow(L[[1]]),3,ncol(L[[1]])),
               list(NULL,NULL,colnames(L[[1]])))
  for(j in 1:3) arr[,j,] <- L[[j]]
  s <- posterior::summarise_draws(posterior::as_draws_array(arr),
    median,~quantile(.x,c(.025,.975)),"rhat","ess_bulk","ess_tail")
  as.data.frame(s)
}
raw <- summ(M); der <- summ(D)
bad <- raw[raw$rhat>1.01|raw$ess_bulk<400|raw$ess_tail<400,]
cat("\nRank-normalised: max R-hat",round(max(raw$rhat),3),
    "; min bulk ESS",round(min(raw$ess_bulk)),
    "; min tail ESS",round(min(raw$ess_tail)),
    "; failing R-hat<=1.01/ESS>=400:",nrow(bad),"/",nrow(raw),"\n")
if(nrow(bad)) print(head(bad[order(-bad$rhat),
  c("variable","median","rhat","ess_bulk","ess_tail")],20),digits=3)
cat("\n=== Effects and period risks (PROVISIONAL) ===\n")
print(der[,c("variable","median","2.5%","97.5%","rhat","ess_bulk")],
      digits=3,row.names=FALSE)

# ---- post-last-capture check: FFBS over a subset of draws ----------------
inp <- readRDS("outputs/joint_sei_survival/joint_PRIVATE_UNRECONCILED_input.rds")
N <- as.integer(inp$N); TT <- as.integer(inp$T)
y <- inp$y; first <- as.integer(inp$first_q); birth <- inp$birth_q
death <- ifelse(is.na(inp$death_q),0L,as.integer(inp$death_q))
hom <- inp$hom_z; male <- as.numeric(inp$sex=="MALE")
em <- inp$log_diagnostic_emission
pid <- (seq_len(TT)-1L)%/%20L+1L; sid <- (seq_len(TT)-1L)%%4L+1L
ndz <- as.integer(Sys.getenv("NDRAW_Z","100"))
step <- ch[[1]]$thin_z%/%ch[[1]]$thin
pick <- do.call(rbind,lapply(1:3,function(j){
  nz <- nrow(ch[[j]]$samples_z)
  if(nrow(ch[[j]]$samples)!=nz*step) fail("z/parameter draws misaligned")
  k <- round(seq(1,nz,length.out=min(ndz,nz)))
  cbind(j,k)
}))
P <- do.call(rbind,lapply(seq_len(nrow(pick)),function(r)
  ch[[pick[r,1]]]$samples[pick[r,2]*step,]))
Zs <- do.call(rbind,lapply(seq_len(nrow(pick)),function(r)
  ch[[pick[r,1]]]$samples_z[pick[r,2],]))
nd <- nrow(P); g <- function(n) P[,n]
B <- g("bGM"); A <- g("h4")*(1-g("f4"))*exp(-4*B); C <- g("h4")*g("f4")
ez <- expm1(B*.25)/(B*.25); rho <- g("rho")
ySE <- P[,sprintf("yearSE[%d]",1:10)]; yEI <- P[,sprintf("yearEI[%d]",1:10)]
per_e <- P[,sprintf("periodEff[%d]",1:10)]
sea_e <- P[,sprintf("seasonEff[%d]",1:4)]
draw_state <- function(W){
  W <- W/rowSums(W); u <- runif(nrow(W))
  1L+(u>W[,1])+(u>W[,1]+W[,2])+(u>W[,1]+W[,2]+W[,3])
}
set.seed(20261120)
onset <- matrix(NA_integer_,nd,N)
for(i in seq_len(N)){
  z1 <- Zs[,sprintf("z[%d, 1]",i)]; z2 <- Zs[,sprintf("z[%d, 2]",i)]
  z4 <- Zs[,sprintf("z[%d, 4]",i)]
  ts <- first[i]:(if(death[i]>0) death[i] else TT); K <- length(ts)
  piE <- plogis(g("entryBeta[1]")+g("entryBeta[2]")*(first[i]-birth[i])/16+
    g("entryBeta[3]")*male[i]+g("entryBeta[4]")*hom[i]+
    g("entryBeta[5]")*(first[i]-1)/40)
  a <- cbind(1-piE,piE*(1-g("entryRho")),piE*g("entryRho"),0)
  Al <- array(0,c(nd,4,K))
  phi <- array(0,c(nd,3,K)); lam <- tau <- matrix(0,nd,K)
  for(k in seq_len(K)){
    t <- ts[k]
    if(t>first[i] && y[i,t]!=-1){
      p <- plogis(g("alphaP")+g("sigmaCapture")*z4+per_e[,pid[t]]+
                  sea_e[,sid[t]])
      if(y[i,t]==1){a[,1:3] <- a[,1:3]*p; a[,4] <- 0
      } else a[,1:3] <- a[,1:3]*(1-p)
    }
    if(y[i,t]==1){
      e <- em[i,t,]; e <- exp(e-max(e))
      a[,1:3] <- sweep(a[,1:3,drop=FALSE],2,e,"*"); a[,4] <- 0
    }
    a <- a/rowSums(a); Al[,,k] <- a
    Bi <- B*exp(g("gammaB")*hom[i]); Ai <- g("h4")*(1-g("f4"))*exp(-4*Bi)
    age <- .25*(t-birth[i]); h <- .25*(Ai*exp(Bi*age)*expm1(Bi*.25)/(Bi*.25)+C)
    lp <- g("betaD[1]")*male[i]+g("betaD[2]")*hom[i]
    B2 <- B*exp((g("gammaB")+g("deltaB"))*hom[i]); A2 <- g("h4")*(1-g("f4"))*exp(-4*B2)
    h2 <- .25*(A2*exp(B2*age)*expm1(B2*.25)/(B2*.25)+C)
    phi[,,k] <- cbind(exp(-h*exp(lp)),exp(-h2*exp(lp+g("betaD[3]"))),
                      exp(-h2*exp(lp+g("betaD[4]"))))
    lam[,k] <- plogis(g("betaSE[1]")+g("betaSE[2]")*age/4+
      g("betaSE[3]")*male[i]+g("betaSE[4]")*hom[i]+
      g("sigmaSE")*z1+ySE[,pid[t]])
    tau[,k] <- plogis(g("betaEI[1]")+g("betaEI[2]")*age/4+
      g("betaEI[3]")*male[i]+g("betaEI[4]")*hom[i]+
      g("load21")*z1+g("sigmaEI")*z2+yEI[,pid[t]])
    if(k==K) break
    f <- phi[,,k]; dy <- rowSums(a[,1:3]*(1-f))
    a <- cbind(a[,1]*f[,1]*(1-lam[,k]),
               a[,1]*f[,1]*lam[,k]+a[,2]*f[,2]*(1-tau[,k]),
               a[,2]*f[,2]*tau[,k]+a[,3]*f[,3],a[,4]+(1-rho)*dy)
    a <- a/rowSums(a)
  }
  W <- Al[,,K]
  if(death[i]>0) W <- cbind(W[,1:3]*(1-phi[,,K]),0)
  s <- draw_state(W); on <- ifelse(s==3L,ts[K],NA_integer_)
  for(k in rev(seq_len(K-1L))){
    al <- Al[,,k]; f <- phi[,,k]; W <- matrix(0,nd,4)
    r <- s==1L; W[r,1] <- al[r,1]*f[r,1]*(1-lam[r,k])
    r <- s==2L; W[r,1] <- al[r,1]*f[r,1]*lam[r,k]
      W[r,2] <- al[r,2]*f[r,2]*(1-tau[r,k])
    r <- s==3L; W[r,2] <- al[r,2]*f[r,2]*tau[r,k]; W[r,3] <- al[r,3]*f[r,3]
    r <- s==4L; W[r,1:3] <- al[r,1:3]*(1-rho[r])*(1-f[r,]); W[r,4] <- al[r,4]
    s <- draw_state(W); on[s==3L] <- ts[k]
  }
  onset[,i] <- on
}
lastlive <- apply(y==1,1,function(v)max(which(v)))
after <- sweep(onset,2,lastlive,">")
share <- function(cols) rowSums(after[,cols,drop=FALSE],na.rm=TRUE)/
  pmax(rowSums(!is.na(onset[,cols,drop=FALSE])),1)
q3 <- function(v) round(quantile(v,c(.5,.025,.975)),3)
chk <- rbind(
  all_badgers=q3(share(seq_len(N))),
  recorded_death=q3(share(which(death>0))),
  no_recorded_death=q3(share(which(death==0))),
  n_with_I_onset=q3(rowSums(!is.na(onset))))
cat("\n=== Share of inferred I onsets after last live capture ===\n")
cat("(",nd,"draws; SEID collaborator's estimate: 0.41-0.43 )\n")
print(chk)

out <- file.path(root,paste0("review_",format(Sys.time(),"%Y%m%dT%H%M%SZ",tz="UTC")))
dir.create(out,recursive=TRUE)
write.csv(raw,file.path(out,"parameters.csv"),row.names=FALSE)
write.csv(der,file.path(out,"effects_and_period_risks.csv"),row.names=FALSE)
write.csv(data.frame(set=rownames(chk),chk,check.names=FALSE),
          file.path(out,"post_last_capture_I_onset.csv"),row.names=FALSE)
write.csv(data.frame(badger=seq_len(N),
  p_I_onset=colMeans(!is.na(onset)),
  p_onset_after_last_capture=colMeans(after&!is.na(after)),
  recorded_death=death>0),
  file.path(out,"PRIVATE_per_badger_onset.csv"),row.names=FALSE)
tp <- intersect(c("betaSE[1]","betaEI[1]","betaSE[4]","betaEI[4]","betaD[2]",
  "betaD[4]","h4","f4","bGM","sigmaSE","sigmaEI","sigmaYearSE",
  "sigmaYearEI","sigmaCapture"),colnames(M[[1]]))
pdf(file.path(out,"traces.pdf"),11,8); par(mfrow=c(3,2))
for(p in tp) matplot(sapply(M,function(m)m[,p]),type="l",lty=1,
                     main=p,ylab=p,xlab="saved draw")
dev.off()
cat("Saved:",out,"\n")
