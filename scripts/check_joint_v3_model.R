# ============================================================================
# Read-only checks on a converged joint run (default full_joint_v3_long).
#  A. Pathway split: change in survival from +1 SD homozygosity, total vs
#     through infection only (acquisition + progression) vs direct only
#     (mortality; includes the v3c senescence term if gammaB is present).
#     Female/male, S at age 1, long-run calendar, frailties integrated.
#  B. Known-death check: inferred state at death for recovered vs unseen
#     deaths (forward-filter backward-sample). The model assumes recovery
#     does not depend on state, so the two should have similar I shares;
#     a much higher I share among unseen deaths means the I-mortality
#     effect is carried mainly by disappearance.
#  C. Posterior predictive check: live captures (badgers without a recorded
#     death) and recovered deaths (all badgers) per 5-year period, observed
#     vs simulated forward from each badger's state at first capture.
# Usage: RUN_NAME=full_joint_v3_long Rscript scripts/check_joint_v3_model.R
# ============================================================================
options(stringsAsFactors=FALSE)
fail <- function(...) stop(paste0(...),call.=FALSE)
run <- Sys.getenv("RUN_NAME","full_joint_v3_long")
root <- file.path("outputs/joint_sei_survival",run)
dirs <- file.path(root,sprintf("chain_%02d",1:3))
ch <- lapply(file.path(dirs,"posterior_samples.rds"),readRDS)
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
  cbind(j,round(seq(1,nz,length.out=min(ndz,nz))))}))
P <- do.call(rbind,lapply(seq_len(nrow(pick)),function(r)
  ch[[pick[r,1]]]$samples[pick[r,2]*step,]))
Zs <- do.call(rbind,lapply(seq_len(nrow(pick)),function(r)
  ch[[pick[r,1]]]$samples_z[pick[r,2],]))
nd <- nrow(P); g <- function(n) P[,n]
gam <- if("gammaB"%in%colnames(P)) g("gammaB") else rep(0,nd)
q3 <- function(v) round(quantile(v,c(.5,.025,.975),na.rm=TRUE),4)
hz <- function(age,hd,ml){   # quarterly cumulative hazard, S state
  B <- g("bGM")*exp(gam*hd); A <- g("h4")*(1-g("f4"))*exp(-4*B)
  .25*(A*exp(B*age)*expm1(B*.25)/(B*.25)+g("h4")*g("f4"))*
    exp(g("betaD[1]")*ml+g("betaD[2]")*hd)}
set.seed(20261125)
out <- file.path(root,paste0("checks_",format(Sys.time(),"%Y%m%dT%H%M%SZ",tz="UTC")))
dir.create(out,recursive=TRUE)

# ---- A. pathway split ------------------------------------------------------
M <- 200L; ages <- c(2,4,6,8,10)
surv <- function(hse,hei,hd,ml){
  z1 <- matrix(rnorm(nd*M),nd); z2 <- matrix(rnorm(nd*M),nd)
  a <- list(S=matrix(1,nd,M),E=matrix(0,nd,M),I=matrix(0,nd,M))
  res <- matrix(NA,nd,length(ages))
  for(q in 4:39){ age <- q/4
    h <- hz(age,hd,ml)
    fS <- exp(-h); fE <- exp(-h*exp(g("betaD[3]"))); fI <- exp(-h*exp(g("betaD[4]")))
    lam <- plogis(g("betaSE[1]")+g("betaSE[2]")*age/4+g("betaSE[3]")*ml+
                  g("betaSE[4]")*hse+g("sigmaSE")*z1)
    tau <- plogis(g("betaEI[1]")+g("betaEI[2]")*age/4+g("betaEI[3]")*ml+
                  g("betaEI[4]")*hei+g("load21")*z1+g("sigmaEI")*z2)
    a <- list(S=a$S*fS*(1-lam),E=a$S*fS*lam+a$E*fE*(1-tau),I=a$E*fE*tau+a$I*fI)
    k <- match((q+1)/4,ages); if(!is.na(k)) res[,k] <- rowMeans(a$S+a$E+a$I)
  }
  res}
pw <- do.call(rbind,lapply(c(0,1),function(ml){
  s0 <- surv(0,0,0,ml); d <- list(total=surv(1,1,1,ml)-s0,
    infection_route=surv(1,1,0,ml)-s0,direct_route=surv(0,0,1,ml)-s0)
  do.call(rbind,lapply(seq_along(ages),function(k) data.frame(
    sex=ifelse(ml==1,"male","female"),age=ages[k],
    survival_from_age1_at_mean_hom=paste(q3(s0[,k]),collapse=" "),
    total_change=paste(q3(d$total[,k]),collapse=" "),
    infection_route=paste(q3(d$infection_route[,k]),collapse=" "),
    direct_route=paste(q3(d$direct_route[,k]),collapse=" "))))}))
cat("\n=== A. Change in survival from age 1, +1 SD homozygosity",
    "(median 2.5% 97.5%) ===\n"); print(pw,row.names=FALSE)
write.csv(pw,file.path(out,"A_pathway_split.csv"),row.names=FALSE)

# ---- B + C. per-badger filter, backward sample, forward simulation --------
B0 <- g("bGM"); rho <- g("rho")
ySE <- P[,sprintf("yearSE[%d]",1:10)]; yEI <- P[,sprintf("yearEI[%d]",1:10)]
per_e <- P[,sprintf("periodEff[%d]",1:10)]; sea_e <- P[,sprintf("seasonEff[%d]",1:4)]
draw <- function(W){ W <- W/rowSums(W); u <- runif(nrow(W))
  1L+(u>W[,1])+(u>W[,1]+W[,2])+(u>W[,1]+W[,2]+W[,3])}
dS <- array(0,c(nd,2,3))   # [draw, recovered/unseen, S/E/I] deaths
per5 <- function(t) (t-1L)%/%20L+1L
obs_cap <- obs_dead <- numeric(10)
rep_cap <- rep_dead <- matrix(0,nd,10)
for(i in seq_len(N)){
  z1 <- Zs[,sprintf("z[%d, 1]",i)]; z2 <- Zs[,sprintf("z[%d, 2]",i)]
  z4 <- Zs[,sprintf("z[%d, 4]",i)]
  ta <- first[i]:TT; KA <- length(ta)
  K <- if(death[i]>0) death[i]-first[i]+1L else KA
  phi <- array(0,c(nd,3,KA)); lam <- tau <- pc <- matrix(0,nd,KA)
  for(k in seq_len(KA)){ t <- ta[k]; age <- .25*(t-birth[i])
    h <- hz(age,hom[i],male[i])
    phi[,,k] <- cbind(exp(-h),exp(-h*exp(g("betaD[3]"))),exp(-h*exp(g("betaD[4]"))))
    lam[,k] <- plogis(g("betaSE[1]")+g("betaSE[2]")*age/4+g("betaSE[3]")*male[i]+
      g("betaSE[4]")*hom[i]+g("sigmaSE")*z1+ySE[,pid[t]])
    tau[,k] <- plogis(g("betaEI[1]")+g("betaEI[2]")*age/4+g("betaEI[3]")*male[i]+
      g("betaEI[4]")*hom[i]+g("load21")*z1+g("sigmaEI")*z2+yEI[,pid[t]])
    pc[,k] <- plogis(g("alphaP")+g("sigmaCapture")*z4+per_e[,pid[t]]+sea_e[,sid[t]])
  }
  piE <- plogis(g("entryBeta[1]")+g("entryBeta[2]")*(first[i]-birth[i])/16+
    g("entryBeta[3]")*male[i]+g("entryBeta[4]")*hom[i]+
    g("entryBeta[5]")*(first[i]-1)/40)
  a <- cbind(1-piE,piE*(1-g("entryRho")),piE*g("entryRho"),0)
  Al <- array(0,c(nd,4,K))
  for(k in seq_len(K)){ t <- ta[k]
    if(k>1 && y[i,t]!=-1){
      if(y[i,t]==1){a[,1:3] <- a[,1:3]*pc[,k]; a[,4] <- 0
      } else a[,1:3] <- a[,1:3]*(1-pc[,k])}
    if(y[i,t]==1){ e <- em[i,t,]; e <- exp(e-max(e))
      a[,1:3] <- sweep(a[,1:3,drop=FALSE],2,e,"*"); a[,4] <- 0}
    a <- a/rowSums(a); Al[,,k] <- a
    if(k==K) break
    f <- phi[,,k]; dy <- rowSums(a[,1:3]*(1-f))
    a <- cbind(a[,1]*f[,1]*(1-lam[,k]),a[,1]*f[,1]*lam[,k]+a[,2]*f[,2]*(1-tau[,k]),
               a[,2]*f[,2]*tau[,k]+a[,3]*f[,3],a[,4]+(1-rho)*dy)
    a <- a/rowSums(a)
  }
  # B: backward sample, record state at death
  W <- Al[,,K]; if(death[i]>0) W <- cbind(W[,1:3]*(1-phi[,,K]),0)
  s <- draw(W); r <- cbind(seq_len(nd),1L,pmin(s,3L))
  if(death[i]>0) dS[r] <- dS[r]+1
  for(k in rev(seq_len(K-1L))){
    al <- Al[,,k]; f <- phi[,,k]; W <- matrix(0,nd,4); nx <- s
    w <- nx==1L; W[w,1] <- al[w,1]*f[w,1]*(1-lam[w,k])
    w <- nx==2L; W[w,1] <- al[w,1]*f[w,1]*lam[w,k]; W[w,2] <- al[w,2]*f[w,2]*(1-tau[w,k])
    w <- nx==3L; W[w,2] <- al[w,2]*f[w,2]*tau[w,k]; W[w,3] <- al[w,3]*f[w,3]
    w <- nx==4L; W[w,1:3] <- al[w,1:3]*(1-rho[w])*(1-f[w,]); W[w,4] <- al[w,4]
    s <- draw(W); u <- nx==4L & s<4L
    if(any(u)) dS[cbind(which(u),2L,s[u])] <- dS[cbind(which(u),2L,s[u])]+1
  }
  # C: observed and simulated captures/deaths by 5-year period
  if(death[i]>0) obs_dead[per5(death[i])] <- obs_dead[per5(death[i])]+1
  if(death[i]==0){ tt <- ta[-1][y[i,ta[-1]]==1]
    for(t in tt) obs_cap[per5(t)] <- obs_cap[per5(t)]+1}
  s <- draw(Al[,,1]); alive <- s<4L
  for(k in seq_len(KA)){ t <- ta[k]
    if(k>1 && death[i]==0 && y[i,t]!=-1){
      c1 <- alive & runif(nd)<pc[,k]; rep_cap[,per5(t)] <- rep_cap[,per5(t)]+c1}
    f <- phi[,,k][cbind(seq_len(nd),pmin(s,3L))]
    die <- alive & runif(nd)>f
    rep_dead[,per5(t)] <- rep_dead[,per5(t)]+(die & runif(nd)<rho)
    u <- runif(nd)
    s <- ifelse(die,4L,ifelse(s==1L & u<lam[,k],2L,ifelse(s==2L & u<tau[,k],3L,s)))
    alive <- s<4L
  }
}
sh <- function(m) m[,3]/pmax(rowSums(m),1)
kd <- rbind(recovered_deaths_share_I=q3(sh(dS[,1,])),
  unseen_deaths_share_I=q3(sh(dS[,2,])),
  recovered_deaths_share_E=q3(dS[,1,2]/pmax(rowSums(dS[,1,]),1)),
  unseen_deaths_share_E=q3(dS[,2,2]/pmax(rowSums(dS[,2,]),1)),
  n_recovered=q3(rowSums(dS[,1,])),n_unseen=q3(rowSums(dS[,2,])))
cat("\n=== B. Inferred state at death: recovered vs unseen ===\n"); print(kd)
write.csv(data.frame(quantity=rownames(kd),kd,check.names=FALSE),
          file.path(out,"B_known_death_check.csv"),row.names=FALSE)
lab <- sprintf("%d-%02d",seq(1976,2021,5),seq(80,125,5)%%100)
ppc <- function(o,r) data.frame(period=lab,observed=o,
  sim_median=apply(r,2,median),sim_2.5=apply(r,2,quantile,.025),
  sim_97.5=apply(r,2,quantile,.975),p_sim_ge_obs=round(colMeans(sweep(r,2,o,">=")),3))
pc_cap <- ppc(obs_cap,rep_cap); pc_dead <- ppc(obs_dead,rep_dead)
cat("\n=== C. PPC: live recaptures, badgers with no recorded death ===\n")
print(pc_cap,row.names=FALSE)
cat("\n=== C. PPC: recovered deaths, all badgers ===\n"); print(pc_dead,row.names=FALSE)
write.csv(pc_cap,file.path(out,"C_ppc_captures.csv"),row.names=FALSE)
write.csv(pc_dead,file.path(out,"C_ppc_recovered_deaths.csv"),row.names=FALSE)
cat("Saved:",out,"\n")
