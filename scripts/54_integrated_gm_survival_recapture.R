# ================================================================
# Stage 54: Integrated Gompertz-Makeham survival + recapture +
#          uncertainty-marginalised time-varying infection.
#
# MODE=compile [default]: tiny compilation/mathematical gate only.
# MODE=full: one substantive MCMC chain, after approving inputs.
#
# Survival is marginalised via a 2-state forward algorithm (alive,
# dead but not recovered); NOT 1603 x 200 sampled latent states.
#
# Infection timing: exact per-individual marginalisation over all
# unique onset quarters and their frequencies among the 500 draws.
# This keeps FULL marginal uncertainty (no mean infection time).
# IMPORTANT: cross-badger shared-column correlations are NOT
# represented in this factorised modular likelihood. Full shared-
# draw dependence requires a separate joint-likelihood treatment.
# ================================================================
options(stringsAsFactors=FALSE)
stopifnot(requireNamespace("nimble", quietly=TRUE))
suppressPackageStartupMessages(library(nimble))
dir.create("outputs/integrated_gm",recursive=TRUE,showWarnings=FALSE)
mode <- Sys.getenv("MODE", "compile")
if(!mode %in% c("compile","full")) stop("MODE must be compile or full")
file <- Sys.getenv("MODEL_INPUT",
                   "data/model_inputs/integrated_gm_inputs_group_proxy_2026.rds")
if(!file.exists(file)) stop("Missing model input: ", file)
x <- readRDS(file)
stopifnot(
  x$T == 200L,
  identical(colnames(x$infection_time),as.character(x$infection_draw_ids)),
  ncol(x$infection_time)==500L,
  all(x$sex %in% c("FEMALE","MALE")),
  all(x$y %in% c(-1L,0L,1L)),
  all(x$first_q >= x$birth_q)
)

# Conditional-on-first-capture mark-recapture likelihood.
# y=-1 missing/unmonitored, y=0 noncapture when monitoring proxy
# supported, y=1 confirmed live; a terminal observed death is a
# separate observation with recovery probability rho.
#
# a,b,c are YEAR-SCALE GM parameters:
# h(age in years) = a exp(b*age)+c, multiplied by effects.
# Exact within-quarter cumulative-hazard integration incorporates
# within-quarter linear-age interactions.
# NIMBLE promotes the y[i,1:T] constant row to NimArr<double> in
# compiled model code, even when supplied as an integer R matrix.
# Both the single-history and marginal-history nimbleFunction
# signatures must therefore accept double(1) for y; observed
# values still remain exactly -1, 0, 1.
gm_cmr_single <- nimbleFunction(
  run=function(
    y=double(1), first=integer(0), death=integer(0),
    birth=integer(0), infection=double(0),
    hom=double(0), male=double(0),
    a=double(0), b=double(0), c1=double(0),
    bh=double(0), bs=double(0), bi=double(0),
    bhi=double(0), bha=double(0), bia=double(0),
    alphaP=double(0), ui=double(0),
    periodEff=double(1), seasonEff=double(1),
    periodIndex=integer(1), seasonIndex=integer(1),
    rho=double(0)) {
    returnType(double(0))
    alive <- 1.0
    deadstate <- 0.0
    logAccum <- 0.0
    nt <- length(y)
    for(t in first:nt) {
      # Do not use the first, ascertainment-defining capture to
      # estimate recapture probability.
      if(t > first) {
        lp <- alphaP + ui +
          periodEff[periodIndex[t]] + seasonEff[seasonIndex[t]]
        p <- 1.0/(1.0+exp(-lp))
        if(y[t] == 1L) {
          alive <- alive*p
          deadstate <- 0.0
        } else if(y[t] == 0L) {
          alive <- alive*(1.0-p)
        }
      }

      age <- (t-birth)*0.25
      infected <- 0.0
      if(infection > 0L & t >= infection) infected <- 1.0

      # Hazard multiplicative effect at current age.
      ageSlope <- bha*hom + bia*infected
      intercept <- bh*hom + bs*male + bi*infected +
        bhi*hom*infected + ageSlope*(age-4.0)

      # Exact integral of (a*exp(b*age)+c1)*exp(intercept)
      # as age increases by 0.25 year during the quarter.
      z1 <- (b+ageSlope)*0.25
      z2 <- ageSlope*0.25
      e1 <- 1.0
      e2 <- 1.0
      if(abs(z1)<1.0e-6) {
        e1 <- 1.0 + z1/2.0 + z1*z1/6.0
      } else {
        e1 <- (exp(z1)-1.0)/z1
      }
      if(abs(z2)<1.0e-6) {
        e2 <- 1.0 + z2/2.0 + z2*z2/6.0
      } else {
        e2 <- (exp(z2)-1.0)/z2
      }
      H <- 0.25*exp(intercept) *
        (a*exp(b*age)*e1 + c1*e2)
      phi <- exp(-H)
      die <- 1.0-phi

      if(t == death) {
        # An animal caught alive and then killed in the same
        # quarter can contribute both events.
        prob <- alive*die*rho
        if(prob <= 0.0) return(-1.0e12)
        return(logAccum+log(prob))
      }

      # If no recovery this quarter, it may be:
      #  - alive at next quarter, or
      #  - dead, but never observed as dead.
      nextAlive <- alive*phi
      nextDead <- deadstate + alive*die*(1.0-rho)
      scale <- nextAlive+nextDead
      if(scale <= 0.0) return(-1.0e12)
      logAccum <- logAccum+log(scale)
      alive <- nextAlive/scale
      deadstate <- nextDead/scale
    }
    return(logAccum)
  }
)

# NIMBLE converts a row-slice of an integer constant matrix to
# NimArr<double> in the compiled calculation graph. The support
# array therefore has type double(1), and its scalar infection onset
# passed to gm_cmr_single has type double(0). Values remain exact
# whole-quarter indices (0, 1, 2, ...) and are never averaged.
# Log-sum-exp over the *empirical* 500-draw infection posterior,
# compressed to its unique onset quarters and exact frequencies.
# This avoids 500 redundant repeats per individual.
gm_cmr_marginal <- nimbleFunction(
  run=function(
    y=double(1), first=integer(0), death=integer(0),
    birth=integer(0),
    support=double(1), weight=double(1), nsupport=integer(0),
    hom=double(0), male=double(0),
    a=double(0), b=double(0), c1=double(0),
    bh=double(0), bs=double(0), bi=double(0),
    bhi=double(0), bha=double(0), bia=double(0),
    alphaP=double(0), ui=double(0),
    periodEff=double(1), seasonEff=double(1),
    periodIndex=integer(1), seasonIndex=integer(1),
    rho=double(0)) {
    returnType(double(0))
    maxll <- -1.0e12
    for(k in 1:nsupport) {
      ll <- gm_cmr_single(
        y,first,death,birth,support[k],
        hom,male,a,b,c1,bh,bs,bi,bhi,bha,bia,
        alphaP,ui,periodEff,seasonEff,
        periodIndex,seasonIndex,rho
      )
      if(ll > maxll) maxll <- ll
    }
    total <- 0.0
    for(k in 1:nsupport) {
      ll <- gm_cmr_single(
        y,first,death,birth,support[k],
        hom,male,a,b,c1,bh,bs,bi,bhi,bha,bia,
        alphaP,ui,periodEff,seasonEff,
        periodIndex,seasonIndex,rho
      )
      total <- total + weight[k]*exp(ll-maxll)
    }
    if(total <= 0.0) return(-1.0e12)
    return(maxll+log(total))
  }
)

# Independent analytic toy check: first capture quarter 1,
# missed capture quarter 2, known death during quarter 2.
toy <- gm_cmr_single(
  as.integer(c(1,0)),1L,2L,1L,0L,
  0,0,0,0,0.08,0,0,0,0,0,0,
  qlogis(0.4),0,c(0),c(0),
  as.integer(c(1,1)),as.integer(c(1,1)),0.3
)
phi <- exp(-0.08*.25)
expected <- phi*(1-0.4)*(1-phi)*0.3
cat("Independent toy-likelihood check:",
    round(exp(toy),12),"vs",round(expected,12),"\n")
if(!is.finite(toy) || abs(exp(toy)-expected)>1e-10)
  stop("Analytic forward-likelihood check FAILED")

# Toy infection onset in second quarter doubles quarter-2 hazard.
toyInf <- gm_cmr_single(
  as.integer(c(1,0)),1L,2L,1L,2L,
  0,0,0,0,0.08,0,0,log(2),0,0,0,
  qlogis(0.4),0,c(0),c(0),
  as.integer(c(1,1)),as.integer(c(1,1)),0.3
)
expectedInf <- phi*(1-0.4)*(1-exp(-0.16*.25))*0.3
cat("Infection-onset hazard check:",
    round(exp(toyInf),12),"vs",round(expectedInf,12),"\n")
if(!is.finite(toyInf) || abs(exp(toyInf)-expectedInf)>1e-10)
  stop("Infection-onset numerical check FAILED")

# Compress the complete infection distribution exactly per animal.
n_total <- x$N
periodIndex <- as.integer((seq_len(x$T)-1L)%/%20L + 1L)
seasonIndex <- as.integer((seq_len(x$T)-1L)%%4L + 1L)
homRaw <- x$hom
hom <- (homRaw-mean(homRaw))/sd(homRaw)
male <- as.numeric(x$sex=="MALE")
death <- as.integer(ifelse(is.na(x$death_q),0L,x$death_q))
freq <- lapply(seq_len(n_total),function(i) {
  t <- table(x$infection_time[i,])
  list(onset=as.integer(names(t)),weight=as.numeric(t)/500)
})
maxSupport <- max(vapply(freq,function(z)length(z$onset),integer(1)))
supportAll <- matrix(0L,n_total,maxSupport)
weightAll <- matrix(0,n_total,maxSupport)
nSupportAll <- integer(n_total)
for(i in seq_len(n_total)) {
  nk <- length(freq[[i]]$onset)
  nSupportAll[i] <- nk
  supportAll[i,seq_len(nk)] <- freq[[i]]$onset
  weightAll[i,seq_len(nk)] <- freq[[i]]$weight
}
if(any(abs(rowSums(weightAll)-1)>1e-12))
  stop("Infection draw compression lost posterior weights")
cat("500 draws compressed to max",maxSupport,
    "distinct onset quarters; median",
    median(nSupportAll),"\n")
cat("IMPORTANT: Infection-time uncertainty is marginalised per individual.\n")
cat("Shared cross-individual posterior draw correlations are not yet included.\n")

if(mode=="compile") {
  deadCandidates <- which(death > 0L)
  aliveCandidates <- which(death == 0L)
  set.seed(20261113)
  sel <- sort(c(
    sample(deadCandidates,min(5L,length(deadCandidates))),
    sample(aliveCandidates,min(12L,length(aliveCandidates)))
  ))
} else {
  sel <- seq_len(n_total)
}
N <- length(sel)
cat("Model mode:",mode,"; individuals:",N,
    "of",n_total,"\n")

# Poisson-zero trick: zeros~Poisson(-log L);
# P(zeros=0)=exp(log L). This is a proper likelihood
# because the forward observation-sequence probability <= 1.
code <- nimbleCode({
  a ~ dexp(10)
  b ~ dexp(5)
  c1 ~ dexp(30)
  bh ~ dnorm(0,sd=0.6)
  bs ~ dnorm(0,sd=0.6)
  bi ~ dnorm(0,sd=0.6)
  bhi ~ dnorm(0,sd=0.6)
  bha ~ dnorm(0,sd=0.3)
  bia ~ dnorm(0,sd=0.3)
  alphaP ~ dnorm(-1.5,sd=1)
  sigmaU ~ dunif(0.01,2)
  rho ~ dbeta(2,3)

  periodEff[1] <- 0
  seasonEff[1] <- 0
  for(j in 2:10) periodEff[j] ~ dnorm(0,sd=0.5)
  for(j in 2:4) seasonEff[j] ~ dnorm(0,sd=0.5)

  for(i in 1:N) {
    ui[i] ~ dnorm(0,sd=sigmaU)
    ll[i] <- gm_cmr_marginal(
      y[i,1:T],first[i],death[i],birth[i],
      support[i,1:K],weight[i,1:K],nsupport[i],
      hom[i],male[i],
      a,b,c1,bh,bs,bi,bhi,bha,bia,
      alphaP,ui[i],
      periodEff[1:10],seasonEff[1:4],
      periodIndex[1:T],seasonIndex[1:T],
      rho
    )
    zeros[i] ~ dpois(-ll[i])
  }
})

constants <- list(
  N=N,T=x$T,K=maxSupport,
  y=matrix(as.integer(x$y[sel,]),nrow=N),
  first=as.integer(x$first_q[sel]),
  death=death[sel],birth=as.integer(x$birth_q[sel]),
  support=matrix(as.integer(supportAll[sel,]),nrow=N),
  weight=matrix(as.numeric(weightAll[sel,]),nrow=N),
  nsupport=as.integer(nSupportAll[sel]),
  hom=as.numeric(hom[sel]),male=as.numeric(male[sel]),
  periodIndex=periodIndex,seasonIndex=seasonIndex
)
data <- list(zeros=rep(0L,N))
inits <- function() list(
  a=0.08,b=0.18,c1=0.015,
  bh=0,bs=0,bi=0,bhi=0,bha=0,bia=0,
  alphaP=-1.5,sigmaU=0.5,rho=0.4,
  periodEff=c(NA,rep(0,9)),
  seasonEff=c(NA,rep(0,3)),
  ui=rep(0,N)
)

cat("Defining model...\n")
model <- nimbleModel(
  code,constants=constants,data=data,inits=inits(),
  check=TRUE,calculate=FALSE
)
logp <- model$calculate()
cat("Initial log probability:",logp,"\n")
if(!is.finite(logp)) stop("Non-finite model log probability")

cat("Compiling integrated forward likelihood...\n")
cm <- tryCatch(
  compileNimble(model),
  error = function(e) {
    message("\n=== NIMBLE C++ COMPILATION DIAGNOSTICS ===")
    try(printErrors(), silent = TRUE)
    stop("Integrated model compilation failed: ",
         conditionMessage(e), call. = FALSE)
  }
)
compiledLogp <- cm$calculate()
cat("Compiled log probability:",compiledLogp,"\n")
if(!is.finite(compiledLogp) ||
   abs(logp-compiledLogp)>1e-5)
  stop("Compiled vs R model likelihood mismatch")

if(mode=="compile") {
  saveRDS(list(
    n=N, ids=x$id[sel],
    model_logp_R=logp,model_logp_compiled=compiledLogp,
    max_support=maxSupport,
    median_support=median(nSupportAll),
    input=file
  ),"outputs/integrated_gm/compile_gate.rds")
  cat("\nCOMPILE GATE PASSED.\n")
  cat("Not a fitted model. Simulation recovery still required.\n")
} else {
  # A full fit is opt-in only, never triggered by the audit.
  niter <- as.integer(Sys.getenv("NITER","12000"))
  nburn <- as.integer(Sys.getenv("NBURN","4000"))
  thin <- as.integer(Sys.getenv("THIN","8"))
  chain <- as.integer(Sys.getenv("CHAIN_ID","1"))
  if(!is.finite(niter) || !is.finite(nburn) ||
     niter <= nburn) stop("Bad MCMC settings")
  cfg <- configureMCMC(
    model,
    monitors=c(
      "a","b","c1","bh","bs","bi","bhi","bha","bia",
      "alphaP","sigmaU","rho","periodEff","seasonEff"
    ),
    thin=thin
  )
  mcmc <- buildMCMC(cfg)
  cmcmc <- tryCatch(
    compileNimble(mcmc, project=model),
    error = function(e) {
      message("\n=== NIMBLE MCMC COMPILATION DIAGNOSTICS ===")
      try(printErrors(), silent = TRUE)
      stop("MCMC compilation failed: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  set.seed(20261113+chain*1009L)
  samples <- runMCMC(
    cmcmc,niter=niter,nburnin=nburn,thin=thin,
    nchains=1,progressBar=TRUE,
    samplesAsCodaMCMC=FALSE
  )
  saveRDS(
    list(samples=samples,chain=chain,n=N,
         input=file,
         input_md5=unname(tools::md5sum(file)),
         assumptions=c(
           "individual posterior infection marginals; joint draw correlation not included",
           "SG-specific monitoring proxy, not independent effort",
           "last observed group carried forward, no movement submodel",
           "shared recovery probability rho; all-cause mortality",
           "historical birth ages, first capture conditioned"
         )),
    sprintf("outputs/integrated_gm/full_chain_%02d.rds",chain)
  )
  cat("Saved full chain",chain,"\n")
}
