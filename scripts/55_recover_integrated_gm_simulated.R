# =====================================================================
# Stage 55: ONE targeted parameter-recovery test for the SAME NIMBLE
#           forward likelihood used by the substantive Stage 54 model.
#
# Data: simulated only. No Supabase access, no source-file changes.
# Synthetic first capture is conditioned on as in the real model.
# Infection timing is drawn from an empirical 500-draw posterior.
#
# Validation scope: a, MLH(hom), infection, infection*age,
# mean detection, and observed-death recovery rho.
#
# Other covariates, true individual detection random effects and
# calendar effects are fixed at their generating values for this
# conditional/oracle recovery gate. THIS IS NOT a complete recovery
# test for sigmaU, period effects or joint cross-badger infection
# dependence. It is a focused validation of the actual likelihood.
# =====================================================================
options(stringsAsFactors=FALSE)
stopifnot(requireNamespace("nimble",quietly=TRUE))
suppressPackageStartupMessages(library(nimble))
source("R/gm_cmr_forward_nimble.R")
dir.create("outputs/integrated_gm",recursive=TRUE,showWarnings=FALSE)
set.seed(20261113)

N <- as.integer(Sys.getenv("NSIM","600"))
T <- 40L
if(!is.finite(N) || N < 200L || N > 3000L)
  stop("NSIM should be between 200 and 3000")
periodIndex <- as.integer((seq_len(T)-1L)%/%20L+1L)
seasonIndex <- as.integer((seq_len(T)-1L)%%4L+1L)
periodEff <- c(0,.16)
seasonEff <- c(0,.12,-.09,.14)
monitored <- rep(TRUE,T)
monitored[seq(6,T,by=7)] <- FALSE

truth <- list(
  a=.085, b=.12, c1=.018,
  bh=.25, bs=.12, bi=.55, bhi=.12, bha=.025, bia=.08,
  alphaP=-.6, sigmaU=.45, rho=.45
)

# Independently implement integrated quarter hazard as a reference.
exprel <- function(z) {
  if(abs(z)<1e-7) 1+z/2+z^2/6+z^3/24 else expm1(z)/z
}
quarter_hazard <- function(t,birth,hom,male,tau) {
  age <- (t-birth)/4
  infected <- as.numeric(tau>0L && t>=tau)
  sl <- truth$bha*hom+truth$bia*infected
  int <- truth$bh*hom+truth$bs*male+truth$bi*infected+
    truth$bhi*hom*infected+sl*(age-4)
  .25*exp(int) * (
    truth$a*exp(truth$b*age)*exprel(.25*(truth$b+sl)) +
    truth$c1*exprel(.25*sl)
  )
}

hom <- rnorm(N)
male <- rbinom(N,1,.5)
ui <- rnorm(N,0,truth$sigmaU)
entry_eligible <- which(monitored & seq_len(T)<=8L)
first <- sample(entry_eligible,N,replace=TRUE)
birth <- first-sample(0:32,N,replace=TRUE)
y <- matrix(-1L,N,T)
death <- integer(N)
true_tau <- integer(N)
infection_time <- matrix(0L,N,500L)
support <- matrix(0L,N,3L)
weight <- matrix(0,N,3L)
nsupport <- integer(N)
never_observed_dead <- integer(N)

for(i in seq_len(N)) {
  center <- if(runif(1)<.28) 0L else sample(6:31,1L)
  if(center==0L) {
    possible <- c(0L,T+6L)
    probs <- c(.88,.12)
  } else {
    possible <- c(center-3L,center,center+3L)
    probs <- c(.20,.60,.20)
  }
  # Draws themselves define the exact uncertainty mixture used by
  # inference: select the latent true tau from these 500 draws.
  infection_time[i,] <- sample(possible,500L,replace=TRUE,prob=probs)
  true_tau[i] <- sample(infection_time[i,],1)
  tab <- table(infection_time[i,])
  k <- length(tab)
  nsupport[i] <- k
  support[i,seq_len(k)] <- as.integer(names(tab))
  weight[i,seq_len(k)] <- as.numeric(tab)/500

  alive <- TRUE
  for(t in first[i]:T) {
    if(monitored[t]) y[i,t] <- 0L
    if(alive) {
      if(t==first[i]) {
        y[i,t] <- 1L # model conditions on initial ascertainment
      } else if(monitored[t]) {
        p <- plogis(
          truth$alphaP+ui[i]+
            periodEff[periodIndex[t]]+seasonEff[seasonIndex[t]]
        )
        y[i,t] <- rbinom(1,1,p)
      }
      die <- rbinom(1,1,-expm1(-quarter_hazard(
        t,birth[i],hom[i],male[i],true_tau[i]
      )))
      if(die==1L) {
        alive <- FALSE
        if(rbinom(1,1,truth$rho)==1L) {
          death[i] <- t
          if(t<T) y[i,(t+1L):T] <- -1L
          break
        } else {
          never_observed_dead[i] <- 1L
        }
      }
    }
  }
}
stopifnot(
  all(y[cbind(seq_len(N),first)]==1L),
  all(y %in% c(-1L,0L,1L)),
  all(rowSums(weight)==1),
  all(infection_time>=0L),
  !anyNA(y)
)
cat("\n=== SYNTHETIC MODEL RECOVERY DESIGN ===\n")
cat("Individuals:",N,"quarters:",T,"\n")
cat("Known death recoveries:",sum(death>0L),"\n")
cat("Deaths not recovered:",sum(never_observed_dead>0L),"\n")
cat("Confirmed-live cells:",sum(y==1L),"\n")
cat("Nondetection cells:",sum(y==0L),"\n")
cat("Missing/after-death cells:",sum(y==-1L),"\n")
cat("Infection-time posterior draws:",ncol(infection_time),"\n")

# Same function(s) as Stage 54. No rewritten likelihood.
# The nuisance parameters and individual ui are fixed to their
# generating values; estimated core parameters are stochastic.
model_code <- nimbleCode({
  a ~ dexp(10)
  bh ~ dnorm(0,sd=.6)
  bi ~ dnorm(0,sd=.6)
  bia ~ dnorm(0,sd=.3)
  alphaP ~ dnorm(-1.5,sd=1)
  rho ~ dbeta(2,3)
  for(i in 1:N) {
    ll[i] <- gm_cmr_marginal(
      y[i,1:T],first[i],death[i],birth[i],
      support[i,1:K],weight[i,1:K],nsupport[i],
      hom[i],male[i],
      a,b,c1,bh,bs,bi,bhi,bha,bia,
      alphaP,ui[i],periodEff[1:2],seasonEff[1:4],
      periodIndex[1:T],seasonIndex[1:T],
      rho
    )
    zeros[i] ~ dpois(-ll[i])
  }
})
constants <- list(
  N=N,T=T,K=3L,
  y=matrix(as.numeric(y),nrow=N),
  first=as.integer(first),birth=as.integer(birth),
  death=as.integer(death),
  support=matrix(as.numeric(support),nrow=N),
  weight=weight,nsupport=as.integer(nsupport),
  hom=as.numeric(hom),male=as.numeric(male),
  b=truth$b,c1=truth$c1,bs=truth$bs,
  bhi=truth$bhi,bha=truth$bha,
  ui=as.numeric(ui),
  periodEff=periodEff,seasonEff=seasonEff,
  periodIndex=periodIndex,seasonIndex=seasonIndex
)
inits <- list(
  a=.075,bh=.1,bi=.3,bia=0,
  alphaP=-.7,rho=.4
)
cat("\nDefining same joint forward likelihood on simulated data...\n")
model <- nimbleModel(
  model_code,constants=constants,
  data=list(zeros=rep(0L,N)),
  inits=inits,check=TRUE,calculate=FALSE
)
lpR <- model$calculate()
if(!is.finite(lpR))stop("Synthetic likelihood non-finite in R")
cat("Compiling synthetic-data likelihood...\n")
cm <- tryCatch(
  compileNimble(model),
  error=function(e){
    try(printErrors(),silent=TRUE)
    stop(e)
  }
)
lpC <- cm$calculate()
cat("Interpreted log probability:",lpR,"\n")
cat("Compiled log probability:",lpC,"\n")
if(!is.finite(lpC)||abs(lpC-lpR)>1e-5)
  stop("Compiled and interpreted synthetic-data likelihood mismatch")

to_model <- function(par) {
  list(
    a=exp(par[1]),bh=par[2],bi=par[3],bia=par[4],
    alphaP=par[5],rho=plogis(par[6])
  )
}
assign_params <- function(par) {
  v <- to_model(par)
  cm$a <<- v$a
  cm$bh <<- v$bh
  cm$bi <<- v$bi
  cm$bia <<- v$bia
  cm$alphaP <<- v$alphaP
  cm$rho <<- v$rho
  invisible(NULL)
}
objective <- function(par) {
  assign_params(par)
  ll <- cm$calculate()
  if(!is.finite(ll)) return(1e15)
  -ll
}
parameters <- c("a","bh","bi","bia","alphaP","rho")
truth_unconstrained <- c(
  log(truth$a),truth$bh,truth$bi,truth$bia,
  truth$alphaP,qlogis(truth$rho)
)
start <- c(log(.075),.1,.3,0,-.7,qlogis(.4))
cat("Synthetic log posterior at generating truth:",
    -objective(truth_unconstrained),"\n")
cat("Optimising six key parameters jointly; nuisance values known...\n")
fit <- optim(
  par=start,fn=objective,method="L-BFGS-B",
  lower=c(log(.008),-1.5,-2,-.6,-3.5,-3.5),
  upper=c(log(.4),1.5,2,.6,2.5,3.5),
  control=list(maxit=180,factr=1e7),
  hessian=TRUE
)
estimate <- unlist(to_model(fit$par))
true <- unlist(truth[parameters])
ses <- rep(NA_real_,length(parameters))
names(ses) <- parameters
if(all(is.finite(fit$hessian))) {
  vc <- tryCatch(solve(fit$hessian),error=function(e) NULL)
  if(!is.null(vc)&&all(diag(vc)>0)) {
    rawSE <- sqrt(diag(vc))
    ses <- c(
      a=estimate["a"]*rawSE[1],
      bh=rawSE[2],bi=rawSE[3],bia=rawSE[4],
      alphaP=rawSE[5],
      rho=estimate["rho"]*(1-estimate["rho"])*rawSE[6]
    )
    names(ses) <- parameters
  }
}
table_out <- data.frame(
  parameter=parameters,
  truth=as.numeric(true),
  estimate=as.numeric(estimate[parameters]),
  approximate_SE=as.numeric(ses),
  abs_error=abs(as.numeric(estimate[parameters])-as.numeric(true)),
  stringsAsFactors=FALSE
)
table_out$error_in_SE <- table_out$abs_error/table_out$approximate_SE
print(table_out,row.names=FALSE,digits=4)
cat("Optimiser convergence:",fit$convergence,"\n")
if(!is.null(fit$message))cat("Optimiser message:",fit$message,"\n")
cat("Synthetic log posterior at mode:",-fit$value,"\n")
cat("Likelihood peaked relative to truth:",
    -fit$value-(-objective(truth_unconstrained)),"\n")
cat("Parameters within 2 approximate SE:",
    sum(table_out$error_in_SE <= 2,na.rm=TRUE),"/",
    length(parameters),"\n")

report <- list(
  table=table_out,truth=truth,fit=fit[c("par","value","convergence","message")],
  n=N,T=T,counts=c(
    recovered_deaths=sum(death>0L),
    unrecovered_deaths=sum(never_observed_dead),
    detections=sum(y==1L),nondetections=sum(y==0L)
  ),
  validations=c(
    interpreted_compiled_match=abs(lpR-lpC)<1e-5,
    optimum_at_least_truth=(fit$value <= objective(truth_unconstrained)+1e-5),
    optimiser_converged=(fit$convergence==0L),
    finite_estimates=all(is.finite(table_out$estimate))
  ),
  limitations=c(
    "Conditional recovery: true individual detection random effects known",
    "Baseline b,c and other nuisance effects fixed at generating truth",
    "One simulated dataset is not a repeated-simulation coverage study",
    "Infection posteriors are per-individual empirical distributions",
    "Does not test shared draw dependence, monitoring misclassification or immigration"
  )
)
saveRDS(report,"outputs/integrated_gm/simulated_recovery_gate.rds")
write.csv(table_out,
  "outputs/integrated_gm/simulated_recovery_parameters.csv",
  row.names=FALSE)
if(fit$convergence!=0L || !all(is.finite(table_out$estimate)))
  stop("SIMULATED RECOVERY OPTIMISATION DID NOT CONVERGE")
cat("\nSIMULATED CONDITIONAL RECOVERY COMPLETED.\n")
cat("Examine truth vs estimates and error_in_SE; do not interpret",
    "a single replicate as proof of statistical identifiability.\n")
cat("This is NOT a full-data MCMC fit.\n")
