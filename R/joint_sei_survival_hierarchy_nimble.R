# ============================================================================
# Bayesian hierarchy for the SINGLE four-state joint SEI + survival model.
#
# This is the exact intended scientific model. MODE=hierarchical_compile
# checks the SAME model with a small representative individual subset.
# MODE=hierarchical_benchmark and eventual full MCMC use the SAME nimbleCode.
#
# Persistent individual latent frailties:
#   S->E:   f[i,1] = sigmaSE*z[i,1]
#   E->I:   f[i,2] = load21*z[i,1] + sigmaEI*z[i,2]
#   Death:  f[i,3] = load31*z[i,1] + load32*z[i,2]+sigmaD*z[i,3]
#   Capture:f[i,4] = sigmaCapture*z[i,4]
# The lower-triangular loadings imply a positive-definite covariance
# among the first three frailties, and permit correlated latent risk.
# Detection heterogeneity remains separate. Priors on loadings
# regularise off-diagonal covariance toward zero.
#
# The fixed test accuracy is already encoded in precomputed assay
# log likelihoods.  A non-capture y=0 indicates global/group
# monitoring proxy, NOT measured individual trap exposure.
#
# Conditional upon an observed first live capture and known genotype.
# ============================================================================

joint_sei_survival_model_code <- nimbleCode({
  # An anchored Gompertz-Makeham hazard: h(age)=a*exp(bGM*age)+c.
  h4 ~ dlnorm(meanlog=log(.25),sdlog=1.0)
  f4 ~ dbeta(2,2)
  bGM ~ dexp(5)
  aGM <- h4*(1.0-f4)*exp(-4.0*bGM)
  cGM <- h4*f4
  gm[1] <- log(aGM)
  gm[2] <- log(bGM)
  gm[3] <- log(cGM)

  # S->E acquisition and E->I progression coefficients:
  # intercept, age/4y, male, homozygosity per 1 full-cohort SD.
  betaSE[1] ~ dnorm(-3.3,sd=1.5)
  betaEI[1] ~ dnorm(-3.0,sd=1.5)
  for(j in 2:4) {
    betaSE[j] ~ dnorm(0,sd=.8)
    betaEI[j] ~ dnorm(0,sd=.8)
  }
  # Survival log-hazard multiplier: male, hom, E vs S, I vs S.
  for(j in 1:4) betaD[j] ~ dnorm(0,sd=.8)

  # Initial infection prevalence (at first observed live capture):
  # intercept, age/4y, male, hom, calendar/10y.
  entryBeta[1] ~ dnorm(-2.2,sd=1.5)
  for(j in 2:5) entryBeta[j] ~ dnorm(0,sd=.8)
  entryRho ~ dbeta(2,2)

  # Smooth, regularized 5-year disease-period pressure.
  # 1976-1980 is reference level; the intercept anchors risk there.
  sigmaYearSE ~ dunif(.01,1.5)
  sigmaYearEI ~ dunif(.01,1.5)
  yearSE[1] <- 0
  yearEI[1] <- 0
  for(j in 2:10) {
    yearSE[j] ~ dnorm(yearSE[j-1],sd=sigmaYearSE)
    yearEI[j] ~ dnorm(yearEI[j-1],sd=sigmaYearEI)
  }
  for(t in 1:T) {
    calSE[t] <- yearSE[yearIndex[t]]
    calEI[t] <- yearEI[yearIndex[t]]
  }

  # Individual heterogeneity, separate from sex/age/genotype fixed effects.
  # Noncentered triangular covariance parameterization.
  sigmaSE ~ dunif(.01,1.5)
  sigmaEI ~ dunif(.01,1.5)
  sigmaD ~ dunif(.01,1.5)
  load21 ~ dnorm(0,sd=.35)
  load31 ~ dnorm(0,sd=.35)
  load32 ~ dnorm(0,sd=.35)
  sigmaCapture ~ dunif(.01,2.0)

  # Recapture detection and dead recovery.
  alphaP ~ dnorm(-1.5,sd=1)
  rho ~ dbeta(2,3)
  periodEff[1] <- 0
  seasonEff[1] <- 0
  for(j in 2:10) periodEff[j] ~ dnorm(0,sd=.4)
  for(j in 2:4) seasonEff[j] ~ dnorm(0,sd=.4)

  for(i in 1:N){
    for(j in 1:4) z[i,j] ~ dnorm(0,sd=1)
    f[i,1] <- sigmaSE*z[i,1]
    f[i,2] <- load21*z[i,1]+sigmaEI*z[i,2]
    f[i,3] <- load31*z[i,1]+load32*z[i,2]+sigmaD*z[i,3]
    f[i,4] <- sigmaCapture*z[i,4]

    entryAge4[i] <- (first[i]-birth[i])*.25/4
    entryCalendar10[i] <- (first[i]-1.0)/40.0
    piEntry[i] <- 1.0/(1.0+exp(-(
      entryBeta[1]+entryBeta[2]*entryAge4[i]+
      entryBeta[3]*male[i]+entryBeta[4]*hom[i]+
      entryBeta[5]*entryCalendar10[i])))
    entry[i,1] <- 1.0-piEntry[i]
    entry[i,2] <- piEntry[i]*(1.0-entryRho)
    entry[i,3] <- piEntry[i]*entryRho

    ll[i] <- joint_sei_survival_ll(
      y[i,1:T],first[i],death[i],birth[i],
      em[i,1:T,1:3],hom[i],male[i],
      entry[i,1:3],
      betaSE[1:4],betaEI[1:4],betaD[1:4],
      gm[1:3],f[i,1:4],
      calSE[1:T],calEI[1:T],
      alphaP,periodEff[1:10],seasonEff[1:4],
      periodIndex[1:T],seasonIndex[1:T],rho)
    # Poisson-zero trick: Pr(0|-ll) = exp(ll), for ll<=0.
    zeros[i] ~ dpois(-ll[i])
  }
})
