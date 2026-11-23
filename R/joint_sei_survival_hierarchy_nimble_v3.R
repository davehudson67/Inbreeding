# ============================================================================
# Joint SEI + survival hierarchy, v3 (branch joint-v3-detcentred).
# v2 (R/joint_sei_survival_hierarchy_nimble_v2.R) plus ONE change:
#  4. Detection period effects centred: periodEff[j]=periodRaw[j]-mean(raw),
#     j=1..10; alphaP is now average (logit) detection over 1976-2025.
#     Why: v2 anchored periodEff[1]=0 (1976-80), giving an alphaP-periodEff
#     ridge (R-hat 1.04, bulk ESS 75-110 for alphaP, periodEff[3:7]) - the
#     same corner-anchor problem v2 change 1 fixed for the calendar RW.
# v2 changes kept: centred calendar RW, half-normal scales, no mortality
# frailty. Season effects stay anchored (not flagged by v2 diagnostics).
# ============================================================================

joint_sei_survival_model_code <- nimbleCode({
  # Anchored Gompertz-Makeham hazard: h(age)=a*exp(bGM*age)+c.
  h4 ~ dlnorm(meanlog=log(.25),sdlog=1.0)
  f4 ~ dbeta(2,2)
  bGM ~ dexp(5)
  aGM <- h4*(1.0-f4)*exp(-4.0*bGM)
  cGM <- h4*f4
  gm[1] <- log(aGM)
  gm[2] <- log(bGM)
  gm[3] <- log(cGM)

  # Intercept (long-run average), age/4y, male, homozygosity (1 SD).
  betaSE[1] ~ dnorm(-3.3,sd=1.5)
  betaEI[1] ~ dnorm(-3.0,sd=1.5)
  for(j in 2:4) {
    betaSE[j] ~ dnorm(0,sd=.8)
    betaEI[j] ~ dnorm(0,sd=.8)
  }
  # Survival log-hazard multiplier: male, hom, E vs S, I vs S.
  for(j in 1:4) betaD[j] ~ dnorm(0,sd=.8)

  # Entry prevalence: intercept, age/4y, male, hom, calendar/10y.
  entryBeta[1] ~ dnorm(-2.2,sd=1.5)
  for(j in 2:5) entryBeta[j] ~ dnorm(0,sd=.8)
  entryRho ~ dbeta(2,2)

  # Change 1 + 2: centred non-centred RW, half-normal step SDs.
  sigmaYearSE ~ T(dnorm(0,sd=1),0,)
  sigmaYearEI ~ T(dnorm(0,sd=1),0,)
  rwSE[1] <- 0
  rwEI[1] <- 0
  for(j in 2:10) {
    zYearSE[j] ~ dnorm(0,sd=1)
    zYearEI[j] ~ dnorm(0,sd=1)
    rwSE[j] <- rwSE[j-1]+sigmaYearSE*zYearSE[j]
    rwEI[j] <- rwEI[j-1]+sigmaYearEI*zYearEI[j]
  }
  rwMeanSE <- mean(rwSE[1:10])
  rwMeanEI <- mean(rwEI[1:10])
  for(j in 1:10) {
    yearSE[j] <- rwSE[j]-rwMeanSE
    yearEI[j] <- rwEI[j]-rwMeanEI
  }
  for(t in 1:T) {
    calSE[t] <- yearSE[yearIndex[t]]
    calEI[t] <- yearEI[yearIndex[t]]
  }

  # Change 2 + 3: half-normal frailty scales; no mortality frailty.
  sigmaSE ~ T(dnorm(0,sd=1),0,)
  sigmaEI ~ T(dnorm(0,sd=1),0,)
  load21 ~ dnorm(0,sd=.35)
  sigmaCapture ~ T(dnorm(0,sd=1),0,)

  # Recapture detection and dead recovery (unchanged).
  alphaP ~ dnorm(-1.5,sd=1)
  rho ~ dbeta(2,3)
  # Change 4: centred detection period effects.
  for(j in 1:10) periodRaw[j] ~ dnorm(0,sd=.4)
  periodMean <- mean(periodRaw[1:10])
  for(j in 1:10) periodEff[j] <- periodRaw[j]-periodMean
  seasonEff[1] <- 0
  for(j in 2:4) seasonEff[j] ~ dnorm(0,sd=.4)

  for(i in 1:N){
    z[i,1] ~ dnorm(0,sd=1)
    z[i,2] ~ dnorm(0,sd=1)
    z[i,4] ~ dnorm(0,sd=1)
    f[i,1] <- sigmaSE*z[i,1]
    f[i,2] <- load21*z[i,1]+sigmaEI*z[i,2]
    f[i,3] <- 0
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
    zeros[i] ~ dpois(-ll[i])
  }
})
