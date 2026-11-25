# ============================================================================
# Joint SEI + survival hierarchy, v3d (homozygosity x infection on
# senescence). v3c plus ONE change:
#  7. Infected (E, I) badgers have their own senescence rate
#     bInf_i = bGM*exp((gammaB+deltaB)*hom_i), same age-4 anchor; S badgers
#     keep b_i = bGM*exp(gammaB*hom_i). deltaB ~ N(0,.5) per SD homozygosity.
#     Why: Hudson et al. 2023 found homozygosity steepened senescence
#     credibly only in badgers positive as cubs (about 1.21 per SD vs 1.07
#     in never-positive, read off Fig. 3); v3c (1.06) could not test this.
#     Here infection is the time-varying latent state, not a cub label.
#  Kernel: R/joint_sei_survival_forward_nimble_infsen.R (takes gmInf).
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
  # Change 6: homozygosity on the senescence rate, anchored at age 4.
  gammaB ~ dnorm(0,sd=.5)
  # Change 7: extra homozygosity effect on senescence when infected.
  deltaB ~ dnorm(0,sd=.5)

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
    bI[i] <- bGM*exp(gammaB*hom[i])
    gmI[i,1] <- log(h4*(1.0-f4))-4.0*bI[i]
    gmI[i,2] <- log(bI[i])
    gmI[i,3] <- gm[3]
    bInf[i] <- bGM*exp((gammaB+deltaB)*hom[i])
    gmInfI[i,1] <- log(h4*(1.0-f4))-4.0*bInf[i]
    gmInfI[i,2] <- log(bInf[i])
    gmInfI[i,3] <- gm[3]
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
      gmI[i,1:3],gmInfI[i,1:3],f[i,1:4],
      calSE[1:T],calEI[1:T],
      alphaP,periodEff[1:10],seasonEff[1:4],
      periodIndex[1:T],seasonIndex[1:T],rho)
    zeros[i] ~ dpois(-ll[i])
  }
})
