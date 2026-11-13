# ============================================================
# Shared compiled GM mark-recapture forward likelihood functions
# Used by scripts/54 and scripts/55; keep a single source of truth.
# ============================================================
# Requires library(nimble) in calling script.

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
    periodIndex=double(1), seasonIndex=double(1),
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
          periodEff[as.integer(periodIndex[t])] +
          seasonEff[as.integer(seasonIndex[t])]
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

# Period/season calendar vectors are also propagated as double
# through compiled constant-slice expressions. Both nested functions
# accept them as double(1), with explicit integer indexing on access.
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
    periodIndex=double(1), seasonIndex=double(1),
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

