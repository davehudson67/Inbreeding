# REMOVAL VARIANT of R/joint_sei_survival_forward_nimble.R (v4 only).
# Only change: death<0 codes a removal (euthanised/killed) at quarter
# -death. The badger is alive at the start of that quarter; the quarter's
# observation counts, then the history is censored (no death, no
# non-recovery terms afterwards).
# ========================================================================
# Single intended joint model: compiled four-state SEI+survival likelihood
#
# One individual, conditional on FIRST observed live capture:
# - live capture and tests occur at quarter start;
# - quarterly S/E/I mortality + S->E/E->I transitions follow;
# - an observed recovered death in the SAME quarter is permitted;
# - an unobserved death stays in the absorbing D state;
# - non-capture y=0 differs from no effort y=-1;
# - fixed per-assay diagnostic log emissions are precomputed in the
#   source-fingerprinted input and are summed for repeated captures.
#
# This is a *likelihood kernel*, NOT an alternative model. The eventual
# hierarchical model supplies shared regression parameters, frailties,
# smoothed calendar baselines and capture/detection random effects.
#
# Requires library(nimble) in caller; compile with compileNimble().
# See docs/one_joint_infection_survival_model_target.md.
# ========================================================================
joint_sei_survival_ll <- nimbleFunction(
  run=function(
    y=double(1),first=integer(0),death=integer(0),birth=integer(0),
    em=double(2),hom=double(0),male=double(0),
    entry=double(1),
    betaSE=double(1),betaEI=double(1),betaD=double(1),
    gm=double(1),frailty=double(1),
    calSE=double(1),calEI=double(1),
    alphaP=double(0),periodEff=double(1),seasonEff=double(1),
    periodIndex=double(1),seasonIndex=double(1),
    rho=double(0)) {
    returnType(double(0))
    # Parameter indices:
    # betaSE/EI: [intercept, age/4, male, hom]
    # betaD: [male, hom, E effect, I effect]
    # gm: [log a, log b, log c], h(age)=a*exp(b*age)+c per year
    # frailty: [S->E, E->I, mortality, capture]
    ns <- length(y)
    aS <- entry[1]
    aE <- entry[2]
    aI <- entry[3]
    aD <- 0.0
    LL <- 0.0
    A <- exp(gm[1])
    B <- exp(gm[2])
    C <- exp(gm[3])
    z <- B*0.25
    ez <- 1.0
    if(abs(z)<1.0e-6) {
      ez <- 1.0+z/2.0+z*z/6.0+z*z*z/24.0
    } else {
      ez <- (exp(z)-1.0)/z
    }
    for(t in first:ns) {
      if(t>first) {
        lp <- alphaP+frailty[4]+
          periodEff[as.integer(periodIndex[t])]+
          seasonEff[as.integer(seasonIndex[t])]
        p <- 1.0/(1.0+exp(-lp))
        if(y[t]==1) {
          aS <- aS*p
          aE <- aE*p
          aI <- aI*p
          aD <- 0.0
        } else if(y[t]==0) {
          aS <- aS*(1.0-p)
          aE <- aE*(1.0-p)
          aI <- aI*(1.0-p)
        }
      }
      # A diagnostic emission is permitted only at a confirmed live
      # capture. First-capture detection probability is conditioned out.
      off <- 0.0
      if(y[t]==1) {
        off <- em[t,1]
        if(em[t,2]>off)off<-em[t,2]
        if(em[t,3]>off)off<-em[t,3]
        aS <- aS*exp(em[t,1]-off)
        aE <- aE*exp(em[t,2]-off)
        aI <- aI*exp(em[t,3]-off)
        aD <- 0.0
      }
      scale <- aS+aE+aI+aD
      if(scale<=0.0)return(-1.0e12)
      LL <- LL+log(scale)+off
      if(death<0) {
        if(t==-death) return(LL)
      }
      aS<-aS/scale
      aE<-aE/scale
      aI<-aI/scale
      aD<-aD/scale

      age <- 0.25*(t-birth)
      h <- 0.25*(A*exp(B*age)*ez+C)
      deathLP <- betaD[1]*male+betaD[2]*hom+frailty[3]
      phiS <- exp(-h*exp(deathLP))
      phiE <- exp(-h*exp(deathLP+betaD[3]))
      phiI <- exp(-h*exp(deathLP+betaD[4]))
      lpSE <- betaSE[1]+betaSE[2]*(age/4.0)+
              betaSE[3]*male+betaSE[4]*hom+
              frailty[1]+calSE[t]
      lpEI <- betaEI[1]+betaEI[2]*(age/4.0)+
              betaEI[3]*male+betaEI[4]*hom+
              frailty[2]+calEI[t]
      lambda <- 1.0/(1.0+exp(-lpSE))
      tau <- 1.0/(1.0+exp(-lpEI))
      dying <- aS*(1.0-phiS)+aE*(1.0-phiE)+
               aI*(1.0-phiI)
      if(t==death) {
        prob <- rho*dying
        if(prob<=0.0)return(-1.0e12)
        # All future occasions after an observed death have y=-1.
        return(LL+log(prob))
      }
      # No recovered death this quarter.
      nS <- aS*phiS*(1.0-lambda)
      nE <- aS*phiS*lambda+aE*phiE*(1.0-tau)
      nI <- aE*phiE*tau+aI*phiI
      nD <- aD+(1.0-rho)*dying
      scale <- nS+nE+nI+nD
      if(scale<=0.0)return(-1.0e12)
      LL <- LL+log(scale)
      aS<-nS/scale
      aE<-nE/scale
      aI<-nI/scale
      aD<-nD/scale
    }
    return(LL)
  }
)
