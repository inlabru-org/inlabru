
#
# Load whales data
#

  data(whales)
  plot(whales)

#
# A simple point process model (LGCP)
#

  # Reload whales data from scratch (if package version is outdated)
  # whales = io_whales.pkgdata.load()

  fpp =  pproc_lgcp(whales,int.scheme=cylindric.integration)
  plot(fpp$INLA$result) # plot INLA result
  plot(fpp) # plot SPDE model result (mean)
  plot(fpp,field="spde",property="sd") # plot SPDE model result (sd)
  plot.dfun_halfnormal(fpp) # Plot resulting detection function
  
  # with sst and explicit parameters
  fpp.sst =  pproc_lgcp(whales,
                 int.scheme=cylindric.integration,
                 int.args = list(truncation=6),
                 covariates=list(sst=yearly.mean.sst),
                 spde.args = list(alpha=2,prior.variance.nominal=10,theta.prior.prec=0.01),
                 formula = y ~ -1  + lat + lon + sst + Intercept +
                   f(distance, model="clinear", range = c(0, Inf),
                     hyper=list(beta=list(initial=0, param=c(0,0.1)))) +
                   f(spde, model=spde.mdl)
          )
  plot(fpp.sst$INLA$result) # plot INLA result
  plot(fpp.sst) # plot SPDE model result



#
# Toy data set point process
#
  
  data(toy1)
  plot(toy1)
  toy.pp =  pproc_lgcp(toy1,int.scheme=grid.integration,
                       int.args = list(n.length=50,n.distance=20,truncation=3),
                       formula = y ~ -1 + Intercept + f(distance, model="clinear", range = c(0, Inf), hyper=list(beta=list(initial=0, param=c(0,0.1)))) + 
                       f(spde, model=spde.mdl))
  plot(toy.pp,rgl=FALSE)
  plot.dfun_halfnormal(toy.pp)
  summary(toy.pp$INLA$result)
  plot(toy.pp$INLA$result)
  
#
# Detection function inference
#

  # Half-normal
  hn = dfun("halfnormal",whales)
  summary(hn)
  plot(hn)
  # Log concave
  lc = dfun("logconcave",whales)
  summary(lc)
  plot(lc)
  # SPDE
  sp = dfun("spde1d",whales)
  plot(sp)

