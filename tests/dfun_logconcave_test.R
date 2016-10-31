#
# Test Data
#
  nSegments = 4
  truncation = 4 # truncation
  ns = 400
  nIntegrate = ns
  
#
# Sampling using the basis we are going to use
#
  source("dfun_logconcave.R")
  coeffs = list(nlbasis=c(0,0,1,1),lbasis=0.0,intercept=1)
  xbasis = seq(0,truncation,0.01)
  basis = dfun_logconcave.basis.lincomb(xbasis,nSegments,truncation,coeffs=coeffs)
  plot(basis)
  samples = sample(xbasis,ns,prob=exp(-basis),replace=TRUE)
  hist(samples,breaks=seq(0,truncation,0.5),xlim=c(0,truncation))
  par(new=TRUE)
  plot(exp(-basis))
  data = data.frame(distance=samples)
  #data = data.frame(distance=seq(0,truncation,0.01))

#
# Run log concave detection funciton estimation
#

  source("dfun_logconcave.R")
  lMod=dfun_logconcave(data,truncation=truncation,
                       nIntegrate=nIntegrate,nSegments=nSegments,clinear=TRUE)
  
#
# Plot results
#   
    
  source("dfun_plot.R")
  plot(lMod,truth=exp(-basis))
  
    
#
# Other test data
#
  
  # Whales
  source("io_whales.R")
  source("dfun_logconcave.R")
  whales = io_whales.pkgdata()
  lcwMod=dfun_logconcave(whales,nSegments=4)
  source("inla_essmod.R")
  source("plot.idfmodel.R")
  plot(lcwMod)
  
  # Whales without linear basis
  lcwMod.nolin=dfun_logconcave(whales,nSegments=4,linbasis=FALSE)
  plot(lcwMod.nolin)
  
#
# Sampling from the whales model
#
  
  x = seq(0,lcwMod$truncation,0.01)
  
  # sample detection functions
  smp = dfun_logconcave.sample.value(lcwMod,x,5)
  smp.matrix = do.call(cbind,smp)
  matplot(x,smp.matrix)
  
  # sample in the log domain
  lsmp = dfun_logconcave.sample.logvalue(lcwMod,x,5)
  lsmp.matrix = do.call(cbind,lsmp)
  matplot(x,lsmp.matrix)
  
  # sample population
  obs.samples = dfun_logconcave.sample.population(mdl,n=2)
  hist(obs.samples[[1]])
