#
# Whales data
#

  source("io_whales.R")
  whales = io_whales.pkgdata()

#
# Infer detection function
#
  source("dfun_spde.R")
  swMod=dfun_spde1d(whales,nIntegrate=20,nKnots=40)

#
# Inspect INLA result
#

  result = swMod$INLA$result
  plot(result)
  lines(result$summary.fitted.values$mean[index])
  
source("inla_essmod.R")
source("dfun_plot.R")
