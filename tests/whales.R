data(whales)

#
# Inference with default settings
#

  wpp = pproc_lgcp(whales)

  summary(wpp)
  plot.dfun_halfnormal(wpp)
  plot(wpp)

  # The actual result of the INLA run can be found here:
  summary(wpp$INLA$result)

  #plot(wpp, rgl = FALSE) # posterior mode of the spde
  #plot(wpp, rgl = FALSE, logscale = FALSE) # posterior mode of the spde

  # plot.marginal(wpp,"Intercept") # bug!


#
# Inference with default settings written out / inserted
#
#
#   wpp = pproc_lgcp(whales,
#             int.points = whales$ips,
#             int.scheme = NULL,
#             int.args = list(truncation=max(detdata(whales)$distance)),
#             dist.trafo = function(dst) { return(-0.5*dst^2) },
#             covariates = list(),
#             formula = y ~ -1
#               + Intercept
#               + f(distance, model="clinear", range = c(0, Inf), hyper = list(beta=list(initial=0, param=c(0,0.1))))
#               + f(spde, model=spde.mdl),
#             spde.args = list(alpha = 2,
#                           prior.variance.nominal = 10,
#                           theta.prior.prec = 0.01),
#             predict = list(),
#             constant = list(),
#             inla.args = list(),
#             sgh.E = 0)

#
# Predicting combined effects while using constrained spde
#

  wpp = pproc_lgcp(whales,
            spde.args = list(alpha = 2,
                          prior.variance.nominal = 10,
                          theta.prior.prec = 0.01,
                          constr = TRUE),
            predict = list(my.tag = list(vars = c("Intercept","spde"))))

  summary(wpp)
  # Plot combined effect
  plot(wpp,field="my.tag",rgl=FALSE,property="mean") # log scale
  plot(wpp,field="my.tag",rgl=FALSE,logscale = FALSE, property="mean") # natural scale


#
# Passing arguments to INLA, e.g.: obtain correlation matrix
#

  wpp = pproc_lgcp(whales,
                   inla.args = list( control.fixed = list(correlation.matrix = TRUE),
                                     control.inla = list(lincomb.derived.correlation.matrix = TRUE)))
  result = wpp$INLA$result
  corrmat = cov2cor(result$misc$lincomb.derived.covariance.matrix)


#
# Changing the formula, e.g.: adding a covariate
#

  wpp = pproc_lgcp(whales, formula = ~ . + lat + lon )
  summary(wpp)

#
# More complicated covariates (centered coordinates)
#

  clonFun = function(x) {return( x$lon - mean(range(whales$ips$lon)))}
  clatFun = function(x) {return( x$lat - mean(range(whales$ips$lat)))}

  wpp = pproc_lgcp(whales,
                   formula = ~ . + clat + clon,
                   covariates = list(clat = clatFun, clon = clonFun))

  summary(wpp)


#
# More covariates (SST)
#
# Note: we have to use the yearly integration points here
#

  data(sst) # SST averaged for each year and interpolated to the mesh
  sstFun1 = function(x) { return(get.tmean(sst,x) - get.mean(sst,x)) }
  sstFun2 = function(x) { return(get.smean(sst,x) - get.mean(sst,x)) }
  sstFun3 = function(x) { return( ( get.value(sst,x) - get.smean(sst,x) ) - (  get.tmean(sst,x) - get.mean(sst,x) )) }

  wpp = pproc_lgcp(whales,
                   int.points = whales$ips.yearly,
                   formula = y ~ -1
                   + Intercept
                   + sst1
                   + sst2
                   + sst3
                   + f(distance, model="clinear", range = c(0, Inf), hyper = list(beta=list(initial=0, param=c(0,0.1))))
                   + f(spde, model=spde.mdl),
                   covariates = list(sst1 = sstFun1, sst2 = sstFun2, sst3 = sstFun3),
                   spde.args = list(alpha = 2,
                                    prior.variance.nominal = 10,
                                    theta.prior.prec = 0.01,
                                    constr = TRUE))


#
# More covariates (SST), with averages performed over survey area and proper constraints to the SPDE
#
# Note: we have to use the yearly integration points here
#

  data(sst)

  load("innerBnd.Rdata")
  library(sp)
  sst.mask <- point.in.polygon(sst$mesh$loc[,1], sst$mesh$loc[,2],
                             innerBnd[,1], innerBnd[,2]) > 0

  sstFun1 = function(x) { return(get.tmean(sst,x) - get.mean(sst,x,weights = sst.mask)) }
  sstFun2 = function(x) { return(get.smean(sst,x, weights = sst.mask) - get.mean(sst,x,weights = sst.mask)) }
  sstFun3 = function(x) { return( ( get.value(sst,x) - get.smean(sst,x,weights = sst.mask) ) - (  get.tmean(sst,x) - get.mean(sst,x,weights = sst.mask) )) }
  clonFun = function(x) { return( x$lon - mean(range(whales$ips$lon))) }
  clatFun = function(x) { return( x$lat - mean(range(whales$ips$lat))) }

  wpp2 <- pproc_lgcp(whales,
                     int.points = whales$ips.yearly,
                     formula = (y ~ -1
                                + Intercept
                                + clon
                                + clat
                                + sst1
                                + sst2
                                + sst3
                                + f(distance, model="clinear",
                                    range = c(0, Inf),
                                    hyper = list(beta=list(initial=0,
                                                 param=c(0,0.1))))
                                + f(spde, model=spde.mdl)),
                     covariates =
                       list(sst1 = sstFun1,
                            sst2 = sstFun2,
                            sst3 = sstFun3,
                            clon = clonFun,
                            clat = clatFun),
                     spde.args = list(alpha = 2,
                     prior.variance.nominal = 10,
                     theta.prior.prec = 0.01,
                     constr = TRUE))

  wpp3 <- pproc_lgcp(whales,
                     int.points = whales$ips.yearly,
                     formula = (y ~ -1
                                + Intercept
                                + clon
                                + clat
                                + sst1
                                + sst2
                                + sst3
                                + f(distance, model="clinear",
                                    range = c(0, Inf),
                                    hyper = list(beta=list(initial=0,
                                                 param=c(0,0.1))))
                                + f(spde, model=spde.mdl)),
                     covariates =
                       list(sst1 = sstFun1,
                            sst2 = sstFun2,
                            sst3 = sstFun3,
                            clon = clonFun,
                            clat = clatFun),
                     spde.args = list(alpha = 2,
                     prior.variance.nominal = 10,
                     theta.prior.prec = 0.01,
                     extraconstr.int =
                       list(A=matrix(sst.mask, 1, length(sst.mask)), e=0)))

#
# Latitude and Longitude as covariate data structures
#
  data(sst)
  lon = sst; lon$values = data.frame(matrix(whales$mesh$loc[,1],nrow =whales$mesh$n,ncol=22))
  colnames(lon$values) = colnames(sst$values)
  lat = sst; lat$values = data.frame(matrix(whales$mesh$loc[,2],nrow =whales$mesh$n,ncol=22))
  colnames(lat$values) = colnames(sst$values)


  lonFun = function(x) { get.value(lon,x) - get.smean(lon,x) }
  latFun = function(x) { get.value(lat,x) - get.smean(lat,x) }


#
# Model with lat,lon,non-centered SST, no SPDE
#

  sstFun.nocenter = function(x) { return(get.value(sst,x)) }

  wpp = pproc_lgcp(whales,
                   int.points = whales$ips.yearly,
                   formula = y ~ -1
                             + lon
                             + lat
                             + sst
                             + Intercept
                             + f(distance, model="clinear", range = c(0, Inf), hyper = list(beta=list(initial=0, param=c(0,0.1)))),
                   covariates = list(sst = sstFun.nocenter))

    summary(wpp)



#
# Full model but without SPDE
#
# Effects: Intercept, clon, clat, sst1, sst2, sst3, constraint distance
#

  sstFun1 = function(x) { return(get.tmean(sst,x) - get.mean(sst,x)) }
  sstFun2 = function(x) { return(get.smean(sst,x) - get.mean(sst,x)) }
  sstFun3 = function(x) { return( ( get.value(sst,x) - get.smean(sst,x) ) - (  get.tmean(sst,x) - get.mean(sst,x) )) }
  clonFun = function(x) { return( x$lon - mean(range(whales$ips$lon))) }
  clatFun = function(x) { return( x$lat - mean(range(whales$ips$lat))) }

  wpp = pproc_lgcp(whales,
                   int.points = whales$ips.yearly,
                   formula = y ~ -1
                   + clon + clat + sst1 + sst2 + sst3 + Intercept
                   + f(distance, model="clinear", range = c(0, Inf), hyper = list(beta=list(initial=0, param=c(0,0.1)))),
                   covariates = list(clon = clonFun, clat = clatFun,
                                     sst1 = sstFun1, sst2 = sstFun2, sst3 = sstFun3))

  summary(wpp)

  # ISSUE: Estimate for distance effect makes no sense
  # BUG 2: Prediction not working, e.g. predict = list(my.tag = list(vars = c("Intercept","clon","clat","sst1","sst2","sst3")))


#
# Full model with SPDE
#
# Effects: Intercept, clon, clat, sst1, sst2, sst3, constraint distance, SPDE
#

  sstFun1 = function(x) { return(get.tmean(sst,x) - get.mean(sst,x)) }
  sstFun2 = function(x) { return(get.smean(sst,x) - get.mean(sst,x)) }
  sstFun3 = function(x) { return( ( get.value(sst,x) - get.smean(sst,x) ) - (  get.tmean(sst,x) - get.mean(sst,x) )) }
  clonFun = function(x) { return( x$lon - mean(range(whales$ips$lon))) }
  clatFun = function(x) { return( x$lat - mean(range(whales$ips$lat))) }

  wpp = pproc_lgcp(whales,
                   int.points = whales$ips.yearly,
                   formula = y ~ -1
                   + clon + clat + sst1 + sst2 + sst3 + Intercept
                   + f(distance, model="clinear", range = c(0, Inf), hyper = list(beta=list(initial=0, param=c(0,0.1))))
                   + f(spde, model=spde.mdl),
                   covariates = list(clon = clonFun, clat = clatFun,
                                     sst1 = sstFun1, sst2 = sstFun2, sst3 = sstFun3))

  summary(wpp)
