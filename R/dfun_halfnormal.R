#' Half-normal detection function fitting using INLA.
#' 
#'
#' @aliases dfun_halfnormal
#' @export 
#' @param data
#' @param data A data set containing a list of sightings.
#' @param truncation The maximum observed distance the model will take into account.
#' @param nIntegrate The number of integration points used to determine the denominator
#' of the Cox process.
#' @param dist.trafo Function used to transform the observed distances into coefficients in
#' the model.
#' @return \code{dfun}, detection function object. 
#' @examples \\dontrun{data(whales) ; fun = dfun_halfnormal(whales$sighting); plot(fun)}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'
dfun_halfnormal = function(data,truncation=max(detdata(data)$distance),dist.trafo=function(x){-0.5*x^2},nIntegrate=100,clinear=TRUE)
{
  
require(INLA)

#
# Abbreviations
#
  d = detdata(data)$distance # the distance data
  tr = truncation # truncate distance
  nInt = nIntegrate # Number of integration points
  nDat = length(d) # Number of of sightings (observations)

#
# Obervations (data and fake integration observations)
#

  # actual observations
  yObs = rep(1,nDat) 
  # integration observations
  yInt = rep(0,nInt)

  # concatenate
  fake_data = c(yInt,yObs)


#
# Integration points
#

  ips = seq(from=0,to=tr,length.out=nInt)


#
# Distance covariates concatenated with fake covariates
#

  # fake distances
  fake_dist = dist.trafo(c(ips,d))



#
# Integration weights (fake E)
# 

  # width of integration element
  hi = truncation/(nInt-1)

  # weights for trapeziodal rule
  wtr = hi*c(0.5,rep(1,nInt-2),0.5)

  # fake integration weights
  fake_E = c(wtr,rep(0,nDat))


#
# Define INLA formula
#

if (clinear) { 
  formula = y ~ 1 + f(distance, model="clinear", range = c(0, Inf),hyper=list(beta=list(initial=0, param=c(0, 0.1))))
}
else { 
  formula = y ~ 1 + distance 
}

#
# Run INLA
#

  result = inla(formula,data=list(y = fake_data,distance=fake_dist), family="poisson", E=fake_E,control.compute=list(config = TRUE)) # , verbose=TRUE,debug=TRUE,keep=TRUE

#
# Return object
#
ret = list(INLA=list(result=result,formula=formula),integration.points=ips,
           int.args = list(truncation=truncation),
           dist.trafo=dist.trafo,
           type="halfnormal",
           data=data,
           clinear=clinear)
class(ret) = c("dfun_halfnormal","dfun")
return (ret)
}

#' Sample a detection function at given distances
#'
#' @aliases sample.value.dfun_halfnormal
#' @export
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'
sample.value.dfun_halfnormal = function(mdl,data=FALSE,n=1){
  return(lapply(dfun_halfnormal.sample.logvalue(mdl,n=n,data=data),exp))
}


#' Sample a log detection function at given distances
#' minla(list(blah=list(1,5,8)),test=1,test2="q")
#' @aliases sample.logvalue.dfun_halfnormal
#' @export
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'
sample.logvalue.dfun_halfnormal = function(mdl,data=FALSE,n=1){
  smp = inla.posterior.sample.structured(mdl$INLA$result,n)
  if (!is.data.frame(data)) {
    sampled.logvalues = lapply(smp,function(x) mdl$dist.trafo(mdl$data$distance)*x$distance)
  } else {
    sampled.logvalues = lapply(smp,function(x) mdl$dist.trafo(data$distance)*x$distance)
  }
  return(t(do.call(rbind,sampled.logvalues)))
}

