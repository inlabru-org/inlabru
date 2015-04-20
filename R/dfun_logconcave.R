#' Log concave detection function fitting using INLA.
#' 
#'
#' @aliases dfun_logconcave
#' @export 
#' @param data
#' @param data A data set containing a list of sightings.
#' @param truncation The maximum observed distance the model will take into account.
#' @param nIntegrate The number of integration points used to determine the denominator
#' of the Cox process.
#' @param nSegments Number of segments the distance domain is separated into, i.e.
#' the number of log concave basis functions modeling the log detection probability.
#' @param clinear Boolean determining wether a linear basis function should be used.
#' @param linbasis Boolean determining if the basis coefficients are to be modeled by
#' the INLA \code{clinear} or \code{linear} model. If \code{clinear=FALSE} the
#' detection function is not constrained to log concavity.
#' @return \code{dfun}, detection function object. 
#' @examples \\dontrun{data(whales) ; fun = dfun_logconcave(whales$sighting); plot(fun)}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

dfun_logconcave = function(data,truncation=max(detdata(data)$distance),nIntegrate=100,segments=5,clinear=TRUE,linbasis=TRUE)
{
  
  require(INLA)
  

#
# Abbreviations
#
d = detdata(data)$distance # the distance data
tr = truncation # truncate distance
nInt = nIntegrate # Number of integration points
nDat = length(d) # Number of of sightings (observations)
nSeg = segments

#
# Obervations (data and fake integration observations)
#

# actual observations
yObs = rep(1,nDat) 
# integration observations
yInt = rep(0,nInt)

# concatenate
if (nInt>0){
  fake_data = c(yInt,yObs)
}
else {
  fake_data = yObs
}


#
# Integration points
#

ips = seq(from=0,to=tr,length.out=nInt)


#
# Basis function evaluated at distance covariates and integration points
#
if (nInt>0){
  bweights = dfun_logconcave.basis.value(c(ips,d),nSeg,truncation)
} else {
  bweights = dfun_logconcave.basis.value(d,nSeg,truncation)
}
#plot(c(ips,d),bweights[[1]])

#
# Integration weights (fake E)
# 

if (nInt>0){
  # width of integration element
  hi = truncation/(nInt-1)
  # weights for trapeziodal rule
  wtr = hi*c(0.5,rep(1,nInt-2),0.5)
  fake_E = c(wtr,rep(0,nDat))
} else {
  fake_E = rep(0,nDat)
}

#
# Beta prior distribution (the default normal is too low for beta=0)
#

  loggamma = "expression:
  a = 1;
  b = 0.001;
  precision = exp(log_precision);
  logdens = log(b^a) - lgamma(a)
  + (a-1)*log_precision - b*precision;
  log_jacobian = log_precision;
  return(logdens + log_jacobian);"
  hyper.new = list(theta = list(prior = loggamma))


#
# Define INLA formula
#
if (clinear) {
 params = ",model = 'clinear',hyper=list(theta=list(prior=loggamma)),range = c(-10,0)" # 
 if (linbasis){
   fml = paste0("y ~ 1 + f(linbasis",params,")+",paste0("f(basis_",1:nSeg,params,")",collapse="+"))
 } else {
   fml = paste0("y ~ 1 +",paste0("f(basis_",1:nSeg,params,")",collapse="+"))
 }
 
} else {
  if (linbasis){
    fml = paste("y ~ 1 + linbasis",paste(paste(rep(" + basis_",nSeg),seq(1:nSeg),sep=""),collapse=''),sep="")  
  } else {
    fml = paste("y ~ 1",paste(paste(rep(" + basis_",nSeg),seq(1:nSeg),sep=""),collapse=''),sep="")  
  }
}
formula = as.formula(fml)


#
# Construct INLA data frame
#

data.INLA = data.frame(bweights)
data.INLA$y = fake_data
if (linbasis) {
  if (nInt>0){
    data.INLA$linbasis = c(ips,d)
  } else {
    data.INLA$linbasis = d
  }
}

#
# Run INLA
#

result = inla(formula,data=data.INLA,
              family="poisson",
              E=fake_E,control.compute=list(config = TRUE)) # , verbose=TRUE,debug=TRUE,keep=TRUE

#
# Return object
#
ret = list(INLA=list(result=result,formula=formula),
           integration.points=ips,
           truncation=truncation,
           int.args = list(truncation=truncation),
           type="logconcave",
           data=data,
           nSegments=nSeg,
           nData.after.filtering = nDat, # Number of data points actually used to perform the inference 
           clinear=clinear,
           linbasis=linbasis)

class(ret) = c("dfun_logconcave","dfun")
return (ret)
}


dfun_logconcave.basis.value = function(d,nSeg,truncation){
  pts = seq(0,truncation,length.out=nSeg+1) # segment boundary points
  wdh = pts[2]-pts[1] # segment width
  lidx = findInterval(d,pts) # for each data point the index of left segment boundary
  offs = d-pts[lidx] # local offsets: (x-x_i) for data point x and left boundary neighbor
  w = list()
  for (j in 1:nSeg) { 
    w[[paste("basis_",j,sep="")]] =  (pmax(0,lidx-j-1)*2*wdh*wdh + (j<lidx)*(2*wdh)*(offs) + (j<lidx)*(wdh^2) + (j==lidx)*offs^2)
    #w[[paste("basis_",j,sep="")]] = w[[paste("basis_",j,sep="")]]/max(w[[paste("basis_",j,sep="")]])
    #mx = (pmax(0,nSeg-j-1)*2*wdh*wdh + (j<nSeg)*(2*wdh)*(wdh) + (j<nSeg)*(wdh^2) + (j==nSeg)*wdh^2)
    #w[[paste("basis_",j,sep="")]] =  w[[paste("basis_",j,sep="")]]/mx
  } 
  return(w)
}

dfun_logconcave.logvalue = function(mdl,data=FALSE,field="mode"){
  if (is.vector(data)) {distances = data}
  else if (is.data.frame(data)) {
    distances = data$distance
  } else {
    distances = mdl$data$distance
  }
  return(dfun_logconcave.basis.lincomb(distances,mdl$nSegments,mdl$truncation,mdl=mdl,field=field))
}

dfun_logconcave.value = function(...){
  return (exp(dfun_logconcave.logvalue(...)))
}

dfun_logconcave.basis.lincomb = function(d,nSeg,truncation,coeffs=FALSE,mdl=FALSE,field="mode"){
  bweights = dfun_logconcave.basis.value(d,nSeg,truncation)
  # If we are not provided with coefficients we extract them from the model
  
  if (!is.list(coeffs)){
    nlbasis = vector()
    if (!mdl$clinear){
      for (i in 1:mdl$nSegments){
        nlbasis[i] = mdl$INLA$result$summary.fixed[paste("basis_",i,sep=""),field]
        if ("linbasis" %in% rownames(mdl$INLA$result$summary.fixed)) {
          lbasis = mdl$INLA$result$summary.fixed["linbasis",field]
        } else{
          lbasis = 0
        }
        coeffs = list(nlbasis=nlbasis,lbasis=lbasis,intercept=mdl$INLA$result$summary.fixed["(Intercept)",field])
      } 
    } else {
      for (i in 1:mdl$nSegments){
        nlbasis[i] = mdl$INLA$result$summary.hyperpar[paste("Beta for basis_",i,sep=""),field]
        if ("Beta for linbasis" %in% rownames(mdl$INLA$result$summary.hyperpar)) {
          lbasis = mdl$INLA$result$summary.hyperpar["Beta for linbasis",field]
        } else {
          lbasis = 0
        }
        coeffs = list(nlbasis=nlbasis,lbasis=lbasis,intercept=mdl$INLA$result$summary.fixed["(Intercept)",field])
      } 
    }
    
  }
  lincomb = rep(0,length(d))
  for (i in 1:nSeg){
    lincomb = lincomb + as.numeric(coeffs$nlbasis[i])*bweights[[i]]
  }
  lincomb = lincomb + as.numeric(coeffs$lbasis)*d
  return (lincomb)
}

#
# Sample values of the log detection function 
#
dfun_logconcave.sample.predictor = function(mdl,n=1,data=FALSE){
  source("inla_essmod.R")
  smp = inla.posterior.sample.structured(mdl$INLA$result,n)
  return(smp)
}

#
# Sample values of the detection function 
#
dfun_logconcave.sample.value = function(mdl,data=FALSE,n=1){
  if (n==1) {return(exp(dfun_logconcave.sample.logvalue(mdl,data=data,n=n)))}
  else {return(lapply(dfun_logconcave.sample.logvalue(mdl,data=data,n=n),exp))}

}


#
# Sample values of the log detection function 
#
dfun_logconcave.sample.logvalue = function(mdl,data=FALSE,n=1){

  if (is.data.frame(data)) { distance = data$distance }
  else if (is.numeric(data)) {distance = data} 
  else { distance = mdl$data$distance }
  structured.samples = inla.posterior.sample.structured(mdl$INLA$result,n)
  smp = list() # List of samples we will return
  for (i in 1:n) {
    ssmp = structured.samples[[i]]
    basis.names = paste0("basis_",1:mdl$nSegments)
    if (mdl$linbasis){ linbasis.coeff=as.numeric(ssmp["linbasis"]) } else {linbasis.coeff=0}
    coeffs = list(lbasis=linbasis.coeff,nlbasis=ssmp[basis.names],intercept=ssmp["(Intercept)"])
    smp[[i]] = dfun_logconcave.basis.lincomb(distance,mdl$nSegments,mdl$truncation,coeffs=coeffs)
  }
  if (n==1) {smp = unlist(smp)}
  return(smp)
   
}


#
# Sample a population
#

dfun_logconcave.sample.population = function(mdl=NULL,n=1,truncation=4,n.seg=4) {
  if (!is.null(mdl)){
    x = seq(0,mdl$truncation,0.1)
    smp = dfun_logconcave.sample.value(mdl,x,n)
    if (n==1) {smp = list(smp)}
    obs.samples = list()
    for (i in 1:length(smp)) {
      obs.samples[[i]] = sample(x,mdl$nData.after.filtering,prob=smp[[i]],replace=TRUE)
    }
    return(obs.samples)
  }
  else {
    coeffs = list(nlbasis=c(0,0,1,1),lbasis=0.0,intercept=1)
    xbasis = seq(0,truncation,0.01)
    basis = dfun_logconcave.basis.lincomb(xbasis,n.seg,truncation,coeffs=coeffs)
    samples = sample(xbasis,ns,prob=exp(-basis),replace=TRUE)
    return(samples)
  }
}
