#' SPDE detection function fitting using INLA.
#' 
#'
#' @aliases dfun_spde1d
#' @export 
#' @param data
#' @param data A data set containing a list of sightings.
#' @param truncation The maximum observed distance the model will take into account.
#' @param nIntegrate The number of integration points used to determine the denominator
#' of the Cox process.
#' @param nKnots Number of knots for the latent SPDE model.
#' @param boundary Left and right boundary condition for the latent SPDE model.
#' @return \code{dfun}, detection function object. 
#' @examples \\dontrun{data(whales) ; fun = dfun_spde1d(whales$sighting); plot(fun)}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

dfun_spde1d = function(data,truncation=max(detdata(data)$distance),nIntegrate=100,nKnots=6,boundary=c("neumann","neumann"))
{
  
  require(INLA)
  
  
  #
  # Abbreviations
  #
  d = data$distance # the distance data
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
  # Integration weights (fake E)
  # 
  

    # width of integration element
    hi = truncation/(nInt-1)
    # weights for trapeziodal rule
    wtr = hi*c(0.5,rep(1,nInt-2),0.5)
  
    fake_E = c(wtr,rep(0,nDat))

  
  #
  # Define mesh
  #
    
    knots = seq(0, truncation, length=nKnots)
    mesh = inla.mesh.1d(knots, interval=c(0, truncation), degree=2, boundary=boundary)
    
  
  #
  # Define A matrices
  #
  
    # integration points
    A.int = inla.spde.make.A(mesh, loc=ips)
    
    # data points
    A.data = inla.spde.make.A(mesh, loc=data$distance)
  
    # Joint A matrix
    A.joint <- rBind(A.int, A.data)

  #
  # SPDE model
  #
    
    # parameters
    sigma0 = 1
    kappa0 = 1e-3
    tau0 = 1/(4*kappa0^3*sigma0^2)^0.5
    
    # the model
    spde = inla.spde2.matern(mesh, constr=FALSE,
                             B.tau = cbind(log(tau0), 1),
                             B.kappa = cbind(log(kappa0), 0),
                             theta.prior.prec = 1e-4)
  
  #
  # Stack
  #
  
    distance.index = inla.spde.make.index("distance", n.spde=spde$n.spde)  
    stack <- inla.stack(data=list(y=fake_data, e=fake_E),
                       A=list(A.joint), tag='pp',
                       effects=list(distance.index))
  
  #
  # Define INLA formula
  #
  
  formula = y ~ -1 + f(distance, model=spde) # + f(rep(1,length(fake_data)),model="linear")
  
  
  #
  # Run INLA
  #

  result = inla(formula,
                family='poisson', data=inla.stack.data(stack),
                control.compute=list(config = TRUE),
                control.predictor=list(A=inla.stack.A(stack),compute=TRUE),
                E=inla.stack.data(stack)$e) # , verbose=TRUE,debug=TRUE,keep=TRUE
  
  
  #
  # Return object
  #
  ret = list(INLA=list(result=result,formula=formula,stack=stack,mesh=mesh),
             integration.points=ips,
             truncation=truncation,
             type="spde1d",
             data=data,
             nKnots=nKnots)
  
  class(ret) = c("dfun_spde1d","dfun")
  return (ret)
}


dfun_spde1d.logvalue = function(mdl,data=FALSE,field="mode"){
  
  if (is.data.frame(data)) { distance = data$distance }
  else if (is.numeric(data)) {distance = data} 
  else { distance = mdl$data$distasnce }

  latent = mdl$INLA$result$summary.random$distance[,field]
  A = inla.spde.make.A(mdl$INLA$mesh, loc=distance) 
  return(as.vector(A%*%latent))
}

dfun_spde1d.value = function(...){
  return (exp(dfun_spde1d.logvalue(...)))
}



#
# Sample values of the detection function 
#
dfun_spde1d.sample.value = function(mdl,data=FALSE,n=1){
  if (n==1) {return(exp(dfun_spde1d.sample.logvalue(mdl,data=data,n=n)))}
  else {return(lapply(dfun_spde1d.sample.logvalue(mdl,data=data,n=n),exp))}
  
}


#
# Sample values of the log detection function 
#
dfun_spde1d.sample.logvalue = function(mdl,data=FALSE,n=1){
  
  
  if (is.data.frame(data)) { distance = data$distance }
  else if (is.numeric(data)) {distance = data} 
  else { distance = mdl$data$distance }
  
  sismp=inla.posterior.sample.structured(mdl$INLA$result,n)
  smp = list()
  
  for (i in 1:n){
    latent = sismp[[i]]$distance
    x = seq(0,mdl$truncation,0.1)
    A = inla.spde.make.A(mdl$INLA$mesh, loc=distance) 
    smp[[i]] = as.vector(A%*%latent)
  }
  if (n==1) {smp = unlist(smp)}
  return(smp) 
}
