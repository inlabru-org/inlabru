#' INLA inference in log Gaussian Cox processes for distance sampling data.
#'
#'
#' @aliases pproc_lgcp
#' @export
#' @param data A distance sampling or SECR data set
#' @param int.points Integration points for the approximation of the LGCP denominator
#' @param int.scheme Function that maps a data set to integration points. Needed if no integration points are provided. 
#' @param int.args List of arguments passed on to int.scheme()
#' @param dist.trafo Function that transforms distance covariate if integration and observation points. Default: -0.5*d^2
#' @param covariates List of functions that add covariates to data frames of integration/observation points
#' @param formula The formula that is used in the call to INLA
#' @param spde.args List of parameters passed on to inla.spde2.matern() when constructing the SPDE
#' @param predict List of linear combinations to predict. Each list entry is a list with two fields. "vars" holds the list of variables to include in the prediction and "loc" describes the locations to predict at.
#' @param constant Function providiong constants to add to the linear predictor
#' @param inla.args List of arguments passed on to INLA. Default: list(family = "poisson"). 
#' @param sgh.E Poisson likelihood exposure parameter for detections (see INLA documentation)
#' @param int.filter Function applied to the data frame of integration points. Use this to filter out particular integration points.
#' @param det.filter Function applied to the data frame of detection points. Use this to filter out particular detection points.
#' @param sgh.y Observation values for sightings. Default: 1
#' @params int.y Observation values for intefration points. Default: 0, later in the alogrithm replaced by integration weights
#' @return \code{pproc_lgcp}, lgcp object.
#' @examples \\dontrun{data(whales) ; pp = pproc_lgcp(whales); plot(pp)}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

pproc_lgcp = function(data,int.points=NULL,
                      int.scheme = NULL,
                      int.args = list(truncation=max(detdata(data)$distance)),
                      dist.trafo = function(dst) { return(-0.5*dst^2) },
                      covariates = list(),
                      formula = y ~ -1 + Intercept + f(distance, model="clinear", range = c(0, Inf), hyper=list(beta=list(initial=0, param=c(0,0.1)))) + f(spde, model=spde.mdl),
                      spde.args = list(alpha=2,prior.variance.nominal=10,theta.prior.prec=0.01),
                      predict = list(),
                      constant = list(),
                      inla.args = list(family = "poisson"),
                      sgh.E = 0,
                      int.filter = NULL,
                      det.filter = NULL,
                      sgh.y = 1,
                      int.y = 0,
                      stack = list()){
  
  #
  # Check if formula actually is a formula update
  #
  
  if (!("y" %in% all.vars(formula))) {
    formula.default =  y ~ -1 + Intercept + f(distance, model="clinear", range = c(0, Inf), hyper=list(beta=list(initial=0, param=c(0,0.1)))) + f(spde, model=spde.mdl)
    formula = update(formula.default,formula)
  }
  
  #
  # Infer fields storing coordinates
  #
  # compatibility
  if (!("mesh.coords" %in% names(data))) { data$mesh.coords = c("lon","lat")}
  
  #
  # Integration and sighting points
  #
  if (is.character(int.scheme)) {
    if (int.scheme == "none") { int.points = NULL }
  } else {
    if (is.null(int.scheme) & is.null(int.points)) { int.scheme = select.integration(data) }
    if (is.null(int.points) ) { 
      if ("ips" %in% names(data)) { int.points = data$ips } # backwards compatibility
      else if ("int.points" %in% names(data)) { int.points = data$int.points }
      else {
        int.points = do.call(int.scheme,c(list(data=data),int.args))
      }    
    }
  }
  sgh.points = detdata(data)
  
  # Filter out points outside the mesh boundary
  sgh.inside = is.inside(data$mesh,sgh.points,mesh.coords=data$mesh.coords)
  if ( any(!sgh.inside) ) {
    warning("pproc_lgcp: Sighting ouside mesh boundary. Discarding.")
    sgh.points = sgh.points[sgh.inside,]
  }
  
  #
  # Apply filters provided by user
  #
  if (!is.null(det.filter)) { sgh.points = det.filter(sgh.points) }
  if (!is.null(int.filter)) { int.points = int.filter(int.points) }
  
  #
  # Count integration and detection points
  #
  
  if ( is.null(sgh.points) ) { sgh.n = 0 } else { sgh.n = nrow(sgh.points) }
  if ( is.null(int.points) ) { int.n = 0 } else { int.n = nrow(int.points) }
  
  #
  # Make SPDE model
  #
  
  spde.mdl = do.call(inla.spde2.matern,c(list(mesh=data$mesh),spde.args)) # n.iid.group=1
  
  #
  # Projection matrices (A)
  #
  
  if ( sgh.n > 0 ) { locmat = inla.spde.make.A(data$mesh, loc=as.matrix(sgh.points[,data$mesh.coords]),repl=1) } else { locmat = NULL}
  if ( int.n > 0 ) { 
    intmat = inla.spde.make.A(data$mesh, loc=as.matrix(int.points[,data$mesh.coords]),repl=1)
    
    # Sanity check. Are all of out points within the mesh?
    if ( any(apply(!(abs(as.matrix(intmat))==0),MARGIN=1,sum)==0) ) {
      warning("Integration points outside of mesh boundary! Discarding. This affects your results!")
      int.points = int.points[!(apply(!(abs(as.matrix(intmat))==0),MARGIN=1,sum)==0),]
      intmat = inla.spde.make.A(data$mesh, loc=as.matrix(int.points[,c("lon","lat")]))
    }
  } else { 
    intmat = NULL
  }
  
  if (int.n ==0 ) { A.pp = locmat }
  else if ( sgh.n == 0 ) { A.pp = intmat }
  else if( int.n > 0 & sgh.n > 0 ) { A.pp = rBind(intmat,locmat) }
  else { stop("The A matrix of integration and sighting points are both NULL.") }
  
  #
  # Observations and expectation (Y and E)
  #
  
  # y of integration points
  if ( is.function(int.y) ) { int.y.values = int.y(int.points) }
  else if ( is.numeric(int.y) ) {
    if ( length(int.y) == 1 && int.n >1 ) { int.y.values = rep(int.y,int.n) }
    else { int.y.values = int.y }
    if ( int.n == 0 ) { int.y.values = NULL }
  } else { 
    stop("unsupported data type of parameter int.y") 
  }
  
  # y of sightings
  if ( is.function(sgh.y) ) { sgh.y.values = sgh.y(sgh.points) }
  else if ( is.numeric(sgh.y) ) {
    if ( length(sgh.y) == 1 && sgh.n >1 ) { sgh.y.values = rep(sgh.y,sgh.n) }
    else { sgh.y.values = sgh.y }
    if ( sgh.n == 0 ) { sgh.y.values = NULL }
  } else { 
    stop("unsupported data type of parameter sgh.y")
  }
  
  # concatenated y
  y.pp = c(int.y.values, sgh.y.values)
  
  # E parameter
  if ( int.n > 0) {
    e.pp = c(int.points[,"weight"],rep(sgh.E,nrow(sgh.points)))
  } else {
    e.pp = rep(sgh.E,nrow(sgh.points))
  }
  
  
  
  
  #
  # Covariates
  #
  if (int.n > 0 ) { int.covar = fetch.covariate(formula,int.points,covariates) } else { int.covar = NULL }
  if (sgh.n > 0 ) { sgh.covar = fetch.covariate(formula,sgh.points,covariates) } else { sgh.covar = NULL }
  covar.data = rbind(int.covar,sgh.covar)
  
  
  #
  # Convert distances
  #
  if ("distance" %in% names(covar.data)) {
    covar.data[,"distance"] = dist.trafo(covar.data[,"distance"])
  }
  
  
  #
  # Constants
  #
  
  if (length(constant)>0){
    fetched.int.const = list()
    fetched.sgh.const = list()
    for (cname in names(constant)){
      if (is.function(constant[[cname]])){
        cfun = constant[[cname]]
        fetched.int.const[[cname]] = cfun(int.points)
        fetched.sgh.const[[cname]] = cfun(sgh.points)
      } 
    } 
    int.const = do.call(cbind,fetched.int.const)
    sgh.const = do.call(cbind,fetched.sgh.const)
    e.pp[1:nrow(int.points)] = e.pp[1:nrow(int.points)] * exp(apply(int.const,1,sum))
  }
  
  
  
  #
  # Stack everything together
  #
  
  stk = inla.stack(data=list(y=y.pp, e = e.pp),
                   A=list(A.pp,1), tag='fPP1',
                   effects=list(list(spde=1:spde.mdl$n.spde), covar.data)
  )
  
  #
  # Add combined predictions on mesh to stack
  #
  
  for (tag in names(predict)){
    tmp.pred = predict[[tag]]
    if (!("loc" %in% names(tmp.pred))) { 
      pred.loc = data$mesh$loc
      colnames(pred.loc) = c(data$mesh.coords,"tmp")
      spde.pred.loc = NULL # as.matrix(pred.loc[,data$mesh.coords])
    }
    else { 
      pred.loc = tmp.pred$loc
      spde.pred.loc = NULL
    }
    
    if ("spde" %in% tmp.pred$vars){
      tmp.stk = inla.stack(data=list(y=rep(NA,nrow(pred.loc))),
                           A=list(inla.spde.make.A(data$mesh, loc=as.matrix(pred.loc[,data$mesh.coords])),1), tag=tag,
                           effects=list(list(spde=1:spde.mdl$n.spde), 
                                        fetch.covariate(tmp.pred$vars,data.frame(pred.loc),covariates)))
    } else {
      tmp.stk = inla.stack(data=list(y=rep(NA,dim(pred.loc)[1])),
                           A=list(1), tag=tag,
                           effects=list(fetch.covariate(tmp.pred$vars,data.frame(pred.loc),covariates)))
    }
    stk = inla.stack(stk,tmp.stk)
  }
  
  #
  # Join internal stack with stack provided by user
  #
  
  if ( !(length(stack)==0) ) {
    y.joint = inla.stack.y(stk,stack)
    e.joint = inla.stack.y(stk,stack)
    stk = inla.stack.join(stk,stack)
    inla.data = inla.stack.data(stk)
    inla.data$y = y.joint
    inla.data$e = e.joint
  } else {
    inla.data = inla.stack.data(stk)
  }
  
  
  #
  # Modify environment of the formula
  #
  
  environment(formula) = environment()
  
  #
  # Some more defaults
  #
  if (!("family" %in% names(inla.args))) { inla.args$family = "poisson" }
  
  #
  # Run INLA
  #
  
  resultFPP = do.call(inla, c(inla.args,
                              list(
                                formula, 
                                data = inla.data,
                                control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
                                E = inla.stack.data(stk)$e,
                                verbose = TRUE)
  ))
  
  #
  # Return object
  #
  
  ret = list(INLA=list(result=resultFPP,formula=formula),
             stack=stk,
             data=data,
             int.args=int.args,
             int.points=int.points,
             sgh.points=sgh.points,
             spde=spde.mdl,
             covariates = covariates,
             dist.trafo = dist.trafo,
             sgh.covar = sgh.covar,
             int.covar = int.covar)
  class(ret) = c("pproc_lgcp","pproc")
  return (ret)
  
}

fetch.covariate = function(formula,points,covariates=NULL){
  if (is.character(formula)){ all.varnames = formula }
  else { all.varnames = all.vars(formula[[3]]) }
  covar.data = list(Intercept=rep(1,nrow(points)))
  for (vname in all.varnames){
    if (vname %in% names(points)) {
      covar.data[[vname]] = points[,vname]
    }
  }
  covar.data = data.frame(do.call(cbind,covar.data))
  
  # Add additional covariates defined by the user
  fetched.covar = list()
  for (cov.name in names(covariates)){
    if (cov.name %in% all.varnames){
      cov.fun = covariates[[cov.name]]
      fetched.covar[[cov.name]] = cov.fun(points)
    }
  }
  fetched.covar = do.call(cbind,fetched.covar)
  if (length(fetched.covar)>0){ covar.data = cbind(covar.data,fetched.covar)}
  
  return(covar.data)
}


#' Simulate a lgcp
#'
#'
#' @aliases simulate.lgcp
#' @export
#' @param pp LGCP point process object
#' @return points
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' @author Daniel Simposon <\email{ZZZ@@bath.ac.uk}>
#'

simulate.pproc = function(mdl,formula=NULL,covariates=mdl$covariates,as.data=FALSE,property=NULL,time=1986) {
  if (is.null(property)) {
    loc = data.frame(mdl$data$mesh$loc[,1:length(mdl$data$mesh.coords)])
    colnames(loc) = mdl$data$mesh.coords
    loc[[mdl$data$time.coords]] = time
    weights = sample.logintensity.pproc(mdl,points=loc,covariates=covariates,formula=formula,n=1)
  } else {
    weights = logintensity.pproc(mdl,points=loc,covariates=covariates,formula=formula,property=property)
  }
  
  pts = simulate_lgcp(mdl$data$mesh,weights)
  colnames(pts) = mdl$data$mesh.coords
  if (as.data){
    data = mdl$data
    data$sighting = pts
    class(data) = c("etpdata","dsdata")
    return(data)
  } else {
    return(pts)
  }
}



#' Sample log intensity of LGCP at given locations
#' 
#'
#'
#' @aliases sample.logintensity.pproc
#' @export
#' @param pp LGCP point process object
#' @param points locations
#' @return log intensity sample
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

sample.logintensity.pproc = function(mdl,points=mdl$data$mesh$loc,covariates=mdl$covariates,formula=mdl$INLA$formula,n=1){

  if (is.character(formula)) { fml.vnames = formula }
  else { fml.vnames = all.vars(terms.formula(formula)[[3]]) }
  
  rate = rep(0,nrow(points))
  vars = variables.inla(mdl$INLA$result)
  covs = fetch.covariate(formula,points,covariates=covariates)
  smps = inla.posterior.sample.structured(mdl$INLA$result,n)
  
  for (s in 1:length(smps)){
    smp = smps[[s]]
    for (v in 1:nrow(vars)){
      varname = rownames(vars)[[v]]
      if (varname %in% fml.vnames){ 
        if (vars[v,"type"] == "random"){
          if (vars[v,"model"]=="SPDE2 model") {
            weights = smp[[varname]]
            if (is.null(colnames(points))){fmlocs = points}
            else { fmlocs = as.matrix(points[,mdl$data$mesh.coords]) }
            A = inla.mesh.project(mdl$data$mesh,loc=fmlocs)$A
            vals = A%*%weights
            rate[,s] = rate[,s] + as.vector(vals)
          }
        } else {
          if (is.null(covs[[varname]])) { rate = rate + smp[[varname]]}
          else { rate = rate + smp[[varname]]*covs[[varname]] }       
        }
      }
    }
  }
  if (length(rate)==0) {
    warning("Log intensity is numeric(0), most likely you forgot to provide sufficient covariates with the argument points.")
  }
  return(rate)
}


#' Log intensity of LGCP at given locations
#' 
#' Currently based on the marginal mean of each latent variable
#'
#'
#' @aliases logintensity.pproc
#' @export
#' @param pp LGCP point process object
#' @param loc locations
#' @return log intensity
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

logintensity.pproc = function(mdl,points=mdl$data$mesh$loc,covariates=mdl$covariates,formula=mdl$INLA$formula,property="mean"){
  if (is.character(formula)) { fml.vnames = formula }
  else { fml.vnames = all.vars(terms.formula(formula)[[3]]) }
  
  rate = rep(0,nrow(points))
  vars = variables.inla(mdl$INLA$result)
  covs = fetch.covariate(formula,points,covariates=covariates)
  
  for (v in 1:nrow(vars)){
    varname = rownames(vars)[[v]]
    if (varname %in% fml.vnames){   
      if ( vars[v,"type"] == "random" ) {
        if (vars[v,"model"]=="SPDE2 model") {
          weights = mdl$INLA$result$summary.ran[[varname]][[property]]
          fmlocs = as.matrix(points[,mdl$data$mesh.coords])
          A = inla.mesh.project(mdl$data$mesh,loc=fmlocs)$A
          vals = A%*%weights
          rate = rate + as.vector(vals)
        }
      } else {
        if (is.null(covs[[varname]])) { rate = rate + vars[v,property]}
        else { rate = rate + vars[v,property]*covs[[varname]] }       
      }
    }
  }
  if (length(rate)==0) {
    warning("Log intensity is numeric(0), most likely you forgot to provide sufficient covariates with the argument points.")
  }
  return(rate)
}



#' Two-stage LGCP
#' 
#'
#'
#' @aliases logintensity.pproc
#' @export
#' @param pp LGCP point process object
#' @param loc locations
#' @return log intensity
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

twostage.lgcp = function(s1.args=NULL,s2.args=NULL,s1.logintensity=NULL,s2.logintensity=NULL,s1conditional.args,n1=1,n2=1,n.iter=1){
  
  #
  # Inital estimate(s) of stage 1
  #
  
  stage1.pp = list()
  stage1.pp[[1]] = do.call(lgcp,args=s1.args)
  
  # DEBUG
  # plot.dfun_halfnormal(stage1.pp[[1]])
  # summary(stage1.pp[[1]])
  
  
  #
  # Stage 2 inference
  #
  stage2.pp = list()
  s1.history=list()
  s2.history=list()
  
  for (k in 1:n.iter){
    
    for (s1 in 1:n1){
      s1logint = function(x) {return(s1.logintensity(stage1.pp[[1]],x))}  
      stage2.pp[[s1]] = do.call(lgcp,args=c(s2.args, list(constant = list(s1logint = s1logint))))
      
      # DEBUG
      # plot(stage2.pp[[1]])
      # summary(stage2.pp[[1]])
      
    }
    # If we iterate, re-estimate stage 1
    if (n.iter>1){
      s2logint = function(x) {return(s2.logintensity(stage2.pp[[1]],x))}  
      stage1.pp[[1]] = do.call(lgcp,args=c(s1conditional.args, list(constant = list(s2logint = s2logint))))
    }
    s2.history[[k]] = variables.inla(stage2.pp[[1]]$INLA$result)
    s1.history[[k]] = variables.inla(stage1.pp[[1]]$INLA$result)
  }
  
  return(list(stage1 = stage1.pp, 
              stage2 = stage2.pp,
              s1.history = s1.history,
              s2.history = s2.history))
}




###########################################################################
## simulate_lgcp.R
## Code for simulating LGCPs on triangulated domains (planar or spherical)
## used throughout Simpson et al. (2011)
## Author: Daniel Simpson
## This code comes with no warranty or guarantee of any kind.
###########################################################################


simulate_lgcp <- function(mesh,weights) {
  
  if (is.null(geometry)) { geometry = mesh$manifold}
  radius.of.earth = 6371
  
  if(geometry == "R2") {
    # print("It's flat!")
    #Construct bounding rectangle
    loc <- mesh$loc
    xmin = min(loc[,1])
    xmax = max(loc[,1])
    ymin = min(loc[,2])
    ymax = max(loc[,2])
    area =(xmax- xmin)*(ymax- ymin)
    
    #Simulate number of points
    lambda_max <- max(weights)
    Npoints <- rpois(n=1,lambda=area*exp(lambda_max))
    
    #if (Npoints > 1e5) {
    #  print(Npoints)
    #  print("You've got to be joking")
    #  return(NULL)
    #} else if (Npoints==0) 
    #{ return(NULL) }
    
    #Simulate uniform points on the bounding rectangle
    x <- runif(n=Npoints, min=xmin, max=xmax)
    y <- runif(n=Npoints, min=ymin, max=ymax)
    
    points <-cbind(x,y)
    
    #Do some thinning
    A <- inla.mesh.project(mesh,points)$A
    
    weights =exp( weights-lambda_max)
    
    pointValues = A%*%weights
    
    keep =which(runif(Npoints) < pointValues)
    
    return(points[keep,])
  }
  
  if (geometry=="S2") {
    print("Ground control to Major Tom")
    area =   4*pi*radius.of.earth^2 #we're on the earth
    
    #Simulate number of points
    lambda_max <- max(weights)
    Npoints <- rpois(n=1,lambda=area*exp(lambda_max))
    print(area*exp(lambda_max))
    print(Npoints)
    
    if (Npoints > 1e5) {
      print(Npoints)
      print("You've got to be joking")
      return(NULL)
    }
    
    #Simulate random points on the sphere
    #this could be faster
    points= c()
    for (i in 1:Npoints)
    {
      tmp = rnorm(3)
      points = rbind(points,tmp/sqrt(t(tmp)%*%tmp))
    }
    
    #Do some thinning
    A <- inla.mesh.project(mesh,points)$A
    
    weights =exp( weights-lambda_max)
    
    pointValues = A%*%weights
    
    keep =which(runif(Npoints) < pointValues)
    
    return(points[keep,])
    
  }
  
  if (geometry == "geo-slice") {
    
    # determine area (slice of earth)
    lon.range = range(mesh$loc[,1])
    factor = abs(diff(lon.range))/360
    area =   factor*4*pi*radius.of.earth^2 #we're on the earth
    print(area)
    
    # Simulate number of points
    lambda_max <- max(weights)
    Npoints <- rpois(n=1,lambda=area*exp(lambda_max))
    
    lon <- runif(n=Npoints,min=lon.range[1],max=lon.range[2])
    lat <- runif(n=Npoints,min=-90,max=90)
    
    # Do some thinning
    A <- inla.mesh.project(mesh,loc=cbind(lon,lat))$A
    weights =exp( weights - lambda_max )
    pointValues = A%*%weights
    
    print(mean(pointValues>0))
    keep =which(runif(Npoints) < pointValues)
    points = data.frame(lon=lon[keep],lat=lat[keep])
    
    print(mean(pointValues>0)*area)
    
    return(points)
    
  }
  
}
