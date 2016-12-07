#' INLA stack constructors for distance sampling
#' 
#' Distance sampling using INLA requires different kinds of \code{inla.stack} structures to be set up. 
#' These often represent different kinds of likelihoods but may also represent parts of the model 
#' that we might want to compute predictions for. iDistance provides the following stack constructors:
#' 
#' \itemize{
#'  \item{\link{detection.stack}: }{Stack for distance sampling detections}
#'  \item{\link{integration.stack}: }{Stack for integration points}
#'  \item{\link{prediction.stack}: }{Stack for predictions}
#'  \item{\link{grpsize.stack}: }{Stack for group size models}
#' }
#' @name stack
NULL


#' Create a stack for animal detections
#'
#' @aliases detection.stack
#' 

detection.stack = function(data,
                           model = NULL,
                           filter = NULL,
                           y = 1,
                           E = 0,
                           tag = "detection.stack"){
  
  # Extract points from data set
  if (class(data) == "SpatialPointsDataFrame" || class(data) == "SpatialPoints" || is.data.frame(data)) { pts = data } 
  else { pts = detdata(data) }
  
  
  # Apply filters provided by user
  if (!is.null(filter)) { pts = filter(pts) }
  
  # Projection matrices (A) and mesh index effects
  A = list.A.model(model, pts)
  eff = list.covariates.model(model, pts)
  idx = list.indices.model(model, pts)
  
  # Observations y
  y.pp = eval.if.function(y, pts)
  
  # Expectation parameter E
  if ( !(length(E) == nrow(as.data.frame(pts))) ) { e.pp = rep(E, nrow(as.data.frame(pts))) } else { e.pp = E } # as.data.frame required for SpatialPoints because they have no nrow function
  
  # Reweight things if dealing with taylor approximated model
  if ( !is.null(model$expr) ) {
    rw = nlinla.reweight(A, eff, model, pts)
    A = rw$A
    # eff = rw$weights
    e.pp = e.pp * exp(rw$const)
    if ( any((e.pp > 0) & (e.pp < 0.00001)) ) { warning("Exposure E is smaller than 0.00001. This may lead to numerical problems and non-convergence. Consider setting scale = 10 or higher when calling lgcp().")}
  }
  
  # Create and return stack
  stk = inla.stack(data = list(y.inla = y.pp, e = e.pp),
                     A = A,
                     tag = tag,
                     effects = c(idx, eff))
}


#' Create a stack for integration points
#'
#' @aliases integration.stack
#' 

integration.stack = function(data = NULL,
                                scheme = NULL,
                                scheme.args = NULL,
                                model = NULL,
                                filter = NULL,
                                y = 0,
                                E = "weight",
                                const = model$const,
                                tag = "integration.stack"){

  # Extract points from data set
  if (is.null(scheme)) { 
    scheme = select.integration(data) 
    pts = do.call(scheme, c(list(data=data), scheme.args))
  } else if (is.function(scheme)){
    pts = do.call(scheme, c(list(data=data), scheme.args))
  } else if (is.data.frame(scheme)) {
    pts = scheme
  } else if (class(scheme) == "SpatialPointsDataFrame") {
    pts = scheme
  }
  
  # Apply filters provided by user
  if (!is.null(filter)) { pts = filter(pts) }
    
  # Projection matrices (A) and mesh index effects
  A = list.A.model(model, pts)
  eff = list.covariates.model(model, pts)
  idx = list.indices.model(model, pts)
  
  # Observations y
  y.pp = eval.if.function(y, pts)
  
  # Expectation parameter E
  e.pp = as.data.frame(pts[,E, drop = FALSE])[,E] # Weird workaround for SpatialPointsDF
  
  # If const is empty list, set it to NULL
  if ( is.list(const) & length(const)==0 ) { const = NULL }
  # If const is a list, turn it into a function
  if ( is.list(const) ) { 
    clist = const
    const = function(loc) {apply(do.call(cbind, lapply(clist, function(cf) { cf(loc)})), MARGIN = 1, sum)}  }
  
  # Adding constants to the predictor turns into re-weighting the integration points
  if ( !is.null(const) ){ 
    const.values = const(pts)
    if ( !(length(const.values) == length(e.pp)) ) { 
      stop("E-vector length does not equal length of vector returned by const function") }
    e.pp = as.numeric(e.pp * exp( const.values )) 
  }
  
  # Reweight things if dealing with taylor approximated model
  if ( !is.null(model$expr) ) {
    rw = nlinla.reweight(A, eff, model, pts)
    A = rw$A
    # eff = rw$weights
    e.pp = e.pp * exp(rw$const) # min(exp(rw$const)
    if ( any((e.pp > 0) & (e.pp < 0.00001)) ) { warning("Exposure E is smaller than 0.00001. This may lead to numerical problems and non-convergence. Consider setting scale = 100 or higher (order of magnitude) when calling lgcp().")}
  }
  
  # Create and return stack
  stk = inla.stack(data = list(y.inla = y.pp, e = e.pp),
                   A = A,
                   tag = tag,
                   effects = c(idx, eff))
}


#' Create a stack for predictions
#'
#' @aliases prediction.stack
#' 

prediction.stack = function(data, model = NULL, loc = NULL, pts = NULL, tag = "prediction.stack"){
  
  if ( !is.null(pts) ) {
    warning("pts is a deprecated parameter. Use loc instead. Setting loc=pts.")
    loc = pts
  }
  pts = loc
  
  # Where to predict
  if (is.null(pts)) { 
    pts = data.frame(data$mesh$loc)[,1:length(data$mesh.coords)]
    colnames(pts) = data$mesh.coords
  }
  
  # Projection matrices (A) and mesh index effects
  A = list.A.model(model, pts)
  eff = list.covariates.model(model, pts)
  idx = list.indices.model(model, pts)
  
  # Number of prediction rows
  if (is.vector(A)) {siz = length(A)} else { siz = dim(A)[1] }
  np = max(ifelse( length(A) > 0, ifelse(is.matrix(A[[1]]),nrow(A[[1]]), siz), 0),
        ifelse( length(eff) > 0, nrow(eff[[1]]), 0))
  
  # Create and return stack
  return(inla.stack(data = list(y.inla = rep(NA,np), e = rep(NA,np)),
                    A = A,
                    tag = tag,
                    effects = c(idx, eff)))
}


#' Create a stack for group size observations
#'
#' @aliases grpsize.stack
#' 

grpsize.stack = function(data,
                           model = NULL,
                           filter = NULL,
                           y = function(x) { log(1+x$grpsize) },
                           tag = "grpsize.stack"){
  
  # Extract points from data set  
  pts = detdata(data)
  
  # Apply filters provided by user
  if (!is.null(filter)) { pts = det.filter(pts) }
  
  # Projection matrices (A) and mesh index effects
  A = list.A.model(model, pts)
  eff = list.covariates.model(model, pts)
  idx = list.indices.model(model, pts)
  
  # Observations y
  y.pp = eval.if.function(y, pts)
  
  # Expectation parameter E
  e.pp = rep(NA, nrow(pts))
  
  # Create and return stack
  return(inla.stack(data = list(y.inla = y.pp, e = e.pp),
                    A = A,
                    tag = tag,
                    effects = c(idx, eff)))
}

#' Create a stack from a model and data
#'
#' @aliases generic.stack
#' @name generic.stack

generic.stack = function(model = NULL,
                    data = NULL,   
                    filter = NULL,
                    y = 0,
                    E = NULL, ...) {

  # Extract points from data set
  if (is.dsdata(data)) { pts = data$effort } else { pts = data }
  
  # Apply filters provided by user
  if (!is.null(filter)) { pts = det.filter(pts) }
  
  # Projection matrices (A) and mesh index effects
  A = list.A.model(model, pts)
  eff = list.covariates.model(model, pts)
  idx = list.indices.model(model, pts)
  
  # Observations y
  y.inla = eval.if.function(y, pts)
  inla.data = list(y.inla = y.inla)
  
  # Expectation parameter E
  if ( !is.null(E) ) { 
    E = eval.if.function(E, pts) 
    inla.data$e = E
  }
  
  # Create and return stack
  stack = inla.stack(data = inla.data,
                     A = A,
                     effects = c(idx, eff), ...)
  return(stack)
}




#####################################
# Helpers
#####################################

eval.if.function = function(fun, x, n = NULL){
  if (is.null(n)) { n = nrow(as.data.frame(x)) }
  if ( is.function(fun) ) { fx = fun(x) }
  else if ( is.numeric(fun) ) {
    if ( length(fun) == 1 && n >1 ) { fx = rep(fun, n) }
    else { fx = fun }
  } else { 
    stop("unsupported data.")
  }
  return(fx)
}
