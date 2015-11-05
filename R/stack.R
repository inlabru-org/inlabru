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
  e.pp = rep(E, nrow(pts))
  
  # Create and return stack
  return(inla.stack(data = list(y.inla = y.pp, e = e.pp),
                     A = A,
                     tag = tag,
                     effects = c(idx, eff),
                    compress = FALSE))
}


#' Create a stack for integration points
#'
#' @aliases integration.stack
#' 

integration.stack = function(data,
                                scheme = NULL,
                                scheme.args = NULL,
                                model = NULL,
                                filter = NULL,
                                y = 0,
                                E = "weight",
                                const = NULL,
                                tag = "integration.stack"){

  # Extract points from data set
  if (is.null(scheme)) { 
    scheme = select.integration(data) 
    pts = do.call(scheme, c(list(data=data), scheme.args))
  } else if (is.function(scheme)){
    pts = do.call(scheme, c(list(data=data), scheme.args))
  } else if (is.data.frame(scheme)) {
    pts = scheme
  }
  
  # Apply filters provided by user
  if (!is.null(filter)) { pts = det.filter(pts) }
    
  # Projection matrices (A) and mesh index effects
  A = list.A.model(model, pts)
  eff = list.covariates.model(model, pts)
  idx = list.indices.model(model, pts)
  
  # Observations y
  y.pp = eval.if.function(y, pts)
  
  # Expectation parameter E
  e.pp = pts[,E]
  
  # Adding constants to the predictor turns into re-weighting the integration points
  if ( !is.null(const) ){ e.pp = e.pp * exp( const(pts) ) }
  names(e.pp) = "e.inla" # INLA will complain if the column name coincides with a vraible name on the RHS 
  
  # Create and return stack
  return(inla.stack(data = list(y.inla = y.pp, e = e.pp),
                    A = A,
                    tag = tag,
                    effects = c(idx, eff)))
}


#' Create a stack for predictions
#'
#' @aliases prediction.stack
#' 

prediction.stack = function(data, model = NULL, pts = NULL, tag = "prediction.stack"){
  
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
  np = max(ifelse( length(A) > 0, ifelse(is.matrix(A[[1]]),nrow(A[[1]]),size(A)), 0),
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


#####################################
# Helpers
#####################################

eval.if.function = function(fun, x, n = NULL){
  if (is.null(n)) { n = nrow(x) }
  if ( is.function(fun) ) { fx = fun(x) }
  else if ( is.numeric(fun) ) {
    if ( length(fun) == 1 && n >1 ) { fx = rep(fun, n) }
    else { fx = fun }
  } else { 
    stop("unsupported data.")
  }
  return(fx)
}
