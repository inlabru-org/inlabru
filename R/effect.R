####################################################################################################
# GENERICS
####################################################################################################
# 
# setClass("effect")
# setGeneric("value", valueClass = "numeric", function(object) { standardGeneric("value") })
# setGeneric("map", valueClass = "numeric", function(object) { standardGeneric("map") })
# setGeneric("amatrix", valueClass = "numeric", function(object) { standardGeneric("amatrix") })
# 
# setMethod("value", signature("effect"), function(object) value.effect(object))
# setMethod("amatrix", signature("effect"), function(object) amatrix.effect(object))
# setMethod("map", signature("effect"), function(object) map.effect(object))
# 
# 

value = function(...){UseMethod("value")}
amatrix = function(...){UseMethod("amatrix")}
map = function(...){UseMethod("map")}
index = function(...){UseMethod("index")}


####################################################################################################
# CONSTRUCTORS
####################################################################################################

effect = function(formula, ...) {
  effs = make.model(formula)$effects
  effs = lapply(effs, function(effect) {
      effect$env = environment(formula)
      class(effect) = "effect"
      effect
    })
  if ( length(effs) == 1) {effs = effs[[1]]}
  
  effs
}


make.effect = function(label,
                       model,
                       map = NULL,
                       group = NULL,
                       values = NULL,
                       n = NULL,
                       A.msk = NULL,
                       ...){
  
  
  if ( inherits(model, "inla.spde") ) { 
    spde.model = model
    model = "spde"
  }
  
  # map.char = as.character(substitute(map))
  # group.char = as.character(substitute(group))
  # A.msk.char = as.character(substitute(A.msk))
  # 
  # Determine effect depending on the type of latent model
  if ( model %in% c("offset") ) {
    
  }
  
  else if ( model %in% c("linear", "clinear") ) {
    xxx = NULL; 
    effect = INLA::f(xxx, ..., group = group, model = model)
  }
  
  else if ( model %in% c("factor") ) {
    
  }
  
  else if ( model %in% c("iid") ) {
    
  }
  
  else if ( model %in% c("seasonal") ) {

  }
  else if ( model %in% c("rw1", "rw2", "ar", "ar1", "ou") ) {

  }
  else if ( model %in% c("spde") ) {

  } else {
    stop(paste0("Effect type '", model, "' not implemented."))
  }
  
  
}

####################################################################################################
# OPERATORS
####################################################################################################


#' Sumarize an effect
#'  
#' @aliases summary.effect
#' @keywords internal
#' @param effect An effect.
#' @param ... ignored.
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

summary.effect = function(object, ...) {
  
  eff = object
  cat(paste0("Label:\t ", eff$label, "\n"))
  #en = ifelse("n" %in% names(eff), eff$n, 1)
  cat(sprintf("model:\t %s\n", eff$model))
  cat(sprintf("map:\t %s [class: %s] \n", deparse(eff$map), class(eff$map)))
  cat(sprintf("f call:\t %s \n", eff$fchar))
  cat("\n")
  
}

#' Evaluate an effect
#' 
#' Calculates a latent effect given some data and the state of the effect's internal random variables.
#' 
#' TODO: Improve speed for iterated calls by making 'mapped' a parameter 
#' 
#' @aliases value.effect
#' @keywords internal
#' @param effect An effect.
#' @param data A \code{data.frame} or Spatial* object of covariates and/or point locations.
#' @param state Either a numeric vector or a list with a numeric entry whose name is equal to the name parameter.
#' @param A A matrix overriding the default projection matrix.
#' @return A numeric vector of effect values
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

value.effect = function(effect, data, state, A = NULL) {
  
  # Extract state if a list of states was provided
  if ( is.list(state) &!is.data.frame(state)) { state = state[[effect$label]] }
  
  # Obtain covariates
  mapped = mapper(effect$map, data, effect, effect$env)
  
  # Obtain A-matrix
  if ( is.null(A) ) { A = amatrix(effect, data) }
  
  # Determine effect depending on the type of latent model
  if ( effect$model %in% c("linear", "clinear") ) {
    if (is.data.frame(mapped)) { mapped = mapped[,1,drop=TRUE] }
    values = A %*% (state * mapped)
  }
  else if ( effect$model %in% c("factor") ) {
    if (is.data.frame(mapped)) { flevel = mapped[,effect$label] } else {flevel = mapped}
    values = A %*% state[flevel]
  } 
  else if ( effect$model %in% c("iid", "seasonal") ) {
    values = A %*% state[mapped]
  }
  else if ( effect$model %in% c("rw1", "rw2", "ar", "ar1", "ou") ) {
    if ( ncol(data.frame(state)) == 2 ) {
      values = A %*% state[, 2, drop = TRUE]
    } else {
      values = A %*% state
    }
    
  }
  else if ( effect$model %in% c("spde2") ) {
    if (is.data.frame(state)) { state = state$value }
    values = A %*% state
  } else {
    stop(paste0("Evaluation of ", effect$model, " not implemented."))
  }
  
  as.vector(values)
}



#' Construct A-matrix
#' 
#' Constructs the A-matrix for a given effect and and some data
#'  
#' @aliases amatrix.effect
#' @export
#' @keywords internal
#' @param effect An effect.
#' @param data A \code{data.frame} or Spatial* object of covariates and/or point locations.
#' @return An A-matrix.
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

amatrix.effect = function(effect, data) {
  
  if ( effect$model %in% c("spde2") ) {
    if ( effect$map == "coordinates" ) {
      if ( is.na(proj4string(data)) | is.null(effect$mesh$crs) ) {
        loc = data
      } else {
        loc = stransform(data, crs = effect$mesh$crs)
      }
      
    } else {
      loc = mapper(effect$map, data, effect)
      if ( !is.matrix(loc) & !inherits(loc,"Spatial") ) loc = as.matrix(loc)
    }
    
    if (is.null(effect$ngroup)) { ng = 1 } else { ng = effect$ngroup }
    if (ng > 1) {
      group = data[[effect$group.char]]
    } else { group = NULL }
    
    A = INLA::inla.spde.make.A(effect$mesh, loc = loc, group = group, n.group = ng)
  } 
  else if ( effect$model %in% c("rw1", "rw2", "ar", "ar1", "ou") ) {
    if ( is.null(effect$values) ) { stop("Parameter 'values' not set for effect '", effect$label, "'.") }
    mesh = INLA::inla.mesh.1d(effect$values)
    loc = mapper(effect$map, data, effect)
    A = INLA::inla.spde.make.A(mesh, loc)
  }
  else {
    if ( is.null(data) ) { A = 1 }
    else { A = Matrix::Diagonal(nrow(data.frame(data))) }
  }

  # Mask columns of A
  if (!is.null(effect$A.msk)) { 
    A = A[, as.logical(effect$A.msk), drop=FALSE]
  }
  
  # Weight rows of A
  # if (!is.null(effect$weights)) {
  #   A = as.matrix(A)*as.vector(effect$weights)
  # }
  
  A
}

#' Obtain covariate
#' 
#' Uses the effet's map value/function to extract covariates.
#'  
#' @aliases map.effect
#' @keywords internal
#' @param effect An effect.
#' @param data A \code{data.frame} or Spatial* object of covariates and/or point locations. If null, return the effect's map.
#' @return An A-matrix.
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

map.effect = function(effect, data) {
  
  mp = mapper(effect$map, data, effect, effect$env)
  
}



#' Obtain indices
#' 
#' Idexes into to the effects
#'  
#' @aliases index.effect
#' @keywords internal
#' @param effect An effect.
#' @param data A \code{data.frame} or Spatial* object of covariates and/or point locations. If null, return the effect's map.
#' @return indices
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

index.effect = function(effect, data) {
  
  if ( effect$model %in% c("spde2") ) {
    if ( "m" %in% names(effect$mesh) ) {
      idx = 1:effect$mesh$m # support inla.mesh.1d models
      # If a is masked, correct number of indices
      if ( !is.null(effect$A.msk) ) { idx = 1:sum(effect$A.msk) }
    } else {
      ng = effect$ngroup
      if ( is.null(ng) ) { ng = 1 }
      
      idx = INLA::inla.spde.make.index(effect$label, 
                                       n.spde = effect$inla.spde$n.spde, 
                                       n.group = ng)
    }
  }
  else {
    idx = map.effect(effect, data)
  }

  idx
  
}



mapper = function(map, points, eff, env = NULL) {
  
  # Evaluate the map with the points as an environment
  emap = tryCatch(eval(map, c(data.frame(points), as.list(env))), error = function(e) {})
  
  # 0) Eval failed. map everything to 1. This happens for automatically added Intercept
  # 1) If we obtain a function, apply the function to the points
  # 2) If we obtain a SpatialGridDataFrame extract the values of that data frame at the point locations using the over() function
  # 3) Else we obtain a vector and return as-is. This happens when map states a column of the data points
  
  if ( is.null(emap) ) { 
    loc = data.frame(x = rep(1, max(1,nrow(as.data.frame(points)))))
    colnames(loc) = deparse(map)
  }
  else if ( is.function(emap) ) { loc = emap(points) }
  else if ( inherits(emap, "SpatialGridDataFrame") | inherits(emap, "SpatialPixelsDataFrame")) {
    if ( length(eff$group.char) == 0 ) {
      loc = over(points, emap)[,1,drop=FALSE]
      colnames(loc) = eff$label
    } else {
      layr = points[[eff$group.char]]
      loc = vector()
      for (l in unique(layr)) { loc[layr == l] = over(points[layr == l,], emap)[,l] }
      loc = data.frame(loc = loc)
      colnames(loc) = eff$label
    }
    
  }
  else if ( eff$label == "offset" && is.numeric(emap) && length(emap)==1 ) { loc = data.frame(offset = rep(emap, nrow(points)))}
  else { loc = emap }
  
  # Check if any of the locations are NA. If we are dealing with SpatialGridDataFrame try
  # to fix that by filling in nearest neighbor values.
  if ( any(is.na(loc)) & ( inherits(emap, "SpatialGridDataFrame") | inherits(emap, "SpatialPixelsDataFrame") )) {
    
    warning(sprintf("Map '%s' has returned NA values. As it is a SpatialGridDataFrame I will try to fix this by filling in values of spatially nearest neighbors. In the future, please design your 'map=' argument as to return non-NA for all points in your model domain/mesh. Note that this can also significantly increase time needed for inference/prediction!",
                    deparse(map), eff$label))
    
    BADpoints = points[as.vector(is.na(loc)),]
    GOODpoints = points[as.vector(!is.na(loc)),]
    dst = rgeos::gDistance(SpatialPoints(GOODpoints),SpatialPoints(BADpoints), byid=T)
    nn = apply(dst, MARGIN = 1, function(row) which.min(row)[[1]])
    loc[is.na(loc)] = loc[!is.na(loc)][nn]
    colnames(loc) = eff$label
  }
  
  # Check for NA values.    
  if ( any(is.na(loc)) ) {
    stop(sprintf("Map '%s' of effect '%s' has returned NA values. Please design your 'map='
                 argument as to return non-NA for all points in your model domain/mesh.",
                 as.character(map)[[1]], eff$label))
  }
  
  loc
}



