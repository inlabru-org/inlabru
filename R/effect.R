####################################################################################################
# GENERICS
####################################################################################################

#' @export
#' @rdname value.effect
value = function(...){UseMethod("value")}
#' @export
#' @rdname amatrix.effect
amatrix = function(...){UseMethod("amatrix")}
#' @export
#' @rdname map.effect
map = function(...){UseMethod("map")}
#' @export
#' @rdname index.effect
index = function(...){UseMethod("index")}


####################################################################################################
# CONSTRUCTORS
####################################################################################################


#' Latent effects
#'  
#' @aliases effect
#' @export
#' @param ... EXPERIMENTAL
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'
effect = function(...){UseMethod("effect")}


#' Latent effects
#'  
#' @aliases effect.formula
#' @export
#' @method effect formula
#' @param formula A formula defining the effect.
#' @param ... EXPERIMENTAL
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'
effect.formula = function(formula, ...) {
  code = code.components(formula)
  parsed = lapply(code, function(x) parse(text=x))
  effects = lapply(parsed, function(x) eval(x, envir = environment(formula)))
  names(effects) = lapply(effects, function(x) x$label)
  
  #if ( length(effects)==1 ) effects = effects[[1]]
  effects
}

#' Latent effects
#'  
#' @aliases effect.character
#' @export
#' @method effect character
#' @param label EXPERIMENTAL
#' @param data EXPERIMENTAL
#' @param model EXPERIMENTAL
#' @param map EXPERIMENTAL
#' @param n EXPERIMENTAL
#' @param season.length EXPERIMENTAL
#' @param group EXPERIMENTAL
#' @param values EXPERIMENTAL
#' @param A.msk EXPERIMENTAL
#' @param ... EXPERIMENTAL
#' 
#' 
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'
effect.character = function(label,
                             data,
                             model,
                             map,
                             n = NULL,
                             season.length = NULL,
                             group = NULL,
                             values = NULL,
                             A.msk = NULL,
                             ...){
  
  # INLA models:
  # itypes = c(linear, iid, mec, meb, rgeneric, rw1, rw2, crw2, seasonal, besag, besag2, bym, bym2, besagproper, 
  #            besagproper2, fgn, fgn2, ar1, ar1c, ar, ou, generic, generic0, generic1, generic2, generic3, spde, 
  #            spde2, spde3, iid1d, iid2d, iid3d, iid4d, iid5d, 2diid, z, rw2d, rw2diid, slm, matern2d, copy,
  #            clinear, sigm, revsigm, log1exp, logdist) 
  #
  # Supported:
  # btypes = c("offset", "factor", "linear", "clinear", "iid", "seasonal", "rw1", "rw2", "ar", "ar1", "ou", "spde")
  
  model.type = model
  if ( inherits(model, "inla.spde") ) { model.type = "spde" }
  
  miss.msg = "Effect '%s' (type '%s') requires argument '%s'. Check out f() for additional information on this argument."
  
  # map.char = as.character(substitute(map))
  group.char = as.character(substitute(group))
  # A.msk.char = as.character(substitute(A.msk))
  
  effect = list(label = label,
                type = model.type,
                map = substitute(map),
                group.char = group.char, # Name of the data column holding the group index
                values = values,
                A.msk = A.msk,
                model = model)
  

  if ( model.type %in% c("offset") ) {
    effect$inla.formula = as.formula(paste0("~ . + offset(offset)"))
  } 
  else if ( model.type %in% c("factor") ) {
    effect$inla.formula = as.formula(paste0("~ . + ", label))
  }
  else {
    
    fcall = sys.call()
    fcall[[1]] = "f"
    fcall[[2]] = ""
    fcall[[2]] = as.symbol(label)
    fcall = fcall[!(names(fcall) %in% c("map","A.msk", "mesh"))]
    # For SPDE models we need a little nasty trick
    if ( model.type %in% c("spde") ) {
      tmp = NA
      fcall$group = as.symbol("tmp") 
    }
    fvals = do.call(INLA::f, as.list(fcall[2:length(fcall)])) # eval(parse(text=deparse(fcall)))
    effect$f = fvals
    # Second part of the trick above
    if ( model.type %in% c("spde") ) { 
      fcall$group = as.symbol(paste0(label, ".group"))
    }
    effect$inla.formula = as.formula(paste0("~ . + ", paste0(deparse(fcall), collapse = "")))
    
    

    if ( model.type %in% c("linear", "clinear") ) {
      effect$mesh = NA
    }
    
    else if ( model.type %in% c("iid") ) {
      if ( missing(n) ) { stop(sprintf(miss.msg, effect$label, model.type, "n")) }
      effect$mesh = INLA::inla.mesh.1d(1:fvals$n)
    }
    else if ( model.type %in% c("seasonal") ) {
      if ( missing(season.length) ) { stop(sprintf(miss.msg, effect$label, model.type, "season.length")) }
      effect$mesh = INLA::inla.mesh.1d(1:fvals$season.length)
    }
    else if ( model.type %in% c("rw1", "rw2", "ar", "ar1", "ou") ) {
      if ( missing(values) ) { stop(sprintf(miss.msg, effect$label, model.type, "values")) }
      effect$mesh = INLA::inla.mesh.1d(sort(unique(fvals$values)))
    }
    else if ( model.type %in% c("spde") ) {
      effect$mesh = model$mesh
      if (is.null(effect$f$ngroup)) { effect$ngroup = 1 } else { effect$ngroup = effect$f$ngroup }
      if (is.null(effect$f$nrep)) { effect$nrep = 1 } else { effect$nrep = effect$f$nrep }
    } 
    else {
      stop(paste0("Effect type '", model.type, "' not implemented."))
    }
  }

  
  class(effect) = c("effect","lists")
  effect
}



#' Convert components to R code
#'  
#' @aliases code.components
#' @keywords internal
#' @param components A \link{formula} describing latent model components.
#' @param fname Chracter setting the name of the function that will interpret the components.
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

code.components = function(components) {
  fname = "effect"
  tms = terms(components)
  codes = attr(tms, "term.labels")
  
  # Check for offset()
  isoff = as.vector(unlist(lapply(rownames(attr(tms, "factors")), function(s) substr(s,1,6)=="offset")))
  if (any(isoff)) {
    codes[[length(codes)+1]] = rownames(attr(tms, "factors"))[isoff]
  }
  
  for (k in 1:length(codes)){
    code = codes[[k]]
    
    # Function syntax or fixed effect?
    ix = regexpr("(", text = code, fixed = TRUE)
    is.offset = FALSE
    if (ix > 0) {
      label = substr(code, 1, ix-1)
      is.fixed = FALSE
      if (label == "offset") {is.offset = TRUE} 
    } else {
      label = code
      is.fixed = TRUE
    }
    
    # Make code
    if ( is.fixed ) {
      codes[[k]] = sprintf("%s(\"%s\", map = %s, model = 'linear')", fname, label, label)
    }
    else if ( is.offset ) {
      codes[[k]] = gsub(paste0(label,"("), paste0(fname,"(\"",label,"\", model = \"offset\", map = "), code, fixed = TRUE)
    }
    else {
      codes[[k]] = gsub(paste0(label,"("), paste0(fname,"(\"",label,"\", "), code, fixed = TRUE)
    }

    
  }
  
  codes
}


####################################################################################################
# OPERATORS
####################################################################################################


#' Summarize an effect
#'  
#' @aliases summary.effect
#' @export
#' @method summary effect
#' @keywords internal
#' @param object An effect.
#' @param ... ignored.
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

summary.effect = function(object, ...) {
  
  eff <- list('Label' = object$label,
              'Type' = object$type,
              'Map' = sprintf("%s [class: %s]",
                              deparse(object$map),
                              class(object$map)),
              'f call' = object$fchar)
  class(eff) <- c('summary.effect', 'list')
  eff
}

#' Print the summary of an effect
#'  
#' @aliases print.summary.effect
#' @export
#' @method print summary.effect
#' @keywords internal
#' @param x A 'summary.effect' object.
#' @author Finn Lindgren <\email{finn.lindgren@@gmail.com}>
#' @rdname summary.effect

print.summary.effect = function(x, ...) {
  for (name in names(x)) {
    # Split TAB character to attempt proper printing in RStudio,
    # but even though this makes a difference on the command line,
    # there's no effect when inside the function. /FL
    cat(name, ":", "\t", x[[name]], "\n", sep="")
  }
  invisible(x)
}



#' Evaluate an effect
#' 
#' Calculates a latent effect given some data and the state of the effect's internal random variables.
#' 
#' TODO: Improve speed for iterated calls by making 'mapped' a parameter 
#' 
#' @aliases value.effect
#' @export
#' @method value effect
#' @keywords internal
#' @param effect An effect.
#' @param data A \code{data.frame} or Spatial* object of covariates and/or point locations.
#' @param state Either a numeric vector or a list with a numeric entry whose name is equal to the name parameter.
#' @param A A matrix overriding the default projection matrix.
#' @param ... Unused.
#' @return A numeric vector of effect values
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'


value.effect = function(effect, data, state, A = NULL, ...) {
  
  # Convenience: extract state if a list of states was provided
  if ( is.list(state) &!is.data.frame(state)) { state = state[[effect$label]] }
  
  # Obtain covariates
  mapped = mapper(effect$map, data, effect, effect$env)
  
  if (is.data.frame(mapped)) { 
    if (effect$label %in% names(mapped)) {
      mapped = mapped[,effect$label,drop=TRUE]
    } else {
      mapped = mapped[,1,drop=TRUE] 
    }
  }
    
  # Make A-matrix (if not provided)
  if ( is.null(A) ) { A = amatrix(effect, data) }
  
  # Determine effect depending on the type of latent model
  if ( effect$type %in% c("linear", "clinear") ) {

    values = A %*% (state * mapped)
  }
  else if ( effect$type %in% c("offset") ) {

    values = A %*% mapped
  }
  else if ( effect$type %in% c("factor") ) {

    values = A %*% state[mapped]
  } 
  else if ( effect$type %in% c("iid", "seasonal") ) {
    
    values = A %*% state[mapped]
    
  }
  else if ( effect$type %in% c("rw1", "rw2", "ar", "ar1", "ou") ) {
    
    values = A %*% state
    
  }
  else if ( effect$type %in% c("spde") ) {

    values = A %*% state
    
  } else {
    stop(paste0("Evaluation of ", effect$type, " not implemented."))
  }
  
  as.vector(values)
}




#' Construct A-matrix
#' 
#' Constructs the A-matrix for a given effect and and some data
#'  
#' @aliases amatrix.effect
#' @export
#' @method amatrix effect
#' @keywords internal
#' @param effect An effect.
#' @param data A \code{data.frame} or Spatial* object of covariates and/or point locations.
#' @param ... Unused.
#' @return An A-matrix.
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

amatrix.effect = function(effect, data, ...) {
  
  if ( effect$type %in% c("spde") ) {
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
    
    
    if (effect$ngroup > 1) {
      group = data[[effect$group.char]]
    } else { group = NULL }
    
    A = INLA::inla.spde.make.A(effect$mesh, 
                               loc = loc, 
                               group = group, 
                               n.group = effect$ngroup, 
                               n.repl = effect$nrepl)
  } 
  else if ( effect$type %in% c("rw1", "rw2", "ar", "ar1", "ou") ) {
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
#' @export
#' @method map effect
#' @keywords internal
#' @param effect An effect.
#' @param data A \code{data.frame} or Spatial* object of covariates and/or point locations. If null, return the effect's map.
#' @param ... Unused.
#' @return An A-matrix.
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

map.effect = function(effect, data, ...) {
  
  mp = mapper(effect$map, data, effect, effect$env)

}



#' Obtain indices
#' 
#' Idexes into to the effects
#'  
#' @aliases index.effect
#' @export
#' @method index effect
#' @keywords internal
#' @param effect An effect.
#' @param data A \code{data.frame} or Spatial* object of covariates and/or point locations. If null, return the effect's map.
#' @param ... Unused.
#' @return a data.frame of indices or list of indices into the effects latent variables
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

index.effect = function(effect, data, ...) {
  
  if ( effect$type %in% c("spde") ) {
    if ( "m" %in% names(effect$mesh) ) {
      # idx = 1:effect$mesh$m # support inla.mesh.1d models
      idx = INLA::inla.spde.make.index(name = effect$label, n.spde = effect$mesh$m)
      # If a is masked, correct number of indices
      if ( !is.null(effect$A.msk) ) { idx = 1:sum(effect$A.msk) }
    } else {
      
      idx = INLA::inla.spde.make.index(effect$label, 
                                       n.spde = effect$model$n.spde, 
                                       n.group = effect$ngroup,
                                       n.repl = effect$nrep)
    }
  }
  else if ( effect$type %in% c("factor") ) {
    idx = map.effect(effect, data)
    if (!is.data.frame(idx)) { idx = data.frame(idx) ; colnames(idx) = effect$label }
    idx[,1] = as.factor(paste0(effect$label, idx[,1]))
  }
  else {
    idx = map.effect(effect, data)
    if (!is.data.frame(idx)) {
      idx = data.frame(idx)
      colnames(idx) = effect$label
    }
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
  
  # Always return as many rows as points has
  if (is.vector(loc) && (length(loc) == 1) && (nrow(as.data.frame(points))>1)) {
    loc = rep(loc, nrow(as.data.frame(points)))
  }
  
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



