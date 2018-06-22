# Internal \link{inlabru} model structure
# 
# See \link{make.model}.
# 
# @name model
#NULL


##########################################################################
# GENERICS
##########################################################################

evaluate = function(...){UseMethod("evaluate")}
list.A = function(...){UseMethod("list.A")}
list.covariates = function(...){UseMethod("list.covariates")}
list.indices = function(...){UseMethod("list.indices")}
list.data = function(...){UseMethod("list.data")}


##########################################################################
# Constructor
##########################################################################

# Create an inlabru \link{model} from a formula
# 
# The \link{inlabru} syntax for model forulae is different from what \link{inla} considers a valid.
# In inla most of the effects are defined by adding an f(...) expression to the formula. 
# In \link{inlabru} the f is replaced by an arbitrary (exception: 'offset') string that will
# determine the label of the effect. For instance
# 
# \code{y ~ f(myspde, ...)}
# 
# is equivalent to
# 
# \code{y ~ myspde(...)}
# 
# A disadvantage of the inla way is that there is no clear separation between the name of the covariate
# and the label of the effect. Furthermore, for some models like SPDE it is much more natural to
# use spatial coordinates as covariates rather than an index into the SPDE vertices. For this purpose
# \link{inlabru} provides the new \code{map} agument. For convenience, the map argument ca be used
# like the first argument of the f function, e.g.
# 
# \code{y ~ f(temperature, model = 'fixed')}
# 
# is equivalent to
# 
# \code{y ~ temperature(map = temperature, model = fixed)}
# as well as
# \code{y ~ temperature(model = fixed)}
# 
# On the other hand, map can also be a function mapping, e.g the \link{coordinates} function of the
# \link{sp} package :
# 
# \code{y ~ mySPDE(map = coordinates, ...)}
#
# Morevover, \code{map} can be any expression that evaluate within your data as an environment.
# For instance, if your data has columns 'a' and 'b', you can create a fixed effect of 'a+b' by
# setting \code{map} in the following way:
# 
# \code{y ~ myEffect(map = sin(a+b))} 
#
#
# @export
# @param fml A formula
# @return A \link{model} object
# 

make.model = function(fml) {
  submodel = list()
  covariates = list()
  
  # Define function for shifting the effect name into the g-call
  shift.names = function(fml) {
    tms = terms(fml)
    labels = attr(tms, "term.labels")

    # Check for offset()
    isoff = as.vector(unlist(lapply(rownames(attr(tms, "factors")), function(s) substr(s,1,6)=="offset")))
    if (any(isoff)) {
      labels[[length(labels)+1]] = rownames(attr(tms, "factors"))[isoff]
    }
    
    for (k in 1:length(labels)){
      lb = labels[[k]]
      
      # Function syntax or fixed effect?
      ix = regexpr("(", text = lb, fixed = TRUE)
      if (ix > 0) {
        ename = substr(lb, 1, ix-1)
        is.fixed = FALSE
      } else {
        ename = lb
        is.fixed = TRUE
      }
      
      # Construct g() call
      if ( is.fixed ) {
        labels[[k]] = sprintf("g(%s, map = %s, model = 'linear')", ename, ename)
      } else {
        if ( !(ename %in% c("g","f"))) {
          labels[[k]] = gsub(paste0(ename,"("), paste0("g(",ename,", "), lb, fixed = TRUE)
        }
        if ( ename %in% c("f") ) {
          labels[[k]] = gsub(paste0(ename,"("), paste0("g("), lb, fixed = TRUE)
        }
        
      }
    }
    
    labels
  }
  

  lbl = shift.names(fml)
  
  for ( k in 1:length(lbl) ) {
    
    lb = lbl[[k]]

    # Call g()
    # We have to add g() to the environment because it is not exported by inlabru and
    # therefore note visible within the parent environment.
    env = environment(fml)
    env$g = g
    ge = eval(parse(text = lb), envir = env)

    # Replace function name by INLA f function
    if ( substr(lb,1,2) == "g(" ) { lb = paste0("f(", substr(lb, 3, nchar(lb)))}

    # Remove extra mesh argument
    lb = gsub("[,][ ]*mesh[ ]*=[^),]*", "", lb)

    # Remove extra map argument
    pat = paste0("", deparse(ge$map))
    lb = gsub(paste0(" ", deparse(ge$map)), "", lb, fixed =TRUE)
    lb = gsub(paste0("=", deparse(ge$map)), "", lb, fixed =TRUE)
    lb = gsub("[,][ ]*map[ ]*=[^),]*", "", lb)
    
    
    # Remove extra A.msk argument
    pat = paste0("", deparse(ge$A.msk))
    lb = gsub(paste0(" ", deparse(ge$A.msk)), "", lb, fixed =TRUE)
    lb = gsub(paste0("=", deparse(ge$A.msk)), "", lb, fixed =TRUE)
    lb = gsub("[,][ ]*A.msk[ ]*=[^),]*", "", lb)

    # For SPDE models replace group=XXX by group=effectname.group
    if (ge$f$model == "spde2") {
      lb = gsub("[,][ ]*group[ ]*=[^),]*", paste0(", group = ",ge$f$label,".group"), lb)
    }
    
    
    lbl[[k]] = lb
    
    smod = c(ge$f, list(mesh = ge$mesh, 
                        A.msk = ge$A.msk,
                        mesh.coords = ge$f$term,
                        map = ge$map,
                        group.char = ge$group.char,
                        inla.spde = ge$model),
                        fchar = lb)
    
    # For copy model extract the mesh from the model that is copied
    if ( !is.null(smod$of) ) {
      smod$mesh = submodel[[smod$of]]$mesh
      smod$n = submodel[[smod$of]]$n
      smod$inla.spde = submodel[[smod$of]]$inla.spde
    }
    
    # Fix label for factor effects
    if (smod$model == "factor") { lbl[[k]] = smod$label ; smod$fchar = smod$label }
    
    # Fix label for offset effect
    # if (smod$label == "offset") lbl[[k]] = paste0("offset(",deparse(smod$map),")")
    if (smod$label == "offset") lbl[[k]] = "offset(offset)"
    
    # Set submodel
    submodel[[ge$f$label]] = smod
    
  }
  
  # Did the old formula have an intercept?
  if ( attr(terms(fml), "intercept") == 0 ) {icpt = "-1"} else { icpt = ""}
  
  # Combine labels into new formula
  new.rhs = as.formula(paste0(". ~ ", paste0(lbl, collapse = " + "), icpt))
  new.fml = update.formula(fml, new.rhs)
  environment(new.fml) = environment(fml)
  
  # Return
  ret = list(effects = submodel, formula = new.fml, in.formula = fml)
  class(ret) = c("model","list")
  ret
}


# A wrapper for the \link{inla} \link{f} function
# 
# @aliases g
# @export
# @param covariate A string defining the label of the INLA effect. If \code{map} is provided this also sets the coariate used as a first argument to the \link{f} call.
# @param map A name, call or function that maps points to effect indices or locations that the model understands.
# @param group A name, call or function that maps the data to groups
# @param model See \link{f} model specifications
# @param mesh An \link{inla.mesh} object required for SPDE models
# @param A.msk A boolean vector for masking A matrix columns. 
# @param ... Arguments passed on to inla \link{f}
# @return A list with mesh, model and the return value of the f-call

g = function(covariate, 
             map = NULL,
             group = NULL,
             model = "linear", 
             mesh = NULL, 
             A.msk = NULL, ...){
 
  label = as.character(substitute(covariate))
  map.char = as.character(substitute(map))
  group.char = as.character(substitute(group))
  A.msk.char = as.character(substitute(A.msk))
  is.copy = !is.null(list(...)$copy)
  
  if ( length(map.char) == 0 ) { map = NULL } else { map = substitute(map) }
  if ( length(A.msk.char) == 0 ) { A.msk = NULL } else { A.msk = A.msk }
  
  # Only call f if we are  not dealing with an offset
  if ( label == "offset" ) { fvals = list(model="offset") }
  else if ( is.character(model) && model == "factor" ) { fvals = list(model="factor") } # , n = list(...)$n
  else { xxx = NULL ; fvals = INLA::f(xxx, ..., group = group, model = model) }
  fvals$label = label
  
  # Default map
  if ( is.null(map) ) { 
    if ( fvals$model == "spde2" && class(mesh) == "inla.mesh") { map = expression(coordinates) }
    else { map = substitute(covariate) }
  }
  
  # Check if we got a mesh argument if the model is an SPDE
  if ( fvals$model == "spde2" & is.null(mesh) ) { mesh = model$mesh } 
  
  # Check if n is present for models that are not fixed effects
  if (is.null(fvals$n) & !(fvals$model == "linear") & !(fvals$model == "offset") & !(fvals$model == "factor") & !is.copy) {
    stop(sprintf("Please provide parameter 'n' for effect '%s'", label))
    }
  
  
  ret = list(mesh = mesh, 
             map = map,
             group.char = group.char,
             A.msk = A.msk,
             model = model, 
             f = fvals)
}


##########################################################################
# OPERATORS ON MODELS
##########################################################################

# Model summary
#
# @aliases summary
# @export
# @param object An \link{inlabru} \link{model}
# @param ... ignored arguments (S3 generic compatibility)
# 
summary.model = function(object, ...) {
  for (label in elabels(object)) {
    eff = effect(object,label)
    en = ifelse("n" %in% names(eff), eff$n, 1)
    cat(paste0(label, ":\n"))
    cat(sprintf("\t model: %s , n = %s \n", eff$model, en))
    cat(sprintf("\t map: %s [class: %s] \n", deparse(eff$map), class(eff$map)))
    cat(sprintf("\t f call: %s \n", eff$fchar))
    
    cat("\n")
    
  }
  cat(paste0("--- FORMULA ---\n\n"))
  print(object$formula)
}


# Retrieve effect labels
#
# @aliases elabels
# @export
# @param model An \link{inlabru} \link{model}
# 
elabels = function(model) { names(model$effects) }


# Retrieve effect properties
#
# @aliases effect
# @export
# @param model An \link{inlabru} \link{model}
# @param label String defining the label of the effect. If NULL, return all effects
# 
effect = function(model, label = NULL) {
  if (is.null(label)) { model$effects } 
  else { model$effects[[label]] }
}


# Retrieve linear predictor expression
#
# @aliases predictor
# @export
# @param model An \link{inlabru} \link{model}
#
default.predictor = function(model) { parse(text=paste0(elabels(model),collapse="+"))}



# List data needed to run INLA
#
# @aliases list.data.model
# @param model An \link{inlabru} \link{model}
# @export
# 

list.data.model = function(model){
  
  # Formula environment as list
  elist = as.list(environment(model$formula))
  
  # Remove previous inla results. For some reason these slow down the next INLA call.
  elist = elist[unlist(lapply(elist, function(x) !inherits(x, "inla")))]
  
  # Remove functions. This can cause problems as well
  elist = elist[unlist(lapply(elist, function(x) !is.function(x)))]
  
  # Remove functions. This can cause problems as well
  elist = elist[unlist(lapply(elist, function(x) !inherits(x, "formula")))]
  
  elist = elist[names(elist) %in% all.vars(model$formula)]
}

# List of covariates effects needed to run INLA
#
# @aliases list.covariates.model
# @param model An \link{inlabru} \link{model}
# @param points A data fram to extract the data from
# 

list.covariates.model = function(model, points){
    
  all.varnames = c(all.vars(model$formula), all.vars(model$expr)) 
  covar.data = list()
    
    # Points may be annotated with covariates
    
    for (vname in all.varnames){
      if (vname %in% names(data.frame(points))) {
        covar.data[[vname]] = data.frame(points)[,vname]
      }
    }
    covar.data = data.frame(do.call(cbind,covar.data))
}

# List of A matrices needed to run INLA
#
# @aliases list.A.model
# @param model An \link{inlabru} \link{model}
# @param points Locations to create the A-matrices for

list.A.model = function(model, points){
  A.lst = list()

  # For each effect create indices
  for ( eff in effect(model) ) {
    if ( is.null(eff$n) ) {
      A.lst[[eff$label]] = Matrix::Diagonal(nrow(data.frame(points)))
    } else if ( is.null(eff$mesh) ) {
      idx = mapper(eff$map, points, eff)
      if (is.data.frame(idx)) idx = idx[,1]
      N <- nrow(data.frame(points))
      A = Matrix::sparseMatrix(i = seq_len(N),
                               j = idx,
                               x = rep(1, N),
                               dims = c(N, eff$n))
      A.lst[[eff$label]] = A
    } else {
      
      if ( eff$map == "coordinates" ) {
        if ( is.na(proj4string(points)) | is.null(eff$mesh$crs) ) {
          loc = points
        } else {
          loc = stransform(points, crs = eff$mesh$crs)
        }
        
      } else {
        loc = mapper(eff$map, points, eff)
        if ( !is.matrix(loc) & !inherits(loc,"Spatial") ) loc = as.matrix(loc)
      }
  
      if (is.null(eff$ngroup)) { ng = 1 } else { ng = eff$ngroup }
      if (ng > 1) {
        group = points[[eff$group.char]]
      } else { group = NULL }
      
      A = INLA::inla.spde.make.A(eff$mesh, loc = loc, group = group, n.group = ng)
      # Mask columns of A
      if (!is.null(eff$A.msk)) { A = A[, as.logical(eff$A.msk), drop=FALSE]}
      # Weights for models with A-matrix are realized in the following way:
      if (!is.null(eff$weights)) {
        w = as.vector(eff$weights)
        A = Matrix::Diagonal(length(w), w) %*% A
      }
      A.lst[[eff$label]] = A

    }
  }
  
  return(A.lst)
}


# List of spde indexing effects needed to run INLA
#
# @aliases list.indices.model
# @param model An \link{inlabru} \link{model}
# @param points A data frame to extract indices from

list.indices.model = function(model, points){
  
  # The list of indices to be returned
  idx = list()

  # For each effect create indices
  for ( eff in effect(model) ) {
    name = eff$label
    
    # Only create indices if the number n of effect variables is known
    if (!is.null(eff$n)){
      
      # Effects with meshes need a special treatment ...
      if ( is.null(eff$mesh) ) { 
        idx[[name]] = seq_len(eff$n) 
      } else {
        if ( "m" %in% names(eff$mesh) ) {
          idx[[name]] = 1:eff$mesh$m # support inla.mesh.1d models
          # If a is masked, correct number of indices
          if (!is.null(eff$A.msk)) { idx[[name]] = 1:sum(eff$A.msk)}
        } else {
          ng = eff$ngroup
          if (is.null(ng)) { ng = 1 }
          idx[[name]] = INLA::inla.spde.make.index(name, n.spde = eff$inla.spde$n.spde, n.group = ng)
          # NOTE: MODELS WITH REPLICATES ARE NOTE IMPLEMENTED YET.
          #      - when using group= or replicate, the f-formula has to be adjusted to read the groups/replicates from the stack
          #      - see make.model() for an example on how to do this for groups
        }
      }
    } else {
      idx[[name]] = mapper(eff$map, points, eff, environment(model$in.formula))
    }
    
  }
  idx
}

# Evaluate or sample from a posterior result given a model and locations
# 
# @aliases evaluate.model evaluate
# @export
# @param model An \link{inlabru} \link{model}
# @param result Posterior of an \link{inla}, \link{bru} or \link{lgcp} run.
# @param points Locations and covariates needed to evaluate the model.
# @param predictor A formula or an expression to be evaluated given the posterior or for each sample thereof. The default (\code{NULL}) returns a \code{data.frame} containing the sampled effects. In case of a formula the right hand side is used for evaluation.
# @param property Property of the model compnents to obtain value from. Default: "mode". Other options are "mean", "0.025quant", "0.975quant", "sd" and "sample". In case of "sample" you will obtain samples from the posterior (see \code{n} parameter).
# @param n Number of samples to draw.
# 
evaluate.model = function(model, 
                          result, 
                          points,
                          predictor = NULL,
                          property = "mode",
                          n = 1) {
  
  data = points # Within the evaluation make points available via the name "data" 
  
  
  if ( inherits(predictor, "formula") ) {
    fml.envir = as.list(environment(predictor))
    predictor = parse(text = as.character(predictor)[length(as.character(predictor))])
  } else { fml.envir = list() }
  
  # Do we otain our values from sampling or from a property of a summary?
  if ( property == "sample") {
      smp = inla.posterior.sample.structured(result, n = n) 
  } else {
      result$model = model
      smp = rep(list(extract.summary(result, property)), n)
  }
  
  # Which effects do we want? Remove effect from model that are not required for the evaluation
  if ( is.null(predictor) ) {
    vars = setdiff(names(smp[[1]]), "Predictor")
  } else {
    effs = effect(model)
    model$effects = effs[intersect(names(effs), all.vars(predictor))]
    vars = intersect(names(smp[[1]]), all.vars(predictor))
  }

  # Obtain A-matrices
  A = list.A.model(model, points)
  
  # Make a function that will apply the A-matrices
  apply.A = function(name, s) {
      mult = as.vector(s[[name]])
      if ( !is.null(A[[name]]) && ( nrow(A[[name]]) > 0) ) {
        if (length(mult) == 1) { as.vector(A[[name]] %*% rep(mult, nrow(A[[name]]))) 
          }
        else { as.vector(A[[name]] %*% s[[name]]) }
      } else { 
        mult 
      }
  }
  

  for ( k in 1:n ) {
    # Discard variables we do not need
    sm = smp[[k]][vars]
    
    # Factor models
    for (label in names(sm)) { 
      if (label %in% names(effect(model)) && effect(model)[[label]]$model == "factor") {
        fc = mapper(effect(model)[[label]]$map, points, effect(model)[[label]], environment(model$in.formula)) 
        sm[[label]] = as.vector(sm[[label]][as.data.frame(fc)[,1]])
      }
    }
    
    # Linear models
    for (label in names(sm)) { 
      if (label %in% names(effect(model)) && effect(model)[[label]]$model == "linear") {
        fc = mapper(effect(model)[[label]]$map, points, effect(model)[[label]], environment(model$in.formula)) 
        if (is.data.frame(fc)) {fc = fc[,1]}
        sm[[label]] = as.vector(sm[[label]] * as.vector(fc))
      }
    }
 
    # Apply A matrices
    for (label in names(sm)) { if (label %in% vars) sm[[label]] = apply.A(label, sm) }
    
    # If no predictor is provided simply return the samples. 
    # Otherwise evaluate predictor with each sample as an environment
    if ( is.null(predictor) ) {
      smp[[k]] = data.frame(sm)
    } else {
      envir = c(sm, as.list(data.frame(points)), fml.envir, as.list(environment(model$formula)))
      smp[[k]] = eval(predictor, envir = envir)
    }
  }
  
  # Return
  if ( property == "sample") { smp } else { smp[[1]] }
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
