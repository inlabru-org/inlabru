#' Models for \link{inlabru} 
#'
#' As with most inference methods inlabru models are set up via model formulae.
#' However, inlabru supports a slightly more advanced syntax which separates the name of an
#' effect from the name of a covariate. For instance, a very simple formula that inlabru
#' can work with is 
#' 
#' \code{y ~ myFixedEffect(temperature)}
#'
#' which defines a latent effect named 'myFixedEffect' with 'tepmerature' as a multiplier.
#' Hence, except for explicitly giving a name to the effect, the above formula is equivalent
#' to
#' 
#' \code{y ~ temperature}
#' 
#' In fact, the latter will also work in inlabru but the name of the effect will be
#' 'temperature'. The advantage of having an explicit name for the effect is that it 
#' allows for a more sophisticated syntax when making predictions.
#' 
#' In general, inlabru formula components are of the form
#' 
#' MyEffectName(map = aFunctionOfMyData, model = 'modelType', ...)
#' 
#' where 
#' 
#' \itemize{
#' \item{\code{MyEffectName} is an aribtrary string defining the name of the effect}
#' \item{\code{map} can be 
#' (A) a name of a covariate, 
#' (B) a function's name that will be called with your data as a parameter
#' (C) an expression that is valid with your data as an environment}
#' \item{model} is one of the \link{inla.models}
#' \item{... are further arguments passed that inla \link{f} understands}
#' }
#' 
#' Internally, inlabru will construct a \code{model} form the formula that you define
#' using the \link{make.model} function. Useful operators on \code{model} objects are:
#' 
#' \itemize{
#'  \item{\link{summary}: }{ Summarize an \link{inlabru} model }
#'  \item{\link{labels}: }{ Character array of effect labels }
#'  \item{\link{effect}: }{ Retrieve one or more effects from a model }
#'  \item{\link{default}: }{ Default linear predictor expression }
#'  \item{\link{list.A}: }{ List A-matrices that will be constructed for the INLA stack }
#'  \item{\link{list.data}: }{ List environmental variables needed to call INLA}
#'  \item{\link{list.indices}: }{ List indices that will be constructed for the INLA stack}
#'  \item{\link{list.effects}: }{ List weights that will be constructed for the INLA stack}
#'  \item{\link{evaluate}: }{ Determine the linear predictor's value at some location }
#'  
#' }
#' 
#' 
#' @name model
NULL


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

#' Create an inlabru \link{model} from a formula
#' 
#' The \link{inlabru} syntax for model forulae is different from what \link{inla} considers a valid.
#' In inla most of the effects are defined by adding an f(...) expression to the formula. 
#' In \link{inlabru} the f is replaced by an arbitrary (exception: 'offset') string that will
#' determine the label of the effect. For instance
#' 
#' \code{y ~ f(myspde, ...)}
#' 
#' is equivalent to
#' 
#' \code{y ~ myspde(...)}
#' 
#' A disadvantage of the inla way is that there is no clear separation between the name of the covariate
#' and the label of the effect. Furthermore, for some models like SPDE it is much more natural to
#' use spatial coordinates as covariates rather than an index into the SPDE vertices. For this purpose
#' \link{inlabru} provides the new \code{map} agument. For convenience, the map argument ca be used
#' like the first argument of the f function, e.g.
#' 
#' \code{y ~ f(temperature, model = 'fixed')}
#' 
#' is equivalent to
#' 
#' \code{y ~ temperature(map = temperature, model = fixed)}
#' as well as
#' \code{y ~ temperature(model = fixed)}
#' 
#' On the other hand, map can also be a function mapping, e.g the \link{coordinates} function of the
#' \link{sp} package :
#' 
#' \code{y ~ mySPDE(map = coordinates, ...)}
#'
#' Morevover, \code{map} can be any expression that evaluate within your data as an environment.
#' For instance, if your data has columns 'a' and 'b', you can create a fixed effect of 'a+b' by
#' setting \code{map} in the following way:
#' 
#' \code{y ~ myEffect(map = sin(a+b))} 
#'
#'
#' @export
#' @param fml A formula
#' @return A \link{model} object
#' 

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
    ge = eval(parse(text = lb), envir = environment(fml))

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

    lbl[[k]] = lb
    
    smod = c(ge$f, list(mesh = ge$mesh, 
                        A.msk = ge$A.msk,
                        mesh.coords = ge$f$term,
                        map = ge$map,
                        inla.spde = ge$model),
                        fchar = lb)
    
    # Fix label for factor effects
    if (smod$model == "factor") { lbl[[k]] = smod$label ; smod$fchar = smod$label }
    
    # Fix label for offset effect
    if (smod$label == "offset") lbl[[k]] = paste0("offset(",as.character(smod$map),")")
    
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


#' A wrapper for the \link{inla} \link{f} function
#' 
#' @aliases g
#' @export
#' @param covariate A string defining the label of the INLA effect. If \code{map} is provided this also sets the coariate used as a first argument to the \link{f} call.
#' @param map A name, call or function that maps points to effect indices or locations that the model understands.
#' @param model See \link{f} model specifications
#' @param mesh An \link{inla.mesh} object required for SPDE models
#' @param A.mask A boolean vector for masking A matrix columns. 
#' @param ... Arguments passed on to inla \link{f}
#' @return A list with mesh, model and the return value of the f-call

g = function(covariate, 
             map = NULL, 
             model = "linear", 
             mesh = NULL, 
             A.msk = NULL, ...){
 
  label = as.character(substitute(covariate))
  map.char = as.character(substitute(map))
  A.msk.char = as.character(substitute(A.msk))
  
  if ( length(map.char) == 0 ) { map = NULL } else { map = substitute(map) }
  if ( length(A.msk.char) == 0 ) { A.msk = NULL } else { A.msk = A.msk }
  
  # Only call f if we are  not dealing with an offset
  if ( label == "offset" ) { fvals = list(model="offset") }
  else if ( is.character(model) && model == "factor" ) { fvals = list(model="factor") } # , n = list(...)$n
  else { fvals = f(xxx, ..., model = model) }
  fvals$label = label
  
  # Default map
  if ( is.null(map) ) { 
    if ( fvals$model == "spde2") { map = expression(coordinates) }
    else { map = substitute(covariate) }
  }
  
  # Check if we got a mesh argument if the model is an SPDE
  if ( fvals$model == "spde2" & is.null(mesh) ) stop(sprintf("Model %s is an SPDE but no mesh was provided", label))
  
  
  ret = list(mesh = mesh, 
             map = map,
             A.msk = A.msk,
             model = model, 
             f = fvals)
}


##########################################################################
# OPERATORS ON MODELS
##########################################################################

#' Model summary
#'
#' @aliases summary.model
#' @export
#' @param model An \link{inlabru} \link{model}
#' 
summary.model = function(model) {
  for (label in elabels(model)) {
    eff = effect(model,label)
    en = ifelse("n" %in% names(eff), eff$n, 1)
    cat(paste0(label, ":\n"))
    cat(sprintf("\t model: %s , n = %s \n", eff$model, en))
    cat(sprintf("\t map: %s [class: %s] \n", deparse(eff$map), class(eff$map)))
    cat(sprintf("\t f call: %s \n", eff$fchar))
    
    cat("\n")
    
  }
  cat(paste0("--- FORMULA ---\n\n"))
  print(model$formula)
}


#' Retrieve effect labels
#'
#' @aliases elabels
#' @export
#' @param model An \link{inlabru} \link{model}
#' 
elabels = function(model) { names(model$effects) }


#' Retrieve effect properties
#'
#' @aliases effect
#' @export
#' @param model An \link{inlabru} \link{model}
#' @param label String defining the label of the effect. If NULL, return all effects
#' 
effect = function(model, label = NULL) {
  if (is.null(label)) { model$effects } 
  else { model$effects[[label]] }
}


#' Retrieve linear predictor expression
#'
#' @aliases predictor
#' @export
#' @param model An \link{inlabru} \link{model}
#'
default.predictor = function(model) { parse(text=paste0(elabels(model),collapse="+"))}



#' List data needed to run INLA
#'
#' @aliases list.data.model
#' @param model An \link{inlabru} \link{model}
#' @export
#' 

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

#' List of covariates effects needed to run INLA
#'
#' @aliases list.covariates.model
#' @param model An \link{inlabru} \link{model}
#' @param points A data fram to extract the data from
#' 

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

#' List of A matrices needed to run INLA
#'
#' @aliases list.A.model
#' @param model An \link{inlabru} \link{model}
#' @param points Locations to create the A-matrices for

list.A.model = function(model, points){
  A.lst = list()

  # For each effect create indices
  for ( eff in effect(model) ) {
    if ( is.null(eff$n) ) {
      A.lst[[eff$label]] = Matrix::Diagonal(nrow(data.frame(points)))
    } else if ( is.null(eff$mesh) ) {
      idx = mapper(eff$map, points, eff)
      A = matrix(0, nrow = nrow(data.frame(points)), ncol=eff$n)
      A[cbind(1:nrow(A), idx)] = 1
      A.lst[[eff$label]] = A
    } else {
      
      loc = mapper(eff$map, points, eff)
      loc = as.matrix(loc) # inla.spde.make.A requires matrix format as input

      if (is.null(eff$ngroup)) { ng = 1 } else { ng = eff$ngroup }
      if (ng > 1) {
        group = as.matrix(eff$group)
      } else { group = NULL }
      A = inla.spde.make.A(eff$mesh, loc = loc, group = group)
      # Mask columns of A
      if (!is.null(eff$A.msk)) { A = A[, as.logical(eff$A.msk), drop=FALSE]}
      # Weights for models with A-matrix are realized in the follwoing way:
      if (!is.null(eff$weights)) {
        w = eff$weights
        A = as.matrix(A)*as.vector(w)
      }
      A.lst[[eff$label]] = A

    }
  }
  
  return(A.lst)
}


#' List of spde indexing effects needed to run INLA
#'
#' @aliases list.indices.model
#' @param model An \link{inlabru} \link{model}
#' @param points A data frame to extract indices from

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
          ng = eff$inla.spde[[name]]$n.group
          if (is.null(ng)) { ng = 1 }
          idx[[name]] = inla.spde.make.index(name, n.spde = eff$inla.spde$n.spde, n.group = ng)
        }
      }
    } else {
      idx[[name]] = mapper(eff$map, points, eff)
    }
    
  }
  idx
}

#' Evaluate or sample from a posterior result given a model and locations
#' 
#' @aliases evaluate.model evaluate
#' @export
#' @param model An \link{inlabru} \link{model}
#' @param result Posterior of an \link{inla}, \link{bru} or \link{lgcp} run.
#' @param points Locations and covariates needed to evaluate the model.
#' @param predictor An expression to be ealuated given the posterior or for each sample thereof. The default (\code{NULL}) returns a \code{data.frame} containing the sampled effects.
#' @param property Property of the model compnents to obtain value from. Default: "mode". Other options are "mean", "0.025quant", "0.975quant", "sd" and "sample". In case of "sample" you will obtain samples from the posterior (see \code{n} parameter).
#' @param n Number of samples to draw.
#' 
evaluate.model = function(model, 
                          result, 
                          points,
                          predictor = NULL,
                          property = "mode",
                          n = 1) {
  
  # Do we otain our values from sampling or from a property of a summary?
  if ( property == "sample") {
      smp = inla.posterior.sample.structured(result, n = n) 
  } else {
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
      if (!is.null(A[[name]])) {
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
    
    # For factor models we need to do a little trick
    for (label in names(sm)) { 
      if (label %in% names(effect(model)) && effect(model)[[label]]$model == "factor") {
        fc = mapper(effect(model)[[label]]$map, points, effect(model)[[label]]) 
        sm[[label]] = as.vector(sm[[label]][fc])
      }
    }
 
    # Apply A matrices
    for (label in names(sm)) { if (label %in% vars) sm[[label]] = apply.A(label, sm) }
    
    # If no predictor is provided simply return the samples. 
    # Otherwise evaluate predictor with each sample as an environment
    if ( is.null(predictor) ) {
      smp[[k]] = data.frame(sm)
    } else {
      smp[[k]] = eval(predictor, envir = c(sm, as.list(data.frame(points))))
    }
  }
  
  # Return
  if ( property == "sample") { smp } else { smp[[1]] }
}


mapper = function(map, points, eff) {
  # If map is a function, return map(points)
  # If map is not a function but a column in points, return points$map
  # Otherwise assume map is an expression that can be evaluated with points as an environment

  if (is.function(map)) { loc = map(points) } 
  else if (class(map) == "call") {
    if ( inherits(eval(map), "SpatialGridDataFrame") ){
      loc = over(points, eval(map))
      if (eff$label %in% names(loc)) { loc = loc[[eff$label]] }
    } else {
      loc = tryCatch( eval(map, data.frame(points)), error = function(e){
        loc = data.frame(points)[[eff$label]]
      })
    }
  } 
  else {
    fetcher = get0(as.character(map))
    if (is.function(fetcher)) { loc = fetcher(points) } 
    else {
      if(as.character(map) %in% names(as.data.frame(points))) {
        loc = as.data.frame(points)[,as.character(map)] 
      } else {
        loc = rep(1, nrow(data.frame(points))) 
      }
    }
  }
}
