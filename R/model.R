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

make.model = function(fml){
  
  # Create effects
  effects = effect(fml)
  
  # Create joint formula that will be used by inla
  formula = y ~ 1
  for (fm in lapply(effects, function(eff) {eff$inla.formula})) {
    formula = update.formula(formula, fm)
  }
  
  mdl = list(effects = effects, formula = formula, in.formula = fml)
  class(mdl) = c("model","list")
  return(mdl)
  
}


make.model.old = function(fml) {
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
    
    # Create an effect structure
    smod$env = env
    class(smod) = "effect"
    
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
  # if (is.null(fvals$n) & 
  #     !(fvals$model == "linear") & 
  #     !(fvals$model == "offset") & 
  #     !(fvals$model == "factor") &
  #     !(fvals$model == "rw1") &
  #     !(fvals$model == "seasonal") & # Obtains n via season.length parameter
  #     !is.copy) 
  #   {
  #   stop(sprintf("Please provide parameter 'n' for effect '%s'", label))
  #   }
  # 
  
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
    eff = effects(object,label)
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
effects = function(model, label = NULL) {
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


list.A.model = function(model, points){ lapply(model$effects, amatrix.effect, points) }
list.indices.model = function(model, points){ lapply(model$effects, index.effect, points) }



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
    effs = effects(model)
    model$effects = effs[intersect(names(effs), all.vars(predictor))]
    vars = intersect(names(smp[[1]]), all.vars(predictor))
  }
  
  # Pre-calculate projection matrices
  As = lapply(model$effects, amatrix, points)

  for ( k in 1:n ) {
    # Discard variables we do not need
    sm = smp[[k]][vars]
    
    # Evaluate effects. Note that the expression may only contain hyper parameters in which case there 
    # are no effects to evaluate.
    enm = intersect(names(sm), names(model$effects))
    
    for (label in enm) {
      sm[[label]] = value(model$effects[[label]], data = points, state = sm[[label]], A = As[[label]]) 
      }
    
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

