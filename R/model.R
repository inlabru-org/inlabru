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

make.model = function(components){
  
  # Automatically add Intercept and -1 to components unless -1 is in components formula
  components = auto.intercept(components)
  
  # Back up environment
  env = environment(components)
  
  # Create effects
  effects = component(components)

  # Create joint formula that will be used by inla
  formula = y ~ -1 
  for (fm in lapply(effects, function(eff) {eff$inla.formula})) {
    formula = update.formula(formula, fm)
  }
  
  # Restore environment
  environment(formula) = env
  
  
  # Make model
  mdl = list(effects = effects, formula = formula)
  class(mdl) = c("model","list")
  return(mdl)
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
    effs = model$effects
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
      if ( is.data.frame(sm[[label]])) { sm[[label]] = sm[[label]]$value } 
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

