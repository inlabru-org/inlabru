#' Create a stack from a model, points and prediction alues
#' 

make.stack = function(points,
                       model,
                       y,
                       E = 1,
                       expr = NULL,
                       tag = "LeStack"){
  
  
  # Projection matrices (A) and mesh index effects
  A = list.A.model(model, points)
  eff = list.covariates.model(model, points)
  idx = list.indices.model(model, points)
  
  # Observations y
  y.pp = eval.if.function(y, points)
  
  # Expectation parameter E
  if ( !(length(E) == nrow(as.data.frame(points))) ) { 
    e.pp = rep(E, nrow(as.data.frame(points))) 
  } else { e.pp = E } 
  
  # Reweight things if dealing with taylor approximated model
  if ( !is.null(expr) ) {
    rw = nlinla.reweight(A, model, points, expr)
    A = rw$A
    e.pp = e.pp * exp(rw$const)
    if ( any((e.pp > 0) & (e.pp < 0.00001)) ) { msg("Exposure E is smaller than 0.00001. This may lead to numerical problems and non-convergence. Consider setting scale = 10 or higher when calling lgcp().")}
  }
  
  # Sort effects and A by names
  effects = c(idx, eff)
  if ("from.points" %in% names(eff)) A = c(A, list(from.points = 1))
  effects = effects[names(A)]

  # Create and return stack
  stk = inla.stack(data = list(y.inla = y.pp, e = e.pp),
                     A = A,
                     tag = tag,
                     effects = effects)
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
