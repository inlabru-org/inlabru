#' Linear predictor models for usage with INLA
#' 
#' This class facilitates the usage of sophisticated modeling approaches using INLA. 
#' A \code{model} has a formula that describes one or more of INLA's \code{f} objects. 
#' It describes how given a data frame of locations or other coordinates translates into
#' evaluating the predictors of the formula. For manually setting up a \code{model} see the
#' constructor \link{make.model}. Useful operators on \code{model} objects are:
#' \itemize{
#'  \item{\link{summary}: }{Summarize a model}
#'  \item{\link{join}: }{Join multiple models}
#'  \item{\link{list.covariates}: }{List covariates of a model}
#'  \item{\link{list.A}: }{List projection matrices of a model}
#'  \item{\link{list.indices}: }{List indices of a model}
#' }
#' Default models are:
#' \itemize{
#' \item Intercept: \link{model.intercept}
#' \item Half-normal detection function: \link{model.halfnormal}
#' \item Log-concave detection function: \link{model.logconcave}
#' \item Spatial SPDE model for animal intensity: \link{model.spde}
#' \item Spatial group size model: \link{model.grpsize}
#' }     
#' @name model
NULL


#####################################
# GENERICS
#####################################

join = function(...){UseMethod("join")}
evaluate = function(...){UseMethod("evaluate")}
list.covariates = function(...){UseMethod("list.covariates")}
list.A = function(...){UseMethod("list.A")}
list.indices = function(...){UseMethod("list.indices")}

#####################################
# OPERATORS ON MODELS
#####################################

#' Model summary
#'
#' @aliases summary.model
#' 

summary.model = function(mdl) {
  cat("--- Model summary --- \n")
  cat(paste0("Model name: ",mdl$name, "\n"))
  cat(paste0("Formula:" , mdl$formula, "\n"))
  cat(paste0("Covariates: ", names(mdl$covariates), "\n"))
}

#' Is it a model?
#'
#' @aliases is.model
#' 

is.model = function(mdl) { return(class(mdl)[[1]] == "model")}


#' Join two or more models by summing up their linear predictors
#'
#' @aliases join.model
#' 

join.model = function(...){
  models = list(...)
  A = list()
  formula = y.inla ~ - 1
  environment(formula) = new.env()
  mesh = list()
  mesh.coords = list()
  covariates = list()
  eval = list()
  
  env.upd = function(env1,env2) 
    {
    nms = names(env2)[which(!(names(env2) %in% names(env1)))]
    lapply(nms,function(nm) {env1[[nm]] = env2[[nm]]})
    return(env1)
    }
  
  for (k in 1:length(models)){
    mdl = models[[k]]
    formula = update.formula(formula,mdl$formula)  
    covariates = c(covariates,mdl$covariates)
    mesh = c(mesh,mdl$mesh)
    mesh.coords = c(mesh.coords,mdl$mesh.coords)
    A = c(A,mdl$A)
    environment(formula) = env.upd(environment(formula), environment(mdl$formula))
    eval = c(eval, mdl$eval)
  }
  return(make.model(
    formula = formula,
    covariates = covariates,
    mesh = mesh,
    mesh.coords = mesh.coords,
    eval = eval
    ))
  
}

#' List of covariates effects needed to run INLA
#'
#' @aliases list.covariates.model
#' 

list.covariates.model = function(mdl, pts){
  if ( length(mdl$covariates)>0 ) {
    formula = mdl$formula
    covariates = mdl$covariates
    all.varnames = all.vars(formula)
    covar.data = list()
    
    # Points may be annotated with covariates
    
    for (vname in all.varnames){
      if (vname %in% names(pts)) {
        covar.data[[vname]] = pts[,vname]
      }
    }
    covar.data = data.frame(do.call(cbind,covar.data))
    
    # Add additional covariates defined by the user
    
    fetched.covar = list()
    for (cov.name in names(covariates)){
      cov.fun = covariates[[cov.name]]
      fetched.covar[[cov.name]] = cov.fun(pts)
    }
    fetched.covar = do.call(cbind,c(fetched.covar))
    
    if (length(fetched.covar)>0 & length(covar.data)>0){ 
      covar.data = cbind(covar.data,fetched.covar)
    } else { 
      covar.data = fetched.covar 
    }
    
    if (is.null(covar.data)) {
      return(list())
    } else { 
      return(list(as.data.frame(covar.data)))
    }
  } else { 
    E = list() 
  }
  return(E)
}

#' List of A matrices needed to run INLA
#'
#' @aliases list.A.model
#' 

list.A.model = function(mdl, points){
  A.lst = list()
  if ( length(mdl$mesh)>0 ) {
    for (k in 1:length(mdl$mesh)) {
      name = names(mdl$mesh)[[k]]
      loc = as.matrix(points[,mdl$mesh.coords[[k]]])
      A.lst[[name]] = inla.spde.make.A(mdl$mesh[[k]], loc)
    }
  }  
  if ( length(mdl$covariates) > 0 ) { A.lst = c(A.lst,1) }
  return(A.lst)
}


#' List of spde indexing effects needed to run INLA
#'
#' @aliases list.indices.model
#' 

list.indices.model = function(mdl, ...){
  idx.lst = list()
  if ( length(mdl$mesh) > 0) {
    for (k in 1:length(mdl$mesh)) {
      name = names(mdl$mesh)[[k]]
      if ( "m" %in% names(mdl$mesh[[k]]) ) {
        idx.lst[[name]] = 1:mdl$mesh[[k]]$m # support inla.mesh.1d models  
      } else {
        idx.lst[[name]] = 1:mdl$mesh[[k]]$n
      }
    }
  }
  return(idx.lst)
}


#' Create a model
#'
#' @aliases make.model
#' 

make.model = function(formula = NULL, name = NULL, mesh = NULL, mesh.coords = list(), covariates = list(), eval = list(), ...){
  mdl = list(
    formula = formula,
    name = name,
    mesh = mesh,
    mesh.coords = mesh.coords,
    covariates = covariates,
    eval = eval,
    args = list(...)
  )
  class(mdl) = c("model","list")
  return(mdl)
}

#' Evaluate model at given locations
#' 
#' Compute an approximation to the linear predictor at given locations and gicen coordinates.
#'
#' @aliases evaluate.model
#' @param model An iDistance \link{model}
#' @param inla.result The result of an \link{inla} run or a sample obtained from \link{inla.posterior.sample.structured}
#' @param loc Locations and covariates needed to evaluate the model. If \code{NULL}, SPDE models will be evaluated at the mesh coordinates.
#' @param property Property of the model compnents to obtain value from. Default: "mode". Other options are "mean", "0.025quant", "0.975quant" and "sd".

evaluate.model = function(model, inla.result, loc = NULL, property = "mode") {
  if (is.list(model$eval)){
    vals = list()
    for (k in 1:length(model$eval)){
      vals[[k]] = model$eval[[k]](inla.result, loc, property = property)
    }
    return(do.call(cbind,vals))
  }
  else {
    return(model$eval(inla.result, loc, property = property))
  }
}


#' Sample from a model
#'
#' @aliases sample.model
#' 

sample.model = function(model, inla.result, n = 1, loc = NULL) {
  samples = inla.posterior.sample.structured(inla.result, n)
  samples = lapply(samples, function(x) { evaluate.model(model, inla.result = x, loc) } )
  return(samples)
}


#####################################
# MODELS: Point process
#####################################

#' Intercept model
#'
#' This model represents a simpel intercept added to the formula, i.e.:
#' 
#' ~ . + Intercept
#'
#' where the effects are 1 for each observation.
#'
#' @aliases model.intercept
#' 

model.intercept = function(data) {
  formula = ~ . -1 + Intercept
  covariates = list( Intercept = function(x) { return(data.frame(Intercept = rep(1, nrow(x)) )) } )
  
  eval = function(inla.result, loc, property = "mode") {
    if (class(inla.result)[[1]] == "inla") { w = inla.result$summary.fixed["Intercept", property]  }
    else if ( is.list(inla.result) ){ w = inla.result$Intercept }
    else if ( is.numeric(inla.result) ){ w = inla.result }
    else { stop("Type of result parameter not supported") }
    val = data.frame(Intercept = rep(w, dim(loc)[1]))
    
    return(val)
  }

  return(make.model(formula = formula,
                    covariates = covariates,
                    eval = eval))
}


#' Spatial SPDE model
#' 
#' Constructs a spatial SPDE model with formula
#' 
#'  ~ . + f(spde, model = spde.mdl)
#'
#' @aliases model.spde
#' @param data Data set to read the names of the mesh coordinates from (mesh.coords)
#' @param mesh The mesh used to construct the SPDE. If not provided, the mesh is read form the data set
#' @param ... Arguments passed on to inla.spde2.matern. If none, the defaults are alpha = 2, prior.variance.nominal = 10, theta.prior.prec = 0.01
#'  

model.spde = function(data, mesh = data$mesh, ...) {
  
  vargs = list(...)
  if ( length(vargs) == 0 ){ 
    spde.args = list(alpha = 2, prior.variance.nominal = 10, theta.prior.prec = 0.01) }
  else { 
    spde.args = vargs 
  }
  
  spde.mdl = do.call(inla.spde2.matern,c(list(mesh=mesh),spde.args)) 
  formula = ~ . + f(spde, model=spde.mdl)
  
  eval = function(inla.result, loc, property = "mode") {
    if (class(inla.result)[[1]] == "inla"){ weights = inla.result$summary.random$spde[,property] }
    else if ( is.list(inla.result) ){ weights = inla.result$spde }
    else if ( is.numeric(inla.result) ){ weights = inla.result }
    else { stop("Type of field parameter not supported") }
    A = inla.spde.make.A(mesh, loc = as.matrix(loc[, data$mesh.coords]))
    val = data.frame(spde = as.vector(A%*%as.vector(weights)))
    
    return(val)
  }
  
  return(make.model(name = "Spatial SPDE model",
                    formula = formula,  
                    mesh = list(spde = mesh), 
                    mesh.coords = list(spde = data$mesh.coords),
                    eval = eval))
}


#####################################
# MODELS: Detection functions
#####################################

#' Half-normal detection function model
#'
#' Construct a half-normal detection function model using the formula
#' 
#'  ~ . + f(nhsd, model = 'clinear', range = c(0, Inf), hyper = list(beta = list(initial = 0, param = c(0,0.1))))
#'
#' Hence, the effect is obtained from a contrained linear INLA model. The covariate "nhsd" for this model is as a 
#' function of the data.frame provided during the inference process:
#'  
#'  nhsd = function(X) { return(-0.5*X$distance^2) },
#'  
#' that is, nhsd is half of the negative squared distance of an observation/integration point.
#'
#' @aliases model.halfnormal
#' @param truncation Distance to truncate at. Currently unused but passed on the the formula environment for later usage.
#' @param constrained If set to false a non-constrained linear effect is used for estimating the detection function. Handy for debugging. 
#' 

model.halfnormal = function(truncation, constrained = TRUE, ...){
  # Formula
  if (!constrained) { formula = ~ . + nhsd }
  else { formula = ~. + f(nhsd, model = 'clinear', range = c(0, Inf), hyper = list(beta = list(initial = 0, param = c(0,0.1))))}
  
  # Covariate
  covariates = list(nhsd = function(X) { return(-0.5*X$distance^2) })
  
  return(make.model(name = "Half-normal detection function",
                    formula = formula, 
                    covariates = covariates))
}


#' Log-concave detection function model
#' 
#' Constructs a log-concave detection function model with formula
#' 
#' ~ . + f(linbasis, model = "clinear", ...) 
#'      + f(basis_1, model = "clinear", ...) 
#'      + ... + f(basis_K, model = "clinear", ...)
#'
#'  whare K is the number of segments the basis lives on.
#'   
#' @aliases model.logconcave
#' @param truncation Maximum distance assumed to be observed..
#' @param segments Number of quadratic basis basis functions
#' @param clinear If set to false, non-constrained linear INLA effects are used for this mode. 
#' This is handy for debugging but the detection function is not necessarily log-concave anymore.
#' @param linbasis Wether to use a linear basis function. 
#' If set to FALSE this implies that that the detection function has zero derivative at the origin. 

model.logconcave = function(truncation, segments = 5, constrained = TRUE, linbasis = TRUE, quadbasis = TRUE, ...){
  
  # Formula
  nSeg = segments
  if (constrained) {
    params = ",model = 'clinear',hyper=list(theta=list(prior=loggamma)),range = c(-10,0)" # 
    if (linbasis){
      if (quadbasis) {
        fml = paste0("~ . + f(linbasis",params,")+",paste0("f(basis_",1:nSeg,params,")",collapse="+"))
      } else {
        fml = paste0("~ . + f(linbasis", params, ")")
      }
    } else {
      fml = paste0("~ . +",paste0("f(basis_",1:nSeg,params,")",collapse="+"))
    }
    
  } else {
    if (linbasis){
        if (quadbasis) {
          fml = paste("~ . + linbasis",paste(paste(rep(" + basis_",nSeg),seq(1:nSeg),sep=""),collapse=''),sep="")  
        } else {
          fml = "~ . + linbasis"
        }
    } else {
      fml = paste("~ . + ",paste(paste(rep(" + basis_",nSeg),seq(1:nSeg),sep=""),collapse=''),sep="")  
    }
  }
  
  loggamma = "expression:
  a = 1;
  b = 0.001;
  precision = exp(log_precision);
  logdens = log(b^a) - lgamma(a)
  + (a-1)*log_precision - b*precision;
  log_jacobian = log_precision;
  return(logdens + log_jacobian);"
  hyper.new = list(theta = list(prior = loggamma))
  
  formula = as.formula(fml)
  
  # Covariate
  covariates = list()
  if (linbasis) { covariates = c(covariates, list(linbasis = function(x) { return(x$distance) }) ) }
  for ( k in 1:nSeg ) {
    covariates[[paste0("basis_",k)]] = function(x) { 
      return( data.frame(dfun_logconcave.basis.value(x$distance,nSeg,truncation))[,paste0("basis_",k)] )
    }
    environment(covariates[[paste0("basis_",k)]]) = new.env(parent=environment())
    environment(covariates[[paste0("basis_",k)]])$k = k
  }
  
  return(make.model(name = "Log-concave detection function",
                    formula = formula, 
                    covariates = covariates))
}

#####################################
# MODELS: Group size
#####################################

#' Spatial SPDE model for group size with formula
#' 
#'  ~ . + f(grps, model = gs.mdl)
#'
#' @aliases model.grpsize
#' @param data Data set to read the names of the mesh coordinates from (mesh.coords)
#' @param mesh The mesh used to construct the SPDE. If not provided, the mesh is read form the data set
#' @param ... Arguments passed on to inla.spde2.matern. If none, the defaults are alpha = 2, prior.variance.nominal = 10, theta.prior.prec = 0.01

model.grpsize = function(data, mesh = data$mesh, ...) {
  
  vargs = list(...)
  if ( length(vargs) == 0 ){ 
    gs.mdl.args = list(alpha = 2, prior.variance.nominal = 10, theta.prior.prec = 0.01) }
  else { 
    gs.mdl.args = vargs 
  }
  
  gs.mdl = do.call(inla.spde2.matern, c(list(mesh=mesh), gs.mdl.args)) 
  formula = ~ . + f(grps, model = gs.mdl)
  covariates = list()
  
  return(make.model(name = "Spatial group size model",
                    formula = formula, 
                    covariates = covariates,
                    mesh = list(grps = mesh),
                    mesh.coords = list(grps = data$mesh.coords)))
}

