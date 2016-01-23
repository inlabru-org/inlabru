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
  cat("--- Effects --- \n")
  for (k in 1:length(mdl$effects)){
    cat(paste0("'", mdl$effect[k],"' from sub-model '", mdl$name[k], "'\n"))
  }
  cat("\n")
  cat("--- Formula --- \n")
  cat(as.character(mdl$formula)[c(2,1,3)])

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
  
  # Check for duplicated effect names
  dup = duplicated(lapply(models, function(x) {x$effects}))
  if (any(dup)) { stop(paste0("Some of your models have identical effect names: ", lapply(models, function(x) {x$effects})[dup])) }
  
  A = list()
  name = character()
  formula = y.inla ~ - 1
  environment(formula) = new.env()
  mesh = list()
  effects = character()
  inla.spde = list()
  mesh.coords = list()
  time.coords = list()
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
    name = c(name, mdl$name)
    formula = update.formula(formula,mdl$formula)  
    covariates = c(covariates,mdl$covariates)
    effects = c(effects, mdl$effects)
    mesh = c(mesh,mdl$mesh)
    inla.spde = c(inla.spde, mdl$inla.spde)
    mesh.coords = c(mesh.coords,mdl$mesh.coords)
    time.coords = c(time.coords,mdl$time.coords)
    A = c(A,mdl$A)
    environment(formula) = env.upd(environment(formula), environment(mdl$formula))
    eval = c(eval, mdl$eval)
  }
  return(make.model(
    name = name,
    formula = formula,
    covariates = covariates,
    effects = effects,
    mesh = mesh,
    inla.spde = inla.spde,
    mesh.coords = mesh.coords,
    time.coords = time.coords,
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
      if (is.null(mdl$mesh[[cov.name]])) {
        cov.fun = covariates[[cov.name]]
        fetched.covar[[cov.name]] = cov.fun(pts)
        # Check if we actually got values back
        if(length(fetched.covar[[cov.name]]) == 0 ) { 
          stop(paste0("Evaluating cvoariate '", cov.name, "' returned no values. Do your points have all the data columns required by the covariate function?"))
        }
      } else {}
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
      if (mdl$inla.spde[[k]]$n.group > 1) {
        group = as.matrix(points[,mdl$time.coords[[k]]])
      } else { group = NULL }
      A = inla.spde.make.A(mdl$mesh[[k]], loc = loc, group = group)
      # Covariates for models with A-matrix are realized in the follwoing way:
      if (name %in% names(mdl$covariates)) {
        w = mdl$covariates[[name]](points)
        if (is.data.frame(w)) { stop("Your covariate function returns a data.frame. This is not allowed, a numeric vector is required.") }
        A = as.matrix(A)*as.vector(w)
      }
      A.lst[[name]] = A
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
        idx.lst = c(idx.lst, list(inla.spde.make.index(name, n.spde = mdl$inla.spde[[k]]$n.spde, n.group = mdl$inla.spde[[k]]$n.group)))
      }
    }
  }
  return(idx.lst)
}


#' Create a model
#'
#' @aliases make.model
#' 

make.model = function(formula = NULL, name = NULL, effects = NULL, mesh = NULL, inla.spde = list(), mesh.coords = list(), time.coords = list(), covariates = list(), eval = list(), ...){
  mdl = list(
    formula = formula,
    name = name,
    effects = effects,
    mesh = mesh,
    inla.spde = inla.spde,
    mesh.coords = mesh.coords,
    time.coords = time.coords,
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

evaluate.model.old = function(model, inla.result, loc = NULL, property = "mode") {
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


evaluate.model = function(model, inla.result, loc, property = "mode", do.sum = TRUE, link = identity) {
  cov = do.call(cbind, list.covariates.model(model, loc))
  Amat = list.A.model(model, loc)
  
  posts = list()
  
  for (k in 1:length(model$effects)){
    name = model$effects[k]
    if (is.null(name)) { name =  names(model$mesh)[[k]] }
    # Either fixed effect or hyper parameter
    if (name %in% rownames(inla.result$summary.fixed)){
      # Fixed effect
      post = inla.result$summary.fixed[name, property] * cov[[name]]
    } else if (paste0("Beta for ",name) %in% rownames(inla.result$summary.hyperpar)) {
      # Hyper parameter
      post = inla.result$summary.hyperpar[paste0("Beta for ",name), property] * cov[[name]]
    } else {
      # SPDE model   
      post = inla.result$summary.random[[name]][,property]
      A = Amat[[name]]
      # A workaround, needed for make.A called with group=1
      if (length(post) == 2*dim(A)[2]) { post = post[1:dim(A)[2]]}
      post = as.vector(A%*%as.vector(post))
    }
    posts[[name]] = post
  }
  ret = do.call(cbind, posts)
  if ( do.sum ) { ret = apply(ret, MARGIN = 1, sum) }
  return(link(ret))
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

model.intercept = function(data, effects = "Intercept") {
  formula = as.formula(paste0("~ . -1 + ", effects))
  covariates = list()
  covariates[[effects]] = function(x) {
    v = rep(1, nrow(x))
    ret = data.frame(v)
    colnames(ret) = effects
    return(ret) 
  }
  
  return(make.model(name = "Basic Intercept", formula = formula, effects = effects, covariates = covariates))
}


#' Spatial SPDE model
#' 
#' Constructs a spatial SPDE model with formula
#' 
#'  ~ . + f(spde, model = spde.mdl, group = spde.group)
#'
#' @aliases model.spde
#' @param data Data set to read the names of the mesh coordinates from (mesh.coords)
#' @param mesh The mesh used to construct the SPDE. If not provided, the mesh is read form the data set
#' @param n.group Number of SPDE model groups (see \link{f}), e.g. for temporal models.
#' @param ... Arguments passed on to inla.spde2.matern. If none, the defaults are alpha = 2, prior.variance.nominal = 10, theta.prior.prec = 0.01
#'  

model.spde = function(data, mesh = data$mesh, n.group = 1, covariates = NULL, ...) {
  
  vargs = list(...)
  if ( length(vargs) == 0 ){ 
    spde.args = list(alpha = 2, prior.variance.nominal = 10, theta.prior.prec = 0.01) }
  else { 
    spde.args = vargs 
  }
  
  spde.mdl = do.call(inla.spde2.matern,c(list(mesh=mesh),spde.args)) 
  spde.mdl$n.group = n.group
  formula = ~ . + f(spde, model=spde.mdl, group = spde.group)
  
  return(make.model(name = "Spatio-temporal SPDE model",
                    formula = formula,
                    effects = "spde",
                    mesh = list(spde = mesh),
                    inla.spde = list(spde = spde.mdl),
                    mesh.coords = list(spde = data$mesh.coords),
                    time.coords = list(spde = data$time.coords)))
}

#####################################
# MODELS: Detection functions
#####################################


#' Create detection function model
#' 
#' A wrapper for convenient construction of different detection functions (see parameter \code{type}).
#' For further help use the help files of the wrapped functions, e.g. \code{model.halfnormal()}
#' 
#' @aliases model.detfun
#' @param type Character setting the type of detection function. Options: "halfnormal", "exponential", "logconcave".

model.detfun = function(type, ...) {
  if (type == "halfnormal") { return(model.halfnormal(...))}
  else if (type == "exponential") { return(model.exponential(...))}
  else if (type == "logconcave") { return(model.logconcave(...))}
  else { stop(paste0("Unknown detection function type: ", type)) }
}


#' Half-normal detection function model
#'
#' Construct a half-normal detection function model using the formula
#' 
#'  ~ . + f(nhsd, model = 'clinear', range = c(0, Inf), hyper = list(beta = list(initial = 0, param = c(0,0.1))))
#'
#' Hence, the effect is obtained from a contrained linear INLA model. The covariate "nhsd" for this model is as a 
#' function of the data.frame provided during the inference process:
#'  
#'  nhsd = function(X) { return(-0.5*X[,colname]^2) },
#'  
#' that is, nhsd is half of the negative squared distance of an observation/integration point.
#'
#' @aliases model.halfnormal
#' @param data \link{dsdata} data set. Used to determine distance truncation (if not provided)
#' @param truncation Distance to truncate at. Currently unused but passed on the the formula environment for later usage.
#' @param constrained If set to false a non-constrained linear effect is used for estimating the detection function. Handy for debugging. 
#' @param colname Effort data column to use for the covariate extraction. Default: "distance"
#' @param effect Character setting the name of the effect in the model. Default: "nhsd" ("negative half squared distance")
#' 

model.halfnormal = function(data = NULL, truncation = NULL, constrained = TRUE, colname = "distance", effect = "nhsd"){
  
  # Extract truncation limit from data (if data is given but limit is not)
  if ( !is.null (data) & is.null(truncation) ) { truncation = max(detdata(data)[, colname]) }
  
  # Formula
  if (!constrained) { formula = as.formula(paste0("~ . + ", effect)) }
  else { formula =  as.formula(paste0(" ~ . + f(", effect, ", model = 'clinear', range = c(0, Inf), hyper = list(beta = list(initial = 0, param = c(0,0.1))))" ))}
  
  # Covariate
  covariates = list()
  covariates[[effect]] = function(X) { return(-0.5*X[,colname]^2) }
  
  return(make.model(name = "Half-normal detection function",
                    formula = formula, 
                    effects = effect,
                    covariates = covariates))
}


#' Expoential detection function model
#'
#' Construct a exponential detection function model using the formula
#' 
#'  ~ . + f(df.exp, model = 'clinear', range = c(0, Inf), hyper = list(beta = list(initial = 0, param = c(0,0.1))))
#'
#' Hence, the effect is obtained from a contrained linear INLA model. The covariate "df.exp" for this model is as a 
#' function of the data.frame provided during the inference process:
#'  
#'  df.exp = function(X) { return(-X[,colname]) },
#'  
#'
#' @aliases model.exponential
#' @param colname Effort data column to use for the covariate extraction. Default: "distance"
#' @param truncation Distance to truncate at. Currently unused but passed on the the formula environment for later usage.
#' @param constrained If set to false a non-constrained linear effect is used for estimating the detection function. Handy for debugging. 
#' 

model.exponential = function(colname = "distance", truncation = NULL,  constrained = TRUE){
  # Formula
  if (!constrained) { formula = ~ . + df.exp }
  else { formula = ~. + f(df.exp, model = 'clinear', range = c(0, Inf), hyper = list(beta = list(initial = 0, param = c(0,0.1))))}
  
  # Covariate
  covariates = list(df.exp = function(X) { return(-X[,colname]) })
  
  return(make.model(name = "Exponential detection function",
                    effects = "df.exp",
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
#' @param colname Effort data column to use for the covariate extraction. Default: "distance"
#' @param truncation Maximum distance assumed to be observed..
#' @param segments Number of quadratic basis basis functions
#' @param clinear If set to false, non-constrained linear INLA effects are used for this model. 
#' This is handy for debugging but the detection function is not necessarily log-concave anymore.
#' @param linbasis Wether to use a linear basis function. 
#' If set to FALSE this implies that that the detection function has zero derivative at the origin. 

model.logconcave = function(colname = "distance", truncation = NULL, segments = 5, constrained = TRUE, linbasis = TRUE, quadbasis = TRUE, ...){
  
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
  if (linbasis) { covariates = c(covariates, list(linbasis = function(x) { return(x[,colname]) }) ) }
  for ( k in 1:nSeg ) {
    covariates[[paste0("basis_",k)]] = function(x) { 
      return( data.frame(dfun_logconcave.basis.value(x[,colname],nSeg,truncation))[,paste0("basis_",k)] )
    }
    environment(covariates[[paste0("basis_",k)]]) = new.env(parent=environment())
    environment(covariates[[paste0("basis_",k)]])$k = k
  }
  
  effects = paste0("basis_", 1:nSeg)
  if (linbasis) { effects = c("linbasis", effects)}
  
  return(make.model(name = "Log-concave detection function",
                    effects = effects,
                    formula = formula, 
                    covariates = covariates))
}

#' Expoential detection function model in two dimensions
#'
#'
#' @aliases model.exponential2d
#' @param colname Effort data column to use for the covariate extraction. Default: c("distance","lgrpsize")
#' @param truncation Distance and log group size at which to truncate. Default: c(6,6)
#' @param constrained If set to false a non-constrained linear effect is used for estimating the detection function. Handy for debugging. 
#' 

model.exponential2d = function(colname = c("distance","lgrpsize"), truncation = c(6,6), constrained = TRUE, ...){
  
  truncation = as.list(truncation)
  names(truncation) = colname
  
  # Formula
  if (!constrained) { formula = ~ . + gdet.beta0 + gdet.beta1 }
  else { 
    formula = ~. + f(gdet.beta0, model = 'clinear', range = c(0, Inf), hyper = list(beta = list(initial = 0, param = c(0,0.1)))) +
      f(gdet.beta1, model = 'clinear', range = c(0, Inf), hyper = list(beta = list(initial = 0, param = c(0,0.1))))
  }
  gdet.beta0.cov = function(X) { return( ( (1 - X[,colname[2]]/truncation[[colname[1]]])) * -X[,colname[1]] ) }
  gdet.beta1.cov = function(X) { return( ( (X$colname[2]) / truncation[[colname[1]]]) * -X[,colname[1]] ) }
  covariates = list(gdet.beta0 = gdet.beta0.cov, gdet.beta1 = gdet.beta1.cov)
  
  
  eval = function(inla.result, loc, property = "mode", mdl = NULL) {
    if (class(inla.result)[[1]] == "inla") {       
      if ( "gdet.beta0" %in% colnames(inla.result$summary.fixed) ) { gdet.beta0 = inla.result$summary.fixed["gdet.beta0", property] }
      else { gdet.beta0 = inla.result$summary.hyperpar["Beta for gdet.beta0", property] }
      if ( "gdet.beta1" %in% colnames(inla.result$summary.fixed) ) { gdet.beta1 = inla.result$summary.fixed["gdet.beta1", property] }
      else { gdet.beta1 = inla.result$summary.hyperpar["Beta for gdet.beta1", property] }
    }
    else if ( is.list(inla.result) ){ 
      gdet.beta0 = inla.result$gdet.beta0
      gdet.beta1 = inla.result$gdet.beta1
    }
    else if ( is.numeric(inla.result) ){ 
      gdet.beta0 = inla.result$gdet.beta0
      gdet.beta1 = inla.result$gdet.beta1
    }
    else { stop("Type of result parameter not supported") }
    val = data.frame(gdet.beta0 = gdet.beta0*gdet.beta0.cov(loc), gdet.beta1 = gdet.beta1*gdet.beta1.cov(loc))    
    return(val)
  }
  
  return(make.model(name = "2D exponential detection function",
                    formula = formula, 
                    covariates = covariates,
                    eval = eval,
                    args = list(truncation = truncation, constrained = constrained)))
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
  gs.mdl$n.group = 1
  formula = ~ . + f(grps, model = gs.mdl)
  covariates = list()
    
  return(make.model(name = "Spatial group size model",
                    formula = formula,
                    effects = "grps",
                    covariates = covariates,
                    mesh = list(grps = mesh),
                    inla.spde = list(grps = gs.mdl),
                    mesh.coords = list(grps = data$mesh.coords),
                    time.coords = list(spde = data$time.coords),
                    eval = eval))
}
