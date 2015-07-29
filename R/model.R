#####################################
# GENERICS
#####################################

join = function(...){UseMethod("join")}
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
  formula = y ~ - 1
  environment(formula) = new.env()
  mesh = list()
  mesh.coords = list()
  covariates = list()
  
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
    
  }
  return(make.model(
    formula = formula,
    covariates = covariates,
    mesh = mesh,
    mesh.coords = mesh.coords,
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

make.model = function(formula = NULL, name = NULL, mesh = NULL, mesh.coords = list(), covariates = list(), ...){
  mdl = list(
    formula = formula,
    name = name,
    mesh = mesh,
    mesh.coords = mesh.coords,
    covariates = covariates,
    args = list(...)
  )
  class(mdl) = c("model","list")
  return(mdl)
}


#####################################
# MODELS: Point process
#####################################

#' A simple intercept model
#'
#' @aliases model.intercept
#' 

model.intercept = function(data) {
  formula = ~ . -1 + Intercept
  covariates = list( Intercept = function(x) { return(0*x[,1]+1) } )
  return(make.model(formula = formula,
                    covariates = covariates))
}


#' A SPDE model
#'
#' @aliases model.spde
#' 

model.spde = function(data, mesh = data$mesh) {
  spde.args = list(alpha=2,prior.variance.nominal=10,theta.prior.prec=0.01)
  spde.mdl = do.call(inla.spde2.matern,c(list(mesh=mesh),spde.args)) 
  formula = ~ . + f(spde, model=spde.mdl)
  return(make.model(name = "Spatial SPDE model",
                    formula = formula,  
                    mesh = list(spde = mesh), 
                    mesh.coords = list(spde = data$mesh.coords)))
}


#####################################
# MODELS: Detection functions
#####################################

#' Half-normal detection function model
#'
#' @aliases halfnormal
#' 

model.halfnormal = function(truncation, constrained = TRUE, ...){
  # Formula
  if (!constrained) { formula = ~ . + nhsd }
  else { formula = ~. + f(nhsd, model = 'clinear', range = c(0, Inf), hyper = list(beta = list(initial = 0, param = c(0,0.1))))}
  
  # Covariate
  covariates = list(nhsd = function(X) { return(-0.5*X$distance^2) })
  
  # Evaluator
  evaluate = function(mdl,x) { stop("Not implemented.") }
  log.evaluate = function(mdl,x) { stop("Not implemented.") }
  
  # Environment
  env = list(...)
  
  return(make.model(name = "Half-normal detection function",
                    formula = formula, 
                    covariates = covariates))
}


#' Log-concave detection function model
#'
#' @aliases model.logconcave
#' 

model.logconcave = function(truncation, segments = 5, clinear = TRUE, linbasis = TRUE, ...){
  
  # Formula
  nSeg = segments
  if (clinear) {
    params = ",model = 'clinear',hyper=list(theta=list(prior=loggamma)),range = c(-10,0)" # 
    if (linbasis){
      fml = paste0("~ . + f(linbasis",params,")+",paste0("f(basis_",1:nSeg,params,")",collapse="+"))
    } else {
      fml = paste0("~ . +",paste0("f(basis_",1:nSeg,params,")",collapse="+"))
    }
    
  } else {
    if (linbasis){
      fml = paste("~ . + linbasis",paste(paste(rep(" + basis_",nSeg),seq(1:nSeg),sep=""),collapse=''),sep="")  
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

#' Spatial SPDE model for groupsize
#'
#' @aliases model.grpsize
#' 

model.grpsize = function(data, mesh = data$mesh) {
  gs.mdl.args = list(alpha=2, prior.variance.nominal = 10, theta.prior.prec = 0.01)
  gs.mdl = do.call(inla.spde2.matern, c(list(mesh=mesh), gs.mdl.args)) 
  formula = ~ . + f(grps, model = gs.mdl)
  covariates = list()
  
  return(make.model(name = "Spatial group size model",
                    formula = formula, 
                    covariates = covariates,
                    mesh = list(grps = mesh),
                    mesh.coords = list(grps = data$mesh.coords)))
}

