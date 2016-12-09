#' iDistance models for INLA
#' 
#' Facilitates the spatial modeling approaches using INLA.
#' A \code{model} has a formula that describes one or more of INLA's \code{f} objects.
#' It describes how given a data frame of locations or other coordinates translates into
#' evaluating the predictors of the formula. For manually setting up a \code{model} see the
#' constructor \link{make.model}. Useful operators on \code{model} objects are:
#' \itemize{
#'  \item{\link{make.model}: }{Create a model}
#'  \item{\link{summary}: }{Summarize a model}
#'  \item{\link{join}: }{Join multiple models}
#'  \item{\link{evaluate}: }{ Determine the linear predictor's value at some location }
#'  \item{\link{list.data}: }{ List environmental variables needed to call INLA}
#' }
#' 
#' 
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
list.data = function(...){UseMethod("list.data")}

#####################################
# OPERATORS ON MODELS
#####################################

#' Model summary
#'
#' @aliases summary.model
#' @export
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

#' Join two or more models by summing up their linear predictors
#'
#' @aliases join.model join
#' @export
#' 

join.model = function(...){
  models = list(...)
  models = models[!unlist(lapply(models, is.null))]
  
  # Check for duplicated effect names
  dup = duplicated(lapply(models, function(x) {x$effects}))
  if (any(dup)) { stop(paste0("Some of your models have identical effect names: ", lapply(models, function(x) {x$effects})[dup])) }
  
  A = list()
  A.msk = list()
  name = character()
  formula = y.inla ~ - 1
  environment(formula) = new.env()
  mesh = list()
  effects = character()
  inla.spde = list()
  mesh.coords = list()
  time.coords = list()
  covariates = list()
  const = list()
  eval = list()
  map = list()
  
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
    A.msk = c(A.msk, mdl$A.msk)
    const = c(const, mdl$const)
    environment(formula) = env.upd(environment(formula), environment(mdl$formula))
    eval = c(eval, mdl$eval)
    map = c(map, mdl$map)
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
    const = const,
    eval = eval,
    A.msk = A.msk,
    iterator = do.call(c, lapply(models, function(m){m$iterator})),
    map = map
  ))
  
}


#' List data needed to run INLA
#'
#' @aliases list.data.model
#' @export
#' 

list.data.model = function(model){
  # THIS IS A WORKAROUND. If a formula is constructed and contains the stack of 
  # a previous run INLA will confuse the data from the new stack with the old stack
  # due to the data = ... as.list(environment(jmdl$formula))) ... parameter.
  # INLA reports the following error: 
  # Error in data.frame(..., check.names = FALSE): arguments imply differing number of rows: 2379, 2739
  
  assign("stack", NULL, envir = environment(model$formula))
  
  # Formula environment as list
  elist = as.list(environment(model$formula))
  
  # Remove previous inla results. For some reason these slow down the next INLA call.
  elist = elist[unlist(lapply(elist, function(x) !inherits(x, "inla")))]
}

#' List of covariates effects needed to run INLA
#'
#' @aliases list.covariates.model
#' 

list.covariates.model = function(mdl, pts){
  if ( length(mdl$covariates)>0 ) {
    formula = mdl$formula
    covariates = mdl$covariates
    all.varnames = setdiff(c(all.vars(formula), all.vars(mdl$expr)), names(mdl$mesh)) # Remove mesh names. Effects with names equal to SPDE models will lead to strange errors.
    covar.data = list()
    
    # Points may be annotated with covariates
    
    for (vname in all.varnames){
      if (vname %in% names(data.frame(pts))) {
        covar.data[[vname]] = data.frame(pts)[,vname]
      }
    }
    covar.data = data.frame(do.call(cbind,covar.data))
    
    # Add additional covariates defined by the user
    
    fetched.covar = list()
    for (cov.name in names(covariates)[!(names(covariates) == "")]){
      if (is.null(mdl$mesh[[cov.name]]) && !is.character(covariates[[cov.name]])) {
        cov.fun = covariates[[cov.name]]
        fetched.covar[[cov.name]] = cov.fun(data.frame(pts))
        # Check if we actually got values back
        if(length(fetched.covar[[cov.name]]) == 0 ) { 
          stop(paste0("Evaluating cvoariate '", cov.name, "' returned no values. Do your points have all the data columns required by the covariate function?"))
        }
      } else {}
    }
    # fetched.covar = do.call(cbind,c(fetched.covar))
    
    
    ret = list()
    ret = c(ret, fetched.covar)
    if (nrow(covar.data) > 0 ) { ret = c(ret, list(from.points = covar.data)) }
    ret
  } else { 
    ret = list()
  }
}

#' List of A matrices needed to run INLA
#'
#' @aliases list.A.model
#' 

list.A.model = function(mdl, points){
  A.lst = list()
  
  mapper = function(map) {
    if (class(map) == "call") { loc = eval(map, data.frame(points)) } 
    else {
      fetcher = get0(as.character(map))
      if (is.function(fetcher)) { loc = fetcher(points) } 
      else { loc = as.data.frame(points)[,as.character(map)] }
    }
  }
  
  for ( name in names(mdl$mesh)) {
    
    # What to use as loc input to inla.spde.make.A ?
    #
    # 1) if mesh.coords are provided, use them to select columns from the points data frame
    # 2) if a map function is provided use this function to map points to locations
    #    a) Function provided as call object
    #    b) Function provided as name object
    #    c) Map function is not really a function but a clumn name
    # 3) Otherwise assume points is a SpatialPoints object hence the locations are given as coordinates(points)
    
    if (!is.null(mdl$mesh.coords[[name]]) && all(mdl$mesh.coords[[name]] %in% names(points))) {
      loc = as.data.frame(points)[,mdl$mesh.coords[[name]]]
    } else if ( name %in% names(mdl$map) ) {
      if (class(mdl$map[[name]]) == "call") {
        loc = eval(mdl$map[[name]], data.frame(points))
      } else {
        fetcher = get0(as.character(mdl$map[[name]]))
        if (is.function(fetcher)) { loc = fetcher(points) } else { loc = as.data.frame(points)[,as.character(mdl$map[[name]])] }
      }
    } else { 
      loc = coordinates(points)
    }
    
    # inla.spde.make.A requires matrix format as input
    loc = as.matrix(loc)
    
    if (!("n.group" %in% names(mdl$inla.spde[[name]]))) { ng = 1 } else { ng = mdl$inla.spde[[name]]$n.group }
    if (ng > 1) {
      group = as.matrix(points[,mdl$time.coords[[name]]])
    } else { group = NULL }
    A = inla.spde.make.A(mdl$mesh[[name]], loc = loc, group = group)
    # Mask columns of A
    if (!is.null(mdl$A.msk[[name]])) { A = A[, as.logical(mdl$A.msk[[name]]), drop=FALSE]}
    # Weights for models with A-matrix are realized in the follwoing way:
    if (name %in% names(mdl$weights)) {
      w = mdl$weights[[name]](points)
      if (is.data.frame(w)) { stop("Your covariate function returns a data.frame. This is not allowed, a numeric vector is required.") }
      A = as.matrix(A)*as.vector(w)
    }
    A.lst[[name]] = A
  }
  
  for ( name in setdiff(mdl$effects, names(mdl$mesh)) ) {
    if ( name %in% names(mdl$map) ) {
      idx = mapper(mdl$map[[name]])
      A = matrix(0, nrow = nrow(points), ncol=6)
      A[cbind(1:nrow(A), idx)] = 1
      A.lst[[name]] = A
    } else {
      A.lst[[name]] = Matrix::Diagonal(nrow(data.frame(points)))
    }
  }
  
  # if ( length(mdl$covariates) > 0 ) { A.lst = c(A.lst,1) }
  
  return(A.lst)
}


#' List of spde indexing effects needed to run INLA
#'
#' @aliases list.indices.model
#' 

list.indices.model = function(mdl, points, ...){
  
  mapper = function(map) {
    if (class(map) == "call") { loc = eval(map, data.frame(points)) } 
    else {
      fetcher = get0(as.character(map))
      if (is.function(fetcher)) { loc = fetcher(points) } 
      else { loc = as.data.frame(points)[,as.character(map)] }
    }
  }
  
  idx.lst = list()
  if ( length(mdl$mesh) > 0) {
    for (k in 1:length(mdl$mesh)) {
      name = names(mdl$mesh)[[k]]
      if ( "m" %in% names(mdl$mesh[[k]]) ) {
        idx.lst[[name]] = 1:mdl$mesh[[k]]$m # support inla.mesh.1d models
        # If a is masked, correct number of indices
        if (!is.null(mdl$A.msk[[name]])) { idx.lst[[name]] = 1:sum(mdl$A.msk[[name]])}
      } else {
        ng = mdl$inla.spde[[name]]$n.group
        if (is.null(ng)) { ng = 1 }
        idx.lst[[name]] = inla.spde.make.index(name, n.spde = mdl$inla.spde[[name]]$n.spde, n.group = ng)
      }
    }
  }
  
  # Non-mesh models with indices
  for ( name in setdiff(names(mdl$map), names(mdl$mesh)) ) {
    idx.lst[[name]] = 1:max(mapper(mdl$map[[name]]))
  }
  
  return(idx.lst)
}


#' Create a model
#'
#' @aliases make.model
#' @export
#' @param formula A formula describing the (INLA) effect
#' @param name A character array describing the model
#' @param effects a list of strings naming the effects used in the formula
#' @param mesh An inla.mesh object (needed for SPDE models)
#' @param inla.spde An inla.spde model, if needed
#' @param mesh.coords Names of the data coordinates used to map from data to locations on the mesh
#' @param time.coords Names of the data coordinates mapping to a temporal aspect of the model
#' @param covariates A named list of covariate functions used to map data to weights or groupings of the effects
#' @param const A function that maps data to values that are constant with respect to the inference process

make.model = function(formula = NULL, name = NULL, effects = NULL, mesh = NULL, inla.spde = list(), mesh.coords = list(), time.coords = list(), covariates = list(), eval = list(), const = list(), ...){
  
  # Some defaults and shortcuts
  if ( is.null(effects) ) { 
    if ( is.list(covariates) ) { effects = names(covariates) }
    else {
      stop("make.model(): Please give your effect a name (via the 'effect' parameter)") 
    }
  }
  
  if ( is.function(covariates) ) { 
    tmp = covariates
    covariates = list()
    covariates[[effects]] = tmp
  }
  
  mdl = c(list(
    formula = formula,
    name = name,
    effects = effects,
    mesh = mesh,
    inla.spde = inla.spde,
    mesh.coords = mesh.coords,
    time.coords = time.coords,
    covariates = covariates,
    const = const,
    eval = eval), list(...))
  
  class(mdl) = c("model","list")
  return(mdl)
}



#' Evaluate model at given locations
#' 
#' Compute an approximation to the linear predictor at given locations and gicen coordinates.
#'
#' @aliases evaluate.model evaluate
#' @export
#' @param model An iDistance \link{model}
#' @param inla.result The result of an \link{inla} run or a sample obtained from \link{inla.posterior.sample.structured}
#' @param loc Locations and covariates needed to evaluate the model. If \code{NULL}, SPDE models will be evaluated at the mesh coordinates.
#' @param property Property of the model compnents to obtain value from. Default: "mode". Other options are "mean", "0.025quant", "0.975quant" and "sd".

evaluate.model = function(model, inla.result, loc, property = "mode", do.sum = TRUE, link = identity, n = 1, predictor = model$expr, use.covariate = TRUE) {
  cov = as.data.frame(do.call(cbind, list.covariates.model(model, loc)))
  Amat = list.A.model(model, loc)
  if ( property == "sample") {
    if ( inherits(inla.result, "inla") ) {
      smp = inla.posterior.sample.structured(inla.result, n = n) 
    } else {
      smp = inla.result
      n = length(smp)
    }
  } 
  posts = list()
  
  for (k in 1:length(model$effects)){
    name = model$effects[k]
    if ( property == "sample") {
      if (name %in% names(model$mesh)) {
        # SPDE model
        A = Amat[[name]]
        post = lapply(smp, function(s) { as.vector(A%*%as.vector(s[[name]])) })
      } else {
        A = Amat[[name]]
        post = lapply(smp, function(s) {
          mult = s[[name]]
          if (length(mult) == 1) { as.vector(A %*% rep(mult, nrow(A))) } 
          else { rowSums(t(t(as.matrix(A)) * mult)) }
        })
      }
    } else {
      if (is.null(name)) { name =  names(model$mesh)[[k]] }
      # Either fixed effect or hyper parameter
      if (name %in% rownames(inla.result$summary.fixed)){
        # Fixed effect
        # Check if the fixed effect is actually a function call
        if ( grepl("[()]", name)) {
          post = inla.result$summary.fixed[name, property] * eval(parse(text = name), as.data.frame(loc))
        } else {
          # Note: if there is no covariate, assume covariate = 1
          if (is.null(cov[[name]]) | !use.covariate) {
            post = rep(inla.result$summary.fixed[name, property], nrow(data.frame(loc)))
          } else {
            post = inla.result$summary.fixed[name, property] * cov[[name]]
          }
        }
        
      } else if (paste0("Beta for ",name) %in% rownames(inla.result$summary.hyperpar)) {
        # Hyper parameter
        post = inla.result$summary.hyperpar[paste0("Beta for ",name), property] * cov[[name]]
      } else {
        post = inla.result$summary.random[[name]][,property]
        if ( name %in% names(Amat) ){ # SPDE model
          A = Amat[[name]]
          # A workaround, needed for make.A called with group=1
          if (length(post) == 2*dim(A)[2]) { post = post[1:dim(A)[2]]}
          post = as.vector(A%*%as.vector(post))
        } else { # other models
          post = post[model$covariates[[name]](loc)]
        }
        
      }
    }
    posts[[name]] = post
  }
  
  if ( property == "sample") {
    if (!is.null(predictor) && do.sum) { 
      ret = do.call(Map, c(list(function(...){apply(cbind(...),MARGIN=1,identity)}), posts))
      ret = lapply(ret, function(r) {eval(predictor, envir = cbind(as.data.frame(t(r)), as.data.frame(loc)))})
    } else {
      ret = do.call(Map, c(list(function(...){apply(cbind(...),MARGIN=1,sum)}), posts))
      if( "const" %in% names(model) & !(length(model$const)==0)) {
        const = colSums(do.call(rbind, lapply(model$const, function(f) { f(loc) })))
        ret = lapply(ret, function(x) { x + const })
      }
    }
    ret = lapply(ret, link)
  } else {
    ret = do.call(cbind, posts)
    if (!is.null(predictor) && do.sum) { ret = eval(predictor, c(posts, as.list(data.frame(loc)))) } 
    else if ( do.sum ) { ret = apply(ret, MARGIN = 1, sum) }
    
    if( "const" %in% names(model) & !(length(model$const)==0)) { 
      ret = ret + colSums(do.call(rbind, lapply(model$const, function(f) { f(loc) })))
    }
    ret = link(ret)
  }
  
  return(ret)
}
