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
#'  \item{\link{sample.value}: }{ Sample the linear predictor value for geiven locations}
#'  \item{\link{sample.points}: }{ Sample spatial points from a log gaussion Cox process model}
#'  \item{\link{update}: }{ Update the iterator of an axproximative model}
#' }
#' 
#' After running inference plotting aspects of models is straight forward using:
#' 
#' \itemize{
#'  \item{\link{plot.detfun}: }{Plot detection functions (\link{model.detfun})}
#'  \item{\link{plot.spatial}: }{Plot a spatial models, e.g. \link{model.spde}}
#'  \item{\link{plot.marginal}: }{Plot marginal prosterior of any model effect}
#' }
#' 
#' Some useful default models are:
#' \itemize{
#' \item Intercept: \link{model.intercept}
#' \item Fixed effects models: \link{model.fixed}
#' \item Detection function models: \link{model.detfun}
#' \item Spatial SPDE model \link{model.spde}
#' \item Spatial group size model: \link{model.grpsize}
#' }     
#' @name model
NULL


#####################################
# GENERICS
#####################################

join = function(...){UseMethod("join")}
evaluate = function(...){UseMethod("evaluate")}
iteration.history  = function(...){UseMethod("iteration.history")}
sample.value = function(...){UseMethod("sample.value")}
sample.points = function(...){UseMethod("sample.points")}
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

#' Is it a model?
#'
#' @aliases is.model
#' 

is.model = function(mdl) { return(class(mdl)[[1]] == "model")}


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
    fetched.covar = do.call(cbind,c(fetched.covar))
    
    if (length(fetched.covar)>0 & length(covar.data)>0){ 
      covar.data = cbind(covar.data,fetched.covar)
    } else if (length(fetched.covar)>0) { 
      covar.data = fetched.covar 
    } else { }
    
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
        # If a is masked, correct number of indices
        if (!is.null(mdl$A.msk[[name]])) { idx.lst[[name]] = 1:sum(mdl$A.msk[[name]])}
      } else {
        ng = mdl$inla.spde[[k]]$n.group
        if (is.null(ng)) { ng = 1 }
        idx.lst = c(idx.lst, list(inla.spde.make.index(name, n.spde = mdl$inla.spde[[k]]$n.spde, n.group = ng)))
      }
    }
    
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


#' Update a model's iterators using an INLA result
#'
#' @aliases update.model update
#' @export 
#' @name update.model
#' @param model A model
#' @param result An inla result
#' 

update.model = function(model, result = NULL){
  for (name in names(model$iterator)) {
    iter = model$iterator[[name]]
    hname = paste0(name,".history")
    if (!is.null(result)) {
      if ( name %in% rownames(result$summary.fixed) ) { 
        assign(name, result$summary.fixed[name,"mode"], envir = iter)
        assign(hname, c(iter[[hname]], result$summary.fixed[name,"mode"]), envir = iter)
        }
      if ( name %in% names(result$summary.random) ) { 
        assign(name, result$summary.random[[name]][,"mode"], envir = iter)
        assign(hname, rbind(iter[[hname]], result$summary.random[[name]][,"mode"]), envir = iter)
        }
      if ( name %in% rownames(result$summary.hyperpar) ) { 
        assign(name, result$summary.hyperpar[name,"mode"], envir = iter)
        assign(hname, c(iter[[hname]], result$summary.hyperpar[name,"mode"]), envir = iter)
        }
      if ( paste0("Beta for ",name) %in% rownames(result$summary.hyperpar) ) {
        assign(name, result$summary.hyperpar[paste0("Beta for ",name),"mode"], envir = iter)
        assign(hname, c(iter[[hname]], result$summary.hyperpar[paste0("Beta for ",name),"mode"]), envir = iter)
      }
    } else {
      iter[[hname]] = iter[[name]]
    }
  }
}

#' Extract the history of the model's iterators
#'
#' @aliases iteration.history
#' @export
#' @param model a \link{model}
#' 
iteration.history.model = function(model, effect){
  return(model$iterator[[effect]][[paste0(effect,".history")]])
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
  cov = do.call(cbind, list.covariates.model(model, loc))
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
        post = lapply(smp, function(s) {
          # Note: if there is no covariate, assume covariate = 1
          if (is.null(cov[[name]]) | !use.covariate) {
            rep(s[[name]],nrow(loc))
            } else { 
              s[[name]] * cov[[name]] } 
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



#' Sample from a model
#'
#' @aliases sample.value.model sample.value
#' @export

sample.value.model = function(model, result = NULL, n = 1, data = NULL, loc = NULL) {
  
  if ( is.null(loc) ) {
    loc = data$mesh$loc[,1:2]; 
    colnames(loc) = data$mesh.coords    
  }
  
  samples = inla.posterior.sample.structured(result, n)
  cov = do.call(cbind, list.covariates.model(model, loc))
  Amat = list.A.model(model, loc)
  ret.samples = list()
  for (s in 1:n) {
    posts = list() 
    for (k in 1:length(model$effects)){
      name = model$effects[k]
      post = samples[[s]][[name]]
      if (name %in% names(Amat)) {
        A = Amat[[name]]
        post = as.vector(A%*%as.vector(post))
      }
      posts[[name]] = post
    }
    ret = do.call(cbind, posts)
    ret = apply(ret, MARGIN = 1, sum)
    ret.samples[[s]] = ret
  }
  ret = do.call(cbind, ret.samples)
  return(ret)
}


#' Sample points from a model
#'
#' @aliases sample.points.model sample.points
#' @export
 
sample.points.model = function(model, result = NULL, data = NULL,  property = "random"){
  loc = data$mesh$loc[,1:2]
  colnames(loc) = data$mesh.coords
  if ( property == "random") {
    weights = sample.value.model(model, result = result, n = 1, loc = loc)
  } else {
    weights = evaluate(model = model, inla.result = result, loc = loc, do.sum = TRUE, property = property)
  }
  pts = sample.lgcp(data$mesh, weights = weights, geometry = data$geometry)
  colnames(pts) = data$mesh.coords
  return(data.frame(pts))
}


#' Estimate animal abundance
#' 
#'
#' @aliases abundance.model abundance
#' @export
#' @param model An iDistance \link{model}
#' @param result The result of an \link{inla} run or a sample obtained from \link{inla.posterior.sample.structured}
#' @param n Number of samples used to calculate abundance

abundance.model = function(model, data, result, n = 1, mask = NULL, property = "sample"){
  if (is.null(mask)) {
    loc = data.frame(data$mesh$loc)
    colnames(loc) = data$mesh.coords
    weights = diag(as.matrix(inla.mesh.fem(data$mesh)$c0))
  } else if ( class(mask) == "inla.mesh") {
    loc = data.frame(mask$loc)
    colnames(loc) = data$mesh.coords
    weights = diag(as.matrix(inla.mesh.fem(mask)$c0))
  }
  
  rates = numeric()
  
  for (k in 1:n) {
    rates[k] = sum( evaluate(model=model, inla.result = result, loc = loc, property = property, link = exp) * weights)
  }
  x = floor(range(rates)[1]):ceiling(range(rates)[2])
  dens = matrix(ncol = n, nrow = length(x))
  for (k in 1:n) {
    dens[,k] = dpois(x, rates[k])
  }
  dens = apply(dens, MARGIN=1, sum)
  return(data.frame(n = x, probability = dens ))
}



#####################################
# MODELS: Point process
#####################################

#' Generic INLA f() wrapper
#'
#'
#'
#' @aliases model.inla
#' 

model.f = function(covariates, ...) {
  effects = names(covariates)[[1]]
  formula = as.formula(paste0("~ . + f(",effects, ", ...)"))
  return(make.model(name = "INLA f", formula = formula, effects = effects, covariates = covariates))
}


#' Intercept model
#'
#' This model represents a simpel intercept added to the formula, i.e.:
#' 
#' ~ . + Intercept
#'
#' where the effects are 1 for each observation.
#'
#' @aliases model.intercept
#' @export
#' 

model.intercept = function(data, effects = "Intercept") {
  formula = as.formula(paste0("~ . -1 + ", effects))
  covariates = list()
  
  covariates[[effects]] = function(x) {
    v = rep(1, nrow(as.data.frame(x)))
    ret = data.frame(v)
    colnames(ret) = effects
    return(ret) 
  }
  
  return(make.model(name = "Basic Intercept", formula = formula, effects = effects, covariates = covariates))
}



#' Fixed effect model
#'
#' A simple fixed effect of an equally named covariate
#' 
#' ~ . + fixed
#'
#'
#' @aliases model.fixed
#' @export
#' 

model.fixed = function(effects = NULL, covariates = list()) {
  if ( is.null(effects) ) { effects = names(covariates) }
  formula = as.formula(paste0("~ . ", do.call(paste0,as.list(paste(" + ", effects, "",sep = "")))))
  for ( k in 1:length(effects) ) {
    effect = effects[k]
    if (!(effect %in% names(covariates))) { covariates[[effect]] = function(x) { x[,effect] } }
    environment(covariates[[effect]]) = new.env()
    assign("effect", effect, envir = environment(covariates[[effect]]))
  }
  return(make.model(name = "Fixed effect", formula = formula, effects = effects, covariates = covariates))
}

#' Spatial SPDE model
#' 
#' Constructs a spatial SPDE model
#' 
#' This model uses the formula
#' 
#'  ~ . + f(spde, model = spde.mdl, group = spde.group, ...)
#'  
#'  where the name 'spde' can be replaced by the \code{effect} parameter. The internal INLA spde model
#'  spde.mdl is constructed using the dotted arguments. If these do not contain a mesh, the mesh of the
#'  data set \code{data} is used. If there are no ... arguments, it sets alpha = 2, 
#'  prior.variance.nominal = 10 and theta.prior.prec = 0.01. 
#'
#' @aliases model.spde
#' @export
#' @param data Data set to read the names of the mesh coordinates from (mesh.coords)
#' @param effect A string defining the name of the SPDE effect
#' @param n.group Number of SPDE model groups (see \link{f}), e.g. for temporal models.
#' @param group A string defining the data column to group by, e.g. "season". Groups must be integer and > 0 
#' @param covariate A function that maps a data frame of points to a scaling factor for the SPDE at that point.
#' @param ... Arguments passed on to inla.spde2.matern. If none, the defaults are alpha = 2, prior.variance.nominal = 10, theta.prior.prec = 0.01
#'  

model.spde = function(data, effect = "spde", n.group = 1, group = data$time.coords, covariate = NULL, ...) {
  
  # Set default parameters for SPDE (if none provided)
  vargs = list(...)
  if ( length(vargs) == 0 ){ spde.args = list(alpha = 2, prior.variance.nominal = 10, theta.prior.prec = 0.01) }
  else { spde.args = vargs }
  
  # If no mesh is provided, the SPDE lives on the data mesh 
  if (!("mesh" %in% names(spde.args))) { spde.args[["mesh"]] = data$mesh }
  
  # Create the SPDE model
  tmp.mdl = do.call(inla.spde2.matern, spde.args)
  tmp.mdl$n.group = n.group
  
  # Make the formula
  formula = as.formula(paste0("~ . + f(",effect,", model = ",effect,".mdl, group = ",effect,".group)"))
  assign(paste0(effect,".mdl"), tmp.mdl, envir = environment(formula))
  
  # Little helper
  lassign = function(name,value) { tmp = list(); tmp[[name]] = value ; return(tmp) }
  
  return(make.model(name = "Spatio-temporal SPDE model",
                    formula = formula,
                    effects = effect,
                    mesh = lassign(effect, spde.args$mesh),
                    covariates = lassign(effect, covariate),
                    inla.spde = lassign(effect, tmp.mdl),
                    mesh.coords = lassign(effect, data$mesh.coords),
                    time.coords = lassign(effect, group)))
}

#####################################
# MODELS: Detection functions
#####################################


#' Create detection function model
#' 
#' A wrapper for convenient construction of different detection functions (see parameter \code{type}).
#' For further help use the help files of the wrapped functions,
#' \itemize{
#'  \item{\link{model.halfnormal}: }{Half-normal model}
#'  \item{\link{model.exponential}: }{Exponential model}
#'  \item{\link{model.logconcave}: }{Log-concave model}
#'  \item{\link{model.hazard}: }{Hazard rate model}
#' }
#' 
#' 
#' @aliases model.detfun
#' @name model.detfun
#' @export
#' @param type Character setting the type of detection function. Options: "halfnormal", "exponential", "logconcave".
#' @param ... Parameters that are passed on to the model.X function
#' @return A detection function \link{model}


model.detfun = function(type, ...) {
  if (type == "halfnormal") { return(model.halfnormal(...))}
  else if (type == "hazard") { return(model.hazard(...))}
  else if (type == "exponential") { return(model.exponential(...))}
  else if (type == "logconcave") { return(model.logconcave(...))}
  else if (type == "smooth") { return(model.smoothdf(...))}
  else { stop(paste0("Unknown detection function type: ", type)) }
}


#' Half-normal detection function
#'
#' Construct a half-normal detection function \link{model}
#
#' The log-likelihood of this model is
#' \deqn{log p(z) = -0.5 * q * z^2}
#' 
#' with precision parameter q and distance z. Internally, it is represented by the \link{formula}
#'  ~ . + f(nhsd, model = 'clinear', range = c(0, Inf), hyper = list(beta = list(initial = 0, param = c(0,0.1))))
#' 
#' Hence, the precision effect 'nhsd' is obtained from a contrained linear INLA model. The covariate used for this
#' model is:
#' 
#'  nhsd = function(X) { return(-0.5*scale*X[,colname]^2) },
#'  
#' that is, nhsd is half of the negative squared distance of an observation/integration point. 
#' The parameter 'scale' is currently set to 'truncation'. If truncation is NULL, scale is set to 1.
#'
#' @aliases model.halfnormal
#' @export
#' @param data \link{dsdata} data set. Used to determine distance truncation (if not provided)
#' @param truncation Distance to truncate at. Currently unused but passed on the the formula environment for later usage.
#' @param constrained If set to false a non-constrained linear effect is used for estimating the detection function. Handy for debugging. 
#' @param colname Effort data column to use for the covariate extraction. Default: "distance"
#' @param effect Character setting the name of the effect in the model. Default: "nhsd" ("negative half squared distance")
#' @return A \link{model} 
#' 

model.halfnormal = function(data = NULL, truncation = NULL, constrained = TRUE, colname = "distance", effect = "nhsd"){
  
  # Extract truncation limit from data (if data is given but limit is not)
  if ( !is.null (data) & is.null(truncation) ) { truncation = max(detdata(data)[, colname]) }
  
  # Formula
  if (!constrained) { formula = as.formula(paste0("~ . + ", effect)) }
  else { formula =  as.formula(paste0(" ~ . + f(", effect, ", model = 'clinear', range = c(0, Inf), hyper = list(beta = list(initial = 0, param = c(0,0.1))))" ))}
  
  # Covariate
  covariates = list()
  covariates[[effect]] = function(X) { return(-(0.5/scale)*X[,colname]^2) }
  assign("colname", colname, envir = environment(covariates[[effect]]))
  
  if (!is.null(truncation)) {
    assign("scale", truncation, envir = environment(covariates[[effect]]))
  } else {
    assign("scale", 1, envir = environment(covariates[[effect]]))
  }
  
  return(make.model(name = "Half-normal detection function",
                    formula = formula, 
                    effects = effect,
                    covariates = covariates))
}

#' Hazard rate detection function 
#' 
#' A taylor approximation to the hazard rate detection function. 
#'
#' This model approximates the detection probability
#'
#' \deqn{p(z) = 1 - exp[ - (z/sigma)^(-b) ]}
#' 
#' with parameters sigma and b and distance z. Note that the parameters are modeled in log-space. 
#' For data dependent sigma the 'conditional' parameter can be used. The respective model is
#' 
#' \deqn{f(z) = 1 - exp[ - (z/(sigma_1 + sigma_2*data)^(-b) ]}
#' 
#' @aliases model.hazard
#' @export
#' @name model.hazard
#' @param conditional Character array denoting the data that sigma depends on (linearly)
#' @param iterator An environment used to pass on the current state of the Taylor approximation
#' @return A \link{model}

model.hazard = function(conditional = NULL, iterator = new.env()){ 
  
  if ( !("haz.logb" %in% names(iterator)) ) {assign("haz.logb", log(1), iterator) }
  
  if ( is.null(conditional) ) {
    if ( !("haz.logsigma" %in% names(iterator)) ) {assign("haz.logsigma", log(1), iterator) }
    lhaz = expression( log(1-exp(-(loc$distance/(exp( haz.logsigma )))^(-exp(haz.logb)))) )
    df.mdl = model.taylor(expr = lhaz, iterator = iterator, effects = c("haz.logsigma","haz.logb"))
    df.mdl$iterator = list(haz.logsigma = iterator, haz.logb = iterator)
  } else {
    if ( !("haz.logsigma1" %in% names(iterator)) ) {assign("haz.logsigma1", log(1), iterator) }
    if ( !("haz.logsigma2" %in% names(iterator)) ) {assign("haz.logsigma2", log(1), iterator) }
    assign("conditional", conditional, iterator)
    lhaz = expression( log(1-exp(-(loc$distance/(exp( - haz.logsigma1 + loc[,conditional]*haz.logsigma2 )))^(-exp(haz.logb)))) )
    df.mdl = model.taylor(expr = lhaz, iterator = iterator, effects = c("haz.logsigma1", "haz.logsigma2", "haz.logb"))
    df.mdl$iterator = list(haz.logsigma1 = iterator, haz.logsigma2 = iterator, haz.logb = iterator)
  }
  
  return(df.mdl)
}



#' Expoential detection function model
#'
#' Construct a model representing an exponential detection function.
#' 
#' The log-likelihood of this model is
#' \deqn{log p(z) = -0.5 * q * z}
#' 
#' with precision parameter q and distance z. Internally, it is represented by the \link{formula}
#' 
#'  ~ . + f(df.exp, model = 'clinear', range = c(0, Inf), hyper = list(beta = list(initial = 0, param = c(0,0.1))))
#'
#' Hence, the effect is obtained from a contrained linear INLA model. The covariate "df.exp" for this model is as a 
#' function of the data.frame provided during the inference process:
#'  
#'  df.exp = function(X) { return(-X[,colname]) },
#'
#' @aliases model.exponential
#' @export
#' @param colname Effort data column to use for the covariate extraction. Default: "distance"
#' @param truncation Distance to truncate at. Currently unused but passed on the the formula environment for later usage.
#' @param constrained If set to false a non-constrained linear effect is used for estimating the detection function. Handy for debugging. 
#' @return A \link{model}

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
#' @export
#' @param colname Effort data column to use for the covariate extraction. Default: "distance"
#' @param truncation Maximum distance assumed to be observed..
#' @param segments Number of quadratic basis basis functions
#' @param clinear If set to false, non-constrained linear INLA effects are used for this model. 
#' This is handy for debugging but the detection function is not necessarily log-concave anymore.
#' @param linbasis Wether to use a linear basis function. If set to FALSE this implies that that the detection function has zero derivative at the origin. 
#' @return A \link{model}

model.logconcave = function(colname = "distance", truncation = NULL, segments = 5, constrained = TRUE, linbasis = TRUE, quadbasis = TRUE, ...){
  
  # Formula
  if (length(segments) == 1) {
    nSeg = segments
  } else {
    nSeg = length(segments)-1
  }
  
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
      return( data.frame(dfun_logconcave.basis.value(d = x[,colname], nSeg = segments, truncation = truncation))[,paste0("basis_",k)] )
    }
    environment(covariates[[paste0("basis_",k)]]) = new.env(parent=environment())
    environment(covariates[[paste0("basis_",k)]])$k = k
    environment(covariates[[paste0("basis_",k)]])$colname = colname
    environment(covariates[[paste0("basis_",k)]])$truncation = truncation
    environment(covariates[[paste0("basis_",k)]])$segments = segments
    environment(covariates[[paste0("basis_",k)]])$nSeg = nSeg
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

#' Smooth 2nd order random walk detection function with boundary constraint
#'
#'
#' @aliases model.smoothdf
#' @export
#' @param knots Knots at which the model operates. Usually a sequence from 0 to truncation distance.
#' @param effect Character defining the internal INLA effect name
#' 

model.smoothdf = function(dset, knots = seq(0, 6, by=1), degree = 2, effect = "smoothdf") {
  covariate = NULL
  mesh <- inla.mesh.1d(knots, degree = degree, boundary=c("neumann", "free"))
  # mesh$m = mesh$n-1
  
  tmp.mdl <- inla.spde2.matern(mesh, B.tau=cbind(0), B.kappa=cbind(-6))
  tmp.mdl$n.group = 1
  tmp.Q <- inla.spde.precision(tmp.mdl, theta=c())[-1,-1,drop=FALSE]
  msk = c(FALSE,rep(TRUE,length(knots)-1))
  
  formula <- as.formula(paste0("~ . + f(",effect,", model = 'generic0', Cmatrix = ", effect,".Q", ")"))
  assign(paste0(effect,".mdl"), tmp.mdl, envir = environment(formula))
  assign(paste0(effect,".Q"), tmp.Q, envir = environment(formula))
  
  # Little helper
  lassign = function(name,value) { tmp = list(); tmp[[name]] = value ; return(tmp) }
  
  mdl = make.model(name = "Smooth detection function",
                   formula = formula,
                   effects = effect,
                   mesh = lassign(effect, mesh),
                   covariates = lassign(effect, covariate),
                   inla.spde = lassign(effect, tmp.mdl),
                   mesh.coords = lassign(effect, "distance"),
                   A.msk = lassign(effect, msk))
}




#' Detectability at zero distance
#'
#'
#' @aliases model.g0fix
#' @export
#' @param g0 Probability of detecting an animal
#' 

model.g0fix = function(g0 = 1, effect = "g0fix", covariate = NULL) {
  
  formula = as.formula(paste0("~ . + f(", effect, ", model = 'clinear', fixed = TRUE, range = c(1, 1+1E-5))"))
  
  covariates = list()
  
  if ( is.null(covariate) ) {
    covariates[[effect]] = function(x) {
      v = rep(log(g0), nrow(as.data.frame(x)))
      ret = data.frame(v)
      colnames(ret) = effect
      return(ret) 
    }  
  } else {
    covariates[[effect]] = covariate
  }

  model = make.model(name = effect,
                    formula = formula, 
                    effects = effect,
                    covariates = covariates)
}

#' Prior distribution for detectability at zero distance
#' 
#' The probability g0 of detecting an animal at zero distance is implicitly modeled 
#' via a latent Gaussian distribution p(theta). The parameters "mean" and "sd" determine 
#' the mean and standard deviation of tau. Tau is then linked to g0 as
#' 
#'  g0 = exp(-exp(tau))
#'
#'
#' @aliases model.g0prior
#' @export
#' @param theta.mean Mean of the latent theta
#' @param theta.sd Standard deviation of the latent theta
#' @param effect Character setting the name of this effect in INLA
#' @param covariate A function returning a weight for exp(theta). Can be used to mask out g0 for particular sightings.
#' 

model.g0prior = function(theta.mean, theta.sd, effect = "g0prior", covariate = NULL) {
  
  hyper = list(theta = list(prior = paste0("expression: mean = ", theta.mean,"; sigma = ", theta.sd,"; 
                                          dens = 1/sqrt(2*pi) * 1/sigma * exp(-0.5*(x-mean)^2/sigma^2); 
                                          logdens = log(dens); 
                                          return(logdens)"),
                            initial = log(2), 
                            fixed = FALSE))
  
  formula = as.formula(paste0("~ . + f(", effect, ", model = 'clinear', range = c(0, Inf), fixed = FALSE, initial = 0, hyper =",effect,".hyper)"))
  
  environment(formula)[[paste0(effect,".hyper")]] = hyper
  
  covariates = list()
  if ( is.null(covariate) ) {
    covariates[[effect]] = function(x) {
      v = rep(-1, nrow(as.data.frame(x)))
      ret = data.frame(v)
      colnames(ret) = effect
      return(ret) 
    }
  } else { covariates[[effect]] = covariate }
  
  model = make.model(name = effect,
                     formula = formula, 
                     effects = effect,
                     covariates = covariates)
}



#' Spatial normal distribution
#'
#' @aliases model.normal
#' @export
#' @param data A \link{dsdata} object
#' @param iterator An environment to store the current state of the approximation
#' @param ... Arguments passed on to inla.spde2.matern. If none, the defaults are alpha = 2, prior.variance.nominal = 10, theta.prior.prec = 0.01

model.normal = function(data, effect, covariate, iterator = new.env(), ...) {
  
  mesh = data$mesh
  tau.name = paste0(effect,".tau")
  assign(effect, rep(mean(covariate(detdata(data))),data$mesh$n), envir = iterator)
  assign(tau.name, 1/sd(covariate(detdata(data))), envir = iterator)
  
  vargs = list(...)
  if ( length(vargs) == 0 ){ 
    gs.mdl.args = list(alpha = 2, prior.variance.nominal = 10, theta.prior.prec = 0.01) }
  else { 
    gs.mdl.args = vargs 
  }
  
  gs.mdl = do.call(inla.spde2.matern, c(list(mesh=mesh), gs.mdl.args))
  gs.mdl$n.group = 1
  formula = as.formula(paste0("~ . + f(", effect, ", model = gs.mdl) +", tau.name))
  
  # loglik = function(loc) { -0.5*exp(tau)*(covariate(loc) - as.numeric(inla.spde.make.A(mesh, as.matrix(loc[, data$mesh.coords])) %*% mode))^2 + 0.5*tau - log(sqrt(2*pi))}
  loglik = eval(parse(text=paste0("function(loc) { -0.5*exp(",tau.name,")*(covariate(loc) - 
                  as.numeric(inla.spde.make.A(mesh, as.matrix(loc[, data$mesh.coords])) %*% ",effect,"))^2 + 
                  0.5*",tau.name," - log(sqrt(2*pi))}")))
  environment(loglik) = new.env(parent = iterator)
  assign("covariate", covariate, envir = environment(loglik))
  assign("mesh", mesh, envir = environment(loglik))
  assign("data", data, envir = environment(loglik))
  
  # d.tau = function(loc) {-0.5*exp(tau)*(covariate(loc) - as.numeric(inla.spde.make.A(mesh = mesh, loc = as.matrix(loc[, data$mesh.coords])) %*% mode))^2 + 0.5 }
  d.tau = eval(parse(text=paste0("function(loc) {-0.5*exp(",tau.name,")*(covariate(loc) - 
            as.numeric(inla.spde.make.A(mesh = mesh, loc = as.matrix(loc[, data$mesh.coords])) %*% ",effect,"))^2 + 0.5 }")))
  environment(d.tau) = new.env(parent = iterator)
  assign("covariate", covariate, envir = environment(d.tau))
  assign("mesh", mesh, envir = environment(d.tau))
  assign("data", data, envir = environment(d.tau))
  
  # d.grps = function(loc) { exp(tau)*(covariate(loc) - as.numeric(inla.spde.make.A(mesh = mesh, loc = as.matrix(loc[, data$mesh.coords])) %*% mode)) }
  d.grps = eval(parse(text=paste0("function(loc) { exp(",tau.name,")*(covariate(loc) - 
            as.numeric(inla.spde.make.A(mesh = mesh, loc = as.matrix(loc[, data$mesh.coords])) %*% ",effect,")) }")))
  
  environment(d.grps) = new.env(parent = iterator)
  assign("covariate", covariate, envir = environment(d.grps))
  assign("mesh", mesh, envir = environment(d.grps))
  assign("data", data, envir = environment(d.grps))
  
  covariates = list()
  covariates[[effect]] = d.grps
  covariates[[paste0(effect,".tau")]] = d.tau  
  
  # const = function(loc) { loglik(loc) - d.tau(loc)*tau -d.grps(loc) * (as.numeric(inla.spde.make.A(mesh = mesh, loc = as.matrix(loc[, data$mesh.coords])) %*% mode))}
  const = eval(parse(text=paste0("function(loc) { loglik(loc) - 
          d.tau(loc)*",tau.name," -d.grps(loc) * (as.numeric(inla.spde.make.A(mesh = mesh, loc = as.matrix(loc[, data$mesh.coords])) %*% ",effect,"))}")))
  
  environment(const) = new.env(parent = iterator)
  assign("loglik", loglik, envir = environment(const))
  assign("d.tau", d.tau, envir = environment(const))
  assign("d.grps", d.grps, envir = environment(const))
  assign("mesh", mesh, envir = environment(const))
  assign("data", data, envir = environment(const))
  
  mklist = function(name, value) {ll=list(); ll[[name]]=value; return(ll)} 
  
  ret = make.model(name = "First order Taylor approximation of a spatial normal model",
                   formula = formula,
                   effects = c(effect, tau.name),
                   covariates = covariates,
                   mesh = mklist(effect, mesh),
                   inla.spde = mklist(effect, gs.mdl),
                   mesh.coords = mklist(effect, data$mesh.coords),
                   time.coords = mklist(effect, data$time.coords),
                   eval = eval, 
                   iterator = c(mklist(effect, iterator), mklist(tau.name, iterator))
                  )
  
  ret$const = const
  return(ret)
}


#' Spatial log group size distribution with scaling
#'
#' @aliases model.loggroupsize
#' @export
#' @param data A \link{dsdata} object
#' @param iterator An environment to store the current state of the approximation
#' @param ... Arguments passed on to inla.spde2.matern. If none, the defaults are alpha = 2, prior.variance.nominal = 10, theta.prior.prec = 0.01

model.loggroupsize = function(data, effect, covariate, iterator = new.env(), ...) {
  
  mesh = data$mesh
  tau.name = paste0(effect,".tau")
  assign(effect, rep(mean(covariate(detdata(data))),data$mesh$n), envir = iterator)
  assign(tau.name, 1/sd(covariate(detdata(data))), envir = iterator)
  
  vargs = list(...)
  if ( length(vargs) == 0 ){ 
    gs.mdl.args = list(alpha = 2, prior.variance.nominal = 10, theta.prior.prec = 0.01) }
  else { 
    gs.mdl.args = vargs 
  }
  
  gs.mdl = do.call(inla.spde2.matern, c(list(mesh=mesh), gs.mdl.args))
  gs.mdl$n.group = 1
  formula = as.formula(paste0("~ . + f(", effect, ", model = gs.mdl) +", tau.name))
  
  # loglik = function(loc) { -0.5*exp(tau)*(covariate(loc) - as.numeric(inla.spde.make.A(mesh, as.matrix(loc[, data$mesh.coords])) %*% mode))^2 + 0.5*tau - log(sqrt(2*pi))}
  loglik = eval(parse(text=paste0("function(loc) { -0.5*exp(",tau.name,")*(covariate(loc) - 
                                  as.numeric(inla.spde.make.A(mesh, as.matrix(loc[, data$mesh.coords])) %*% ",effect,"))^2 -
                                  as.numeric(inla.spde.make.A(mesh, as.matrix(loc[, data$mesh.coords])) %*% ",effect,") -
                                  0.5 * exp(-",tau.name,") +
                                  0.5*",tau.name," - log(sqrt(2*pi))}")))
  environment(loglik) = new.env(parent = iterator)
  assign("covariate", covariate, envir = environment(loglik))
  assign("mesh", mesh, envir = environment(loglik))
  assign("data", data, envir = environment(loglik))
  
  # d.tau = function(loc) {-0.5*exp(tau)*(covariate(loc) - as.numeric(inla.spde.make.A(mesh = mesh, loc = as.matrix(loc[, data$mesh.coords])) %*% mode))^2 + 0.5 }
  d.tau = eval(parse(text=paste0("function(loc) {-0.5*exp(",tau.name,")*(covariate(loc) - 
                                 as.numeric(inla.spde.make.A(mesh = mesh, loc = as.matrix(loc[, data$mesh.coords])) %*% ",effect,"))^2 + 0.5 + 0.5 * exp(-",tau.name,")}")))
  environment(d.tau) = new.env(parent = iterator)
  assign("covariate", covariate, envir = environment(d.tau))
  assign("mesh", mesh, envir = environment(d.tau))
  assign("data", data, envir = environment(d.tau))
  
  # d.grps = function(loc) { exp(tau)*(covariate(loc) - as.numeric(inla.spde.make.A(mesh = mesh, loc = as.matrix(loc[, data$mesh.coords])) %*% mode)) }
  d.grps = eval(parse(text=paste0("function(loc) { exp(",tau.name,")*(covariate(loc) - 
                                  as.numeric(inla.spde.make.A(mesh = mesh, loc = as.matrix(loc[, data$mesh.coords])) %*% ",effect,")) - 1 }")))
  
  environment(d.grps) = new.env(parent = iterator)
  assign("covariate", covariate, envir = environment(d.grps))
  assign("mesh", mesh, envir = environment(d.grps))
  assign("data", data, envir = environment(d.grps))
  
  covariates = list()
  covariates[[effect]] = d.grps
  covariates[[paste0(effect,".tau")]] = d.tau  
  
  # const = function(loc) { loglik(loc) - d.tau(loc)*tau -d.grps(loc) * (as.numeric(inla.spde.make.A(mesh = mesh, loc = as.matrix(loc[, data$mesh.coords])) %*% mode))}
  const = eval(parse(text=paste0("function(loc) { loglik(loc) - 
          d.tau(loc)*",tau.name," -d.grps(loc) * (as.numeric(inla.spde.make.A(mesh = mesh, loc = as.matrix(loc[, data$mesh.coords])) %*% ",effect,"))}")))
  
  environment(const) = new.env(parent = iterator)
  assign("loglik", loglik, envir = environment(const))
  assign("d.tau", d.tau, envir = environment(const))
  assign("d.grps", d.grps, envir = environment(const))
  assign("mesh", mesh, envir = environment(const))
  assign("data", data, envir = environment(const))
  
  mklist = function(name, value) {ll=list(); ll[[name]]=value; return(ll)} 
  
  ret = make.model(name = "First order Taylor approximation of a spatial normal model",
                   formula = formula,
                   effects = c(effect, tau.name),
                   covariates = covariates,
                   mesh = mklist(effect, mesh),
                   inla.spde = mklist(effect, gs.mdl),
                   mesh.coords = mklist(effect, data$mesh.coords),
                   time.coords = mklist(effect, data$time.coords),
                   eval = eval, 
                   iterator = c(mklist(effect, iterator), mklist(tau.name, iterator))
  )
  
  ret$const = const
  return(ret)
  }


#####################################
# MODELS: Group size
#####################################

#' Spatial SPDE model for log group size
#'
#' @aliases model.loggroupsize
#' @export
#' @param data A \link{dsdata} object
#' @param iterator An environment to store the current state of the approximation
#' @param ... Arguments passed on to inla.spde2.matern. If none, the defaults are alpha = 2, prior.variance.nominal = 10, theta.prior.prec = 0.01
# 
# model.loggroupsize = function(data, iterator = new.env(), ...) {
#   
#   mesh = data$mesh
#   rep(mean(detdata(dset)$lgrpsize),dset$mesh$n)
#   assign("grps", rep(mean(detdata(data)$lgrpsize),data$mesh$n), envir = iterator)
#   assign("tau", 1/sd(detdata(data)$lgrpsize), envir = iterator)
#   
#   vargs = list(...)
#   if ( length(vargs) == 0 ){ 
#     gs.mdl.args = list(alpha = 2, prior.variance.nominal = 10, theta.prior.prec = 0.01) }
#   else { 
#     gs.mdl.args = vargs 
#   }
#   
#   gs.mdl = do.call(inla.spde2.matern, c(list(mesh=mesh), gs.mdl.args))
#   gs.mdl$n.group = 1
#   formula = ~ . + f(grps, model = gs.mdl) + tau
#   
#   loglik = function(loc) { -0.5*exp(tau)*(loc$lgrpsize - as.numeric(inla.spde.make.A(mesh, as.matrix(loc[, data$mesh.coords])) %*% grps))^2 + 0.5*tau - log(sqrt(2*pi))}
#   environment(loglik) = new.env(parent = iterator)
#   assign("mesh", mesh, envir = environment(loglik))
#   assign("data", data, envir = environment(loglik))
#   
#   d.tau = function(loc) {-0.5*exp(tau)*(loc$lgrpsize - as.numeric(inla.spde.make.A(mesh = mesh, loc = as.matrix(loc[, data$mesh.coords])) %*% grps))^2 + 0.5 }
#   environment(d.tau) = new.env(parent = iterator)
#   assign("mesh", mesh, envir = environment(d.tau))
#   assign("data", data, envir = environment(d.tau))
#   
#   d.grps = function(loc) { exp(tau)*(loc$lgrpsize - as.numeric(inla.spde.make.A(mesh = mesh, loc = as.matrix(loc[, data$mesh.coords])) %*% grps)) }
#   environment(d.grps) = new.env(parent = iterator)
#   assign("mesh", mesh, envir = environment(d.grps))
#   assign("data", data, envir = environment(d.grps))
#   
#   covariates = list()
#   covariates[["grps"]] = d.grps
#   covariates[["tau"]] = d.tau  
#   const = function(loc) { loglik(loc) - d.tau(loc)*tau -d.grps(loc) * (as.numeric(inla.spde.make.A(mesh = mesh, loc = as.matrix(loc[, data$mesh.coords])) %*% grps))}
#   environment(const) = new.env(parent = iterator)
#   assign("loglik", loglik, envir = environment(const))
#   assign("d.tau", d.tau, envir = environment(const))
#   assign("d.grps", d.grps, envir = environment(const))
#   assign("mesh", mesh, envir = environment(const))
#   assign("data", data, envir = environment(const))
#   
#   ret = make.model(name = "First order Taylor expansion of spatial log group size model",
#                    formula = formula,
#                    effects = c("grps", "tau"),
#                    covariates = covariates,
#                    mesh = list(grps = mesh),
#                    inla.spde = list(grps = gs.mdl),
#                    mesh.coords = list(grps = data$mesh.coords),
#                    time.coords = list(spde = data$time.coords),
#                    eval = eval, iterator = list(grps = iterator, tau = iterator)
#   )
#   
#   ret$const = const
#   return(ret)
# }
# 


#' Model for Taylor approximation in a single variable 
#' 
#' See \link{model.taylor} for the real deal. This function is kept for means of debugging.
#'
#' @aliases model.taylor1d

model.taylor1d = function(expr = NULL, effects = NULL, initial = NULL, result = NULL) {
  
  effect = effects
  #formula = as.formula(paste0("~. +  f(", effect ,", model = 'clinear', range = c(0, 10), hyper = list(beta = list(initial = 0, param = c(0,0.1))))"))
  formula = as.formula(paste0("~. + ", effect))
  covariates = list()
    
  covariates[[effect]] = function(loc) {
    myenv <- new.env()
    assign(effect, initial[[effect]], envir = myenv)
    assign("loc", loc, envir = myenv)
    nderiv = numericDeriv(expr[[1]], c(effect), myenv)
    gradient = attr(nderiv,"grad")[,1]
    ret = data.frame(gradient)
    colnames(ret) = effect
    return(ret)
  } 
  
  const = function(loc) {
    myenv <- new.env()
    assign(effect, initial[[effect]], envir = myenv)
    assign("loc", loc, envir = myenv)
    nderiv = numericDeriv(expr[[1]], c(effect), myenv)
    gradient = attr(nderiv,"grad")[,1]
    ret = nderiv - gradient * initial[[effect]]
    return(ret) 
  }

  mdl = make.model(name = paste0("Taylor approximation in '",effect, "' of ",as.character(expr)), 
                   formula = formula, effects = effects, covariates = covariates)
  
  mdl$const = const
  return(mdl)
}



#' Model for Taylor approximations
#' 
#'
#' @aliases model.taylor
#' @param expr Expression approximated
#' @param effect List of charact arrays giving the variables that the expression is approximated in
#' @param iterator An environment giving the (inital) values of the effects
#' @param formula By default, the Variables of the Taylor approximation are modeled using INLA fixed effects. This can be changed by overwriting the formula of the model.

model.taylor = function(expr = NULL, effects = NULL, iterator = NULL, formula = NULL) {
  
  if (is.null(formula)) {
    formula = as.formula(paste0("~. + ", do.call(paste, c(sep = " +", as.list(effects)))   ))
  }
  
  if ( !is.environment(iterator) ) { iterator = as.environment(iterator) }
  
  covariates = list()
  grads = list()
  
  for (k in 1:length(effects)) { 
    effect = effects[k]
    covariates[[effect]] = function(loc) { 
      myenv = new.env(parent = iterator)
      assign("loc", loc, envir = myenv)
      nderiv = numericDeriv(expr[[1]], c(effect), myenv)
      gradient = attr(nderiv,"grad")[,1]
      ret = data.frame(gradient)
      colnames(ret) = effect
      return(as.vector(ret[,1]))
    }
    # clone environment
    environment(covariates[[effect]]) = new.env()
    assign("effect", effect, envir = environment(covariates[[effect]]))
    
    
    grads[[effect]] = function(loc) {
      myenv = new.env(parent = iterator)
      assign("loc", loc, envir = myenv)
      nderiv = numericDeriv(expr[[1]], c(effect), myenv)
      gradient = attr(nderiv,"grad")[,1]
      if (paste0(effect,".projector") %in% names(iterator)) {
        ret = gradient * iterator[[paste0(effect,".projector")]](loc)
      } else { ret = gradient * iterator[[effect]] }
      return(ret) 
    }
    environment(grads[[effect]]) = new.env()
    assign("effect", effect, envir = environment(grads[[effect]]))
    
  }
  
  value = function(loc){
    myenv = new.env(parent = iterator)
    assign(effect, iterator[[effect]], envir = myenv)
    assign("loc", loc, envir = myenv)
    nderiv = numericDeriv(expr[[1]], c(effect), myenv)
    return(nderiv)
  }
  
  
  const = function(loc) { 
    return( value(loc) - apply(do.call(cbind, lapply(grads, function(f) {f(loc)})), MARGIN=1, sum)  )
  }
  environment(const) = new.env()
  # assign("grads", grads, envir = environment(const))
  assign("covariates", covariates, envir = environment(const))
  assign("value", value, envir = environment(const))
  
  mdl = make.model(name = paste0("Taylor approximation in '",effect, "' of ",as.character(expr)), 
                   formula = formula, effects = effects, covariates = covariates)
  
  mdl$grad = grads
  mdl$const = const
  mdl$value = value
  mdl$iterator = lapply(effects, function(eff) { iterator })
  names(mdl$iterator) = effects
  return(mdl)
}


