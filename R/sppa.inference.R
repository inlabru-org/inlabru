#' Normal regression using INLA
#'
#' @aliases gauss
#' @export
#' @param points A SpatialPoints[DataFrame] object
#' @param model Typically a formula or a \link{model} describing the components of the LGCP density. If NULL, an intercept and a spatial SPDE component are used
#' @param mesh An inla.mesh object modelling the domain. If NULL, the mesh is constructed from a non-convex hull of the points provided
#' @return An \link{inla} object

gauss = function(points, model = NULL, predictor = NULL, mesh = NULL, y = ~ k, ...) {
  
  if ( is.null(mesh) ) { mesh = default.mesh(points) }
  if ( is.null(model) ) { model = default.model(mesh) }
  if ( class(model)[[1]] == "formula" ) { 
    base.model = model.intercept()
    model = as.model.formula(model)
    if (!is.null(model)) { model = join(base.model, model) } else { model = join(base.model) }
  }
  if ( !is.null(predictor) ) { model$expr = expr ; stop("Not implemented: gaussian likelihood and non-linear predictor. missing: use const-parameter of INLA instead of E")}
  
  stk = function(points, model) { detection.stack(points, model = model, y = get_all_vars(y, points)[,1], E = 1) }
  
  result = iinla(points, model, stk, family = "gaussian", ...)
  result$mesh = mesh
  return(result)
}

#' Poisson regression using INLA
#' 
#' This function provides an easy-to-use interface to Poisson regression using INLA, in particular for spatial count data. 
#' The \code{points} parameter is used to provide the regression data, e.g. counts, exposures and covariates to regress on.
#' The \code{model} parameter, typically a \link{formula}, defines two aspects of the regression. On the left hand side
#' the (known) counts as well as the exposures are stated. The right hand side defines the log linear predictor of the 
#' regression. For non-linear regression the \code{predictor} parameter can be emplyed to override the linear structure
#' implied by the formula. Note that the \link{INLA} \code{family} argument can be overwritten and thereby other INLA 
#' likelihoods like the binomial and zero-inflated count models are accessible.
#'
#' @aliases poiss
#' @export
#' @param points A data frame or SpatialPointsDataFrame object
#' @param model Typically a formula or a \link{model} describing the linear regression predictor. If NULL, an intercept and a spatial SPDE component are used
#' @param predictor If NULL, the linear combination defined by the model/formula is used as a predictor for the counts. If a (possibly non-linear) expression is provided the respective Taylor approximation is used as a predictor. Multiple runs if INLA are then required for a better approximation of the posterior.
#' @param mesh An inla.mesh object modelling s spatial domain. If NULL, the mesh is constructed from a non-convex hull of the points provided
#' @return An \link{iinla} object

poiss = function(points, model = NULL, predictor = NULL, mesh = NULL, family = "poisson", ...) {
  
  if ( is.null(mesh) ) { mesh = default.mesh(points) }
  if ( is.null(model) ) { 
    model = default.model(mesh)
    model$formula = update.formula(model$formula, coordinates ~ .)
  }
  
  if ( class(model)[[1]] == "formula" ) {
    # Check if right hand side was provided
    if (as.character(model)[length(as.character(model))] == ".") {
      fml = model
      model = join.model(default.model(mesh))
      model$formula = update.formula(model$formula, fml)
      
    } else {
      fml = model
      if (attr(terms(fml), "intercept") == 1) { base.model = model.intercept() } else { base.model = NULL }
      more.model = as.model.formula(model)
      lhs = update.formula(fml, . ~ 0)
      model = join.model(more.model, base.model)
      rhs = reformulate(attr(terms(model$formula), "term.labels"), intercept = FALSE)
      model$formula = update.formula(lhs, rhs)
    }
  }
  if ( !is.null(predictor) ) { model$expr = predictor }
  
  yE = get_all_vars(update.formula(model$formula , . ~ 1), data = points)
  y = yE[,1]
  E = yE[,2]
  # y = as.numeric(get_all_vars(y, points)[,1])
  # E = get_all_vars(E, points)
  stk = function(points, model) { detection.stack(points, model = model, y = y, E = E) }
  
  result = iinla(points, model, stk, family = family, ...)
  result$mesh = mesh
  result$sppa$method = "poiss"
  if ( inherits(points, "SpatialPoints") ) {result$sppa$coordnames = coordnames(points)}
  return(result)
}

#' Log Gaussian Cox process (LGCP) inference using INLA
#' 
#' This function performs inference on a LGCP observed via points residing possibly multiple dimensions. 
#' These dimensions are defined via the left hand side of the formula provided via the model parameter.
#' The left hand side determines the intensity function that is assumed to drive the LGCP. This may include
#' effects that lead to a thinning (filtering) of the point process. By default, the log intensity is assumed
#' to be a linear combination of the effects defined by the formula's RHS. More sofisticated models, e.g.
#' non-linear thinning, can be achieved by using the predictor argument. The latter requires multiple runs
#' of INLA for improving the required approximation of the predictor (see \link{iinla}). In many applications
#' the LGCP is only observed through subsets of the dimensions the process is living in. For example, spatial
#' point realizations may only be known in sub-areas of the modeled space. These observed subsets of the LGCP
#' domain are called samplers and can be provided via the respective parameter. If samplers is NULL it is
#' assumed that all of the LGCP's dimensions have been observed completely. 
#' 
#'
#' @aliases lgcp
#' @export
#' @param points A data frame or SpatialPoints[DataFrame] object
#' @param samplers A data frame or Spatial[Points/Lines/Polygons]DataFrame objects
#' @param model Typically a formula or a \link{model} describing the components of the LGCP density. If NULL, an intercept and a spatial SPDE component are used
#' @param predictor If NULL, the linear combination defined by the model/formula is used as a predictor for the point location intensity. If a (possibly non-linear) expression is provided the respective Taylor approximation is used as a predictor. Multiple runs if INLA are then required for a better approximation of the posterior.
#' @param mesh An inla.mesh object modelling a spatial domain. If NULL and spatial data is provied the mesh is constructed from a non-convex hull of the points
#' @param scale If provided as a scalar then rescale the exposure parameter of the Poisson likelihood. This will influence your model's intercept but can help with numerical instabilities, e.g. by setting scale to a large value like 10000
#' @param ... Arguments passed on to \link{iinla}
#' @return An \link{iinla} object

lgcp = function(points, samplers = NULL, model = NULL, predictor = NULL, mesh = NULL, scale = NULL, ...) {
  
  if ( is.null(mesh) ) { mesh = default.mesh(points) }
  if ( is.null(model) ) { 
    model = default.model(mesh)
    model$formula = update.formula(model$formula, coordinates ~ .)
    }
  
  if ( class(model)[[1]] == "formula" ) {
    fml = model
    if (attr(terms(fml), "intercept") == 1) { base.model = model.intercept() } else { base.model = NULL }
    more.model = as.model.formula(model)
    lhs = update.formula(fml, . ~ 0)
    model = join.model(more.model, base.model)
    rhs = reformulate(attr(terms(model$formula), "term.labels"), intercept = FALSE)
    model$formula = update.formula(lhs, rhs)
  }
  
  if ( !is.null(predictor) ) { model$expr = predictor }
  
  # Create integration points
  icfg = iconfig(samplers, points, model)
  ips = ipoints(samplers, icfg)
  
  # If scale is not NULL, rescale integration weights
  if ( !is.null(scale) ) { ips$weight = scale * ips$weight }

  # Stack
  stk = function(points, model) { 
    inla.stack(detection.stack(points, model = model), 
               integration.stack(scheme = ips, model = model)) 
  }
  
  result = iinla(points, model, stk, family = "poisson", ...)
  result$mesh = mesh
  result$ips = ips
  result$iconfig = icfg
  result$sppa$method = "lgcp"
  result$sppa$model = model
  if ( inherits(points, "SpatialPoints") ) {result$sppa$coordnames = coordnames(points)}
  return(result)
}

#' Generate a simple default mesh
#'
#' @aliases default.mesh
#' @export
#' @param spObject A Spatial* object
#' @param max.edge A parameter passed on to \link{inla.mesh.2d} which controls the granularity of the mesh. If NULL, 1/20 of the domain size is used.
#' @return An \code{inla.mesh} object

default.mesh = function(spObject, max.edge = NULL){
  if (class(spObject) == "SpatialPointsDataFrame") {
    x = c(bbox(spObject)[1,1], bbox(spObject)[1,2], bbox(spObject)[1,2], bbox(spObject)[1,1])
    y = c(bbox(spObject)[2,1], bbox(spObject)[2,1], bbox(spObject)[2,2], bbox(spObject)[2,2])
    # bnd = inla.mesh.segment(loc = cbind(x,y))
    # mesh = inla.mesh.2d(interior = bnd, max.edge = diff(bbox(spObject)[1,])/10)
    if ( is.null(max.edge) ) { max.edge = max.edge = diff(bbox(spObject)[1,])/20 }
    hull = inla.nonconvex.hull(points = coordinates(spObject))
    mesh = inla.mesh.2d(boundary = hull, max.edge = max.edge)
  } else {
    NULL
  }
}


#' Generate a default model with a spatial SPDE component and an intercept
#'
#' @aliases default.model
#' @export
#' @param mesh An inla.mesh object
#' @return A \link{model} object

default.model = function(mesh) {
  model = join(model.spde(list(mesh = mesh)), model.intercept())
}

#' A wrapper for the \link{inla} \link{f} function
#' 
#' g makes the mesh an explicit argument and catches the provided model for later usage.
#' See \link{f} for details on model specifications.
#'
#' @aliases g
#' @export
#' @param mesh An inla.mesh
#' @param model See \link{f} model specifications
#' @param ... Passed on to inla \link{f}
#' @return A list with mesh, model and the return value of the f-call

g = function(mesh, model, ...) { c(list(mesh = mesh, model = model), f(model = model, ...)) }


#' Turn a formula into an iDistance \link{model}
#' 
#' To be used with formulae that use the \link{g} function as wrapper for inla's \link{f} function
#' Warning: this functions mainly exists for reasons of back-compatibilty. Avoid using it at all
#' costs.
#'
#' @aliases as.model.formula
#' @export
#' @param fml A formula
#' @return A \link{model} object

as.model.formula = function(fml) {
  
  tms = terms(fml)
  lbl = attr(tms, "term.labels")
  if ( length(lbl)> 0 ) {
    gidx = which(substr(lbl,1,1) == "g")
    others = which(!substr(lbl,1,1) == "g")
    base.fml.char = paste0("~. +", paste0("", lbl[others], collapse = " +"))
    mesh = list()
    mesh.coords = list()
    # mesh.map = list()
    inla.models = list()
    covariates = list()
    effects = lbl
    
    # Select g-terms
    
    for ( k in gidx ) {
      lb = lbl[[k]]
      # Extract mesh and spde model
      ge = eval(parse(text = lb), envir = environment(fml))
      mesh[[ge$term]] = ge$mesh
      mesh.coords[[ge$term]] = ge$term
      # mesh.map[[ge$term]] = function(x) { x[,ge$term,drop=FALSE]}
      inla.models[[ge$term]] = ge$model
      effects[[gidx]] = ge$term
      # Replace function name by INLA f function
      lb = gsub("g\\(","f(",lb)
      # Remove extra mesh argument
      lb = gsub("[,][ ]*mesh[ ]*=[^),]*", "", lb)
      lbl[[k]] = lb
    }
    for ( k in others ) {
      gpd = getParseData(parse(text=lbl[k]))
      if (gpd[1,"token"] == "SYMBOL" && gpd[1,"text"] %in% names(environment(fml))) {
        covariates[[lbl[k]]] = function(x) {
          v = rep(1, nrow(as.data.frame(x)))
          ret = data.frame(v)
          colnames(ret) = effect
          return(ret) 
        }
      environment(covariates[[lbl[k]]]) = new.env()
      assign("effect", lbl[k], envir = environment(covariates[[lbl[k]]]))
      }
    }
    if ( (length(gidx) > 0) || length(others) > 0 ) {
      new.fml = as.formula(paste0("~.+", paste0("",paste0("", lbl, collapse = " + "))))
      # Add left hand side of fml
      new.fml = update.formula(update.formula(fml, ~ 1),new.fml)
      # Add environment
      environment(new.fml) = environment(fml)
      # Make model
      mdl = make.model(formula = new.fml, name = "", mesh = mesh, effects = effects, covariates = covariates, 
                       inla.spde = inla.models, mesh.coords = mesh.coords, time.coords = NULL)  
    } else { return(NULL) }
  } else { return(NULL)}
}


#' Predictions based on log Gaussian Cox processes
#' 
#' @aliases predict.lgcp 
#' @export
#' @param result An object obtained by ralling lgcp()
#' @return Predicted values

predict.lgcp = function(result,  predictor = result$model$expr, model = NULL, samplers = NULL, property = "sample", n = 100, postproc = summarize) {
  
  # Is the predictor is either supplied as an expression?
  if ( !(substitute(predictor)[[1]] == "expression" ) )  { predictor = as.expression(substitute(predictor)) }
  
  # Determine which dimensions NOT to integrate over. Via the model paramter the user can specify
  # either characters (dimension names), a formula where the left hand side determines the dimensions
  # or a model where model$ formula determines the dimensions that are not integrated over.
  if ( is.null(model) ) {
    dims = NULL
  } else if ( is.character(model) ) {
    dims = model
  } else if ( class(model)[1] == "formula" ) {
    dims = all.vars(update.formula(model, . ~ 0))
  } else {
    dims = all.vars(update.formula(model$formula, . ~ 0))
  }
  
  # Determine all dimensions of the process
  pdims = names(result$iconfig)
  
  # Determine dimensions to intagrate over
  idims = setdiff(pdims, dims)
  
  # Generate points for dimensions to integrate over
  wips = ipoints(samplers, result$iconfig[idims])
  
  if ( is.null(dims) ) { 
    pts = wips
  } else {
    # Generate points for dimensions that we are returning. 
    # Remove weight and then merge with dims to integrate over
    rips = ipoints(NULL, result$iconfig[dims])
    rips = rips[,setdiff(names(rips),"weight"),drop=FALSE]
    pts = merge(rips[,setdiff(names(rips),"weight"),drop=FALSE], wips)
    if ("coordinates" %in% dims ) { coordinates(pts) = coordnames(rips) }
  }
  
  # Evaluate the model for these points
  vals = evaluate.model(result$sppa$model, result, pts, property = property, do.sum = TRUE, link = identity, n = n, predictor = predictor)
  
  # If we sampled, summarize
  if ( is.list(vals) ) { vals = do.call(cbind, vals) }
  
  # Weighting
  vals = vals * pts$weight
  
  # Sum up!
  if ( is.null(dims) ) {
    integral = data.frame(colSums(vals)) ; colnames(integral) = "integral"
    if ( !is.null(postproc) ) { integral = postproc(t(integral)) }
  } else {
    # If we are integrating over space we have to turn the coordinates into data that the by() function understands
    if ("coordinates" %in% dims ) { by.coords = cbind(coordinates(pts), pts@data[,c(setdiff(dims, "coordinates"))]) } 
    else { by.coords = pts[,dims]}
    
    integral = do.call(rbind, by(vals, by.coords, colSums))
    if ( !is.null(postproc) ) { integral = postproc(integral, x = as.data.frame(rips)) }
    if ("coordinates" %in% dims ) { coordinates(integral) = coordnames(rips) }
  }
  
  # return
  integral
}





#' Iterated INLA
#' 
#' This is a wrapper for iterated runs of \link{inla}. Before each run the \code{stackmaker} function is used to
#' set up the \link{inla.stack} for the next iteration. For this purpose \code{stackmaker} is called given the
#' \code{data} and \code{model} arguments. The \code{data} argument is the usual data provided to \link{inla}
#' while \link{model} provides more information than just the usual inla formula. 
#' 
#' @aliases iinla
#' @export
#' @param data A data.frame
#' @param model A \link{model} object
#' @param stackmaker A function creating a stack from data and a model
#' @param n Number of \link{inla} iterations
#' @param iinla.verbose If TRUE, be verbose (use verbose=TRUE to make INLA verbose)
#' @param ... Arguments passed on to \link{inla}
#' @return An \link{inla} object


iinla = function(data, model, stackmaker, n = NULL, result = NULL, 
                 iinla.verbose = iinla.getOption("iinla.verbose"), 
                 control.inla = iinla.getOption("control.inla"), 
                 control.compute = iinla.getOption("control.compute"), ...){
  
  if ( iinla.verbose ) { cat("iinla(): Let's go.\n") }
  
  # Default number of maximum iterations
  if ( !is.null(model$expr) && is.null(n) ) { n = 10 } else { if (is.null(n)) {n = 1} }
  
  # Track variables?
  track = list()
  
  # Inital stack
  stk = stackmaker(data, model)
  
  k = 1
  interrupt = FALSE
  
  while ( (k <= n) & !interrupt ) {
    iargs = list(...) # Arguments passed on to INLA
    
    # When running multiple times propagate theta
    if ( k>1 ) {
      iargs[["control.mode"]] = list(restart = TRUE, theta = result$mode$theta)
    }
    
    # Verbose
    if ( iinla.verbose ) { cat(paste0("iinla() iteration"),k,"[ max:", n,"].") }
    
    # Return previous result if inla crashes, e.g. when connection to server is lost 
    if ( k > 1 ) { old.result = result } 
    
    result <- tryCatch( do.call(inla, c(list(formula = update.formula(model$formula, y.inla ~ .),
                                             data = c(inla.stack.data(stk), list.data(model)),
                                             control.predictor = list( A = inla.stack.A(stk), compute = TRUE),
                                             E = inla.stack.data(stk)$e,
                                             control.inla = control.inla,
                                             control.compute = control.compute), iargs)), 
                        error = function(e) { 
                          if (k == 1) { stop(e) }
                          else { 
                            cat(paste0("INLA crashed during iteration ",k,". It is likely that there is a convergence problem or the connection to the server was lost (if computing remotely)."))
                            return(old.result)
                          }
                        }
    )
    if ( iinla.verbose ) { cat(" Done. ") }
    
    # Update model
    update.model(model, result)
    model$result = result
    track[[k]] = cbind(effect = rownames(result$summary.fixed), iteration = k, result$summary.fixed)
    if ( n > 1 & k < n) { stk = stackmaker(data, model) }
    
    # Stopping criterion
    if ( k>1 ){
      max.dev = 0.01
      dev = do.call(c, lapply(by(do.call(rbind, track), as.factor(do.call(rbind, track)$effect), identity), function(X) { abs(X$mean[k-1] - X$mean[k])/X$sd[k] }))
      cat(paste0("Max deviation from previous: ", signif(100*max(dev),3),"% of SD [stop if: <",100*max.dev,"%]\n"))
      interrupt = all( dev < max.dev)
      if (interrupt) {cat("Convergence criterion met, stopping INLA iteration.")}
    } else {
      cat("\n")
    }
    k = k+1
  }
  result$stack = stk
  result$model = model
  result$iinla$track = do.call(rbind, track)
  class(result) = c("iinla", "inla", "list")
  return(result)
}
