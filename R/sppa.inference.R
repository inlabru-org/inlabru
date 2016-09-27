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
#' @aliases poiss
#' @export
#' @param points A SpatialPoints[DataFrame] object
#' @param model Typically a formula or a \link{model} describing the components of the LGCP density. If NULL, an intercept and a spatial SPDE component are used
#' @param mesh An inla.mesh object modelling the domain. If NULL, the mesh is constructed from a non-convex hull of the points provided
#' @return An \link{inla} object

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

#' Log Gaussian Cox process models using INLA
#'
#' @aliases lgcp
#' @export
#' @param points A SpatialPoints[DataFrame] object
#' @param samplers A Spatial[Points/Lines/Polygons]DataFrame objects
#' @param model Typically a formula or a \link{model} describing the components of the LGCP density. If NULL, an intercept and a spatial SPDE component are used
#' @param mesh An inla.mesh object modelling the domain. If NULL, the mesh is constructed from a non-convex hull of the points provided
#' @param ... Arguments passed on to iinla
#' @return An \link{inla} object

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

#' Generate a mesh from a hull of points
#'
#' @aliases default.mesh
#' @param spObject A Spatial* object
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
#' @param mesh An inla.mesh object
#' @return A \link{model} object

default.model = function(mesh) {
  model = join(model.spde(list(mesh = mesh)), model.intercept())
}

#' A wrapper for inla f() function 
#' 
#' g makes the mesh an explicit argument and catches the provided model for later usage
#'
#' @aliases g
#' @param mesh An inla.mesh
#' @param model See \link{f} model specifications
#' @param ... Passed on to inla \link{f}
#' @return A list with mesh, model and the return value of the f-call

g = function(mesh, model, ...) { c(list(mesh = mesh, model = model), f(model = model, ...)) }


#' Turn a formula into an iDistance \link{model}
#' 
#' To be used with formulae that use the \link{g} function as wrapper for inla's \link{f} function
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


#' Iterated INLA
#' 
#' This is a wrapper for iterated calls to \link{inla}.
#' Before each call the stackmaker function is used to set up the stack.
#' 
#' @aliases iinla
#' @export
#' @param data A data.frame
#' @param model A \link{model} object
#' @param stackmaker A function creating a stack from data and a model
#' @param n Number of \link{inla} iterations
#' @param idst.verbose If TRUE, be verbose (use verbose=TRUE to make INLA verbose)
#' @param ... Arguments passed on to \link{inla}
#' @return An \link{inla} object


iinla = function(data, model, stackmaker, n = 1, idst.verbose = FALSE, result = NULL, ...){
  # Inital stack
  stk = stackmaker(data, model)
  
  k = 1
  
  for ( k in 1:n ) {
    iargs = list(...) # Arguments passed on to INLA
    
    # When running multiple times propagate theta
    if ( k>1 ) {
      iargs[["control.mode"]] = list(restart = TRUE, theta = result$mode$theta)
    }
    
    # Verbose
    if ( idst.verbose ) { cat(paste0("Iteration: "),k, " ...") }
    
    # Return previous result if inla crashes, e.g. when connection to server is lost 
    if ( k > 1 ) { old.result = result } 
    
    result <- tryCatch( do.call(inla, c(list(formula = update.formula(model$formula, y.inla ~ .),
                                             data = c(inla.stack.data(stk), list.data(model)),
                                             control.predictor = list( A = inla.stack.A(stk), compute = TRUE),
                                             E = inla.stack.data(stk)$e), iargs)), 
                        error = function(e) { 
                          if (k == 1) { stop(e) }
                          else { 
                            cat(paste0("INLA crashed during iteration ",k,". It is likely that there is a convergence problem or the connection to the server was lost (if computing remotely)."))
                            return(old.result)
                          }
                        }
    )
    if ( idst.verbose ) { cat("done.\n") }
    
    # Update model
    update.model(model, result)
    model$result = result
    if ( n > 1 & k < n) { stk = stackmaker(data, model) }
  }
  result$stack = stk
  result$model = model
  class(result) = c("iinla", "inla", "list")
  return(result)
}
