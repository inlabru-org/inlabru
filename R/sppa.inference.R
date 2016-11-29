#' Spatial model fitting using INLA
#'
#' @aliases bru
#' @export
#' @param points A SpatialPoints[DataFrame] object
#' @param model a formula describing the model components
#' @param mesh An inla.mesh object modelling the domain. If NULL, the mesh is constructed from a non-convex hull of the points provided
#' @return An \link{inla} object

bru = function(points,
               model = ~ spde(model = inla.spde2.matern(mesh), map = coordinates, mesh = mesh) + Intercept - 1, 
               predictor = y ~ spde + Intercept, 
               mesh = NULL, 
               run = TRUE,
               E = 1,
               family = "gaussian", ...) {
  
  # Construct default mesh
  if ( is.null(mesh) ) { mesh = default.mesh(points) }
  
  # Turn model formula into internal bru model
  model = as.model.formula(model, data.frame(points))
  
  # Set model$expr as RHS of predictor 
  pred.rhs = parse(text = as.character(predictor)[3])
  model$expr = pred.rhs
  
  # Extract y as from left hand side of the predictor
  pred.lhs = all.vars(update(predictor, .~0))
  y = as.data.frame(points)[,pred.lhs]
  
  # Create the stack
  stk = function(points, model) { detection.stack(points, model = model, y = y, E = E) }
  
  # Run iterated INLA
  if ( run ) { result = iinla(points, model, stk, family = family, ...) } 
  else { result = list() }
  
  # Create result object
  result$sppa$expr = pred.rhs
  result$sppa$method = "bru"
  result$sppa$model = model
  result$sppa$points = points
  result$sppa$stack = stk
  result$sppa$family = family
  if ( inherits(points, "SpatialPoints") ) {result$sppa$coordnames = coordnames(points)}
  class(result) = c("bru", class(result))
  return(result)
}


#' Multi-likelihood models
#' 
#' @details 
#' This function realizes model fitting using multiple likelihoods. Single likelihood can
#' be constructed by running the \link{bru} command using the parameter \code{run=FALSE}. 
#' The resulting object can be one of many passed on to \code{multibru} via the
#' \code{brus} parameter.
#'
#' @aliases multibru
#' @export
#' @param bru A list of \code{bru} objects (each obtained by a call of \link{bru}) 
#' @param model a formula describing the model components
#' @param ... further arguments passed on to iinla
#' @return An \link{inla} object

multibru = function(brus, model = NULL, ...) {
  
  # Default model: take it from first bru
  if ( is.null(model) ) { 
    model = brus[[1]]$sppa$model
  } else {
    model = as.model.formula(model, data.frame(brus[[1]]$sppa$points))
  }
  
  # Create joint stackmaker
  stk = function(xx, mdl) { 
    do.call(inla.stack.mjoin, lapply(brus, function(br) {br$sppa$stack(br$sppa$points, br$sppa$model)}))
  }
  
  #' Family
  family = as.character(lapply(brus, function(br) br$sppa$family))
  
  #' Run INLA
  result = iinla(NULL, model = model, stackmaker = stk, family = family, ...)
  
  # Create result object
  result$sppa$expr = lapply(brus, function(br) br$sppa$model$expr)
  result$sppa$method = "multibru"
  result$sppa$model = model
  result$sppa$points = brus[[1]]$sppa$points
  result$sppa$stack = stk
  result$sppa$family = family
  
  class(result) = c("bru", class(result))
  return(result)
}



#' Poisson regression using INLA
poiss = function(...) {
  stop("poiss() has been removed from the inlabru package. Please use bru() with family='poisson' instead. 
       Note: Instead of providing the exposure via the formula use the bru() E-parameter.")
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
#' @param append A list of functions which are evaluated for the \code{points} and the constructed integration points. The Returned values are appended to the respective data frame of points/integration points.
#' @param ... Arguments passed on to \link{iinla}
#' @return An \link{iinla} object

lgcp = function(points, 
                samplers = NULL, 
                model = ~ spde(model = inla.spde2.matern(mesh), map = coordinates, mesh = mesh) + Intercept - 1, 
                predictor = coordinates ~ spde + Intercept, 
                mesh = NULL, 
                scale = NULL,
                append = NULL, 
                ...) {
  
  if ( is.null(mesh) ) { mesh = default.mesh(points) }
  
  # Backwards compatibility: 
  # - If model has left hand side use it as left hand side of predictor
  # - If predictor is default use sum of effects in model as predictor RHS
  if (length(as.character(model)) == 3 ) { #( length(as.character(mdl)) == 3 )
    message("Note: You are using the old syntax where the left hand side of model defines the dimensions of the LGCP. Please use the LHS of predictor in the future.")
    prd.lhs = as.formula(paste0(paste(  all.vars(update.formula(model, .~0))  , collapse = " + "), "~."))
    if ( inherits(predictor, "expression") ) { predictor = as.formula(paste0("~", as.character(predictor))) }
    predictor = update.formula(predictor, prd.lhs)
    tmp = as.model.formula(model, data.frame(points))
    # if (!(toString(model) == toString(. ~ g(spde, model = inla.spde2.matern(mesh), map = coordinates, mesh = mesh) + Intercept - 1))){
    #   prd.rhs = as.formula(paste0(".~ ", paste(tmp$effects, collapse = " + ")))
    #   predictor = update.formula(predictor, prd.rhs)
    # } 
  }
  
  # Turn model formula into internal bru model
  model = as.model.formula(model, data.frame(points))
  
  # Set model$expr as RHS of predictor
  if ( inherits(predictor, "formula") ) {
    model$expr = parse(text = as.character(predictor)[3])
  }
  
  # Figure out LGCP dimensions. This is the left hand side of the predictor
  model$dim.names = all.vars(update(predictor, .~0))
    
  # Create integration points
  icfg = iconfig(samplers, points, model, mesh = mesh)
  ips = ipoints(samplers, icfg)
  
  # Apend covariates to points and integration points
  if (!is.null(append)) {
    append.tool = function(points, append) { 
      for ( k in 1:length(append) ) { 
        points[[names(append)[[k]]]] = append[[k]](points) 
      } 
      points 
    }
    ips = append.tool(ips, append)
    points = append.tool(points, append)
  }
  
  # If scale is not NULL, rescale integration weights
  if ( !is.null(scale) ) { ips$weight = scale * ips$weight }

  # Stack
  stk = function(points, model) { 
    inla.stack(detection.stack(points, model = model), 
               integration.stack(scheme = ips, model = model)) 
  }
  
  result = iinla(points, model, stk, family = "poisson", ...)
  result$mesh = mesh
  result$sppa$mesh = mesh
  result$ips = ips
  result$iconfig = icfg
  result$sppa$method = "lgcp"
  result$sppa$model = model
  result$sppa$points = points
  if ( inherits(points, "SpatialPoints") ) {result$sppa$coordnames = coordnames(points)}
  class(result) = c("lgcp",class(result))
  return(result)
}

#' Summarize a LGCP object
#'
#' @aliases summary.lgcp
#' @export
#' @param result A result object obtained from a lgcp() run
#' 
summary.lgcp = function(result) {
  cat("### LGCP Summary: #################################################################################\n")
  
  cat("\n Dimensions: \n")
  icfg = result$iconfig
  invisible(lapply(names(icfg), function(nm) {
    cat(paste0("  ",nm, " [",icfg[[nm]]$class,"]",
               ": ",
               "n=",icfg[[nm]]$n.points,
               ", min=",icfg[[nm]]$min,
               ", max=",icfg[[nm]]$max,
               ", cardinality=",signif(icfg[[nm]]$max-icfg[[nm]]$min),
               "\n"))
  }))
  
  cat("\n Effects: \n")
  cat(paste0("  ",paste(result$names.fixed, collapse = ", "),", ", paste(names(result$summary.random), collapse = ","), "\n"))
  
  cat("\n Hyper parameters: \n")
  cat(paste0("  ", paste(rownames(result$summary.hyperpar), collapse = ", "), "\n"))
  
  cat("\n Criteria: \n")
  cat(paste0("  Watanabe-Akaike information criterion (WAIC): \t", sprintf("%1.3e", result$waic$waic), " (Effective params: ",sprintf("%1.3e", result$waic$p.eff),")\n"))
  cat(paste0("  Deviance Information Criterion (DIC): \t\t", sprintf("%1.3e", result$dic$dic), " (Effective params: ",sprintf("%1.3e", result$dic$p.eff),")\n"))
}


#' Generate a simple default mesh
#'
#' @aliases default.mesh
#' @export
#' @param spObject A Spatial* object
#' @param max.edge A parameter passed on to \link{inla.mesh.2d} which controls the granularity of the mesh. If NULL, 1/20 of the domain size is used.
#' @return An \code{inla.mesh} object

default.mesh = function(spObject, max.edge = NULL, convex = -0.15){
  if (inherits(spObject, "SpatialPoints")) {
    x = c(bbox(spObject)[1,1], bbox(spObject)[1,2], bbox(spObject)[1,2], bbox(spObject)[1,1])
    y = c(bbox(spObject)[2,1], bbox(spObject)[2,1], bbox(spObject)[2,2], bbox(spObject)[2,2])
    # bnd = inla.mesh.segment(loc = cbind(x,y))
    # mesh = inla.mesh.2d(interior = bnd, max.edge = diff(bbox(spObject)[1,])/10)
    if ( is.null(max.edge) ) { max.edge = max.edge = diff(bbox(spObject)[1,])/20 }
    hull = inla.nonconvex.hull(points = coordinates(spObject), convex = convex)
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

g = function(..., map, A.msk, mesh, model) {
  if (as.character(substitute(map))[[1]] == "" ) { map = NULL } else { map = substitute(map) }
  if (as.character(substitute(A.msk))[[1]] == "" ) { A.msk = NULL } else { A.msk = A.msk }
  ret = list(mesh = mesh, 
             map = map,
             A.msk = A.msk,
             model = model, 
             f = f(..., model = model))
}

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

as.model.formula = function(fml, data) {
  
  # Define function for shifting the effect name into the g-call
  shift.names = function(fml) {
    tms = terms(fml)
    labels = attr(tms, "term.labels")
    effects = character()
    for (k in 1:length(labels)){
      lb = labels[[k]]
      gpd = getParseData(parse(text=lb))
      # Determine the name of the effect
      ename = getParseText(gpd, id=1)
      is.fixed = (gpd[1,"token"] == "SYMBOL")
      if ( !(ename %in% c("g","f","offset")) & !is.fixed ) {
        effects[[k]] = ename
        # Replace effect name by g(ename
        labels[[k]] = gsub(paste0(ename,"("), paste0("g(",ename,", "), lb, fixed = TRUE)
      }
    }
    labels
  }

  lbl = shift.names(fml)
  if ( length(lbl)> 0 ) {
    gidx = which(substr(lbl,1,2) == "g(")
    others = which(!substr(lbl,1,2) == "g(")
    base.fml.char = paste0("~. +", paste0("", lbl[others], collapse = " +"))
    mesh = list()
    mesh.coords = list()
    map = list()
    A.msk = list()
    inla.models = list()
    covariates = list()
    effects = lbl
    
    # Select g-terms
    
    for ( k in gidx ) {
      lb = lbl[[k]]
      # Extract mesh and spde model
      ge = eval(parse(text = lb), envir = environment(fml))
      name = ge$f$term
      effects[[k]] = name
      mesh[[name]] = ge$mesh
      mesh.coords[[name]] = ge$f$term
      inla.models[[name]] = ge$model
      map[[name]] = ge$map
      A.msk[[name]] = ge$A.msk
      
      # Replace function name by INLA f function
      lb = gsub("g\\(","f(",lb)
      
      # Remove extra mesh argument
      lb = gsub("[,][ ]*mesh[ ]*=[^),]*", "", lb)
      
      # Remove extra map argument
      pat = paste0("", deparse(ge$map))
      lb = gsub(deparse(ge$map), "", lb, fixed =TRUE)
      lb = gsub("[,][ ]*map[ ]*=[^),]*", "", lb)
      
      # Remove extra A.msk argument
      pat = paste0("", deparse(ge$A.msk))
      lb = gsub(deparse(ge$map), "", lb, fixed =TRUE)
      lb = gsub("[,][ ]*A.msk[ ]*=[^),]*", "", lb)
    
      lbl[[k]] = lb
      
      
    }
    for ( k in others ) {
      gpd = getParseData(parse(text=lbl[k]))
      if (gpd[1,"token"] == "SYMBOL") {
        if (!(gpd[1,"text"] %in% c(names(environment(fml)), names(data)))) { environment(fml)[[gpd[1,"text"]]] = 0 }
        if ( gpd[1,"text"] %in% names(data) ) {
          # covariates[[lbl[k]]] = function(x) {x[,effect,drop=FALSE]}
          # environment(covariates[[lbl[k]]]) = new.env()
          # assign("effect", lbl[k], envir = environment(covariates[[lbl[k]]]))
        } else {
          covariates[[lbl[k]]] = function(x) {
            v = rep(1, nrow(as.data.frame(x)))
            ret = data.frame(v)
            colnames(ret) = effect
            return(ret)
          }
          environment(covariates[[lbl[k]]]) = new.env()
          assign("effect", lbl[k], envir = environment(covariates[[lbl[k]]]))
        }
      
      } else if (gpd[1,"token"] == "expr") {
        old.label = lbl[k]
        lbl[k] = paste0(gpd[2,"text"],".effect")
        effects[k] = paste0(gpd[2,"text"],".effect")
        covariates[[lbl[k]]] = function(...) {do.call(function(...) {eval(parse(text=old.label), envir = list(...))}, as.list(...))}
      }
    }
    if ( (length(gidx) > 0) || length(others) > 0 ) {
      new.fml = as.formula(paste0("~.+", paste0("",paste0("", lbl, collapse = " + "))))
      # Add left hand side of fml
      new.fml = update.formula(update.formula(fml, ~ 1),new.fml)
      # If input formula had no intercept, add -1 to new formula
      if ( !(attr(terms(update.formula(~ DOT.FILL, fml)), "intercept") == 1) ) { new.fml = update.formula(new.fml, . ~ .-1) }
      # Add environment
      environment(new.fml) = environment(fml)
      # Make model
      mdl = make.model(formula = new.fml, name = "", mesh = mesh, effects = effects, covariates = covariates, 
                       inla.spde = inla.models, mesh.coords = mesh.coords, time.coords = NULL)
      mdl$map = map
      mdl$A.msk = A.msk
      mdl
    } else { return(NULL) }
  } else { return(NULL)}
}


#' Predictions for iinla objects
#' 
#' Currently only a shortcut to \link{predict.lgcp}
#'  
#' @aliases predict.iinla
#' @export

predict.iinla = function(...) { predict.lgcp(...) }


#' Predictions based on log Gaussian Cox processes
#' 
#' Takes a fitted LGCP produced by \link{lgcp} and predicts the expression stated by the right hand side of \code{predictor} for
#' a given set of \code{points}. A typical task is to predict values of the model for a sequence of points with consecutive values x_i, 
#' where x is one of the dimensions that the LGCP lives in. For convenience, this can be done automatically by leaving 
#' \code{points} set to \code{NULL} and stating a left hand side of the \code{predictor} formula, e.g. \code{x ~ f(x)}. The values
#' x_i for which the right hand side is evaluated are then either taken from the predictor's environment or, if not present, from
#' the integration points of the LGCP. If not points are provided and the right hand side is empty as well, \code{predict} assumes
#' that the predicted expression will evaluate so a single values for which the summary statistics are computed. Otherwise the
#' summary statistics are computed for each of the x_i. 
#' 
#' A common task is to look at statistics that result from integrating the 
#' right hand side expression over one or more dimensions, say y and z. This can be achieved by setting the \code{integrate}
#' parameter. The summary statistics are then computed for the values resulting from the integration. If only a subset of
#' the space should be integrated over the \code{samplers} argument can be used. Usually this is a Spatial* object (see \link{sp}).
#' 
#' @aliases predict.lgcp 
#' @export
#' @param result An lgcp object obtained by calling lgcp()
#' @param predictor A formula defining what to predict. The left hand side defines the dimension to predict over. The right hand side ought to be an expression that is evaluated for a given sample of the LGCP.
#' @param points Locations at which to predict. If NULL, the integration points of the LGCP are used
#' @param integrate A character array defining the dimensions to integrate the predicted expression over (e.g. "coordinates")
#' @param samplers A data.frame or Spatial* object that defines subsets of the dimensions to integrate over.
#' @param property By default ("sample") predictions are made using samples from the LGCP posterior. For debugging purposes this can also be set to the usual INLA summary statistics, i.e. "mean", "mode", "median", "sd" etc.
#' @param n Number of samples used to compute the predictor's statistics. The default is a rather low number, chose more for higher precision.
#' @param dens An expression defining an element of a micture density to be computed based on the samples generated by \code{predict}
#' @param discrete Set this to TRUE if \code{dens} is only defined for integer values
#' @param mcerr Monte Carlo error at which the sampling chain is stopped
#' @return Predicted values

predict.lgcp = function(result, 
                        predictor, 
                        points = NULL, 
                        integrate = NULL, 
                        samplers = NULL,  
                        property = "sample", 
                        n = 250,  
                        dens = NULL,
                        discrete = FALSE,
                        mcerr = 5e-4,
                        verbose = FALSE)
  {
  
  # Extract target dimension from predictor (given as right hand side of formula)
  dims = setdiff(all.vars(update(predictor, .~0)), ".")
  pchar = as.character(predictor)
  predictor.rhs = pchar[-1]
  predictor = parse(text = predictor.rhs)
  
  # Alternatively, take dims from points
  if ( !is.null(points) ) {
    dims = colnames(points)
    icfg.points = points
  } else {
    icfg.points = result$sppa$points
  }
  
  
  # Dimensions we are integrating over
  idims = integrate
  
  # Determine all dimensions of the process (pdims)
  pdims = names(result$iconfig)
  
  # Collect some information that will be useful for plotting
  misc = list(predictor = predictor.rhs, dims = dims, idims = idims, pdims = pdims)
  type = "1d"
  
  # Generate points for dimensions to integrate over
  wicfg = iconfig(NULL, result$sppa$points, result$model, idims)
  wips = ipoints(samplers, wicfg[idims])
  
  if ( length(dims) == 0 ) {
    pts = wips
    type = "full" 
  } else {
    # If no points for return dimensions were supplied we generate them
    if (is.null(points)) {
      icfg = iconfig(NULL, result$sppa$points, result$model, dims, mesh = result$sppa$mesh)
      rips = ipoints(NULL, icfg[dims])
      rips = rips[,setdiff(names(rips),"weight"),drop=FALSE]
      # Sort by value in dimension to plot over. Prevents scrambles prediction plots.
      if (!(dims[1] == "coordinates")) {rips = rips[sort(rips[,dims], index.return = TRUE)$ix,,drop=FALSE]}
    } else {
      rips = points
    }
    # Merge in integrations points if we are integrating
    if ( !is.null(wips) ) {
      pts = merge(rips, wips, by = NULL)
      if (!is.data.frame(pts)) {coordinates(pts) = coordnames(rips)}
    } else {
      pts = rips
      pts$weight = 1
      # Generate ifcg to set attribs of return value
      icfg = iconfig(NULL, icfg.points, result$model, dims, mesh = result$sppa$mesh)
    }
    
    # if ("coordinates" %in% dims ) { coordinates(pts) = coordnames(rips) }
  }
  
  sample.fun = function(n) {
    # Evaluate the model for these points
    vals = evaluate.model(result$sppa$model, result, pts, property = property, do.sum = TRUE, link = identity, n = n, predictor = predictor, use.covariate = FALSE)
    
    # If we sampled, summarize
    if ( is.list(vals) ) { vals = do.call(cbind, vals) }
    
    if ( !is.null(integrate) ) {
      # Weighting
      vals = vals * pts$weight
      
      # Sum up!
      if ( length(dims) == 0 ) {
        integral = as.list(colSums(vals))
        integral = lapply(integral, function(s) {list(integral=s)})
      } else {
        # If we are integrating over space we have to turn the coordinates into data that the by() function understands
        if ("coordinates" %in% dims ) { by.coords = cbind(1:nrow(rips), pts@data[,c(setdiff(dims, "coordinates"))]) } 
        else { by.coords = pts[,dims]}
        integral = do.call(rbind, by(vals, by.coords, colSums))
        # integral = summarize(integral, x = as.data.frame(rips))
        # if ("coordinates" %in% dims ) { coordinates(integral) = coordnames(rips) }
      }
    } else { integral = vals}

    integral
  }
  
  # If we are calculating a univariate density, do it properly. For multivariate predictions we currently only
  # compute rough estimates of the mean and the default inla quantiles
  if ( length(dims) == 0 ){
    
    # Pre-sample for bandwidth selection and intal interval x
    pre.smp = sample.fun(n)

    if (is.null(dens)) {
      component = function(x, smp) { approxfun(density(smp$integral, kernel = "triangular", 
                                         bw = bw.nrd(as.vector(unlist(pre.smp)))), rule = 0)(x) }
    } else {
      component = function(x,smp) eval(dens, c(list(x=x),smp))
    }
    
    mixer = function(x, smp) { apply(do.call(rbind, lapply(smp, function(s) component(x,s))), MARGIN = 2, mean)}
    prd = montecarlo.posterior(dfun = mixer, sfun = sample.fun, samples = pre.smp, mcerr = mcerr, n = n, discrete = discrete, verbose = verbose)
    
    class(prd) = c("prediction", class(prd))
    attr(prd, "type") = "full"
    attr(prd,"samples") = data.frame(integral = as.vector(unlist(prd$samples)))
    attr(prd, "misc") = misc

  } 
  
  # This is the case when predicting in more than one dimensions.
  # Calculating the exact posterior in this case is not implemented yet
  else {
  
    smp = sample.fun(n)
    prd = summarize(smp, x = as.data.frame(rips))
    if ("coordinates" %in% dims ) { coordinates(prd) = coordnames(rips) }
    
    if (inherits(prd, "SpatialPointsDataFrame")){
      type = "spatial"
      misc$p4s = icfg$coordinates$p4s
      misc$cnames = icfg$coordinates$cnames
      misc$mesh = icfg$coordinates$mesh
      prd = as.data.frame(prd)
    }
    
    attr(prd, "samples") = smp
    attr(prd, "total.weight") = sum(pts$weight)
    attr(prd, "type") = type
    attr(prd, "misc") = misc
    class(prd) = c("prediction",class(prd))
  }
  
  prd
  
  }


#' Monte Carlo method for estimating aposterior
#'  
#' @aliases montecarlo.posterior
#' @export
#' @param dfun A function returning a density for given x
#' @param sfun A function providing samples from a posterior
#' @param x Inital domain on which to perform the estimation. This will be adjusted as more samples are generated.
#' @param samples An initial set of samples. Not required but will be used to estimate the inital domain \code{x} if \code{x} is \code{NULL}
#' @param mcerr Monte Carlo error at which to stop the chain
#' @param n Inital number of samples. This will be doubled for each iteration.
#' @param discrete St this to \code{TRUE} if the density is only defined for integer \code{x}
#' @param verbose Be verbose?

montecarlo.posterior = function(dfun, sfun, x = NULL, samples = NULL, mcerr = 1e-4, n = 100, discrete = FALSE, verbose = FALSE) {

  xmaker = function(hpd) {
    mid = (hpd[2]+hpd[1])/2
    rg = (hpd[2]-hpd[1])/2
    x = seq(mid-2*rg, mid+2*rg, length.out = 256)
  }
  
  # Inital samples
  if ( is.null(samples) ) { samples = sfun(n) }
  
  # Inital HPD
  if ( is.null(x) ) { x = xmaker(range(as.vector(unlist(samples)))) }
  
  # Round x if needed
  if (discrete) x = round(x)

  # First density estimate
  lest = dfun(x, samples) 
  
  
  converged = FALSE
  while ( !converged ) {
    
    # Compute last HPD interval
    xnew = xmaker(inla.hpdmarginal(0.95, list(x=x, y=lest)))
    
    # Map last estimate to the HPD interval
    if (discrete) xnew = round(xnew)
    lest = inla.dmarginal(xnew, list(x=x, y=lest))  
    x = xnew
    
    # Sample new density
    n = 2 * n
    samples = sfun(n)
    est = dfun(x, samples)
    
    # Compute Monte Carlo error
    err = sd(est/sum(est)-lest/sum(lest))
    
    # Plot new density estimate versus old one (debugging)
    if ( verbose ) {
      cat(paste0("hpd:", min(x)," ",max(x), ", err = ", err, ", n = ",n, "\n")) 
      plot(x, lest, type = "l") ; lines(x, est, type = "l", col = "red")
    }
    
    # Convergence?
    if ( err < mcerr ) { converged = TRUE } 
    else { lest =  0.5*(est + lest) }
  }
  list(x = x, y = est, samples = samples, mcerr = err)
}  


#' Summarize and annotate data
#' 
#' @aliases summarize
#' @export
#' @param data A data.frame
#' @param x Annotations for the resulting summary
#' @return A data.frame with summary statistics
summarize = function(data, x = NULL, gg = FALSE) {
  if ( is.list(data) ) { data = do.call(cbind, data) }
  if ( gg ) {
    smy = rbind(
      data.frame(y = apply(data, MARGIN = 1, mean, na.rm = TRUE), property = "mean"),
      data.frame(y = apply(data, MARGIN = 1, sd, na.rm = TRUE), property = "sd"),
      data.frame(y = apply(data, MARGIN = 1, quantile, 0.025, na.rm = TRUE), property = "lq"),
      data.frame(y = apply(data, MARGIN = 1, quantile, 0.5, na.rm = TRUE), property = "median"),
      data.frame(y = apply(data, MARGIN = 1, quantile, 0.975, na.rm = TRUE), property = "uq"))
  }
  else { 
    smy = data.frame(
      apply(data, MARGIN = 1, mean, na.rm = TRUE),
      apply(data, MARGIN = 1, sd, na.rm = TRUE),
      t(apply(data, MARGIN = 1, quantile, prob = c(0.025, 0.5, 0.975), na.rm = TRUE)))
    colnames(smy) = c("mean", "sd", "lq", "median","uq")
  }
  if ( !is.null(x) ) { smy = cbind(x, smy)}
  return(smy)
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
                                             data = c(inla.stack.mdata(stk), list.data(model)),
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
