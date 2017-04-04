#' @title Convenient model fitting using (iterated) INLA
#'
#' @description This method is a wrapper for \link{inla} and provides multiple enhancements. 
#' (1) For spatial data and models, \code{bru} will construct the required projection matrices automatically
#' (2) Multiple likelihoods can be employed in a convenient way
#' (3) Non-linear predictors are approximated numerically and fitting happens via iterated INLA calls.
#'
#' @aliases bru
#' @export
#' 
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 
#' @param components a formula describing the latent components
#' @param family Character defining one of the likelihoods supported by \link{like}. Alternatively, an object cunstructed by \link{like}.
#' @param data A data.frame or SpatialPoints[DataFrame] object
#' @param ... Additional likelihoods, each constructed by a calling \link{like}
#' @param mesh An \code{inla.mesh} object for spatial models without SPDE components. Mostly used for successive spatial predictions.
#' @param run If TRUE, run inference. Otherwise only return configuration needed to run inference.
#' @param max.iter maximum number of \link{iinla} iterations
#' @param offset the usual \link{inla} offset. If a nonlinear formula is used, the resulting Taylor approximation constant will be added to this automatically.
#' @param result An \code{inla} object returned from previous calls of \link{inla}, \link{bru} or \link{lgcp}. This will be used as a starting point for further improvement of the approximate posterior.
#' 
#' @return A \link{bru} object (inherits from iinla and \link{inla})
#' 
#' @example 
#' 
#' input.df <- data.frame(x=cos(1:10))
#' input.df <- within(input.df, y <- 5 + 2*cos(1:10) + rnorm(10, mean=0, sd=0.1))
#' fit.newbru <- bru(y ~ x, "gaussian", input.df)
#' summary(fit.newbru)
#' 

bru = function(components = y ~ Intercept,
               family = NULL,
               data = NULL,
               ...,
               mesh = NULL, 
               run = TRUE,
               max.iter = 10,
               offset = 0,
               result = NULL, 
               options = bru.options(control.compute = list(config = TRUE, dic = TRUE, waic = TRUE))) {
  
  # Automatically add Intercept and -1 to components unless -1 is in components formula
  components = auto.intercept(components)
  
  # Turn model components into internal bru model
  bru.model = make.model(components)
  
  # The `family` parameter can be either a string or a likelihood constructed
  # by like(). In the former case constrcut a proper likelihood using like() and
  # merge it with the list of likelihood provided via `...`.
  
  lhoods = list()
  if ( inherits(family, "lhood") & inherits(data, "lhood") ) {
    lhoods = c(lhoods, list(default = family, lh2 = data)) ; data = NULL ; family = NULL
  } else if ( inherits(family, "lhood") & !inherits(data, "lhood") ) {
    lhoods = c(list(default = family), lhoods) ; family = NULL 
  } else if ( !inherits(family, "lhood") & inherits(data, "lhood") ) {
    lhoods = c(list(default = like(family), lh2 = data), lhoods); data = NULL 
  } else {
    if( !is.null(family) ) { lhoods = c(list(default = like(family, data = data))); family = NULL }
  }
  lhoods = c(lhoods, list(...))
  
  
  # Handle each likelihood
  for (k in 1:length(lhoods)) {
    
    lh = lhoods[[k]]
    
    # Check if likelihood has data attached to it. If not, attach the 'data' argument or break if not available
    if ( is.null(lh$data) ) {
      if (is.null(data)) {stop(sprintf("Likelihood %s has not data attached to it and no data was supplied to bru() either.", names(lhoods)[[k]]))}
      lh$data = data
    }
    
    # Check the likelihood's response. If NULL, set the component's left hand side as response
    if ( is.null(lh$response) ) lh$response = all.vars(update(components, .~0))
    
    # Extract responses y for each likelihood and its data set
    lh$y = as.data.frame(lh$data)[, lh$response]
    
    # Create stackmaker for each likelihood (each stackmaker needs its own environment!)
    env = new.env() ; env$lh = lh
    
    lh$stackmaker = function(model, result) { 
      make.stack(points = lh$data, model = model, expr = lh$expr, y = lh$y, E = lh$E, result = result) 
    }
    environment(lh$stackmaker) = env
      
      
    lhoods[[k]] = lh
  }
  
  # Create joint stackmaker
  stk = function(xx, mdl, result) {
    do.call(inla.stack.mjoin, lapply(lhoods, function(lh) { lh$stackmaker(bru.model, result) }))
  }
  
  # Set max interations to 1 if all likelihood formulae are linear 
  if (all(sapply(lhoods, function(lh) lh$linear))) { max.iter = 1 }
  
  # Extract the family of each likelihood
  family = sapply(1:length(lhoods), function(k) lhoods[[k]]$family)
  
  # Run iterated INLA
  if ( run ) { result = do.call(iinla, c(list(data, bru.model, stk, family = family, n = max.iter, offset = offset, result = result), options))} 
  else { result = list() }
  
  
  ## Create result object ## 
  result$sppa$method = "bru"
  result$sppa$lhoods = lhoods
  result$sppa$model = bru.model
  result$sppa$mesh = mesh
  class(result) = c("bru", class(result))
  return(result)
}


#' Likelihood construction for usage with \link{bru}
#'
#' @aliases like
#' @export
#' 
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 
#' @param family A character identifying a valid \link{inla} likelihood. Alternatively 'cp' for Cox processes.
#' @param formula a \link{formula} where the right hand side expression defines the predictor used in the optimization.
#' @param data Likelihood-specific data
#' @param E Poisson exposure parameter (if family = 'poisson') passed on to \link{inla} 
#' @param samplers Integration domain for 'cp' family (not implemented!)
#' 
like = function(family, formula = . ~ ., data = NULL, E = 1, samplers = NULL) {
  
  # Does the likelihood formula imply a linear predictor?
  linear = as.character(formula)[2] == "."
  
  # If not linear, set predictor expression according to the formula's RHS
  if ( !linear ) { expr = parse(text = as.character(formula)[length(as.character(formula))]) }
  else { expr = NULL }
  
  response = all.vars(update(formula, .~0))
  if (response == ".") response = NULL
  
  lh = list(family = family, 
         formula = formula, 
         data = data, 
         E = E, 
         samplers = samplers, 
         linear = linear,
         expr = expr,
         response = response)
  
  class(lh) = c("lhood","list")
  
  # Return likelihood
  lh
}


#' Additional \link{bru} options
#'
#' @aliases bru.options
#' @export
#' 
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 
#' @param ...
#' 
bru.options = function(...) {
  list(...)
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
    model = make.model(model)
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
#' @param n maximum number of \link{iinla} iterations
#' @param offset the usual \link{inla} offset. If a nonline predictor is used, the resulting Taylor approximation constant will be added to this automatically.
#' @param result An \code{iinla} object returned from previous calls of \code{bru}/\code{lgcp}. This will be used as a starting point for further improvement of the approximate posterior.
#' @param ... Arguments passed on to \link{iinla}
#' @return An \link{iinla} object

lgcp = function(points, 
                samplers = NULL, 
                model = coordinates ~ spde(model = inla.spde2.matern(mesh), map = coordinates, mesh = mesh) + Intercept - 1, 
                predictor = . ~ ., 
                mesh = NULL,
                run = TRUE,
                scale = NULL,
                append = NULL,
                n = 10,
                offset = 0,
                result = NULL,
                weights.as.offset = FALSE,
                ...) {
  
  # Automatically complete moel and predictor
  ac = autocomplete(model, predictor, points, mesh)
  
  # If automatic completion detected linearity, expr to NULL
  if ( ac$linear ) { ac$expr = NULL ; n = 1 }
  
  # Turn model formula into internal bru model
  model = make.model(ac$model)

  # Create integration points
  model$dim.names = ac$lhs # workaround :(
  icfg = iconfig(samplers, points, model, mesh = ac$mesh)
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
  if ( weights.as.offset ) {
    stk = function(points, model, result) {
      inla.stack(make.stack(points = points, model = model, expr = ac$expr, y = 1, E = 0, result = result),
                 make.stack(points = ips, model = model, expr = ac$expr, y = 0, E = 1, offset = log(ips$weight), result = result))
    }
  } else {
    stk = function(points, model, result) {
      inla.stack(make.stack(points = points, model = model, expr = ac$expr, y = 1, E = 0, result = result),
                 make.stack(points = ips, model = model, expr = ac$expr, y = 0, E = ips$weight, offset = 0, result = result))
    }
  }

  # Run iterated INLA
  if ( run ) { result = iinla(points, model, stk, family = "poisson", n, offset = offset, result = result, ...) } 
  else { result = list() }
  
  ########## Create result object
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
  cat("### LGCP Summary #################################################################################\n\n")
  
  cat(paste0("Predictor: log(lambda) = ", as.character(result$model$expr),"\n"))
  
  cat("\n--- Points & Samplers ----\n\n")
  cat(paste0("Number of points: ", nrow(result$sppa$points)), "\n")
  if ( inherits(result$sppa$points,"Spatial") ) {
    cat(paste0("Coordinate names: ", paste0(coordnames(result$sppa$points), collapse = ", ")), "\n")
    cat(paste0("Coordinate system: ", proj4string(result$sppa$points), "\n"))
  }
  
  cat(paste0("Total integration weight: ", sum(result$ips$weight)), "\n")

  cat("\n--- Dimensions -----------\n\n")
  icfg = result$iconfig
  invisible(lapply(names(icfg), function(nm) {
    cat(paste0("  ",nm, " [",icfg[[nm]]$class,"]",
               ": ",
               "n = ",icfg[[nm]]$n.points,
               ", min = ",icfg[[nm]]$min,
               ", max = ",icfg[[nm]]$max,
               ", cardinality = ",signif(icfg[[nm]]$max-icfg[[nm]]$min),
               "\n"))
  }))
  
  summary.bru(result)
}

#' Summarize a bru object
#'
#' @aliases summary.bru
#' @export
#' @param result A result object obtained from a bru() run
#' 
summary.bru = function(result) {
  
  cat("\n--- Likelihoods ----------------------------------------------------------------------------------\n\n")
  for ( k in 1:length(result$sppa$lhoods) ) {
    lh = result$sppa$lhoods[[k]]
    cat(sprintf("Name: '%s', family: '%s', data class: '%s', \t formula: '%s' \n", names(result$sppa$lhoods)[[k]], lh$family, class(lh$data),deparse(lh$formula)))
  }
  
  #rownames(df) = names(result$sppa$lhoods)
  #print(df)
  
  cat("\n--- Criteria -------------------------------------------------------------------------------------\n\n")
  cat(paste0("Watanabe-Akaike information criterion (WAIC): \t", sprintf("%1.3e", result$waic$waic),"\n"))
  cat(paste0("Deviance Information Criterion (DIC): \t\t", sprintf("%1.3e", result$dic$dic),"\n"))
  
  cat("\n--- Fixed effects -------------------------------------------------------------------------------- \n\n")
  
  if ( nrow(result$summary.fixed)>0 ) {
  fe = result$summary.fixed
  fe$kld=NULL
  fe$signif = sign(fe[,"0.025quant"]) == sign(fe[,"0.975quant"])
  print(fe)
  } else { cat("None.\n") }
  
  cat("\n--- Random effects ------------------------------------------------------------------------------- \n\n")
  for ( nm in names(result$summary.random) ){
    sm = result$summary.random[[nm]]
    cat(paste0(nm,": "))
    cat(paste0("mean = [", signif(range(sm$mean)[1])," : ",signif(range(sm$mean)[2]), "]"))
    cat(paste0(", quantiles = [", signif(range(sm[,c(4,6)])[1])," : ",signif(range(c(4,6))[2]), "]"))
    if (nm %in% names(result$model$mesh)) {
      cat(paste0(", area = ", signif(sum(diag(as.matrix(inla.mesh.fem(result$model$mesh[[nm]])$c0))))))
    }
    cat("\n")
  }
  if ( length(names(result$summary.random)) == 0 ) {cat("None.\n")}
  
  cat("\n--- All hyper parameters (internal representation) ----------------------------------------------- \n\n")
  # cat(paste0("  ", paste(rownames(result$summary.hyperpar), collapse = ", "), "\n"))
  print(result$summary.hyperpar)
  
  
  marginal.summary = function(m, name) {
    df = data.frame(param = name,
                    mean = inla.emarginal(identity, m))
    df$var = inla.emarginal(function(x) {(x-df$mean)^2}, m)
    df$sd = sqrt(df$var)
    df[c("lq","median","uq")] = inla.qmarginal(c(0.025, 0.5, 0.975), m)
    df
  }
  
  cat("\n")
  for (nm in names(result$sppa$model$effects)) {
    eff = result$sppa$model$effects[[nm]]
    if (!is.null(eff$mesh)){
      hyp = inla.spde.result(result, nm, eff$inla.spde)
      cat(sprintf("\n--- Field '%s' transformed hyper parameters ---\n", nm))
      df = rbind(marginal.summary(hyp$marginals.range.nominal$range.nominal.1, "range"), 
                 marginal.summary(hyp$marginals.log.range.nominal$range.nominal.1, "log.range"), 
                 marginal.summary(hyp$marginals.variance.nominal$variance.nominal.1, "variance"),
                 marginal.summary(hyp$marginals.log.variance.nominal$variance.nominal.1, "log.variance"))
      print(df)
    }
  }
  

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


#' Predictions for iinla objects
#' 
#' Currently only a shortcut to \link{predict.lgcp}
#'  
#' @aliases predict.iinla
#' @export

predict.iinla = function(...) { predict.lgcp(...) }


#' Predictions based on bru
#' 
#' @aliases predict.bru
#' @export
#' @param result An object obtained by calling \link{bru})
#' @param data A data.frame or SpatialPointsDataFrame of covariates needed for the prediction
#' @param formula A formula determining which effects to predict and how to combine them
#' 
#' @return Predicted values

predict.bru = function(result,
                       data,
                       formula = NULL,
                       n.samples = 100)
{
  if ( is.null(formula) ) {
    lh.name = "default"
    if ( is.null(result$sppa$lhoods[[lh.name]]$expr) ) { 
      expr = parse(text = paste0(names(result$sppa$model$effects), collapse = "+")) }
    else {
      expr = parse(text = as.character(result$sppa$lhoods[[lh.name]]$formula)[3])
    }
  } else {
    expr = parse(text = as.character(formula)[length(as.character(formula))])
  }
  

  vals = evaluate.model(model = result$sppa$model, result = result, points = data, 
                        property = "sample", n = n.samples, predictor = expr)

  # Summarize
  vals = summarize(vals, x = data)
  
  vals
  
}







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
                        mcerr = 0.01,
                        verbose = FALSE)
  {
  
  # Extract target dimension from predictor (given as left hand side of formula)
  lhs.dims = setdiff(all.vars(update(predictor, .~0)), ".")
  pchar = as.character(predictor)
  predictor.rhs = pchar[length(pchar)]
  predictor = parse(text = predictor.rhs)
  
  # Dimensions we are integrating over
  idims = integrate
  

  # Alternatively, take dims from points
  if ( !is.null(points) ) {
    dims = names(points)
    if ( inherits(points, "Spatial") ) { dims = c(dims, "coordinates") }
    if ( !inherits(points, "SpatialPointsDataFrame") & (nrow(points) == 2)) { 
      points = SpatialPointsDataFrame(points, data = data.frame(weight = rep(1, nrow(coordinates(points)))))}
    icfg.points = points
  } else {
    icfg.points = result$sppa$points
  }
  
  
  # Dimensions we are integrating over
  idims = integrate
  
  # Determine all dimensions of the process (pdims)
  pdims = names(result$iconfig)
  
  # Check for ambiguous dimension definition
  if (any(idims %in% dims)) { stop("The left hand side of the predictor contains dimensions you want to integrate over. Remove them.") }
  
  # Collect some information that will be useful for plotting
  misc = list(predictor = predictor.rhs, dims = dims, idims = idims, pdims = pdims)
  type = "1d"
  
  # Generate points for dimensions to integrate over
  wicfg = iconfig(NULL, result$sppa$points, result$model, idims, mesh = result$sppa$mesh)
  wips = ipoints(samplers, wicfg[idims])
  
  # Determine all dimensions of the process (pdims)
  # pdims = names(result$iconfig)
  
  cat(sprintf("---- Generated points ----- \n"))
  
  
  ## Generate points for dimensions to integrate over
  ## Skip dimensions that are given via points
  gen.idims = setdiff(idims, names(points))
  if ( length(gen.idims) > 0 ) {
    wicfg = iconfig(NULL, result$sppa$points, result$model, gen.idims)
    wips = ipoints(samplers, wicfg[gen.idims])
  } else { wips = NULL }
  cat(sprintf("wips dims : %s \n", paste0(names(wips), collapse = ", ")))

  
  ## Generate points from predictor LHS
  ## Skip dimensions that we are integrating over or that are given via points
  gen.rdims = setdiff(lhs.dims, c(names(wips), names(points)))
  if ( length(gen.rdims) > 0 ) {
    icfg = iconfig(NULL, result$sppa$points, result$model, gen.rdims, mesh = result$sppa$mesh)
    rips = ipoints(NULL, icfg[gen.rdims])
    rips = rips[,setdiff(names(rips),"weight"),drop=FALSE]
    # Sort by value in dimension to plot over. Prevents scrambles prediction plots.
    if (!(gen.rdims[1] == "coordinates") & (length(gen.rdims)==1)) {rips = rips[sort(rips[,gen.rdims], index.return = TRUE)$ix,,drop=FALSE]}
  } else { rips = NULL }
  cat(sprintf("rips dims : %s \n", paste0(names(rips), collapse = ", ")))
  
  ## Provided points
  cat(sprintf("points : %s \n", paste0(names(points), collapse = ", ")))
  
  ## MERGE 
  all = list()
  if (!is.null(points)) { all = c(all, list(points))}
  if (!is.null(wips)) { all = c(all, list(wips))}
  if (!is.null(rips)) { all = c(all, list(rips))}
  if( length(all)>1 ) { pts = do.call(merge, all)} else { pts = all[[1]] }
  cat(sprintf("pts dims : %s \n", paste0(names(pts), collapse = ", ")))
  
  
  remain = list()
  if (!is.null(points)) { remain = c(remain, list(points))}
  if (!is.null(rips)) { remain = c(remain, list(rips))}
  if( length(remain)>1 ) { remain = do.call(merge, remain)} 
  else if (length(remain) ==1 ) { remain = remain[[1]] } 
  else { remain = NULL }
  dims = names(remain)
  if ( inherits(remain, "SpatialPoints") ) dims = c("coordinates", dims)
  cat(sprintf("remain dims : %s \n", paste0(dims, collapse = ", ")))
  
  ## Some additional info for plotting
  misc = list(predictor = predictor.rhs, dims = dims, idims = idims, pdims = names(pts))
  
  ## Define the sampling function 
  
  sample.fun = function(n) {
    vals = evaluate.model(model = result$sppa$model, result = result, points = pts, property = property, n = n, predictor = predictor)
    
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
        if ("coordinates" %in% dims ) { by.coords = cbind(1:nrow(remain), pts@data[,c(setdiff(dims, "coordinates"))]) } 
        else { by.coords = pts[,dims]}
        integral = do.call(rbind, by(vals, by.coords, colSums))
      }
    } else { integral = vals }

    integral
  }
  
  # If we are calculating a univariate density, do it properly. For multivariate predictions we currently only
  # compute rough estimates of the mean and the default inla quantiles
  if ( length(dims) == 0 ){
    
    # Pre-sample for bandwidth selection and intal interval x
    pre.smp = sample.fun(n)

    if (is.null(dens)) {
      component = function(x, smp) { approxfun(density(smp[[1]], kernel = "triangular", 
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
    type = "1d"
    smp = sample.fun(n)
    prd = summarize(smp, x = as.data.frame(remain))
    
    if ("coordinates" %in% dims ) { 
      type = "spatial"
      coordinates(prd) = coordnames(rips) 
      misc$p4s = icfg$coordinates$p4s
      misc$cnames = icfg$coordinates$cnames
      misc$mesh = icfg$coordinates$mesh
      prd = as.data.frame(prd)
      }
    
    attr(prd, "samples") = smp
    attr(prd, "total.weight") = sum(pts$weight)
    attr(prd, "type") = type
    attr(prd, "misc") = misc
    class(prd) = c("prediction", class(prd))
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

montecarlo.posterior = function(dfun, sfun, x = NULL, samples = NULL, mcerr = 0.01, n = 100, discrete = FALSE, verbose = FALSE) {

  xmaker = function(hpd) {
    mid = (hpd[2]+hpd[1])/2
    rg = (hpd[2]-hpd[1])/2
    x = seq(mid-1.2*rg, mid+1.2*rg, length.out = 256)
  }
  xmaker2 = function(hpd) {
    x = seq(hpd[1], hpd[2], length.out = 256)
  }
  
  inital.xmaker = function(smp) {
    mid = median(smp)
    rg = (quantile(smp,0.975)-quantile(smp,0.25))/2
    x = seq(mid-3*rg, mid+3*rg, length.out = 256)
  }
  
  # Inital samples
  if ( is.null(samples) ) { samples = sfun(n) }
  
  # Inital HPD
  if ( is.null(x) ) { x = inital.xmaker(as.vector(unlist(samples))) }
  
  # Round x if needed
  if (discrete) x = unique(round(x))

  # First density estimate
  lest = dfun(x, samples) 
  
  
  converged = FALSE
  while ( !converged ) {
    
    # Compute last HPD interval
    xnew = xmaker2(inla.hpdmarginal(0.999, list(x=x, y=lest)))
    
    # Map last estimate to the HPD interval
    if (discrete) xnew = unique(round(xnew))
    lest = inla.dmarginal(xnew, list(x=x, y=lest))  
    x = xnew
    
    # Sample new density
    n = 2 * n
    samples = sfun(n)
    est = dfun(x, samples)
    
    # Compute Monte Carlo error
    # err = sd(est/sum(est)-lest/sum(lest))
    err = max( ( (est-lest) / max(lest) )^2 )
    
    # Plot new density estimate versus old one (debugging)
    if ( verbose ) {
      cat(paste0("hpd:", min(x)," ",max(x), ", err = ", err, ", n = ",n, "\n")) 
      plot(x, lest, type = "l") ; lines(x, est, type = "l", col = "red")
    }
    
    # Convergence?
    if ( err < mcerr ) { converged = TRUE } 
    else { 
      lest =  0.5*(est + lest) 
    }
  }

  marg = list(x = x, y = est, samples = samples, mcerr = err)
  
  # Append some important statistics
  marg$quantiles = inla.qmarginal(c(0.025, 0.5, 0.975),marg)
  marg$mean = inla.emarginal(identity, marg) 
  marg$sd = sqrt(inla.emarginal(function(x) x^2, marg) - marg$mean^2) 
  marg$cv = marg$sd/marg$mean
  marg$mce = err
  
  marg
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
    smy$cv = smy$sd/smy$mean
    smy$var = smy$sd^2
  }
  if ( !is.null(x) ) { 
    if ( inherits(x, "Spatial") ) {
      if ( class(x) == "SpatialPoints" ) { smy = SpatialPointsDataFrame(x, data = smy) }
      else if ( class(x) == "SpatialPixels" ) { smy = SpatialPixelsDataFrame(x, data = smy) }
      else { x@data = cbind(x@data, smy) ; smy = x }
    }
    else {smy = cbind(x, smy) }
  }
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


iinla = function(data, model, stackmaker, n = 10, result = NULL, 
                 iinla.verbose = iinla.getOption("iinla.verbose"), 
                 control.inla = iinla.getOption("control.inla"), 
                 control.compute = iinla.getOption("control.compute"),
                 offset = NULL, ...){
  
  # # Default number of maximum iterations
  # if ( !is.null(model$expr) && is.null(n) ) { n = 10 } else { if (is.null(n)) {n = 1} }
  
  # Track variables?
  track = list()
  
  # Set old result
  old.result = result
  
  # Inital stack
  stk = stackmaker(data, model, result)
  
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
    result = NULL
    
    icall = expression(result <- tryCatch( do.call(inla, c(list(formula = update.formula(model$formula, y.inla ~ .),
                                                   data = c(inla.stack.mdata(stk), list.data(model)),
                                                   control.predictor = list( A = inla.stack.A(stk), compute = TRUE),
                                                   E = inla.stack.data(stk)$e,
                                                   offset = inla.stack.data(stk)$bru.offset + offset,
                                                   control.inla = control.inla,
                                                   control.compute = control.compute), iargs)), 
                              error = warning
                            )
                       )
    eval(icall)
    
    if ( is.character(result) ) { stop(paste0("INLA returned message: ", result)) }
    
    n.retry = 0
    max.retry = 10
    while ( ( is.null(result) | length(result) == 5 ) & ( n.retry <= max.retry ) ) {
      msg("INLA crashed or returned NULL. Waiting for 60 seconds and trying again.")
      Sys.sleep(60)
      eval(icall)
      n.retry = n.retry + 1
    } 
    if ( ( is.null(result) | length(result) == 5 ) ) { 
      msg(sprintf("The computation failed %d times. Giving up and returning last successfully obtained result.", n.retry-1))
      return(old.result)
    }
    
    if ( iinla.verbose ) { cat(" Done. ") }
    
    # Extract values tracked for estimating convergence
    if ( n > 1 & k < n) track[[k]] = cbind(effect = rownames(result$summary.fixed), iteration = k, result$summary.fixed)
    
    # Update stack given current result
    if ( n > 1 & k < n) { stk = stackmaker(data, model, result) }
    
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



#' Automatically complete model and predictor:
#' - Set RHS and LHS
#' - Add intercept
#' - Add mesh if needed
#' - substitute . with default SPDE model
#'
autocomplete = function(model, predictor, points, mesh) {
  
  
  # Automatically insert default SPDE model
  #env = environment(model)
  #model = update.formula(~ spde(model = inla.spde2.matern(mesh), map = coordinates, mesh = mesh), model)
  #environment(model) = env
  
  # Automatically add intercept (unless '-Intercept' is in the formula)
  env = environment(model)
  if (attr(terms(model),"intercept")) {
    if (!(length(grep("-[ ]*Intercept", as.character(model)[[length(as.character(model))]]))>0)) {
      model = update.formula(model, . ~ . + Intercept-1)
    } else {
      model = update.formula(model, . ~ . -1)
    }
    
  } 
  environment(model) = env
  
  # If predictor RHS has "." fill in effect names of model
  apred = auto.pred(model,predictor)
  
  # If RHS has "." we know that the predictor is linear
  linear = !(predictor == apred)
  
  # Set predictor to automatically generated one
  predictor = apred
  
  # Set model$expr as RHS of predictor 
  expr = parse(text = as.character(predictor)[3])
  
  # Extract y as from left hand side of the predictor
  lhs = all.vars(update(predictor, .~0))
  
  # Add mesh to model environment if needed
  if ( "mesh" %in% all.vars(model) & !("mesh" %in% names(environment(model)))) {
    if ( is.null(mesh) ) { mesh = default.mesh(points) }
    environment(model)$mesh = mesh
  }
  
  
  return(list(model = model, 
              predictor = predictor,
              linear = linear,
              lhs = lhs,
              expr = expr,
              mesh = mesh))
  
}


#'
#' A helper function used in autocomplete()
#'
auto.pred = function(model, predictor) {
  
  # Predictor LHS
  auto.pred.lhs = function(predictor, model){
    lhs = as.character(predictor)[2]
    if (lhs == ".") { lhs = as.character(model)[2] }
    lhs
  }
  
  # Predictor RHS
  auto.pred.rhs = function(model) {
    imodel = make.model(model)
    paste0(elabels(imodel), collapse = " + ")
  }
  
  # Predictor RHS is dot?
  if ( as.character(predictor)[3] == "." ) {
    predictor = as.formula(paste0(auto.pred.lhs(predictor,model)," ~ ",  paste0(auto.pred.rhs(model))))
  }
  predictor
}

auto.intercept = function(components) {
  env = environment(components)
  
  if (attr(terms(components),"intercept")) {
    if (!(length(grep("-[ ]*Intercept", as.character(components)[[length(as.character(components))]]))>0)) {
      components = update.formula(components, . ~ . + Intercept-1)
    } else {
      components = update.formula(components, . ~ . -1)
    }
    
  } 
  environment(components) = env
  components
}
