# GENERICS
generate = function(...){UseMethod("generate")}

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
#' @param options See \link{bru.options} 
#' @return A \link{bru} object
#' 
#' @examples
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
               options = list()) {
  
  # Update default options
  options = do.call(bru.options, options)
  
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
    if( !is.null(family) ) { lhoods = c(list(default = like(family, data = data, E = options$E))); family = NULL }
  }
  lhoods = c(lhoods, list(...))
  
  
  # If a likelihood was provided without data/response, update according to bru's 
  # arguments `data` and LHS of components
  
  for (k in 1:length(lhoods)) {
    
    lh = lhoods[[k]]
    
    # Check if likelihood has data attached to it. If not, attach the 'data' argument or break if not available
    if ( is.null(lh$data) ) {
      if (is.null(data)) {stop(sprintf("Likelihood %s has not data attached to it and no data was supplied to bru() either.", names(lhoods)[[k]]))}
      lh$data = data
    }
    if ( is.null(lh$components) ) { lh$components = components }
    if ( is.null(lh$response) ) { lh$response = all.vars(update(components, .~0)) }
    
    lhoods[[k]] = lh
  }
  
  # Create joint stackmaker
  stk = function(xx, mdl, result) {
    do.call(inla.stack.mjoin, lapply(lhoods, function(lh) { stackmaker.like(lh)(bru.model, result) }))
  }
  
  # Set max interations to 1 if all likelihood formulae are linear 
  if (all(sapply(lhoods, function(lh) lh$linear))) { options$max.iter = 1 }
  
  # Extract the family of each likelihood
  family = sapply(1:length(lhoods), function(k) lhoods[[k]]$inla.family)
  
  # Run iterated INLA
  if ( options$run ) { result = do.call(iinla, list(data, bru.model, stk, family = family, n = options$max.iter, offset = options$offset, result = options$result, inla.options = options$inla.options))} 
  else { result = list() }
  
  ## Create result object ## 
  result$sppa$method = "bru"
  result$sppa$lhoods = lhoods
  result$sppa$model = bru.model
  result$sppa$mesh = options$mesh
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
#' @param components Components
#' @param mesh An inla.mesh object
#' @param E Exposure parameter for family = 'poisson' passed on to \link{inla}. Special case if family is 'cp': rescale all integration weights by E.
#' @param samplers Integration domain for 'cp' family
#' @param ips Integration points for 'cp' family. Overrides \code{samplers}
#' @param domain Named list of domain definitions
#' 
like = function(family, formula = . ~ ., data = NULL, components = NULL, mesh = NULL, E = 1, samplers = NULL, ips = NULL, domain = NULL) {
  
  # Some defaults
  inla.family = family
  
  # Does the likelihood formula imply a linear predictor?
  linear = as.character(formula)[length(as.character(formula))] == "."
  
  # If not linear, set predictor expression according to the formula's RHS
  if ( !linear ) { expr = parse(text = as.character(formula)[length(as.character(formula))]) }
  else { expr = NULL }
  
  #' Set response name
  response = all.vars(update(formula, .~0))
  if (response[1] == ".") response = NULL
  
  
  #' More on special bru likelihoods
  if ( family == "cp" ) {
    if ( is.null(data) ) { stop("You called like() with family='cp' but no 'data' argument was supplied.") }
    #if ( is.null(samplers) ) { stop("You called like() with family='cp' but no 'samplers' argument was supplied.") }
    bru.model = make.model(components)
    if (as.character(formula)[2] == ".") { 
      bru.model$dim.names = all.vars(update(components, .~0)) } 
    else { 
      bru.model$dim.names = all.vars(update(formula, .~0))
    }
    
    if ( is.null(ips) ) {
      icfg = iconfig(samplers, data, bru.model, mesh = mesh, domain = domain)
      ips = ipmaker(samplers, icfg)
    }
    
    inla.family = "poisson"
  }
  
  # Calculate data ranges
  drange = lapply(names(data), function(nm) {  if(is.numeric(data[[nm]])) {range(data[[nm]])} else {NULL} } )
  names(drange) = names(data)
  if ( inherits(data, "Spatial") ) drange[["coordinates"]] = mesh

  
  # The likelihood object that will be returned
  
  lh = list(family = family, 
         formula = formula, 
         data = data, 
         E = E, 
         samplers = samplers, 
         linear = linear,
         expr = expr,
         response = response,
         inla.family = inla.family,
         ips = ips,
         domain = domain,
         drange = drange)
  
  class(lh) = c("lhood","list")
  
  # Return likelihood
  lh
}

stackmaker.like = function(lhood) {
  
  env = new.env() ; env$lhood = lhood
  
  # Special inlabru likelihoods
  if (lhood$family == "cp"){
    sm = function(model, result) {
      inla.stack(make.stack(points = lhood$data, model = model, expr = lhood$expr, y = 1, E = 0, result = result),
                 make.stack(points = lhood$ips, model = model, expr = lhood$expr, y = 0, E = lhood$E * lhood$ips$weight, offset = 0, result = result))
    }
    
  } else { 
    
    sm = function(model, result) { 
      make.stack(points = lhood$data, model = model, expr = lhood$expr, y = as.data.frame(lhood$data)[,lhood$response], E = lhood$E, result = result) 
    }
    
  }
  environment(sm) = env
  sm
}


#' Additional \link{bru} options
#'
#' @aliases bru.options
#' @export
#' 
#' @param mesh An \code{inla.mesh} object for spatial models without SPDE components. Mostly used for successive spatial predictions.
#' @param run If TRUE, run inference. Otherwise only return configuration needed to run inference.
#' @param max.iter maximum number of inla iterations
#' @param offset the usual \link{inla} offset. If a nonlinear formula is used, the resulting Taylor approximation constant will be added to this automatically.
#' @param result An \code{inla} object returned from previous calls of \link{inla}, \link{bru} or \link{lgcp}. This will be used as a starting point for further improvement of the approximate posterior.
#' @param E \link{inla} exposure parameter
#' @param control.compute INLA option, See \link{control.compute}
#' @param control.inla INLA option, See \link{control.inla}
#' @param ... Additional options passed on to \link{inla}
#' 
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 

bru.options = function(mesh = NULL, 
                       run = TRUE,
                       max.iter = 10,
                       offset = 0,
                       result = NULL, 
                       E = 1,
                       control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
                       control.inla = iinla.getOption("control.inla"),
                       ... )
{
  
  args <- as.list(environment())
  args$control.compute = NULL
  args$control.inla = NULL
  args$inla.options = list(...)
  args$inla.options$control.compute = control.compute
  args$inla.options$control.inla = control.inla
  
  args
}


#' Log Gaussian Cox process (LGCP) inference using INLA
#' 
#' This function performs inference on a LGCP observed via points residing possibly multiple dimensions. 
#' These dimensions are defined via the left hand side of the formula provided via the model parameter.
#' The left hand side determines the intensity function that is assumed to drive the LGCP. This may include
#' effects that lead to a thinning (filtering) of the point process. By default, the log intensity is assumed
#' to be a linear combination of the effects defined by the formula's RHS. More sofisticated models, e.g.
#' non-linear thinning, can be achieved by using the predictor argument. The latter requires multiple runs
#' of INLA for improving the required approximation of the predictor. In many applications
#' the LGCP is only observed through subsets of the dimensions the process is living in. For example, spatial
#' point realizations may only be known in sub-areas of the modeled space. These observed subsets of the LGCP
#' domain are called samplers and can be provided via the respective parameter. If samplers is NULL it is
#' assumed that all of the LGCP's dimensions have been observed completely. 
#' 
#'
#' @aliases lgcp
#' @export
#' @param components A formula describing the latent components
#' @param data A data frame or SpatialPoints[DataFrame] object
#' @param samplers A data frame or Spatial[Points/Lines/Polygons]DataFrame objects
#' @param domain Named list of domain definitions
#' @param ips Integration points (overrides \code{samplers})
#' @param formula If NULL, the linear combination implied by the \code{components} is used as a predictor for the point location intensity. If a (possibly non-linear) expression is provided the respective Taylor approximation is used as a predictor. Multiple runs if INLA are then required for a better approximation of the posterior.
#' @param E Single numeric used rescale all integration weights by a fixed factor 
#' @param options See \link{bru.options}
#' @return An \link{bru} object

lgcp = function(components,
                data,
                samplers = NULL,
                domain = NULL,
                ips = NULL,
                formula = . ~ .,
                E = 1,
                options = list()) {
  
  lik = like("cp", formula = formula, data = data, samplers = samplers, components = components, E = E, ips = ips, domain = domain)
  result = bru(components, lik, options = options)
  
}


# Summarize a LGCP object
#
# @aliases summary.lgcp
# @export
# @param object A result object obtained from a lgcp() run
# @param ... ignored arguments (S3 generic compatibility)
 
summary.lgcp = function(object, ...) {
  
  result = object
  
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


#' Summarize a \link{bru} object
#'
#' @aliases summary.bru
#' @export
#' @param object An object obtained from a \link{bru} or \link{lgcp} call
# 

summary.bru = function(object, ...) {
  
  cat("\n--- Likelihoods ----------------------------------------------------------------------------------\n\n")
  for ( k in 1:length(object$sppa$lhoods) ) {
    lh = object$sppa$lhoods[[k]]
    cat(sprintf("Name: '%s', family: '%s', data class: '%s', \t formula: '%s' \n", names(object$sppa$lhoods)[[k]], lh$family, class(lh$data),deparse(lh$formula)))
  }
  
  #rownames(df) = names(object$sppa$lhoods)
  #print(df)
  
  cat("\n--- Criteria -------------------------------------------------------------------------------------\n\n")
  cat(paste0("Watanabe-Akaike information criterion (WAIC): \t", sprintf("%1.3e", object$waic$waic),"\n"))
  cat(paste0("Deviance Information Criterion (DIC): \t\t", sprintf("%1.3e", object$dic$dic),"\n"))
  
  cat("\n--- Fixed effects -------------------------------------------------------------------------------- \n\n")
  
  if ( nrow(object$summary.fixed)>0 ) {
  fe = object$summary.fixed
  fe$kld=NULL
  fe$signif = sign(fe[,"0.025quant"]) == sign(fe[,"0.975quant"])
  print(fe)
  } else { cat("None.\n") }
  
  cat("\n--- Random effects ------------------------------------------------------------------------------- \n\n")
  for ( nm in names(object$summary.random) ){
    sm = object$summary.random[[nm]]
    cat(paste0(nm,": "))
    cat(paste0("mean = [", signif(range(sm$mean)[1])," : ",signif(range(sm$mean)[2]), "]"))
    cat(paste0(", quantiles = [", signif(range(sm[,c(4,6)])[1])," : ",signif(range(c(4,6))[2]), "]"))
    if (nm %in% names(object$model$mesh)) {
      cat(paste0(", area = ", signif(sum(diag(as.matrix(inla.mesh.fem(object$model$mesh[[nm]])$c0))))))
    }
    cat("\n")
  }
  if ( length(names(object$summary.random)) == 0 ) {cat("None.\n")}
  
  cat("\n--- All hyper parameters (internal representation) ----------------------------------------------- \n\n")
  # cat(paste0("  ", paste(rownames(object$summary.hyperpar), collapse = ", "), "\n"))
  print(object$summary.hyperpar)
  
  
  marginal.summary = function(m, name) {
    df = data.frame(param = name,
                    mean = inla.emarginal(identity, m))
    df$var = inla.emarginal(function(x) {(x-df$mean)^2}, m)
    df$sd = sqrt(df$var)
    df[c("lq","median","uq")] = inla.qmarginal(c(0.025, 0.5, 0.975), m)
    df
  }
  
  cat("\n")
  for (nm in names(object$sppa$model$effects)) {
    eff = object$sppa$model$effects[[nm]]
    if (eff$model == "spde2"){
      hyp = inla.spde.result(object, nm, eff$inla.spde)
      cat(sprintf("\n--- Field '%s' transformed hyper parameters ---\n", nm))
      df = rbind(marginal.summary(hyp$marginals.range.nominal$range.nominal.1, "range"), 
                 marginal.summary(hyp$marginals.log.range.nominal$range.nominal.1, "log.range"), 
                 marginal.summary(hyp$marginals.variance.nominal$variance.nominal.1, "variance"),
                 marginal.summary(hyp$marginals.log.variance.nominal$variance.nominal.1, "log.variance"))
      print(df)
    }
  }
  

}

#' Predictions based on bru
#' 
#' @aliases predict.bru
#' @export
#' @param object An object obtained by calling \link{bru})
#' @param data A data.frame or SpatialPointsDataFrame of covariates needed for the prediction
#' @param formula A formula determining which effects to predict and how to combine them
#' @param n.samples Integer setting the number of samples to draw in order to calculate the posterior statistics. The default is rather low but provides a quick approximate result.
#' @param ... ignored arguments (S3 generic compatibility)
#' 
#' @return Predicted values

predict.bru = function(object,
                       data = NULL,
                       formula = NULL,
                       n.samples = 100, ...)
{
  
  # Convert data into list, data.frame or a Spatial object if not provided as such
  if ( is.character(data) ) { data = as.list(setNames(data, data)) }
  else if ( inherits(data, "inla.mesh") ) { data = vertices(data) }
  
  vals = generate.bru(object, data, formula = formula, n.samples = n.samples)

  # Summarize
  
  if (is.data.frame(vals[[1]])) {
    vals.names = names(vals[[1]])
    covar = intersect(vals.names, names(data))
    estim = setdiff(vals.names, covar)
    smy = list()
    
    for ( nm in estim ) {
        smy[[nm]] = summarize(lapply(vals, function(v) v[[nm]]), x = vals[[1]][,covar,drop=FALSE])
    }
    vals = smy
    is.annot = sapply(names(vals), function(v) all(vals[[v]]$sd==0))
    annot = do.call(cbind, lapply(vals[is.annot], function(v) v[,1]))
    vals = vals[!is.annot]
    if ( !is.null(annot) ) {
      vals = lapply(vals, function(v) cbind(data.frame(annot), v))
    }
    
    
    if(length(vals)==1) vals = vals[[1]]
    
  } else {
    # if ( nrow(vals[[1]]) == nrow(data) ) { add.x = vals[[1]][,covar,drop=FALSE] } else { add.x = NULL }
    vals = summarize(vals, x = data)
  }

  if (!inherits(vals, "Spatial")) class(vals) = c("prediction",class(vals))
  vals
  
}

#' Sampling based on bru posteriors
#' 
#' @aliases generate.bru
#' @export
#' @param object An object obtained by calling \link{bru})
#' @param data A data.frame or SpatialPointsDataFrame of covariates needed for the prediction
#' @param formula A formula determining which effects to predict and how to combine them
#' @param n.samples Integer setting the number of samples to draw in order to calculate the posterior statistics. The default is rather low but provides a quick approximate result.
#' 
#' @return Predicted values

generate.bru = function(object,
                       data,
                       formula = NULL,
                       n.samples = 100)
{
  # Convert data into list, data.frame or a Spatial object if not provided as such
  if ( is.character(data) ) { data = as.list(setNames(data, data)) }
  else if ( inherits(data, "inla.mesh") ) { data = vertices(data) }

  # If data is provided as list, generate data automatically for each dimension stated in this list
  if ( class(data)[1] == "list" ) {
    lhs.names = names(data)
    add.pts = lapply(lhs.names, function(nm) { ipoints(object$sppa$lhoods$default$drange[[nm]], name = nm) })
    data = do.call(cprod, add.pts)
  }

  # Turn formula into an expression
  if ( is.null(formula) ) { formula = object$sppa$lhoods[["default"]]$formula }
  
  vals = evaluate.model(model = object$sppa$model, result = object, points = data, 
                        property = "sample", n = n.samples, predictor = formula)

}



# Monte Carlo method for estimating aposterior
#  
# @aliases montecarlo.posterior
# @export
# @param dfun A function returning a density for given x
# @param sfun A function providing samples from a posterior
# @param x Inital domain on which to perform the estimation. This will be adjusted as more samples are generated.
# @param samples An initial set of samples. Not required but will be used to estimate the inital domain \code{x} if \code{x} is \code{NULL}
# @param mcerr Monte Carlo error at which to stop the chain
# @param n Inital number of samples. This will be doubled for each iteration.
# @param discrete St this to \code{TRUE} if the density is only defined for integer \code{x}
# @param verbose Be verbose?

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
      # plot(x, lest, type = "l") ; lines(x, est, type = "l", col = "red")
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
#' @param data A list of samples, each either numeric or a \code{data.frame}
#' @param x A \code{data.frame} of data columns that should be added to the summary data frame
#' @param cbind.only If TRUE, only \code{cbind} the samples and return a matrix where each column is a sample
#' @return A \code{data.frame} or Spatial[Points/Pixels]DataFrame with summary statistics

summarize = function(data, x = NULL, cbind.only = FALSE) {
  if ( is.list(data) ) { data = do.call(cbind, data) }
  if ( cbind.only ) {
    smy = data.frame(data)
    colnames(smy) = paste0("sample.",1:ncol(smy))
  } else {
    smy = data.frame(
      apply(data, MARGIN = 1, mean, na.rm = TRUE),
      apply(data, MARGIN = 1, sd, na.rm = TRUE),
      t(apply(data, MARGIN = 1, quantile, prob = c(0.025, 0.5, 0.975), na.rm = TRUE)),
      apply(data, MARGIN = 1, min, na.rm = TRUE),
      apply(data, MARGIN = 1, max, na.rm = TRUE))
    colnames(smy) = c("mean", "sd", "q0.025", "median","q0.975", "smin", "smax")
    smy$cv = smy$sd/smy$mean
    smy$var = smy$sd^2
  }
  if ( !is.null(x) ) {
    if ( inherits(x, "Spatial") ) {
      if ( nrow( coordinates(x)) == nrow(smy) ) {
        if ( class(x) == "SpatialPoints" ) { smy = SpatialPointsDataFrame(x, data = smy) }
        else if ( class(x) == "SpatialPixels" ) { smy = SpatialPixelsDataFrame(x, data = smy) }
        else { x@data = cbind(x@data, smy) ; smy = x }
      }
    }
    else {
      if ( (nrow(smy) == nrow(x)) ) { smy = cbind(x, smy) }
    }
  }
  return(smy)
}




# Iterated INLA
# 
# This is a wrapper for iterated runs of \link{inla}. Before each run the \code{stackmaker} function is used to
# set up the \link{inla.stack} for the next iteration. For this purpose \code{stackmaker} is called given the
# \code{data} and \code{model} arguments. The \code{data} argument is the usual data provided to \link{inla}
# while \link{model} provides more information than just the usual inla formula. 
# 
# @aliases iinla
# @export
# @param data A data.frame
# @param model A \link{model} object
# @param stackmaker A function creating a stack from data and a model
# @param n Number of \link{inla} iterations
# @param iinla.verbose If TRUE, be verbose (use verbose=TRUE to make INLA verbose)
# @param ... Arguments passed on to \link{inla}
# @return An \link{inla} object


iinla = function(data, model, stackmaker, n = 10, result = NULL, 
                 family,
                 iinla.verbose = iinla.getOption("iinla.verbose"), 
                 offset = NULL, inla.options){
  
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
    
    # When running multiple times propagate theta
    if ( k>1 ) {
      inla.options[["control.mode"]] = list(restart = TRUE, theta = result$mode$theta)
    }
    
    # Verbose
    if ( iinla.verbose ) { cat(paste0("iinla() iteration"),k,"[ max:", n,"].") }
    
    # Return previous result if inla crashes, e.g. when connection to server is lost 
    if ( k > 1 ) { old.result = result } 
    result = NULL
    
    icall = expression(result <- tryCatch( do.call(inla, c(list(formula = update.formula(model$formula, y.inla ~ .),
                                                   data = c(inla.stack.mdata(stk), list.data(model)),
                                                   family = family,
                                                   control.predictor = list( A = inla.stack.A(stk), compute = TRUE),
                                                   E = inla.stack.data(stk)$e,
                                                   offset = inla.stack.data(stk)$bru.offset + offset),
                                                   inla.options)), 
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
    if ( n > 1 & k <= n) track[[k]] = cbind(effect = rownames(result$summary.fixed), iteration = k, result$summary.fixed)
    
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
