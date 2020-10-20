#' Generate samples from fitted bru models
#'
#' @description
#'
#' Generic function for sampling for fitted models. The function invokes particular methods
#' which depend on the class of the first argument.
#'
#' @name generate
#' @export
#' @family sample generators
#' @param object a fitted model.
#' @param ... additional arguments affecting the samples produced.
#' @return The form of the value returned by gg depends on the class of its argument. See the documentation of the particular methods for details of what is produced by that method.
#' @example inst/examples/generate.bru.R

generate <- function(object, ...) {
  UseMethod("generate")
}


#' @title Convenient model fitting using (iterated) INLA
#'
#' @description This method is a wrapper for \code{INLA::inla} and provides
#'   multiple enhancements.
#'
#'   \itemize{
#'     \item
#'       Easy usage of spatial covariates and automatic construction
#'       of inla projection matrices for (spatial) SPDE models. This feature is
#'       accessible via the \code{components} parameter. Practical examples on how
#'       to use spatial data by means of the components parameter can also be found
#'       by looking at the \link{lgcp} function's documentation.
#'     \item
#'       Constructing multiple likelihoods is straight forward. See \link{like} for more
#'       information on how to provide additional likelihoods to \code{bru} using
#'       the \code{...} parameter list.
#'     \item
#'       Support for non-linear predictors. See example below.
#'     \item
#'       Log Gaussian Cox process (LGCP) inference is
#'       available by using the \code{cp} family or (even easier) by using the
#'       \link{lgcp} function.
#'   }
#' @aliases bru
#' @export
#'
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'
#' @param components A \link{formula}-like specification of latent components.
#'   Also used to define a default linear additive predictor.  See
#'   \link{component} for details.
#' @param ... Likelihoods, each constructed by a calling \code{\link{like}}, or named
#'   parameters that can be passed to a single \code{\link{like}} call.
#' @param options A list of name and value pairs that are either interpretable
#'   by \link{bru.options} or valid inla parameters. See \link{bru.options} for
#'   caveats about some of the \code{control.*} inla options.
#'
#' @return bru returns an object of class "bru". A \code{bru} object inherits
#'   from \code{INLA::inla} (see the inla documentation for its properties) and
#'   adds additional information stored in the \code{sppa} field.
#'
#' @example inst/examples/bru.R
#'

bru <- function(components = y ~ Intercept,
                ...,
                options = list()) {
  requireINLA()

  # Update default options
  options <- do.call(bru.options, options)

  lhoods <- list(...)
  dot_is_lhood <- vapply(lhoods, function(lh) inherits(lh, "lhood"), TRUE)
  if (any(dot_is_lhood)) {
    if (!all(dot_is_lhood)) {
      stop("Cannot mix like() parameters with 'lhood' objects.")
    }
  } else {
    if (is.null(lhoods[["formula"]])) {
      lhoods[["formula"]] <- . ~ .
    }
    lhoods[["formula"]] <- auto_copy_response(lhoods[["formula"]], components)
    lhoods <- list(do.call(like, c(lhoods, list(options = options))))
  }

  if (length(lhoods) == 0) {
    stop("No response likelihood models provided.")
  }

  # Turn model components into internal bru model
  bru.model <- make.model(components, lhoods)
  
  
  # Set max iterations to 1 if all likelihood formulae are linear
  if (all(vapply(lhoods, function(lh) lh$linear, TRUE))) {
    options$max.iter <- 1
  }

  # Run iterated INLA
  if (options$run) {
    result <- iinla(
      model = bru.model,
      lhoods = lhoods,
      n = options$max.iter,
      offset = options$offset,
      result = options$result,
      inla.options = options$inla.options
    )
  } else {
    result <- list()
  }

  ## Create result object ##
  result$sppa$method <- "bru"
  result$sppa$lhoods <- lhoods
  result$sppa$model <- bru.model
  result$sppa$mesh <- options$mesh
  class(result) <- c("bru", class(result))
  return(result)
}


#' Likelihood construction for usage with \link{bru}
#'
#' @aliases like
#' @export
#'
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'
#' @param family A string identifying a valid \code{INLA::inla} likelihood family. The default is
#'   \code{gaussian} with identity link. In addition to the likelihoods provided
#'   by inla (see \code{inla.models()$likelihood}) inlabru supports fitting Cox
#'   processes via \code{family = "cp"}.
#'   As an alternative to \code{\link{bru}}, the \code{\link{lgcp}} function provides
#'   a convenient interface to fitting Cox processes.
#' @param formula a \link{formula} where the right hand side expression defines
#'   the predictor used in the optimization.
#' @param data Likelihood-specific data, as a \code{data.frame} or \code{SpatialPoints[DataFrame]}
#'   object.
#' @param components Components.
#' @param mesh An inla.mesh object.
#' @param E Exposure parameter for family = 'poisson' passed on to
#'   \code{INLA::inla}. Special case if family is 'cp': rescale all integration
#'   weights by E. Default taken from \code{options$E}.
#' @param Ntrials A vector containing the number of trials for the 'binomial'
#'  likelihood. Default value is rep(1, n.data).
#'  Default taken from \code{options$Ntrials}.
#' @param samplers Integration domain for 'cp' family.
#' @param ips Integration points for 'cp' family. Overrides \code{samplers}.
#' @param domain Named list of domain definitions.
#' @param include Character vector of component labels that are needed by the
#'   predictor expression; Default: NULL (include all components that are not
#'   explicitly excluded)
#' @param exclude Character vector of component labels that are not used by the
#'   predictor expression. The exclusion list is applied to the list
#'   as determined by the `include` parameter; Default: NULL (do not remove
#'   any components from the inclusion list)
#'   on the `include` parameter). 
#' @param options list of global options overriding \link{bru.options}
#'
#' @return A likelihood configuration which can be used to parameterize \link{bru}.
#'
#' @example inst/examples/like.R

like <- function(family, formula = . ~ ., data = NULL,
                 mesh = NULL, E = NULL, Ntrials = NULL,
                 samplers = NULL, ips = NULL, domain = NULL,
                 include = NULL, exclude = NULL,
                 options = list()) {
  options <- do.call(bru.options, options)

  # Some defaults
  inla.family <- family

  # Does the likelihood formula imply a linear predictor?
  linear <- as.character(formula)[length(as.character(formula))] == "."

  # If not linear, set predictor expression according to the formula's RHS
  if (!linear) {
    expr <- parse(text = as.character(formula)[length(as.character(formula))])
  } else {
    expr <- NULL
  }

  # Set response name
  response <- all.vars(update(formula, . ~ 0))
  if (response[1] == ".") {
    stop("Missing response variable names")
  }

  if (is.null(E)) {
    E <- options$E
  }
  if (is.null(Ntrials)) {
    Ntrials <- options$Ntrials
  }

  # More on special bru likelihoods
  if (family == "cp") {
    if (is.null(data)) {
      stop("You called like() with family='cp' but no 'data' argument was supplied.")
    }
    # if ( is.null(samplers) ) { stop("You called like() with family='cp' but no 'samplers' argument was supplied.") }

    if (is.null(ips)) {
      # TODO: split ipmaker into one domain extractor (or not; should force
      # explicit domain specification!), and one integration point constructor
#      message("ipmaker can't really work here unless it doesn't need the predictor information")
#      warning("ipmaker can't really work here unless it doesn't need the predictor information")
      ips <- ipmaker(
        samplers,
        domain = domain,
        dnames = response,
        data = NULL,
        model = NULL,
        int.args = options$int.args
      )
    }

    inla.family <- "poisson"
  }

  # Calculate data ranges
  drange <- lapply(names(data), function(nm) {
    if (is.numeric(data[[nm]])) {
      range(data[[nm]])
    } else {
      NULL
    }
  })
  names(drange) <- names(data)
  if (inherits(data, "Spatial")) drange[["coordinates"]] <- mesh


  # The likelihood object that will be returned

  lh <- list(
    family = family,
    formula = formula,
    data = data,
    E = E,
    Ntrials = Ntrials,
    samplers = samplers,
    linear = linear,
    expr = expr,
    response = response,
    inla.family = inla.family,
    ips = ips,
    domain = domain,
    drange = drange,
    include_components = include,
    exclude_components = exclude
  )

  class(lh) <- c("lhood", "list")

  # Return likelihood
  lh
}

single_stackmaker <- function(model, lhood, result) {
  if (lhood$family == "cp") {
    INLA::inla.stack(
      make.stack(
        data = lhood$data, model = model, expr = lhood$expr, y = 1,
        E = 0, result = result
      ),
      make.stack(
        data = lhood$ips, model = model, expr = lhood$expr, y = 0,
        E = lhood$E * lhood$ips$weight, offset = 0, result = result
      ),
      # Make sure components with zero derivative are kept:
      remove.unused = FALSE
    )
  } else {
    make.stack(
      data = lhood$data, model = model, expr = lhood$expr,
      y = as.data.frame(lhood$data)[, lhood$response],
      E = lhood$E, Ntrials = lhood$Ntrials, result = result
    )
  }
}

joint_stackmaker <- function(model, lhoods, result) {
  do.call(
    inla.stack.mjoin,
    c(
      lapply(
        lhoods,
        function(lh) {
          single_stackmaker(model, lh, result)
        }
      ),
      # Make sure components with zero derivative are kept:
      remove.unused = FALSE
    )
  )
}





#' Additional \link{bru} options
#'
#' @aliases bru.options
#' @export
#'
#' @param mesh An \code{inla.mesh} object for spatial models without SPDE
#'   components. Mostly used for successive spatial predictions.
#' @param run If TRUE, run inference. Otherwise only return configuration needed
#'   to run inference.
#' @param max.iter maximum number of inla iterations
#' @param offset the usual \code{INLA::inla} offset. If a nonlinear formula is
#'   used, the resulting Taylor approximation constant will be added to this
#'   automatically.
#' @param result An \code{inla} object returned from previous calls of
#'   \code{INLA::inla}, \link{bru} or \link{lgcp}. This will be used as a
#'   starting point for further improvement of the approximate posterior.
#' @param E \code{INLA::inla} 'poisson' likelihood exposure parameter
#' @param Ntrials \code{INLA::inla} 'binomial' likelihood parameter
#' @param control.compute INLA option, See \code{INLA::control.compute}
#' @param control.inla INLA option, See \code{INLA::control.inla}
#' @param control.fixed INLA option, See \code{INLA::control.fixed}. Warning:
#'   due to how inlabru currently constructs the \code{inla()}, the \code{mean},
#'   \code{prec}, \code{mean.intercept}, and \code{prec.intercept} will have no
#'   effect. Until a more elegant alterative has been implemented, use explicit
#'   \code{mean.linear} and \code{prec.linear} specifications in each
#'   \code{model="linear"} component instead.
#' @param int.args List of arguments passed all the way to the integration method
#' \code{ipoints} and \code{int.polygon} for 'cp' family models
#' @param ... Additional options passed on to \code{INLA::inla}
#'
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}> and Finn Lindgren \email{finn.lindgren@@gmail.com}
#'
#' @examples
#'
#' \donttest{
#'
#' # Generate default bru options
#' opts <- bru.options()
#'
#' # Print them:
#' opts
#' }
#'
bru.options <- function(mesh = NULL,
                        run = TRUE,
                        max.iter = 10,
                        offset = 0,
                        result = NULL,
                        E = 1,
                        Ntrials = 1,
                        control.compute = list(),
                        control.inla = list(),
                        control.fixed = list(),
                        int.args = list(method = "stable",nsub = NULL),
                        ...) {
  control.compute <-
    override_config_defaults(
      control.compute,
      iinla.getOption("control.compute")
    )
  control.inla <-
    override_config_defaults(
      control.inla,
      iinla.getOption("control.inla")
    )
  control.fixed <-
    override_config_defaults(
      control.fixed,
      iinla.getOption("control.fixed")
    )

  args <- as.list(environment())
  args$control.compute <- NULL
  args$control.inla <- NULL
  args$control.fixed <- NULL
  args$inla.options <- list(...)
  args$inla.options$control.compute <- control.compute
  args$inla.options$control.inla <- control.inla
  args$inla.options$control.fixed <- control.fixed

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
#' @examples
#'
#' \donttest{
#'
#' # Load the Gorilla data
#' data(gorillas, package = "inlabru")
#'
#' # Plot the Gorilla nests, the mesh and the survey boundary
#' ggplot() +
#'   gg(gorillas$mesh) +
#'   gg(gorillas$nests) +
#'   gg(gorillas$boundary) +
#'   coord_fixed()
#'
#' if (require("INLA", quietly = TRUE)) {
#'   # Define SPDE prior
#'   matern <- inla.spde2.pcmatern(gorillas$mesh,
#'     prior.sigma = c(0.1, 0.01),
#'     prior.range = c(5, 0.01)
#'   )
#'
#'   # Define domain of the LGCP as well as the model components (spatial SPDE effect and Intercept)
#'   cmp <- coordinates ~ mySmooth(map = coordinates, model = matern) + Intercept
#'
#'   # Fit the model
#'   fit <- lgcp(cmp, gorillas$nests, samplers = gorillas$boundary)
#'
#'   # Predict the spatial intensity surface
#'   lambda <- predict(fit, pixels(gorillas$mesh), ~ exp(mySmooth + Intercept))
#'
#'   # Plot the intensity
#'   ggplot() +
#'     gg(lambda) +
#'     gg(gorillas$mesh) +
#'     gg(gorillas$nests) +
#'     gg(gorillas$boundary) +
#'     coord_fixed()
#' }
#' }
#'
lgcp <- function(components,
                 data,
                 samplers = NULL,
                 domain = NULL,
                 ips = NULL,
                 formula = . ~ .,
                 E = NULL,
                 options = list()) {
#  # If formula response missing, copy from components
#  formula <- auto_copy_response(formula, components)
#  lik <- like(family = "cp",
#    formula = formula, data = data, samplers = samplers,
#    E = E, ips = ips, domain = domain,
#    options = options
#  )
  bru(components,
      family = "cp",
      formula = formula, data = data, samplers = samplers,
      E = E, ips = ips, domain = domain,
      options = options)
}


# Summarize a LGCP object
#
# @aliases summary.lgcp
# @export
# @param object A result object obtained from a lgcp() run
# @param ... ignored arguments (S3 generic compatibility)

summary.lgcp <- function(object, ...) {
  result <- object

  cat("### LGCP Summary #################################################################################\n\n")

  cat(paste0("Predictor: log(lambda) = ", as.character(result$model$expr), "\n"))

  cat("\n--- Points & Samplers ----\n\n")
  cat(paste0("Number of points: ", nrow(result$sppa$points)), "\n")
  if (inherits(result$sppa$points, "Spatial")) {
    cat(paste0("Coordinate names: ", paste0(coordnames(result$sppa$points), collapse = ", ")), "\n")
    cat(paste0("Coordinate system: ", proj4string(result$sppa$points), "\n"))
  }

  cat(paste0("Total integration weight: ", sum(result$ips$weight)), "\n")

  cat("\n--- Dimensions -----------\n\n")
  icfg <- result$iconfig
  invisible(lapply(names(icfg), function(nm) {
    cat(paste0(
      "  ", nm, " [", icfg[[nm]]$class, "]",
      ": ",
      "n = ", icfg[[nm]]$n.points,
      ", min = ", icfg[[nm]]$min,
      ", max = ", icfg[[nm]]$max,
      ", cardinality = ", signif(icfg[[nm]]$max - icfg[[nm]]$min),
      "\n"
    ))
  }))

  summary.bru(result)
}


#' Summary for a \link{bru} fit
#'
#' Takes a fitted bru object produced by bru() or lgcp() and creates various summaries from it.
#'
#' @aliases summary.bru
#' @export
#' @method summary bru
#' @param object An object obtained from a \link{bru} or \link{lgcp} call
#' @param ... ignored arguments
#' @example inst/examples/bru.R
#'

summary.bru <- function(object, ...) {
  cat("\n--- Likelihoods ----------------------------------------------------------------------------------\n\n")
  for (k in 1:length(object$sppa$lhoods)) {
    lh <- object$sppa$lhoods[[k]]
    cat(sprintf("Name: '%s', family: '%s', data class: '%s', \t formula: '%s' \n", names(object$sppa$lhoods)[[k]], lh$family, class(lh$data), deparse(lh$formula)))
  }

  # rownames(df) = names(object$sppa$lhoods)
  # print(df)

  cat("\n--- Criteria -------------------------------------------------------------------------------------\n\n")
  cat(paste0("Watanabe-Akaike information criterion (WAIC): \t", sprintf("%1.3e", object$waic$waic), "\n"))
  cat(paste0("Deviance Information Criterion (DIC): \t\t", sprintf("%1.3e", object$dic$dic), "\n"))

  cat("\n--- Fixed effects -------------------------------------------------------------------------------- \n\n")

  if (nrow(object$summary.fixed) > 0) {
    fe <- object$summary.fixed
    fe$kld <- NULL
    fe$signif <- sign(fe[, "0.025quant"]) == sign(fe[, "0.975quant"])
    print(fe)
  } else {
    cat("None.\n")
  }

  cat("\n--- Random effects ------------------------------------------------------------------------------- \n\n")
  for (nm in names(object$summary.random)) {
    sm <- object$summary.random[[nm]]
    cat(paste0(nm, " ranges: "))
    cat(paste0(
      "mean = [",
      signif(range(sm$mean)[1]), ", ",
      signif(range(sm$mean)[2]), "]"
    ))
    cat(paste0(
      ", sd = [",
      signif(range(sm$sd)[1]), ", ",
      signif(range(sm$sd)[2]), "]"
    ))
    cat(paste0(
      ", quantiles = [",
      signif(range(sm[, c(4, 6)])[1]), " : ",
      signif(range(sm[, c(4, 6)])[2]), "]"
    ))
    if (nm %in% names(object$model$mesh)) {
      cat(paste0(
        ", and area = ",
        signif(sum(Matrix::diag(INLA::inla.mesh.fem(object$model$mesh[[nm]])$c0)))
      ))
    }
    cat("\n")
  }
  if (length(names(object$summary.random)) == 0) {
    cat("None.\n")
  }

  cat("\n--- All hyper parameters (internal representation) ----------------------------------------------- \n\n")
  # cat(paste0("  ", paste(rownames(object$summary.hyperpar), collapse = ", "), "\n"))
  print(object$summary.hyperpar)


  marginal.summary <- function(m, name) {
    df <- data.frame(
      param = name,
      mean = INLA::inla.emarginal(identity, m)
    )
    df$var <- INLA::inla.emarginal(function(x) {
      (x - df$mean)^2
    }, m)
    df$sd <- sqrt(df$var)
    df[c("lq", "median", "uq")] <- INLA::inla.qmarginal(c(0.025, 0.5, 0.975), m)
    df
  }

  cat("\n")
  for (nm in names(object$sppa$model$effects)) {
    eff <- object$sppa$model$effects[[nm]]
    if (eff$model == "spde2") {
      hyp <- INLA::inla.spde.result(object, nm, eff$inla.spde)
      cat(sprintf("\n--- Field '%s' transformed hyper parameters ---\n", nm))
      df <- rbind(
        marginal.summary(hyp$marginals.range.nominal$range.nominal.1, "range"),
        marginal.summary(hyp$marginals.log.range.nominal$range.nominal.1, "log.range"),
        marginal.summary(hyp$marginals.variance.nominal$variance.nominal.1, "variance"),
        marginal.summary(hyp$marginals.log.variance.nominal$variance.nominal.1, "log.variance")
      )
      print(df)
    }
  }
}

#' Prediction from fitted bru model
#'
#' Takes a fitted \code{bru} object produced by the function \link{bru}() and produces predictions given
#' a new set of values for the model covariates or the original values used for the model fit. The
#' predictions can be based on any R expression that is valid given these values/covariates and the joint
#' posterior of the estimated random effects.
#'
#' Mean value predictions are accompanied by the standard errors, upper and lower 2.5% quantiles, the
#' median, variance, coefficient of variation as well as the variance and minimum and maximum sample
#' value drawn in course of estimating the statistics.
#'
#' Internally, this method calls \link{generate.bru} in order to draw samples from the model.
#'
#' @aliases predict.bru
#' @export
#' @param object An object obtained by calling \link{bru} or \link{lgcp}.
#' @param data A data.frame or SpatialPointsDataFrame of covariates needed for the prediction.
#' @param formula A formula determining which effects to predict and how to combine them.
#' @param n.samples Integer setting the number of samples to draw in order to calculate the posterior statistics. The default is rather low but provides a quick approximate result.
#' @param seed Random number generator seed passed on to \code{inla.posterior.sample}
#' @param num.threads Specification of desired number of threads for parallel
#' computations. Default NULL, leaves it up to INLA.
#' When seed != 0, overridden to "1:1"
#' @param ... ignored arguments (S3 generic compatibility).
#'
#' @return a data.frame or Spatial* object with predicted mean values and other summary statistics attached.
#' @example inst/examples/predict.bru.R

predict.bru <- function(object,
                        data = NULL,
                        formula = NULL,
                        n.samples = 100,
                        seed = 0L,
                        num.threads = NULL,
                        ...) {

  # Convert data into list, data.frame or a Spatial object if not provided as such
  if (is.character(data)) {
    data <- as.list(setNames(data, data))
  }
  else if (inherits(data, "inla.mesh")) {
    data <- vertices(data)
  }

  vals <- generate.bru(object,
                       data,
                       formula = formula,
                       n.samples = n.samples,
                       seed = seed,
                       num.threads = num.threads)

  # Summarize

  if (is.data.frame(vals[[1]])) {
    vals.names <- names(vals[[1]])
    covar <- intersect(vals.names, names(data))
    estim <- setdiff(vals.names, covar)
    smy <- list()

    for (nm in estim) {
      smy[[nm]] <- summarize(lapply(vals, function(v) v[[nm]]), x = vals[[1]][, covar, drop = FALSE])
    }
    vals <- smy
    is.annot <- vapply(names(vals), function(v) all(vals[[v]]$sd == 0), TRUE)
    annot <- do.call(cbind, lapply(vals[is.annot], function(v) v[, 1]))
    vals <- vals[!is.annot]
    if (!is.null(annot)) {
      vals <- lapply(vals, function(v) cbind(data.frame(annot), v))
    }


    if (length(vals) == 1) vals <- vals[[1]]
  } else {
    # if ( nrow(vals[[1]]) == nrow(data) ) { add.x = vals[[1]][,covar,drop=FALSE] } else { add.x = NULL }
    vals <- summarize(vals, x = data)
  }

  if (!inherits(vals, "Spatial")) class(vals) <- c("prediction", class(vals))
  vals
}

#' Sampling based on bru posteriors
#'
#' @description
#' Takes a fitted \code{bru} object produced by the function \link{bru}() and produces samples given
#' a new set of values for the model covariates or the original values used for the model fit. The
#' samples can be based on any R expression that is valid given these values/covariates and the joint
#' posterior of the estimated random effects.
#'
#' @aliases generate.bru
#' @export
#' @family sample generators
#' @param object A \code{bru} object obtained by calling \link{bru}.
#' @param data A data.frame or SpatialPointsDataFrame of covariates needed for sampling.
#' @param formula A formula determining which effects to sample from and how to combine them analytically.
#' @param n.samples Integer setting the number of samples to draw in order to calculate the posterior statistics.
#'                  The default is rather low but provides a quick approximate result.
#' @param seed Random number generator seed passed on to \code{inla.posterior.sample}
#' @param num.threads Specification of desired number of threads for parallel
#' computations. Default NULL, leaves it up to INLA.
#' When seed != 0, overridden to "1:1"
#' @param ... additional, unused arguments.
#'
#' @return List of generated samples
#' @seealso \link{predict.bru}
#' @example inst/examples/generate.bru.R
#' @rdname generate

generate.bru <- function(object,
                         data,
                         formula = NULL,
                         n.samples = 100,
                         seed = 0L,
                         num.threads = NULL,
                         ...) {
  # Convert data into list, data.frame or a Spatial object if not provided as such
  if (is.character(data)) {
    data <- as.list(setNames(data, data))
  }
  else if (inherits(data, "inla.mesh")) {
    data <- vertices(data)
  }

  # If data is provided as list, generate data automatically for each dimension stated in this list
  if (class(data)[1] == "list") {
    lhs.names <- names(data)
    add.pts <- lapply(lhs.names, function(nm) {
      ipoints(object$sppa$lhoods$default$drange[[nm]], name = nm)
    })
    data <- do.call(cprod, add.pts)
  }

  # Turn formula into an expression
  if (is.null(formula)) {
    formula <- object$sppa$lhoods[["default"]]$formula
  }

  vals <- evaluate.model(
    model = object$sppa$model, result = object, data = data,
    property = "sample", n = n.samples, predictor = formula, seed = seed
  )
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

montecarlo.posterior <- function(dfun, sfun, x = NULL, samples = NULL, mcerr = 0.01, n = 100, discrete = FALSE, verbose = FALSE) {
  xmaker <- function(hpd) {
    mid <- (hpd[2] + hpd[1]) / 2
    rg <- (hpd[2] - hpd[1]) / 2
    x <- seq(mid - 1.2 * rg, mid + 1.2 * rg, length.out = 256)
  }
  xmaker2 <- function(hpd) {
    x <- seq(hpd[1], hpd[2], length.out = 256)
  }

  inital.xmaker <- function(smp) {
    mid <- median(smp)
    rg <- (quantile(smp, 0.975) - quantile(smp, 0.25)) / 2
    x <- seq(mid - 3 * rg, mid + 3 * rg, length.out = 256)
  }

  # Inital samples
  if (is.null(samples)) {
    samples <- sfun(n)
  }

  # Inital HPD
  if (is.null(x)) {
    x <- inital.xmaker(as.vector(unlist(samples)))
  }

  # Round x if needed
  if (discrete) x <- unique(round(x))

  # First density estimate
  lest <- dfun(x, samples)


  converged <- FALSE
  while (!converged) {

    # Compute last HPD interval
    xnew <- xmaker2(INLA::inla.hpdmarginal(0.999, list(x = x, y = lest)))

    # Map last estimate to the HPD interval
    if (discrete) xnew <- unique(round(xnew))
    lest <- INLA::inla.dmarginal(xnew, list(x = x, y = lest))
    x <- xnew

    # Sample new density
    n <- 2 * n
    samples <- sfun(n)
    est <- dfun(x, samples)

    # Compute Monte Carlo error
    # err = sd(est/sum(est)-lest/sum(lest))
    err <- max(((est - lest) / max(lest))^2)

    # Plot new density estimate versus old one (debugging)
    if (verbose) {
      cat(paste0("hpd:", min(x), " ", max(x), ", err = ", err, ", n = ", n, "\n"))
      # plot(x, lest, type = "l") ; lines(x, est, type = "l", col = "red")
    }

    # Convergence?
    if (err < mcerr) {
      converged <- TRUE
    }
    else {
      lest <- 0.5 * (est + lest)
    }
  }

  marg <- list(x = x, y = est, samples = samples, mcerr = err)

  # Append some important statistics
  marg$quantiles <- INLA::inla.qmarginal(c(0.025, 0.5, 0.975), marg)
  marg$mean <- INLA::inla.emarginal(identity, marg)
  marg$sd <- sqrt(INLA::inla.emarginal(function(x) x^2, marg) - marg$mean^2)
  marg$cv <- marg$sd / marg$mean
  marg$mce <- err

  marg
}


# Summarize and annotate data
#
# @aliases summarize
# @export
# @param data A list of samples, each either numeric or a \code{data.frame}
# @param x A \code{data.frame} of data columns that should be added to the summary data frame
# @param cbind.only If TRUE, only \code{cbind} the samples and return a matrix where each column is a sample
# @return A \code{data.frame} or Spatial[Points/Pixels]DataFrame with summary statistics

summarize <- function(data, x = NULL, cbind.only = FALSE) {
  if (is.list(data)) {
    data <- do.call(cbind, data)
  }
  if (cbind.only) {
    smy <- data.frame(data)
    colnames(smy) <- paste0("sample.", 1:ncol(smy))
  } else {
    smy <- data.frame(
      apply(data, MARGIN = 1, mean, na.rm = TRUE),
      apply(data, MARGIN = 1, sd, na.rm = TRUE),
      t(apply(data, MARGIN = 1, quantile, prob = c(0.025, 0.5, 0.975), na.rm = TRUE)),
      apply(data, MARGIN = 1, min, na.rm = TRUE),
      apply(data, MARGIN = 1, max, na.rm = TRUE)
    )
    colnames(smy) <- c("mean", "sd", "q0.025", "median", "q0.975", "smin", "smax")
    smy$cv <- smy$sd / smy$mean
    smy$var <- smy$sd^2
  }
  if (!is.null(x)) {
    if (inherits(x, "Spatial")) {
      if (nrow(coordinates(x)) == nrow(smy)) {
        if (class(x) == "SpatialPoints") {
          smy <- SpatialPointsDataFrame(x, data = smy)
        }
        else if (class(x) == "SpatialPixels") {
          smy <- SpatialPixelsDataFrame(x, data = smy)
        }
        else {
          x@data <- cbind(x@data, smy)
          smy <- x
        }
      }
    }
    else {
      if ((nrow(smy) == nrow(x))) {
        smy <- cbind(x, smy)
      }
    }
  }
  return(smy)
}




#' Iterated INLA
#'
#' This is an internal wrapper for iterated runs of \code{INLA::inla}. Before each run the
#' \code{joint_stackmaker} function is used to set up the \code{INLA::inla.stack} for
#' the next iteration.
#'
#' @aliases iinla
#' @export
#' @param data A data.frame
#' @param model A \link{bru_model} object
#' @param lhoods A list of likelihood objects from \code{\link{like}}
#' @param n Number of \code{INLA::inla} iterations
#' @param result A previous inla result, to be used as starting point
#' @param iinla.verbose If TRUE, be verbose (use verbose=TRUE to make INLA verbose)
#' @param offset An additional predictor offset
#' @param inla.options A list of further arguments passed on to \code{INLA::inla}
#'   Some \code{control.*} arguments are partially overridden, according to
#'   internal inlabru method requirements.
#' @return An \code{INLA::inla} object
#' @keywords internal


iinla <- function(model, lhoods, n = 10, result = NULL,
                  iinla.verbose = inlabru:::iinla.getOption("iinla.verbose"),
                  offset = NULL, inla.options) {

  # # Default number of maximum iterations
  # if ( !is.null(model$expr) && is.null(n) ) { n = 10 } else { if (is.null(n)) {n = 1} }

  # Track variables?
  track <- list()

  # Set old result
  old.result <- result

  # Initialise required local options
  if (is.null(inla.options$control.mode)) {
    inla.options$control.mode <- list()
  }
  if (is.null(inla.options$control.predictor)) {
    inla.options$control.predictor <- list()
  }
  inla.options$control.predictor$compute <- TRUE

  # Extract the family of each likelihood
  family <- vapply(seq_along(lhoods), function(k) lhoods[[k]]$inla.family, "family")
  
  # Inital stack
  stk <- joint_stackmaker(model, lhoods, result)
  stk.data <- INLA::inla.stack.data(stk)
  inla.options$control.predictor$A <- INLA::inla.stack.A(stk)

  k <- 1
  interrupt <- FALSE

  while ((k <= n) & !interrupt) {

    # When running multiple times propagate theta
    if (k > 1) {
      inla.options[["control.mode"]]$restart <- TRUE
      inla.options[["control.mode"]]$theta <- result$mode$theta
    }

    # Verbose
    if (iinla.verbose) {
      cat(paste0("iinla() iteration"), k, "[ max:", n, "].")
    }

    # Return previous result if inla crashes, e.g. when connection to server is lost
    if (k > 1) {
      old.result <- result
    }
    result <- NULL

    inla.formula <- update.formula(model$formula, BRU.response ~ .)
    inla.data <-
      c(stk.data,
        do.call(c, c(lapply(model$effects,
                            function(xx) as.list(xx$env_extra)),
                     use.names = FALSE)))
#        list.data(model$formula))
    icall <- expression(
      result <- tryCatch(
        do.call(
          inla,
          c(
            list(
              formula = inla.formula,
              data = inla.data,
              family = family,
              E = stk.data[["BRU.E"]],
              Ntrials = stk.data[["BRU.Ntrials"]],
              offset = stk.data[["BRU.offset"]] + offset
            ),
            inla.options
          ),
          envir = environment(model$effects)
        ),
        error = warning
      )
    )
    eval(icall)

    if (is.character(result)) {
      stop(paste0("INLA returned message: ", result))
    }

    n.retry <- 0
    max.retry <- 10
    while ((is.null(result) | length(result) == 5) & (n.retry <= max.retry)) {
      msg("INLA crashed or returned NULL. Waiting for 60 seconds and trying again.")
      Sys.sleep(60)
      eval(icall)
      n.retry <- n.retry + 1
    }
    if ((is.null(result) | length(result) == 5)) {
      msg(sprintf("The computation failed %d times. Giving up and returning last successfully obtained result.", n.retry - 1))
      return(old.result)
    }

    if (iinla.verbose) {
      cat(" Done. ")
    }

    # Extract values tracked for estimating convergence
    if (n > 1 & k <= n) {
      # Note: The number of fixed effets may be zero, and strong
      # non-linearities that don't necessarily affect the fixed
      # effects may appear in the random effects, so we need to
      # track all of them.
      track[[k]] <- data.frame(effect = NULL, iteration = NULL, mean = NULL, sd = NULL)
      if (!is.null(result$summary.fixed) && (nrow(result$summary.fixed) > 0)) {
        track[[k]] <-
          data.frame(
            effect = rownames(result$summary.fixed),
            iteration = rep(k, nrow(result$summary.fixed)),
            mean = result$summary.fixed[, "mean"],
            sd = result$summary.fixed[, "sd"]
          )
      }
      if (!is.null(result$summary.random) &&
        (length(result$summary.random) > 0)) {
        ## This wastes memory temporarily.
        joined.random <- do.call(rbind, result$summary.random)
        track[[k]] <-
          rbind(
            track[[k]],
            data.frame(
              effect = rownames(joined.random),
              iteration = rep(k, nrow(joined.random)),
              mean = joined.random[, "mean"],
              sd = joined.random[, "sd"]
            )
          )
      }
    }

    # Update stack given current result
    if (n > 1 & k < n) {
      stk <- joint_stackmaker(model, lhoods, result)
      stk.data <- INLA::inla.stack.data(stk)
      inla.options$control.predictor$A <- INLA::inla.stack.A(stk)
    }

    # Stopping criterion
    if (k > 1) {
      ## TODO: Make max.dev a configurable option.
      max.dev <- 0.01
      dev <- abs(track[[k - 1]]$mean - track[[k]]$mean) / track[[k]]$sd
      ## do.call(c, lapply(by(do.call(rbind, track),
      ##                           as.factor(do.call(rbind, track)$effect),
      ##                           identity),
      ##                        function(X) { abs(X$mean[k-1] - X$mean[k])/X$sd[k] }))
      if (iinla.verbose) {
        cat(paste0("Max deviation from previous: ", signif(100 * max(dev), 3), "% of SD [stop if: <", 100 * max.dev, "%]\n"))
      }
      interrupt <- all(dev < max.dev)
      if (interrupt && iinla.verbose) {
        cat("Convergence criterion met, stopping INLA iteration.")
      }
    } else {
      cat("\n")
    }
    k <- k + 1
  }
  result$stack <- stk
  result$model <- model
  result$iinla$track <- do.call(rbind, track)
  class(result) <- c("iinla", "inla", "list")
  return(result)
}


auto.intercept <- function(components) {
  env <- environment(components)
  
  tm <- terms(components)
  # Check for -1/+0 and +/- Intercept/NULL
  # Required modification:
  # IVar ITerm AutoI Update
  #  T    T     1     -1
  #  T    T     0     -1
  #  T    F     1     -IVar-1
  #  T    F     0     -IVar-1
  #  F    T     1     Impossible
  #  F    T     0     Impossible
  #  F    F     1     +Intercept-1
  #  F    F     0     -1
  
  # Convert list(var1,var2) call to vector of variable names.
  # ([-1] removes "list"!)
  var_names <- as.character(attr(tm, "variables"))[-1]
  # Locate Intercept, if present
  inter_var <- grep("^Intercept[\\($]*", var_names)
  inter_term <- grep("^Intercept[\\($]*", attr(tm, "term.labels"))
  if (length(inter_var) > 0) {
    if (length(inter_term) > 0) {
      components <- update.formula(components, . ~ . + 1)
    } else {
      components <- update.formula(
        components,
        as.formula(paste0(
          ". ~ . - ",
          var_names[inter_var],
          " +1"
        ))
      )
    }
  } else if (attr(tm, "intercept")) {
    components <- update.formula(components, . ~ . + Intercept + 1)
  } else {
    components <- update.formula(components, . ~ . + 1)
  }
  environment(components) <- env
  components
}


auto_linear_formula <- function(formula, components) {
  if (as.character(formula)[
    length(as.character(formula))] != ".") {
    return(formula)
  }
  
  # No RHS; construct from components

  components <- auto.intercept(components)
  tm <- terms(components)
  
  var_names <- attr(tm, "term.labels")
  # Trim component options
  var_names <- gsub(pattern = "^([^\\(]*)\\(.*",
                    replacement = "\\1",
                    x = var_names)
  if (length(attr(tm, "offset")) > 0) {
    var_names <- c(var_names,
                   as.character(attr(tm, "variables"))[-1][attr(tm, "offset")])
  }
  formula <- update.formula(components,
                            paste0(". ~ ",
                                   paste0(var_names, collapse = " + ")))
  environment(formula) <- environment(components)
  formula
}


auto_copy_response <- function(formula, components) {
  if ((length(as.character(formula)) == 3) &&
      (as.character(formula)[2] != ".")) {
    # Already has response; do nothing
    return(formula)
  }
  env <- environment(formula)
  tm_cm <- terms(components)
  LHS <- as.formula(paste0(
    as.character(attr(tm_cm, "variables"))[-1][attr(tm_cm, "response")],
    " ~ ."))
  if (as.character(formula)[
    length(as.character(formula))] == ".") {
    # No formula RHS
    formula <- LHS
  } else {
    tm_fm <- terms(formula)
    if ((attr(tm_fm, "response") == 0) ||
        (as.character(attr(tm_fm, "variables"))[-1][attr(tm_fm, "response")] == "."))
    {
      formula <- update(formula, LHS)
    }
  }
  environment(formula) <- env
  formula
}


# Returns a formula's environemnt as a list. Removes all variable that are of type
# inla, function or formula. Also removes all variables that are not variables of the formula.
list.data <- function(formula) {

  # Formula environment as list
  elist <- as.list(environment(formula))

  # Remove previous inla results. For some reason these slow down the next INLA call.
  elist <- elist[unlist(lapply(elist, function(x) !inherits(x, "inla")))]

#  # Remove functions. This can cause problems as well.
#  elist <- elist[unlist(lapply(elist, function(x) !is.function(x)))]

#  # Remove formulae. This can cause problems as well.
#  elist <- elist[unlist(lapply(elist, function(x) !inherits(x, "formula")))]

  # The formula expression is too general for this to be reliable:
  #  # Keep only purse formula variables
  #  elist <- elist[names(elist) %in% all.vars(formula)]
  
  elist
}
