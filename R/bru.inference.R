#' Generate samples from fitted bru models
#'
#' @description
#'
#' Generic function for sampling for fitted models. The function invokes
#' particular methods which depend on the class of the first argument.
#'
#' @name generate
#' @export
#' @family sample generators
#' @param object a fitted model.
#' @param ... additional arguments affecting the samples produced.
#' @return The form of the value returned by gg depends on the class of its
#' argument. See the documentation of the particular methods for details of
#' what is produced by that method.
#' @example inst/examples/generate.bru.R

generate <- function(object, ...) {
  UseMethod("generate")
}

bru_check_object_bru <- function(object) {
  if (is.null(object[["bru_info"]])) {
    if (is.null(object[["sppa"]])) {
      stop("bru object contains neither current `bru_info` or old `sppa` information")
    }
    object[["bru_info"]] <- object[["sppa"]]
    object[["sppa"]] <- NULL
    old <- TRUE
  } else {
    old <- FALSE
  }
  object[["bru_info"]] <- bru_info_upgrade(object[["bru_info"]], old = old)
  object
}

bru_info_upgrade <- function(object, old = FALSE) {
  if (!is.list(object)) {
    stop("bru_info part of the object can't be converted to `bru_info`; not a list")
  }
  if (!inherits(object, "bru_info")) {
    old <- TRUE
    class(object) <- c("bru_info", "list")
  }
  ver <- getNamespaceVersion("inlabru")
  if (is.null(object[["inlabru_version"]])) {
    object[["inlabru_version"]] <- "0.0.0"
    old <- TRUE
  }
  old_ver <- object[["inlabru_version"]]
  if (old || (utils::compareVersion(ver, old_ver) > 0)) {
    warning(
      "Old bru_info object version ",
      old_ver,
      " detected. Attempting upgrade."
    )

    message("Detected bru_info version ", old_ver)
    if (utils::compareVersion("2.1.14.901", old_ver) > 0) {
      message("Upgrading bru_info to 2.1.14.901")
      # Ensure INLA_version stored
      if (is.null(object[["INLA_version"]])) {
        object[["INLA_version"]] <- "0.0.0"
      }
      # Check for component$mapper as a bru_mapper_multi
      for (k in seq_along(object[["model"]][["effects"]])) {
        cmp <- object[["model"]][["effects"]][[k]]
        cmp[["mapper"]] <-
          bru_mapper_multi(list(
            main = cmp$main$mapper,
            group = cmp$group$mapper,
            replicate = cmp$replicate$mapper
          ))
        object[["model"]][["effects"]][[k]] <- cmp
      }
      object[["inlabru_version"]] <- "2.1.14.901"
    }
    object[["inlabru_version"]] <- ver
  }
  object
}


get_component_type <- function(x, part, components) {
  if (is.null(x[["copy"]])) {
    x[[part]][["type"]]
  } else {
    components[[x[["copy"]]]][[part]][["type"]]
  }
}


#' Methods for bru_info objects
#' @export
#' @method summary bru_info
#' @param object Object to operate on
#' @param \dots Arguments passed on to other methods
#' @rdname bru_info
summary.bru_info <- function(object, ...) {
  result <- list(
    inlabru_version = object[["inlabru_version"]],
    INLA_version = object[["INLA_version"]],
    components =
      lapply(
        object[["model"]][["effects"]],
        function(x) {
          list(
            label = x[["label"]],
            copy_of = x[["copy"]],
            main_type = get_component_type(
              x, "main",
              object[["model"]][["effects"]]
            ),
            group_type = get_component_type(
              x, "group",
              object[["model"]][["effects"]]
            ),
            replicate_type = get_component_type(
              x, "replicate",
              object[["model"]][["effects"]]
            )
          )
        }
      ),
    lhoods =
      lapply(
        object[["lhoods"]],
        function(x) {
          list(
            family = x[["family"]],
            data_class = class(x[["data"]]),
            predictor = deparse(x[["formula"]])
          )
        }
      )
  )
  class(result) <- c("summary_bru_info", "list")
  result
}
#' Summary for bru_info objects
#' @export
#' @param x A `summary_bru_info` object to be printed
#' @rdname bru_info
print.summary_bru_info <- function(x, ...) {
  cat(paste0("inlabru version: ", x$inlabru_version, "\n"))
  cat(paste0("INLA version: ", x$INLA_version, "\n"))
  cat(paste0("Components:\n"))
  for (cmp in x$components) {
    if (!is.null(cmp$copy_of)) {
      cat(sprintf(
        "  %s: Copy of '%s' (types main='%s', group='%s', replicate='%s)\n",
        cmp$label,
        cmp$copy_of,
        cmp$main_type,
        cmp$group_type,
        cmp$replicate_type
      ))
    } else {
      cat(sprintf(
        "  %s: Model types main='%s', group='%s', replicate='%s'\n",
        cmp$label,
        cmp$main_type,
        cmp$group_type,
        cmp$replicate_type
      ))
    }
  }
  cat(paste0("Likelihoods:\n"))
  for (lh in x$lhoods) {
    cat(sprintf(
      "  Family: '%s'\n    Data class: %s\n    Predictor: %s\n",
      lh$family,
      paste0("'", lh$data_class, "'", collapse = ", "),
      if (length(lh$predictor) > 1) {
        paste0(
          "\n        ",
          paste0(
            lh$predictor,
            collapse = "\n        "
          )
        )
      } else {
        lh$predictor
      }
    ))
  }
  invisible(x)
}

#' @export
#' @rdname bru_info
bru_info <- function(...) {
  UseMethod("bru_info")
}

#' @export
#' @param method character; The type of estimation method used
#' @param inlabru_version character; inlabru package version. Default: NULL, for
#' automatically detecting the version
#' @param INLA_version character; INLA package version. Default: NULL, for
#' automatically detecting the version
#' @rdname bru_info
bru_info.character <- function(method,
                               ...,
                               inlabru_version = NULL,
                               INLA_version = NULL) {
  if (is.null(inlabru_version)) {
    inlabru_version <-
      tryCatch(
        {
          getNamespaceVersion("inlabru")
        },
        error = function(e) {
          "unknown"
        }
      )
  }
  if (is.null(INLA_version)) {
    INLA_version <-
      tryCatch(
        {
          getNamespaceVersion("INLA")
        },
        error = function(e) {
          "unknown"
        }
      )
  }
  object <- list(
    method = method,
    ...,
    inlabru_version = inlabru_version,
    INLA_version = INLA_version
  )
  class(object) <- c("bru_info", "list")
  object
}

#' @export
#' @rdname bru_info
bru_info.bru <- function(object, ...) {
  object <- bru_check_object_bru(object)
  object[["bru_info"]]
}


#' @title Convenient model fitting using (iterated) INLA
#'
#' @description This method is a wrapper for `INLA::inla` and provides
#'   multiple enhancements.
#'
#'   \itemize{
#'     \item
#'       Easy usage of spatial covariates and automatic construction
#'       of inla projection matrices for (spatial) SPDE models. This feature is
#'       accessible via the `components` parameter. Practical examples on how
#'       to use spatial data by means of the components parameter can also be found
#'       by looking at the [lgcp] function's documentation.
#'     \item
#'       Constructing multiple likelihoods is straight forward. See [like] for more
#'       information on how to provide additional likelihoods to `bru` using
#'       the `...` parameter list.
#'     \item
#'       Support for non-linear predictors. See example below.
#'     \item
#'       Log Gaussian Cox process (LGCP) inference is
#'       available by using the `cp` family or (even easier) by using the
#'       [lgcp] function.
#'   }
#' @aliases bru
#' @export
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'
#' @param components A `formula`-like specification of latent components.
#'   Also used to define a default linear additive predictor.  See
#'   [component()] for details.
#' @param ... Likelihoods, each constructed by a calling [like()], or named
#'   parameters that can be passed to a single [like()] call.
#' @param options A [bru_options] options object or a list of options passed
#' on to [bru_options()]
#'
#' @return bru returns an object of class "bru". A `bru` object inherits
#'   from `INLA::inla` (see the inla documentation for its properties) and
#'   adds additional information stored in the `bru_info` field.
#'
#' @example inst/examples/bru.R
#'

bru <- function(components = ~ Intercept(1),
                ...,
                options = list()) {
  stopifnot(bru_safe_inla())

  # Update default options
  options <- bru_call_options(options)

  lhoods <- list(...)
  dot_is_lhood <- vapply(
    lhoods,
    function(lh) inherits(lh, "bru_like"),
    TRUE
  )
  dot_is_lhood_list <- vapply(
    lhoods,
    function(lh) inherits(lh, "bru_like_list"),
    TRUE
  )
  if (any(dot_is_lhood | dot_is_lhood_list)) {
    if (!all(dot_is_lhood | dot_is_lhood_list)) {
      stop(paste0(
        "Cannot mix like() parameters with 'bru_like' and `bru_like_list` objects.",
        "\n  Check if the argument(s) ",
        paste0("'", names(lhoods)[!(dot_is_lhood | dot_is_lhood_list)], "'", collapse = ", "),
        " were meant to be given to a call to 'like()',\n  or in the 'options' list argument instead."
      ))
    }
  } else {
    if (is.null(lhoods[["formula"]])) {
      lhoods[["formula"]] <- . ~ .
    }
    if (inherits(components, "formula")) {
      lhoods[["formula"]] <- auto_response(
        lhoods[["formula"]],
        extract_response(components)
      )
    }
    lhoods <- list(do.call(like, c(lhoods, list(options = options))))
    dot_is_lhood <- TRUE
    dot_is_lhood_list <- FALSE
  }
  if (any(dot_is_lhood_list)) {
    lhoods <- like_list(c(
      lhoods[dot_is_lhood],
      do.call(c, lhoods[dot_is_lhood_list])
    ))
  } else {
    lhoods <- like_list(lhoods)
  }

  if (length(lhoods) == 0) {
    stop("No response likelihood models provided.")
  }

  # Turn input into a list of components (from existing list, or a special formula)
  components <- component_list(components)

  # Turn model components into internal bru model
  bru.model <- bru_model(components, lhoods)

  # Set max iterations to 1 if all likelihood formulae are linear
  if (all(vapply(lhoods, function(lh) lh$linear, TRUE))) {
    options$bru_max_iter <- 1
  }

  info <- bru_info(
    method = "bru",
    model = bru.model,
    lhoods = lhoods,
    options = options
  )

  # Run iterated INLA
  if (options$bru_run) {
    result <- iinla(
      model = info[["model"]],
      lhoods = info[["lhoods"]],
      options = info[["options"]]
    )
  } else {
    result <- list()
  }

  # Add bru information to the result
  result$bru_info <- info
  class(result) <- c("bru", class(result))
  return(result)
}


#' @details
#' * `bru_rerun` Continue the optimisation from a previously computed estimate.
#' @param result A previous estimation object of class `bru`
#'
#' @rdname bru
#' @export
bru_rerun <- function(result, options = list()) {
  stopifnot(inherits(result, "bru"))
  info <- result[["bru_info"]]
  info[["options"]] <- bru_call_options(
    bru_options(
      info[["options"]],
      as.bru_options(options)
    )
  )

  result <- iinla(
    model = info[["model"]],
    lhoods = info[["lhoods"]],
    initial = result,
    options = info[["options"]]
  )

  # Add bru information to the result
  result$bru_info <- info
  class(result) <- c("bru", class(result))
  return(result)
}


#' Parse inclusion of component labels in a predictor expression
#' @param thenames Set of labels to restrict
#' @param include Character vector of component labels that are needed by the
#'   predictor expression; Default: NULL (include all components that are not
#'   explicitly excluded)
#' @param exclude Character vector of component labels that are not used by the
#'   predictor expression. The exclusion list is applied to the list
#'   as determined by the `include` parameter; Default: NULL (do not remove
#'   any components from the inclusion list)
#' @keywords internal
parse_inclusion <- function(thenames, include = NULL, exclude = NULL) {
  if (!is.null(include)) {
    intersect(thenames, setdiff(include, exclude))
  } else {
    setdiff(thenames, exclude)
  }
}



#' Likelihood construction for usage with [bru()]
#'
#' @aliases like
#' @export
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'
#' @param formula a `formula` where the right hand side is a general R
#'   expression defines the predictor used in the model.
#' @param family A string identifying a valid `INLA::inla` likelihood family.
#' The default is
#'   `gaussian` with identity link. In addition to the likelihoods provided
#'   by inla (see `names(INLA::inla.models()$likelihood)`)
#'   inlabru supports fitting latent Gaussian Cox
#'   processes via `family = "cp"`.
#'   As an alternative to [bru()], the [lgcp()] function provides
#'   a convenient interface to fitting Cox processes.
#' @param data Likelihood-specific data, as a `data.frame` or
#' `SpatialPoints[DataFrame]`
#'   object.
#' @param mesh An inla.mesh object.
#' @param E Exposure parameter for family = 'poisson' passed on to
#'   `INLA::inla`. Special case if family is 'cp': rescale all integration
#'   weights by E. Default taken from `options$E`.
#' @param Ntrials A vector containing the number of trials for the 'binomial'
#'  likelihood. Default value is rep(1, n.data).
#'  Default taken from `options$Ntrials`.
#' @param samplers Integration domain for 'cp' family.
#' @param ips Integration points for 'cp' family. Overrides `samplers`.
#' @param domain Named list of domain definitions.
#' @param include Character vector of component labels that are needed by the
#'   predictor expression; Default: NULL (include all components that are not
#'   explicitly excluded)
#' @param exclude Character vector of component labels that are not used by the
#'   predictor expression. The exclusion list is applied to the list
#'   as determined by the `include` parameter; Default: NULL (do not remove
#'   any components from the inclusion list)
#' @param allow_latent logical. If `TRUE`, the latent state of each component is
#' directly available to the predictor expression, with a `_latent` suffix.
#' This also makes evaluator functions with suffix `_eval` available, taking
#' parameters `main`, `group`, and `replicate`, taking values for where to
#' evaluate the component effect that are different than those defined in the
#' component definition itself.
#' @param allow_combine logical; If `TRUE`, the predictor expression may
#' involve several rows of the input data to influence the same row.
#' (TODO: review what's needed to allow the result to also be of a different size)
#' @param options A [bru_options] options object or a list of options passed
#' on to [bru_options()]
#'
#' @return A likelihood configuration which can be used to parameterize [bru()].
#'
#' @example inst/examples/like.R

like <- function(formula = . ~ ., family = "gaussian", data = NULL,
                 mesh = NULL, E = NULL, Ntrials = NULL,
                 samplers = NULL, ips = NULL, domain = NULL,
                 include = NULL, exclude = NULL,
                 allow_latent = FALSE, allow_combine = FALSE,
                 options = list()) {
  options <- bru_call_options(options)

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

    if (is.null(ips)) {
      ips <- ipmaker(
        samplers = samplers,
        domain = domain,
        dnames = response,
        int.args = options[["bru_int_args"]]
      )
    }

    if (length(E) > 1) {
      warning("Exposure/effort parameter E should be a scalar for likelihood 'cp'.")
    }

    ips_is_Spatial <- inherits(ips, "Spatial")
    if (ips_is_Spatial) {
      ips_coordnames <- sp::coordnames(ips)
      ips_crs <- fm_sp_get_crs(ips)
    }
    data_is_Spatial <- inherits(data, "Spatial")
    if (data_is_Spatial) {
      data_coordnames <- sp::coordnames(data)
      data_crs <- fm_sp_get_crs(data)
      if (ips_is_Spatial) {
        new_coordnames <- data_coordnames[seq_len(min(
          length(ips_coordnames),
          length(data_coordnames)
        ))]
        new_coordnames[new_coordnames %in% ""] <-
          paste0(
            "BRU_dummy_coordinate_",
            seq_along(new_coordnames)
          )[new_coordnames %in% ""]
        ips_coordnames <- paste0(
          "BRU_dummy_coordinate_",
          seq_along(ips_coordnames)
        )
        data_coordnames <- paste(
          "BRU_dummy_coordinate_",
          seq_along(data_coordnames)
        )
        ips_coordnames[seq_along(new_coordnames)] <- new_coordnames
        data_coordnames[seq_along(new_coordnames)] <- new_coordnames
        coordnames(ips) <- ips_coordnames
        coordnames(data) <- data_coordnames

        # TODO: check that the crs info is the same
      }
    }
    data <- as.data.frame(data)
    ips <- as.data.frame(ips)
    dim_names <- intersect(names(data), names(ips))
    data <- rbind(
      cbind(data[dim_names], BRU_E = 0, BRU_response_cp = 1),
      cbind(ips[dim_names], BRU_E = E * ips[["weight"]], BRU_response_cp = 0)
    )
    if (ips_is_Spatial) {
      non_coordnames <- setdiff(names(data), data_coordnames)
      data <- sp::SpatialPointsDataFrame(
        coords = data[new_coordnames],
        data = data[non_coordnames],
        proj4string = data_crs,
        match.ID = FALSE
      )
    }

    response <- "BRU_response_cp"
    inla.family <- "poisson"
    E <- data[["BRU_E"]]
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
    domain = domain,
    drange = drange,
    include_components = include,
    exclude_components = exclude,
    allow_latent = allow_latent,
    allow_combine = allow_combine
  )

  class(lh) <- c("bru_like", "list")

  # Return likelihood
  lh
}


#' @details
#' * `like_list`: Combine a `bru_like` likelihoods
#' into a `bru_like_list` object
#' @param \dots For `like_list.bru_like`, one or more `bru_like` objects
#' @export
#' @rdname like
like_list <- function(...) {
  UseMethod("like_list")
}

#' @details
#' * `like_list.list`: Combine a list of `bru_like` likelihoods
#' into a `bru_like_list` object
#' @param object A list of `bru_like` objects
#' @param envir An optional environment for the new `bru_like_list` object
#' @export
#' @rdname like
like_list.list <- function(object, envir = NULL, ...) {
  if (is.null(envir)) {
    envir <- environment(object)
  }
  if (any(vapply(object, function(x) !inherits(x, "bru_like"), TRUE))) {
    stop("All list elements must be of class 'bru_like'.")
  }

  class(object) <- c("bru_like_list", "list")
  environment(object) <- envir
  object
}

#' @details
#' * `like_list.bru_like`: Combine several `bru_like` likelihoods
#' into a `bru_like_list` object
#' @export
#' @rdname like
like_list.bru_like <- function(..., envir = NULL) {
  like_list(list(...), envir = envir)
}

#' @export
#' @param x `bru_like_list` object from which to extract element(s)
#' @param i indices specifying elements to extract
#' @rdname like
`[.bru_like_list` <- function(x, i) {
  env <- environment(x)
  object <- NextMethod()
  class(object) <- c("bru_like_list", "list")
  environment(object) <- env
  object
}


bru_like_expr <- function(lhood, components) {
  if (!is.null(lhood[["expr"]])) {
    lhood[["expr"]]
  } else {
    included <-
      parse_inclusion(
        names(components),
        include = lhood[["include_components"]],
        exclude = lhood[["exclude_components"]]
      )
    parse(text = paste0(included, collapse = " + "))
  }
}






single_stackmaker <- function(model, lhood, state) {
  make.stack(
    data = lhood$data, model = model, expr = lhood$expr,
    y = lhood$data[[lhood$response]],
    E = lhood$E, Ntrials = lhood$Ntrials,
    state = state,
    include = lhood[["include_components"]],
    exclude = lhood[["exclude_components"]]
  )
}

joint_stackmaker <- function(model, lhoods, state) {
  do.call(
    inla.stack.mjoin,
    c(
      lapply(
        lhoods,
        function(lh) {
          single_stackmaker(model, lh, state)
        }
      ),
      # Make sure components with zero derivative are kept:
      remove.unused = FALSE
    )
  )
}






#' Log Gaussian Cox process (LGCP) inference using INLA
#'
#' This function performs inference on a LGCP observed via points residing
#' possibly multiple dimensions. These dimensions are defined via the left
#' hand side of the formula provided via the model parameter.
#' The left hand side determines the intensity function that is assumed to
#' drive the LGCP. This may include effects that lead to a thinning (filtering)
#' of the point process. By default, the log intensity is assumed to be a linear
#' combination of the effects defined by the formula's RHS. More sofisticated
#' models, e.g. non-linear thinning, can be achieved by using the predictor
#' argument. The latter requires multiple runs of INLA for improving the
#' required approximation of the predictor. In many applications the LGCP is
#' only observed through subsets of the dimensions the process is living in.
#' For example, spatial point realizations may only be known in sub-areas of
#' the modelled space. These observed subsets of the LGCP
#' domain are called samplers and can be provided via the respective parameter.
#' If samplers is NULL it is assumed that all of the LGCP's dimensions have
#' been observed completely.
#'
#'
#' @aliases lgcp
#' @export
#' @param components A formula describing the latent components
#' @param data A data frame or `SpatialPoints(DataFrame)` object
#' @param samplers A data frame or `Spatial[Points/Lines/Polygons]DataFrame` objects
#' @param domain Named list of domain definitions
#' @param ips Integration points (overrides `samplers`)
#' @param formula If NULL, the linear combination implied by the `components`
#' is used as a predictor for the point location intensity. If a (possibly
#' non-linear) expression is provided the respective Taylor approximation is
#' used as a predictor. Multiple runs if INLA are then required for a better
#' approximation of the posterior.
#' @param E Single numeric used rescale all integration weights by a fixed factor
#' @param options See [bru_options_set()]
#' @return An [bru()] object
#' @examples
#'
#' \donttest{
#' if (bru_safe_inla()) {
#'
#'   # Load the Gorilla data
#'   data(gorillas, package = "inlabru")
#'
#'   # Plot the Gorilla nests, the mesh and the survey boundary
#'   ggplot() +
#'     gg(gorillas$mesh) +
#'     gg(gorillas$nests) +
#'     gg(gorillas$boundary) +
#'     coord_fixed()
#'
#'   # Define SPDE prior
#'   matern <- INLA::inla.spde2.pcmatern(gorillas$mesh,
#'     prior.sigma = c(0.1, 0.01),
#'     prior.range = c(0.01, 0.01)
#'   )
#'
#'   # Define domain of the LGCP as well as the model components (spatial SPDE
#'   # effect and Intercept)
#'   cmp <- coordinates ~ mySmooth(map = coordinates, model = matern) + Intercept
#'
#'   # Fit the model (with int.strategy="eb" to make the example take less time)
#'   fit <- lgcp(cmp, gorillas$nests,
#'     samplers = gorillas$boundary,
#'     domain = list(coordinates = gorillas$mesh),
#'     options = list(control.inla = list(int.strategy = "eb"))
#'   )
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
  # If formula response missing, copy from components (for formula input)
  if (inherits(components, "formula")) {
    response <- extract_response(components)
    formula <- auto_response(formula, response)
  }
  lik <- like(
    family = "cp",
    formula = formula, data = data, samplers = samplers,
    E = E, ips = ips, domain = domain,
    options = options
  )
  bru(components, lik, options = options)
}




expand_to_dataframe <- function(x, data = NULL) {
  if (is.null(data)) {
    data <- data.frame(matrix(nrow = NROW(x), ncol = 0))
  }
  only_x <- setdiff(names(x), names(data))
  if (length(only_x) < length(names(x))) {
    x <- x[!(names(x) %in% names(data))]
  }
  if (inherits(x, "SpatialPixels") &&
    !inherits(x, "SpatialPixelsDataFrame")) {
    result <- sp::SpatialPixelsDataFrame(x, data = data)
  } else if (inherits(x, "SpatialGrid") &&
    !inherits(x, "SpatialGridDataFrame")) {
    result <- sp::SpatialGridDataFrame(x, data = data)
  } else if (inherits(x, "SpatialLines") &&
    !inherits(x, "SpatialLinesDataFrame")) {
    result <- sp::SpatialLinesDataFrame(x, data = data)
  } else if (inherits(x, "SpatialPolygons") &&
    !inherits(x, "SpatialPolygonsDataFrame")) {
    result <- sp::SpatialPolygonsDataFrame(x, data = data)
  } else if (inherits(x, "SpatialPoints") &&
    !inherits(x, "SpatialPointsDataFrame")) {
    # Other classes inherit from SpatialPoints, so need to be handled first
    result <- sp::SpatialPointsDataFrame(x, data = data)
  } else if (inherits(x, "Spatial")) {
    result <- sp::cbind.Spatial(x, data)
  } else {
    result <- cbind(x, data)
  }
  result
}


#
#' Prediction from fitted bru model
#'
#' Takes a fitted `bru` object produced by the function [bru()] and produces
#' predictions given a new set of values for the model covariates or the
#' original values used for the model fit. The predictions can be based on any
#' R expression that is valid given these values/covariates and the joint
#' posterior of the estimated random effects.
#'
#' Mean value predictions are accompanied by the standard errors, upper and
#' lower 2.5% quantiles, the
#' median, variance, coefficient of variation as well as the variance and
#' minimum and maximum sample
#' value drawn in course of estimating the statistics.
#'
#' Internally, this method calls [generate.bru()] in order to draw samples from
#' the model.
#'
#' @aliases predict.bru
#' @export
#' @param object An object obtained by calling [bru()] or [lgcp()].
#' @param data A data.frame or SpatialPointsDataFrame of covariates needed for
#' the prediction.
#' @param formula A formula defining an R expression to evaluate for each generated
#' sample. If `NULL`, the latent and hyperparameter states are generated
#' as named list elements. In addition to the component names (that give the effect
#' of each component evaluated for the input data), the suffix `_latent` can be used
#' to directly access the latent state for a component, and the suffix `_eval`
#' can be used to evaluate a component at other input values than the expressions
#' defined in the component definition itself, e.g. `field_eval(cbind(x,y))` for a
#' component defined with `field(coordinates, ...)`. For "iid" models with
#' `mapper = bru_mapper_index(n)`, `rnorm()` is used to generate new realisations for
#' indices greater than `n`.
#' @param n.samples Integer setting the number of samples to draw in order to
#' calculate the posterior statistics. The default is rather low but provides
#' a quick approximate result.
#' @param seed Random number generator seed passed on to `inla.posterior.sample`
#' @param num.threads Specification of desired number of threads for parallel
#' computations. Default NULL, leaves it up to INLA.
#' When seed != 0, overridden to "1:1"
#' @param include Character vector of component labels that are needed by the
#'   predictor expression; Default: NULL (include all components that are not
#'   explicitly excluded)
#' @param exclude Character vector of component labels that are not used by the
#'   predictor expression. The exclusion list is applied to the list
#'   as determined by the `include` parameter; Default: NULL (do not remove
#'   any components from the inclusion list)
#' @param drop logical; If `keep=FALSE`, `data` is a `Spatial*DataFrame`, and the
#' prediciton summary has the same number of rows as `data`, then the output is
#' a `Spatial*DataFrame` object. Default `FALSE`.
#' @param \dots Additional arguments passed on to `inla.posterior.sample`
#'
#' @return a data.frame or Spatial* object with predicted mean values and other
#' summary statistics attached.
#' @example inst/examples/predict.bru.R

predict.bru <- function(object,
                        data = NULL,
                        formula = NULL,
                        n.samples = 100,
                        seed = 0L,
                        num.threads = NULL,
                        include = NULL,
                        exclude = NULL,
                        drop = FALSE,
                        ...) {
  object <- bru_check_object_bru(object)

  # Convert data into list, data.frame or a Spatial object if not provided as such
  if (is.character(data)) {
    data <- as.list(setNames(data, data))
  }
  else if (inherits(data, "inla.mesh")) {
    data <- vertices(data)
  }
  else if (inherits(data, "formula")) {
    stop("Formula supplied as data to predict.bru(). Please check your argument order/names.")
  }

  vals <- generate.bru(object,
    data = data,
    formula = formula,
    n.samples = n.samples,
    seed = seed,
    num.threads = num.threads,
    include = include,
    exclude = exclude,
    ...
  )

  # Summarise

  data <- expand_to_dataframe(data)
  if (is.data.frame(vals[[1]])) {
    vals.names <- names(vals[[1]])
    covar <- intersect(vals.names, names(data))
    estim <- setdiff(vals.names, covar)
    smy <- list()

    for (nm in estim) {
      smy[[nm]] <-
        bru_summarise(
          lapply(
            vals,
            function(v) v[[nm]]
          ),
          x = vals[[1]][, covar, drop = FALSE]
        )
    }
    is.annot <- vapply(names(smy), function(v) all(smy[[v]]$sd == 0), TRUE)
    annot <- do.call(cbind, lapply(smy[is.annot], function(v) v[, 1]))
    smy <- smy[!is.annot]
    if (!is.null(annot)) {
      smy <- lapply(smy, function(v) cbind(data.frame(annot), v))
    }

    if (!drop) {
      smy <- lapply(
        smy,
        function(tmp) {
          if (NROW(data) == NROW(tmp)) {
            expand_to_dataframe(data, tmp)
          } else {
            tmp
          }
        }
      )
    }

    if (length(smy) == 1) smy <- smy[[1]]
  } else if (is.list(vals[[1]])) {
    vals.names <- names(vals[[1]])
    if (any(vals.names == "")) {
      warning("Some generated list elements are unnamed")
    }
    smy <- list()
    for (nm in vals.names) {
      tmp <-
        bru_summarise(
          lapply(
            vals,
            function(v) v[[nm]]
          )
        )
      if (!drop &&
        (NROW(data) == NROW(tmp))) {
        smy[[nm]] <- expand_to_dataframe(data, tmp)
      } else {
        smy[[nm]] <- tmp
      }
    }
  } else {
    tmp <- bru_summarise(vals)
    if (!drop &&
      (NROW(data) == NROW(tmp))) {
      smy <- expand_to_dataframe(data, tmp)
    } else {
      smy <- tmp
    }
  }

  if (!inherits(smy, "Spatial")) class(smy) <- c("prediction", class(smy))
  smy
}

#' Sampling based on bru posteriors
#'
#' @description
#' Takes a fitted `bru` object produced by the function [bru()] and produces
#' samples given a new set of values for the model covariates or the original
#' values used for the model fit. The samples can be based on any R expression
#' that is valid given these values/covariates and the joint
#' posterior of the estimated random effects.
#'
#' @aliases generate.bru
#' @export
#' @family sample generators
#' @param object A `bru` object obtained by calling [bru()].
#' @param data A data.frame or SpatialPointsDataFrame of covariates needed for
#' sampling.
#' @param formula A formula defining an R expression to evaluate for each generated
#' sample. If `NULL`, the latent and hyperparameter states are returned
#' as named list elements. In addition to the component names (that give the effect
#' of each component evaluated for the input data), the suffix `_latent` can be used
#' to directly access the latent state for a component, and the suffix `_eval`
#' can be used to evaluate a component at other input values than the expressions
#' defined in the component definition itself, e.g. `field_eval(cbind(x,y))` for a
#' component defined with `field(coordinates, ...)`. For "iid" models with
#' `mapper = bru_mapper_index(n)`, `rnorm()` is used to generate new realisations for
#' indices greater than `n`.
#' @param n.samples Integer setting the number of samples to draw in order to
#' calculate the posterior statistics.
#' The default, 100, is rather low but provides a quick approximate result.
#' @param seed Random number generator seed passed on to `INLA::inla.posterior.sample`
#' @param num.threads Specification of desired number of threads for parallel
#' computations. Default NULL, leaves it up to INLA.
#' When seed != 0, overridden to "1:1"
#' @param include Character vector of component labels that are needed by the
#'   predictor expression; Default: NULL (include all components that are not
#'   explicitly excluded)
#' @param exclude Character vector of component labels that are not used by the
#'   predictor expression. The exclusion list is applied to the list
#'   as determined by the `include` parameter; Default: NULL (do not remove
#'   any components from the inclusion list)
#' @param ... additional, unused arguments.
#'
#' @return List of generated samples
#' @seealso [predict.bru]
#' @example inst/examples/generate.bru.R
#' @rdname generate

generate.bru <- function(object,
                         data = NULL,
                         formula = NULL,
                         n.samples = 100,
                         seed = 0L,
                         num.threads = NULL,
                         include = NULL,
                         exclude = NULL,
                         ...) {
  object <- bru_check_object_bru(object)

  # Convert data into list, data.frame or a Spatial object if not provided as such
  if (is.character(data)) {
    data <- as.list(setNames(data, data))
  }
  else if (inherits(data, "inla.mesh")) {
    data <- vertices(data)
  }
  else if (inherits(data, "formula")) {
    stop("Formula supplied as data to generate.bru(). Please check your argument order/names.")
  }

  # If data is provided as list, generate data automatically for each dimension
  # stated in this list
  if (class(data)[1] == "list") {
    # Todo: check if this feature works at all.
    # TODO: add method ipoints.list to handle this;
    # ipoints(list(coordinates=mesh, etc)) and remove this implicit code from
    # generate()
    warning(paste0(
      "Attempting to convert data list into gridded data.\n",
      "This probably doesn't work.\n",
      "Please contact the package developers if you use this feature."
    ))
    lhs.names <- names(data)
    add.pts <- lapply(lhs.names, function(nm) {
      ipoints(object$bru_info$lhoods$default$drange[[nm]], name = nm)
    })
    data <- do.call(cprod, add.pts)
  }

  state <- evaluate_state(
    object$bru_info$model,
    result = object,
    property = "sample",
    n = n.samples,
    seed = seed,
    num.threads = num.threads,
    ...
  )
  if (is.null(formula)) {
    state
  } else {
    # TODO: clarify the output format, and use the format parameter
    vals <- evaluate_model(
      model = object$bru_info$model,
      state = state,
      data = data,
      predictor = formula,
      include = include,
      exclude = exclude
    )
    vals
  }
}



# Monte Carlo method for estimating aposterior
#
# @aliases montecarlo.posterior
# @export
# @param dfun A function returning a density for given x
# @param sfun A function providing samples from a posterior
# @param x Inital domain on which to perform the estimation. This will be
# adjusted as more samples are generated.
# @param samples An initial set of samples. Not required but will be used to
# estimate the inital domain \code{x} if \code{x} is \code{NULL}
# @param mcerr Monte Carlo error at which to stop the chain
# @param n Inital number of samples. This will be doubled for each iteration.
# @param discrete Set this to \code{TRUE} if the density is only defined for integer \code{x}
# @param verbose Be verbose?

montecarlo.posterior <- function(dfun, sfun, x = NULL, samples = NULL,
                                 mcerr = 0.01, n = 100, discrete = FALSE, verbose = FALSE) {
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
      lest <- (est + lest) / 2
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


# Summarise and annotate data
#
# @export
# @param data A list of samples, each either numeric or a \code{data.frame}
# @param x A \code{data.frame} of data columns that should be added to the summary data frame
# @param cbind.only If TRUE, only \code{cbind} the samples and return a matrix where each column is a sample
# @return A \code{data.frame} or Spatial[Points/Pixels]DataFrame with summary statistics

bru_summarise <- function(data, x = NULL, cbind.only = FALSE) {
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
    smy <- expand_to_dataframe(x, smy)
  }
  return(smy)
}





lin_predictor <- function(lin, state) {
  do.call(
    c,
    lapply(
      lin,
      function(x) {
        as.vector(
          x$offset + Matrix::rowSums(do.call(
            cbind,
            lapply(names(x$A), function(xx) {
              x[["A"]][[xx]] %*% state[[xx]]
            })
          ))
        )
      }
    )
  )
}
nonlin_predictor <- function(model, lhoods, state, A) {
  do.call(
    c,
    lapply(
      lhoods,
      function(lh) {
        as.vector(
          evaluate_model(
            model = model,
            data = lh[["data"]],
            state = list(state),
            A = A,
            predictor = bru_like_expr(lh, model[["effects"]]),
            format = "matrix"
          )
        )
      }
    )
  )
}
scale_state <- function(state0, state1, scaling_factor) {
  new_state <- lapply(
    names(state0),
    function(idx) {
      state0[[idx]] + (state1[[idx]] - state0[[idx]]) * scaling_factor
    }
  )
  names(new_state) <- names(state0)
  new_state
}

line_search_optimisation_target <- function(x, param) {
  (x - 1)^2 * param[1] + 2 * (x - 1) * x^2 * param[2] + x^4 * param[3]
}



bru_line_search <- function(model,
                            lhoods,
                            lin,
                            state0,
                            state,
                            A,
                            options) {
  if (length(options$bru_method$search) == 0) {
    return(
      list(
        active = FALSE,
        step_scaling = 1,
        state = state
      )
    )
  }

  fact <- options$bru_method$factor
  pred_norm <- function(delta) {
    if (any(!is.finite(delta))) {
      Inf
    } else {
      sum(delta^2)^0.5
    }
  }

  state1 <- state
  lin_pred0 <- lin_predictor(lin, state0)
  lin_pred1 <- lin_predictor(lin, state1)
  nonlin_pred <- nonlin_predictor(model, lhoods, state1, A)
  step_scaling <- 1

  norm01 <- pred_norm(lin_pred1 - lin_pred0)
  norm1 <- pred_norm(nonlin_pred - lin_pred1)

  do_finite <-
    any(c("all", "finite") %in% options$bru_method$search)
  do_contract <-
    any(c("all", "contract") %in% options$bru_method$search)
  do_expand <-
    any(c("all", "expand") %in% options$bru_method$search)
  do_optimise <-
    any(c("all", "optimise") %in% options$bru_method$search)

  finite_active <- 0
  contract_active <- 0
  expand_active <- 0

  if (do_finite || do_contract) {
    nonfin <- any(!is.finite(nonlin_pred))
    norm0 <- pred_norm(nonlin_pred - lin_pred0)

    while ((do_finite && nonfin) ||
      (do_contract && (norm0 > norm01 * fact))) {
      if (do_finite && nonfin) {
        finite_active <- finite_active - 1
      } else {
        contract_active <- contract_active - 1
      }
      step_scaling <- step_scaling / fact
      state <- scale_state(state0, state1, step_scaling)
      nonlin_pred <- nonlin_predictor(model, lhoods, state, A)
      nonfin <- any(!is.finite(nonlin_pred))
      norm0 <- pred_norm(nonlin_pred - lin_pred0)

      bru_log_message(
        paste0(
          "iinla: Step rescaling: ",
          signif(100 * step_scaling, 3),
          "%, Contract"
        ),
        verbose = options$bru_verbose,
        verbose_store = options$bru_verbose_store,
        verbosity = 3
      )

      if (step_scaling < 2^-10) {
        break
      }
    }
  }

  if (do_expand &&
    finite_active == 0 &&
    contract_active == 0) {
    overstep <- FALSE
    norm0 <- pred_norm(nonlin_pred - lin_pred0)
    norm1 <- pred_norm(nonlin_pred - lin_pred1)
    while (!overstep && (norm0 < norm01)) {
      expand_active <- expand_active + 1
      step_scaling <- step_scaling * fact
      state <- scale_state(state0, state1, step_scaling)

      nonlin_pred <- nonlin_predictor(model, lhoods, state, A)
      norm1_prev <- norm1
      norm0 <- pred_norm(nonlin_pred - lin_pred0)
      norm1 <- pred_norm(nonlin_pred - lin_pred1)
      overstep <-
        (norm1 > norm1_prev) || !is.finite(norm1)

      bru_log_message(
        paste0(
          "iinla: Step rescaling: ",
          signif(100 * step_scaling, 3),
          "%, Expand"
        ),
        verbose = options$bru_verbose,
        verbose_store = options$bru_verbose_store,
        verbosity = 3
      )

      if (step_scaling > 2^10) {
        break
      }
    }
    if (((step_scaling > 1) && overstep) || !is.finite(norm1)) {
      expand_active <- expand_active - 1
      step_scaling <- step_scaling / fact
      state <- scale_state(state0, state1, step_scaling)
      nonlin_pred <- nonlin_predictor(model, lhoods, state, A)
      norm1 <- pred_norm(nonlin_pred - lin_pred1)

      bru_log_message(
        paste0(
          "iinla: Step rescaling: ",
          signif(100 * step_scaling, 3),
          "%, Overstep"
        ),
        verbose = options$bru_verbose,
        verbose_store = options$bru_verbose_store,
        verbosity = 3
      )
    }
  }

  if (do_optimise) {
    lin_pred <- lin_predictor(lin, state)
    delta_lin <- (lin_pred1 - lin_pred0)
    delta_nonlin <- (nonlin_pred - lin_pred) / step_scaling^2
    alpha <-
      optimise(
        line_search_optimisation_target,
        step_scaling * c(1 / fact, fact),
        param = c(
          sum(delta_lin^2),
          sum(delta_lin * delta_nonlin),
          sum(delta_nonlin^2)
        )
      )

    step_scaling_opt <- alpha$minimum
    state_opt <- scale_state(state0, state1, step_scaling_opt)
    nonlin_pred_opt <- nonlin_predictor(model, lhoods, state_opt, A)
    norm1_opt <- pred_norm(nonlin_pred_opt - lin_pred1)

    if (norm1_opt < norm1) {
      step_scaling <- step_scaling_opt
      state <- state_opt
      nonlin_pred <- nonlin_pred_opt
      norm1 <- norm1_opt

      bru_log_message(
        paste0(
          "iinla: Step rescaling: ",
          signif(100 * step_scaling, 4),
          "%, Optimisation"
        ),
        verbose = options$bru_verbose,
        verbose_store = options$bru_verbose_store,
        verbosity = 3
      )
    }
  }

  active <-
    (finite_active + contract_active + expand_active != 0) ||
      (abs(step_scaling - 1) > 0.001)

  if (active && (options$bru_verbose <= 2)) {
    bru_log_message(
      paste0(
        "iinla: Step rescaling: ",
        signif(100 * step_scaling, 3),
        "%"
      ),
      verbose = options$bru_verbose,
      verbose_store = options$bru_verbose_store,
      verbosity = 2
    )
  }

  list(
    active = active,
    step_scaling = step_scaling,
    state = state
  )
}


latent_names <- function(state) {
  nm <- lapply(
    names(state),
    function(x) {
      if (length(state[[x]]) == 1) {
        x
      } else {
        paste0(x, ".", seq_len(length(state[[x]])))
      }
    }
  )
  names(nm) <- names(state)
  nm
}
tidy_tracker <- function(state, value_name = "value") {
  df <- data.frame(
    effect = latent_names(state),
    value = unlist(state)
  )
  names(df) <- c("effect", value_name)
  df
}


#' Iterated INLA
#'
#' This is an internal wrapper for iterated runs of `INLA::inla`.
#' For nonlinear models, a linearisation is done with
#' `bru_compute_linearisation`, with a line search method between each
#' iteration. The `INLA::inla.stack` information is setup by [bru_make_stack()].
#'
#' @aliases iinla
#' @export
#' @param data A data.frame
#' @param model A [bru_model] object
#' @param lhoods A list of likelihood objects from [like()]
#' @param initial A previous `bru` result or a list of named latent variable
#' initial states (missing elements are set to zero), to be used as starting
#' point, or `NULL`. If non-null, overrides `options$bru_initial`
#' @param options A `bru_options` object.
#' @return An `iinla` object that inherits from `INLA::inla`, with an
#' added field `bru_iinla` with elements
#' \describe{
#' \item{log}{The diagnostic log messages produced by the run}
#' \item{states}{The list of linearisation points, one for each inla run}
#' \item{inla_stack}{The `inla.stack` object from the final inla run}
#' \item{track}{A list of convergence tracking vectors}
#' }
#' If an inla run is aborted by an error, the returned object also contains
#' an element `error` with the error object.
#' @keywords internal


iinla <- function(model, lhoods, initial = NULL, options) {
  inla.options <- bru_options_inla(options)

  initial_log_length <- length(bru_log_get())
  # Local utility method for collecting information object:
  collect_misc_info <- function(...) {
    list(
      log = if (initial_log_length <= 0) {
        bru_log_get()
      } else {
        bru_log_get()[-seq_len(initial_log_length)]
      },
      states = states,
      inla_stack = stk,
      track = if (is.null(original_track) ||
        setequal(names(original_track), names(track[[1]]))) {
        do.call(rbind, c(list(original_track), track))
      } else {
        track <- do.call(rbind, track)
        original_names <- names(original_track)
        new_names <- names(track)
        for (nn in setdiff(new_names, original_names)) {
          original_track[[nn]] <- NA
        }
        for (nn in setdiff(original_names, new_names)) {
          track[[nn]] <- NA
        }
        rbind(original_track, track)
      },
      ...
    )
  }

  # Initialise required local options
  inla.options <- modifyList(
    inla.options,
    list(
      control.mode = list(),
      control.predictor = list(compute = TRUE)
    )
  )

  # Extract the family of each likelihood
  family <- vapply(
    seq_along(lhoods),
    function(k) lhoods[[k]]$inla.family,
    "family"
  )

  initial <-
    if (is.null(initial)) {
      options[["bru_initial"]]
    } else {
      initial
    }
  result <- NULL
  if (is.null(initial) || inherits(initial, "bru")) {
    if (!is.null(initial[["bru_iinla"]][["states"]])) {
      states <- initial[["bru_iinla"]][["states"]]
      # The last linearisation point will be used again:
      states <- c(states, states[length(states)])
    } else {
      # Set old result
      states <- evaluate_state(model,
        lhoods = lhoods,
        result = initial,
        property = "joint_mode"
      )
    }
    if (inherits(initial, "bru")) {
      result <- initial
    }
  } else if (is.list(initial)) {
    state <- initial[intersect(names(model[["effects"]]), names(initial))]
    for (lab in names(model[["effects"]])) {
      if (is.null(state[[lab]])) {
        state[[lab]] <- rep(0, ibm_n(model[["effects"]][[lab]][["mapper"]]))
      } else if (length(state[[lab]]) == 1) {
        state[[lab]] <-
          rep(
            state[[lab]],
            ibm_n(model[["effects"]][[lab]][["mapper"]])
          )
      }
    }
    states <- list(state)
    result <- NULL
  } else {
    stop("Unknown previous result information class")
  }
  old.result <- result

  # Track variables
  track <- list()
  if (is.null(old.result[["bru_iinla"]][["track"]])) {
    original_track <- NULL
    track_size <- 0
  } else {
    original_track <- old.result[["bru_iinla"]][["track"]]
    track_size <- max(original_track[["iteration"]])
  }


  do_line_search <- (length(options[["bru_method"]][["search"]]) > 0)
  if (do_line_search || !identical(options$bru_method$taylor, "legacy")) {
    A <- evaluate_A(model, lhoods)
    lin <- bru_compute_linearisation(
      model,
      lhoods = lhoods,
      state = states[[length(states)]],
      A = A
    )
  }
  # Initial stack
  if (identical(options$bru_method$taylor, "legacy")) {
    stk <- joint_stackmaker(
      model,
      lhoods,
      state = list(states[[length(states)]])
    )
  } else {
    idx <- evaluate_index(model, lhoods)
    stk <- bru_make_stack(lhoods, lin, idx)
  }
  stk.data <- INLA::inla.stack.data(stk)
  inla.options$control.predictor$A <- INLA::inla.stack.A(stk)

  k <- 1
  interrupt <- FALSE
  line_search <- list(
    active = FALSE,
    step_scaling = 1,
    state = states[[length(states)]]
  )

  if ((!is.null(result) && !is.null(result$mode))) {
    previous_x <- result$mode$x
  } else {
    previous_x <- 0 # TODO: construct from the initial linearisation state instead
  }

  do_final_integration <- (options$bru_max_iter == 1)
  do_final_theta_no_restart <- FALSE
  while (!interrupt) {
    if ((k >= options$bru_max_iter) && !do_final_integration) {
      do_final_integration <- TRUE
      bru_log_message(
        "iinla: Maximum iterations reached, running final INLA integration.",
        verbose = options$bru_verbose,
        verbose_store = options$bru_verbose_store
      )
    }

    # When running multiple times propagate theta
    if ((k > 1) || (!is.null(result) && !is.null(result$mode))) {
      inla.options[["control.mode"]]$restart <- TRUE
      inla.options[["control.mode"]]$theta <- result$mode$theta
      inla.options[["control.mode"]]$x <-
        previous_x + (result$mode$x - previous_x) * line_search[["step_scaling"]]
    }
    if ((!is.null(result) && !is.null(result$mode))) {
      previous_x <- result$mode$x
    }

    bru_log_message(
      paste0("iinla: Iteration ", k, " [max:", options$bru_max_iter, "]"),
      verbose = options$bru_verbose,
      verbose_store = options$bru_verbose_store
    )

    # Return previous result if inla crashes, e.g. when connection to server is lost
    old.result <- result
    result <- NULL

    inla.formula <- update.formula(model$formula, BRU.response ~ .)
    inla.data <-
      c(
        stk.data,
        do.call(c, c(lapply(
          model$effects,
          function(xx) as.list(xx$env_extra)
        ),
        use.names = FALSE
        ))
      )
    #        list.data(model$formula))

    inla.options.merged <-
      modifyList(
        inla.options,
        list(
          formula = inla.formula,
          data = inla.data,
          family = family,
          E = stk.data[["BRU.E"]],
          Ntrials = stk.data[["BRU.Ntrials"]],
          offset = stk.data[["BRU.offset"]]
        )
      )
    if (do_final_integration) {
      if (do_final_theta_no_restart) {
        # Compute the minimal amount required
        inla.options.merged <-
          modifyList(
            inla.options.merged,
            list(
              control.mode = list(restart = FALSE)
            )
          )
      }
    } else {
      # Compute the minimal amount required
      # Note: configs = TRUE only because latent indexing into mode$x is
      #   defined in configs$contents
      inla.options.merged <-
        modifyList(
          inla.options.merged,
          list(
            control.inla = list(int.strategy = "eb"),
            control.compute = list(
              config = TRUE,
              dic = FALSE,
              waic = FALSE
            ),
            control.predictor = list(compute = FALSE)
          )
        )
    }

    result <- try_callstack(
      do.call(
        INLA::inla,
        inla.options.merged,
        envir = environment(model$effects)
      )
    )

    if (inherits(result, "try-error")) {
      bru_log_message(
        paste0("iinla: Problem in inla: ", result),
        verbose = FALSE,
        verbose_store = options$bru_verbose_store
      )
      warning(
        paste0("iinla: Problem in inla: ", result),
        immediate. = TRUE
      )
      bru_log_message(
        paste0("iinla: Giving up and returning last successfully obtained result for diagnostic purposes."),
        verbose = TRUE,
        verbose_store = options$bru_verbose_store
      )
      if (is.null(old.result)) {
        old.result <- list()
        class(old.result) <- c("iinla", "list")
      } else if (!inherits(old.result, "iinla")) {
        class(old.result) <- c("iinla", class(old.result))
      }

      old.result[["bru_iinla"]] <- collect_misc_info()
      old.result$error <- result
      return(old.result)
    }

    # Extract values tracked for estimating convergence
    if (options$bru_max_iter > 1 & k <= options$bru_max_iter) {
      # Note: The number of fixed effects may be zero, and strong
      # non-linearities that don't necessarily affect the fixed
      # effects may appear in the random effects, so we need to
      # track all of them.
      result_mode <- evaluate_state(model, result, property = "joint_mode")[[1]]
      result_sd <- evaluate_state(model, result,
        property = "sd",
        internal_hyperpar = TRUE
      )[[1]]
      result_names <- latent_names(result_mode)
      track_df <- list()
      for (label in names(result_mode)) {
        track_df[[label]] <-
          data.frame(
            effect = result_names[[label]],
            iteration = track_size + k,
            mode = result_mode[[label]],
            sd = result_sd[[label]]
          )
      }
      track[[k]] <- do.call(rbind, track_df)
    }

    # Only update the linearisation state after the non-final "eb" iterations:
    if (!do_final_integration) {
      # Update stack given current result
      state0 <- states[[length(states)]]
      state <- evaluate_state(model,
        result = result,
        property = "joint_mode"
      )[[1]]
      if ((options$bru_max_iter > 1)) {
        if (do_line_search) {
          line_search <- bru_line_search(
            model = model,
            lhoods = lhoods,
            lin = lin,
            state0 = state0,
            state = state,
            A = A,
            options = options
          )
          state <- line_search[["state"]]
        }
        if (do_line_search || !identical(options$bru_method$taylor, "legacy")) {
          lin <- bru_compute_linearisation(
            model,
            lhoods = lhoods,
            state = state,
            A = A
          )
        }
        if (identical(options$bru_method$taylor, "legacy")) {
          stk <- joint_stackmaker(model, lhoods, state = list(state))
        } else {
          stk <- bru_make_stack(lhoods, lin, idx)
        }
        stk.data <- INLA::inla.stack.data(stk)
        inla.options$control.predictor$A <- INLA::inla.stack.A(stk)
      }
      # Store the state
      states <- c(states, list(state))
    }

    # Stopping criteria
    if (do_final_integration || (k >= options$bru_max_iter)) {
      interrupt <- TRUE
    } else if (!interrupt && (k > 1)) {
      max.dev <- options$bru_method$stop_at_max_rel_deviation
      dev <- abs(track[[k - 1]]$mode - track[[k]]$mode) / track[[k]]$sd
      ## do.call(c, lapply(by(do.call(rbind, track),
      ##                           as.factor(do.call(rbind, track)$effect),
      ##                           identity),
      ##                        function(X) { abs(X$mean[k-1] - X$mean[k])/X$sd[k] }))
      bru_log_message(
        paste0(
          "iinla: Max deviation from previous: ",
          signif(100 * max(dev), 3),
          "% of SD [stop if: <", 100 * max.dev, "%]"
        ),
        verbose = options$bru_verbose,
        verbose_store = options$bru_verbose_store,
        verbosity = 1
      )
      do_final_integration <- all(dev < max.dev) && (!line_search[["active"]])
      if (do_final_integration) {
        do_final_theta_no_restart <- TRUE
        bru_log_message(
          "iinla: Convergence criterion met, running final INLA integration with known theta mode.",
          verbose = options$bru_verbose,
          verbose_store = options$bru_verbose_store
        )
      }
    }
    k <- k + 1
  }

  result[["bru_iinla"]] <- collect_misc_info()
  class(result) <- c("iinla", class(result))
  return(result)
}


auto_intercept <- function(components) {
  if (!inherits(components, "formula")) {
    if (inherits(components, "component_list")) {
      return(components)
    }
    stop("components must be a formula to auto-add an intercept component")
  }
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
          " -1"
        ))
      )
    }
  } else if (attr(tm, "intercept")) {
    components <- update.formula(components, . ~ . + Intercept(1) + 1)
  } else {
    components <- update.formula(components, . ~ . + 1)
  }
  environment(components) <- env
  components
}


# If missing RHS, set to the sum of the component names
auto_additive_formula <- function(formula, components) {
  if (as.character(formula)[length(as.character(formula))] != ".") {
    return(formula)
  }
  stopifnot(inherits(components, "component_list"))

  env <- environment(formula)
  formula <- update.formula(
    formula,
    paste0(
      ". ~ ",
      paste0(names(component), collapse = " + ")
    )
  )
  environment(formula) <- env
  formula
}

# Extract the LHS of a formula, as response ~ .
extract_response <- function(formula) {
  stopifnot(inherits(formula, "formula"))
  if (length(as.character(formula)) == 3) {
    as.formula(paste0(as.character(formula)[2], " ~ ."))
  } else {
    . ~ .
  }
}
# If response missing, add one
auto_response <- function(formula, response) {
  if (!inherits(formula, "formula")) {
    return(formula)
  }
  if ((length(as.character(formula)) == 3) &&
    (as.character(formula)[2] != ".")) {
    # Already has response; do nothing
    return(formula)
  }
  env <- environment(formula)
  if (as.character(formula)[length(as.character(formula))] == ".") {
    # No formula RHS
    formula <- response
  } else {
    formula <- update.formula(formula, response)
  }
  environment(formula) <- env
  formula
}


# Returns a formula's environment as a list. Removes all variable that are of type
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






# Summary methods ----

#  Summarise a LGCP object
#
# @export
# @param object A result object obtained from a lgcp() run
# @param ... ignored arguments (S3 generic compatibility)

summary.lgcp <- function(object, ...) {
  result <- object
  warning("The summary.lgcp() method probably doesn't work with the current devel inlabru version!")

  cat("### LGCP Summary #################################################################################\n\n")

  cat(paste0("Predictor: log(lambda) = ", as.character(result$model$expr), "\n"))

  cat("\n--- Points & Samplers ----\n\n")
  cat(paste0("Number of points: ", nrow(result$bru_info$points)), "\n")
  if (inherits(result$bru_info$points, "Spatial")) {
    cat(paste0("Coordinate names: ", paste0(coordnames(result$bru_info$points), collapse = ", ")), "\n")
    cat(paste0("Coordinate system: ", proj4string(result$bru_info$points), "\n"))
  }

  cat(paste0("Total integration mass, E*weight: ", sum(result$bru_info$lhoods[[1]]$E)), "\n")

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


#' Summary for an inlabru fit
#'
#' Takes a fitted `bru` object produced by [bru()] or [lgcp()] and creates
#' various summaries from it.
#'
#' @export
#' @method summary bru
#' @param object An object obtained from a [bru()] or [lgcp()] call
#' @param \dots ignored arguments
#' @example inst/examples/bru.R
#'

summary.bru <- function(object, ...) {
  object <- bru_check_object_bru(object)

  result <- list(
    bru_info = summary(object[["bru_info"]])
  )

  result$WAIC <- object[["waic"]][["waic"]]
  result$DIC <- object[["dic"]][["dic"]]

  if (inherits(object, "inla")) {
    result$inla <- NextMethod("summary", object)
    result[["inla"]][["call"]] <- NULL
  }

  class(result) <- c("summary_bru", "list")
  return(result)

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
  for (nm in names(object$bru_info$model$effects)) {
    eff <- object$bru_info$model$effects[[nm]]
    if (identical(eff[["main"]][["type"]], "spde")) {
      hyp <- INLA::inla.spde.result(object, nm, eff$main$model)
      cat(sprintf("\n--- Field '%s' transformed hyper parameters ---\n", nm))
      df <- rbind(
        marginal.summary(hyp$marginals.range.nominal$range.nominal.1, "range"),
        marginal.summary(hyp$marginals.variance.nominal$variance.nominal.1, "variance"),
        marginal.summary(hyp$marginals.variance.nominal$variance.nominal.1, "variance"),
      )
      print(df)
    }
  }
  message(
    "The current summary.bru(...) method is outdated and will be replaced.\n",
    "Until then, you may prefer the output from INLA:::summary.inla(...) as an alternative."
  )
  class(result) <- c("summary_bru", "list")
  result
}


#' @export
#' @param x An `summary_bru2` object
#' @rdname summary.bru

print.summary_bru <- function(x, ...) {
  print(x$bru_info)
  print(x$inla)
  invisible(x)
}







#' @describeIn inlabru-deprecated Old summary for an inlabru fit.
#'
#' Takes a fitted `bru` object produced by [bru()] or [lgcp()] and creates
#' various summaries from it.
#'
#' @export
#' @param object An object obtained from a [bru()] or [lgcp()] call
#' @param \dots arguments passed on to other methods or ignored

summary_bru <- function(object, ...) {
  .Deprecated(new = "summary")

  cat("\n--- Likelihoods ----------------------------------------------------------------------------------\n\n")
  for (k in 1:length(object$bru_info$lhoods)) {
    lh <- object$bru_info$lhoods[[k]]
    cat(sprintf("Name: '%s', family: '%s', data class: '%s', \t formula: '%s' \n", names(object$bru_info$lhoods)[[k]], lh$family, class(lh$data), deparse(lh$formula)))
  }

  # rownames(df) = names(object$bru_info$lhoods)
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
    if (inherits(object$model$effects[[nm]]$main$mapper, "inla.mesh")) {
      cat(paste0(
        ", and area = ",
        signif(sum(INLA::inla.mesh.fem(object$model$effects[[nm]]$main$mapper)$va))
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
  for (nm in names(object$bru_info$model$effects)) {
    eff <- object$bru_info$model$effects[[nm]]
    if (identical(eff[["main"]][["type"]], "spde")) {
      hyp <- INLA::inla.spde.result(object, nm, eff$main$model)
      cat(sprintf("\n--- Field '%s' transformed hyper parameters ---\n", nm))
      df <- rbind(
        marginal.summary(hyp$marginals.range.nominal$range.nominal.1, "range"),
        marginal.summary(hyp$marginals.variance.nominal$variance.nominal.1, "variance"),
        marginal.summary(hyp$marginals.variance.nominal$variance.nominal.1, "variance"),
      )
      print(df)
    }
  }
}
