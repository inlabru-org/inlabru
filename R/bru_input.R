#' Obtain component inputs
#'
#' @export
#' @rdname bru_input
#' @param \dots Passed on to sub-methods.
#' @seealso [summary.bru_input()], [bru_comp()]
bru_input <- function(...) {
  UseMethod("bru_input")
}
#' @describeIn bru_input Create a `bru_input` object.
#' @export
#' @examples
#' inp <- bru_input(expression(x), "LABEL")
#' bru_input(inp, data.frame(x = 1:3))
bru_input.default <- function(input,
                              label = NULL,
                              layer = NULL,
                              selector = NULL,
                              ...) {
  inp <- structure(
    list(
      input = input,
      label = label,
      layer = layer,
      selector = selector
    ),
    class = "bru_input"
  )
  inp
}

#' @section Simple covariates and the input parameters:
#'
#' It is not unusual for a random effect act on a transformation of a covariate.
#' In other frameworks this would mean that the transformed covariate would have
#' to be calculated in advance and added to the data frame that is usually
#' provided via the `data` parameter. inlabru provides the option to do this
#' transformation automatically. For instance, one might be interested in the
#' effect of a covariate \eqn{x^2}. In inla and other frameworks this would
#' require to add a column `xsquared` to the input data frame and use the
#' formula
#'
#' \itemize{\item{`formula = y ~ f(xsquared, model = "linear")`,}}
#'
#' In inlabru this can be achieved in several ways of using the `main` parameter
#' (`map` in version 2.1.13 and earlier), which does not need to be named.
#'
#' \itemize{
#' \item{`components = y ~ psi(x^2, model = "linear")`}
#' \item{`components = y ~ psi(main = x^2, model = "linear")`}
#' \item{`components = y ~ psi(mySquareFun(x), model = "linear")`,}
#' \item{`components = y ~ psi(myOtherSquareFun, model = "linear")`,}
#'
#' }
#'
#' In the first example inlabru will interpret the map parameter as an
#' expression to be evaluated within the data provided. Since \eqn{x} is a known
#' covariate it will know how to calculate it. The second example is an
#' expression as well but it uses a function called `mySquareFun`. This function
#' is defined by user but has to be accessible within the work space when
#' setting up the components. The third example provides the function
#' `myOtherSquareFun`. In this case, inlabru will call the function as
#' `myOtherSquareFun(.data.)`, where `.data.` is the data provided via the
#' [bru_obs()] `data` parameter. The function needs to know what parts of the
#' data to use to construct the needed output. For example,
#' ```
#' myOtherSquareFun <- function(data) {
#'   data[ ,"x"]^2
#' }
#' ```
#'
#' @section Spatial Covariates:
#'
#' When fitting spatial models it is common to work with covariates that depend
#' on space, e.g. sea surface temperature or elevation. Although it is
#' straightforward to add this data to the input data frame or write a covariate
#' function like in the previous section there is an even more convenient way in
#' inlabru. Spatial covariates are often stored as `SpatialPixelsDataFrame`,
#' `SpatialPixelsDataFrame` or `RasterLayer` objects. These can be provided
#' directly via the input expressions if they are supported by [eval_spatial()],
#' and the [bru_obs()] data is an `sf` or `SpatialPointsDataFrame` object.
#' `inlabru` will then automatically evaluate and/or interpolate the covariate
#' at your data locations when using code like
#' ```
#' components = y ~ psi(mySpatialPixels, model = "linear")
#' ```
#' For more precise control, use the the `layer` and `selector` arguments (see
#' [bru_comp()]), or call [eval_spatial()] directly, e.g.:
#' ```
#' components = y ~ psi(eval_spatial(mySpatialPixels, where = .data.),
#'                      model = "linear")
#' ```
#'
#' @section Coordinates:
#'
#' A common spatial modelling component when using inla are SPDE models. An
#' important feature of inlabru is that it will automatically calculate the so
#' called A-matrix (a component model matrix) which maps SPDE values at the mesh
#' vertices to values at the data locations. For this purpose, the input can be
#' set to `coordinates`, which is the `sp` package function that extracts point
#' coordinates from the `SpatialPointsDataFrame` that was provided as input to
#' [bru_obs()]. The code for this would look as follows:
#' ```
#' components = y ~ field(coordinates, model = inla.spde2.matern(...))
#' ```
#' Since `coordinates` is a function from the `sp` package, this results in
#' evaluation of `sp::coordinates(.data.)`, which loses any CRS information
#' from the data object.
#'
#' For `sf` data with a geometry column (by default named `geometry`), use
#' ```
#' components = y ~ field(geometry, model = inla.spde2.matern(...))
#' ```
#' Since the CRS information is part of the geometry column of the `sf` object,
#' this retains CRS information, so this is more robust, and allows the model
#' to be built on a different CRS than the observation data.
#'
#'
#' @export
#' @keywords internal
#' @param component A [bru_comp] object.
#' @param data A `data.frame`, `tibble`, `sf`, `list`, or `Spatial*` object of
#' covariates and/or point locations.
#' @return * `bru_input(bru_comp)`: A list of mapper input values, formatted
#'   for the full component mapper (of type [bm_pipe])
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}, Finn Lindgren
#'   \email{finn.lindgren@@gmail.com}
#' @rdname bru_input
bru_input.bru_comp <- function(component,
                               data,
                               ...) {
  bru_log_message(
    paste0("bru_input.bru_comp(", component$label, ")"),
    verbosity = 4
  )

  if (is.null(component[["mapper"]])) {
    part_names <- c("main", "group", "replicate")
  } else {
    stopifnot(inherits(component[["mapper"]], c("bm_pipe", "bru_mapper_pipe")))

    # The names should be a subset of main, group, replicate
    part_names <- ibm_names(component[["mapper"]][["mappers"]][[1]])
  }
  mapper_val <- list()
  for (part in part_names) {
    mapper_val[[part]] <-
      bru_input(
        component[[part]]$input,
        data,
        env = component$env,
        label = part,
        ...
      )
  }
  if (is.null(component[["weights"]])) {
    scale_val <- NULL
  } else {
    scale_val <-
      bru_input(
        component[["weights"]],
        data,
        env = component$env,
        label = paste(component$label, "(weights)"),
        ...
      )
  }
  # Any potential length mismatches must be handled by bm_multi and
  # bm_collect, since e.g. 'main' might take a list() input object.
  list(mapper = mapper_val, scale = scale_val)
}

#' @export
#' @rdname bru_input
bru_input.bru_model <- function(model, lhoods, ...) {
  bru_input(lhoods, components = model[["effects"]])
}

#' @param components A [bru_comp_list].
#' @export
#' @rdname bru_input
bru_input.bru_comp_list <-
  function(components,
           data,
           ...) {
    bru_log_message("bru_input(bru_comp_list)", verbosity = 4)
    lapply(components, function(x) bru_input(x, data = data, ...))
  }

#' @describeIn bru_input Computes the component inputs for included components
#' for each model likelihood
#'
#' @param lhoods A [bru_obs_list] object
#' @export
bru_input.bru_obs_list <- function(lhoods, components, ...) {
  bru_log_message(
    "Evaluate component inputs for each observation model",
    verbosity = 3L
  )
  lapply(
    lhoods,
    function(lh) {
      included <- bru_used(lh)[["effect"]]

      bru_input(
        components[included],
        data = lh[["data"]],
        ...
      )
    }
  )
}


bru_input_layer <- function(layer, selector = NULL, envir, enclos,
                            label,
                            e_input) {
  input_layer <- tryCatch(
    eval(layer, envir = envir, enclos = enclos),
    error = function(e) {
      e
    }
  )
  if (inherits(input_layer, "error")) {
    stop(paste0(
      "Failed to evaluate 'layer' input '",
      paste0(deparse(layer), collapse = "\n"),
      "' for '",
      paste0(label, ":layer"),
      "'."
    ))
  }
  if (is.null(input_layer) && is.null(selector)) {
    if (label %in% names(e_input)) {
      input_layer <- label
    }
  }
  input_layer
}


#' @describeIn bru_input Attempts to evaluate a component input (e.g. `main`,
#' `group`, `replicate`, or `weight`), and process the results:
#' 1. If `eval()` failed, return NULL or map everything to 1
#'    (see the `null.on.fail` argument). This should normally not
#'    happen, unless the component use logic is incorrect,
#'    (e.g. via `used` or `include`/`exclude`)
#'    leading to missing columns for a certain likelihood in a
#'    multi-`bru_obs()` model.
#' 2. If we obtain a function, apply the function to the data object
#' 3. If we obtain an object supported by [eval_spatial()], extract the values
#'    of that data frame at the point locations
#' 4. Else we obtain a vector and return as-is. This happens when input
#'    references a column of the data points, or some other complete expression
#'
#' @export
bru_input.bru_input <- function(input, data, env = NULL,
                                null.on.fail = FALSE, ...) {
  bru_log_message(
    paste0("bru_input(bru_input) for (", input$label, ")"),
    verbosity = 5
  )
  # Evaluate the map with the data in an environment
  enclos <-
    if (is.null(env)) {
      parent.frame()
    } else {
      env
    }
  envir <- new.env(parent = enclos)
  if (is.list(data)) {
    for (nm in names(data)) {
      assign(nm, data[[nm]], envir = envir)
    }
  } else {
    if (inherits(data, "Spatial")) {
      data_df <- as.data.frame(data)
    } else {
      data_df <- tibble::as_tibble(data)
    }
    for (nm in names(data_df)) {
      assign(nm, data_df[[nm]], envir = envir)
    }
  }
  assign(".data.", data, envir = envir)

  e_input <- tryCatch(
    eval(input$input, envir = envir, enclos = enclos),
    error = function(e) {
      e
    }
  )

  # ## Need to handle varying input lengths.
  # ## auto-expansion of scalars needs to happen elsewhere, where it's needed,
  # ## e.g. when using component A-matrices to construct a default linear model
  # ## A matrix.
  #  n <- nrow(as.data.frame(data))

  handle_problems <- function(val) {
    if (is.null(val)) {
      return(NULL)
    }
    if (inherits(e_input, "simpleError")) {
      if (null.on.fail) {
        return(NULL)
      }

      val <- 1
      input_string <- paste0(deparse(input$input), collapse = "\n")
      if (identical(input_string, "coordinates")) {
        warning(
          paste0(
            "The input evaluation '",
            input_string,
            "' for '", input$label,
            "' failed. Perhaps you need to load the 'sp' package ",
            "with 'library(sp)'?",
            " Attempting 'sp::coordinates'."
          ),
          immediate. = TRUE
        )
        input$input <- expression(sp::coordinates)
        return(bru_input(
          input,
          data = data,
          env = env,
          null.on.fail = null.on.fail,
          ...
        ))
      } else {
        stop(
          paste0(
            "The input evaluation '",
            input_string,
            "' for '", input$label,
            "' failed. Perhaps the data object doesn't contain ",
            "the needed variables?"
          )
        )
      }
    }
    val
  }

  e_input <- handle_problems(e_input)
  if (is.null(e_input)) {
    return(NULL)
  }

  if (is.function(e_input)) {
    # Allow but detect failures:
    val <- tryCatch(
      e_input(data),
      error = function(e) {
        e
      }
    )

    val <- handle_problems(val)
    if (is.null(val)) {
      return(NULL)
    }

    if (identical(as.character(input$input), "coordinates")) {
      tryCatch(
        expr = {
          # Return SpatialPoints instead of a matrix
          val <- as.data.frame(val)
          sp::coordinates(val) <- seq_len(ncol(val))
          # Allow proj4string failures:
          data_crs <- tryCatch(fm_CRS(data),
            error = function(e) {
            }
          )
          if (!fm_crs_is_null(data_crs)) {
            sp::proj4string(val) <- data_crs
          }
          val
        },
        error = function(e) {
          NULL
        }
      )
    }
  } else if (inherits(e_input, "formula")) {
    # Allow but detect failures:
    val <- tryCatch(
      MatrixModels::model.Matrix(e_input, data = data, sparse = TRUE),
      error = function(e) {
        e
      }
    )
    val <- handle_problems(val)
    if (is.null(val)) {
      return(NULL)
    }
    val <- as(val, "Matrix")
  } else if (inherits(
    e_input,
    gsub(
      pattern = "^eval_spatial\\.([^*]*)\\*?",
      replacement = "\\1",
      x = format(utils::.S3methods("eval_spatial"))
    )
  )) {
    input_layer <-
      bru_input_layer(
        layer = input[["layer"]],
        selector = input[["selector"]],
        envir = envir,
        enclos = enclos,
        label = input[["label"]],
        e_input = e_input
      )
    layer <- extract_layer(
      data,
      input_layer,
      input[["selector"]]
    )
    check_layer(e_input, data, layer)
    val <- eval_spatial(
      e_input,
      data,
      layer = layer,
      selector = NULL
    )
    #  } else if ((input$label %in% c("offset", "const")) &&
    #    is.numeric(e_input) &&
    #    (length(e_input) == 1)) {
    #    val <- rep(e_input, n)
  } else {
    val <- e_input
  }

  # ## Need to allow different sizes; K %*% effect needs to have same length as
  # response, but the input and effect effect itself doesn't.
  #  # Always return as many rows as data has
  #  if (is.vector(val) && (length(val) == 1) && (n > 1)) {
  #    val <- rep(val, n)
  #  }

  # Check if any of the locations are NA. If we are dealing with
  # SpatialGridDataFrame try to fix that by filling in nearest neighbour values.
  # # TODO: Check how to deal with this fully in the case of multilikelihood
  # models
  # Answer: should respect the lhood "include/exclude" info for the component
  # list
  if ((inherits(e_input, c(
    "SpatialGridDataFrame",
    "SpatialPixelsDataFrame",
    "SpatRaster"
  ))) &&
    any(is.na(as.data.frame(val)))) {
    msg <- paste0(
      "Model input '",
      paste0(deparse(input$input), collapse = "\n"),
      "' for '", input$label,
      "' returned some NA values.\n",
      "Attempting to fill in spatially by nearest available value.\n",
      "To avoid this basic covariate imputation, supply complete data."
    )
    bru_log_warn(msg)

    val <- bru_fill_missing(
      data = e_input, where = data, values = val,
      layer = layer, selector = NULL
    )
  }

  # Check for NA values.
  # TODO: update 2022-11-04:
  # iid models now by default are given _index mappers, which treat NAs
  # as zero.  An option to turn on a warning could be useful for checking
  # models that aren't supposed to have NA inputs.
  #
  #  if (any(is.na(unlist(as.vector(val))))) {
  #    # TODO: remove this check and make sure NAs are handled properly
  #    # elsewhere,
  #    # if possible. Problem: treating NA as "no effect" can have negative side
  #    # effects. For spatial covariates, can be handled by infill, but a
  #    # general
  #    # solution doesn't exist, so not sure how to deal with this.
  #    stop(sprintf(
  #      "Input '%s' of component '%s' has returned NA values. Please design
  #       your argument as to return non-NA for all points in your model
  #       domain/mesh.",
  #      as.character(input$input)[[1]], label
  #    ))
  #  }

  val
}

#' @export
#' @describeIn bru_input `r lifecycle::badge("deprecated")` from `2.12.0.9023`.
#' Use `bru_input()` instead.
input_eval <- function(...) {
  lifecycle::deprecate_warn(
    "2.12.0.9023",
    "input_eval()",
    "bru_input()"
  )
  bru_input(...)
}

#' @describeIn bru_input `r lifecycle::badge("deprecated")` since version
#'   `2.12.0.9023`. Use [bru_input()] instead. Computes the component inputs for
#'   included components for each observation model.
#' @keywords internal
evaluate_inputs <- function(...) {
  lifecycle::deprecate_warn(
    when = "2.12.9023",
    "evaluate_inputs()",
    "bru_input()"
  )
  bru_input(...)
}


#' Summarise component inputs
#'
#' @keywords internal
#' @param object Object to be summarised.
#' @param ... Passed on to other summary methods.
#' @param verbose logical; If `TRUE`, includes more details of the
#' component definitions. When `FALSE`, only show basic component
#' definition information.  Default `TRUE`.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [bru_input()], [bru_comp()]
#' @export
#' @method format bru_input
#' @param label.override character; If not `NULL`, use this label instead of
#' the object's label.
#' @param type character; if non-NULL, added to the output'; `label =
#'   type(input)`.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname summary.bru_input
format.bru_input <- function(x, verbose = TRUE, ..., label.override = NULL,
                             type = NULL) {
  inp <- x[["input"]]
  if (!is.null(label.override)) {
    lab <- label.override
  } else {
    lab <- x[["label"]]
  }
  if (is.null(inp)) {
    text <- ""
  } else if (is.null(type)) {
    text <-
      paste0(
        lab, " = ",
        paste0(deparse(inp), collapse = "\n")
      )
  } else {
    text <-
      paste0(
        lab, " = ",
        type, "(",
        paste0(deparse(inp), collapse = "\n"),
        ")"
      )
  }
  text
}

#' @export
#' @method summary bru_input
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname summary.bru_input
summary.bru_input <- function(object,
                              verbose = TRUE,
                              ...,
                              label.override = NULL) {
  res <- structure(
    list(
      text = format(object,
        verbose = verbose,
        ...,
        label.override = label.override
      )
    ),
    class = "summary_bru_input"
  )
  return(res)
}

#' @export
#' @param x Object to be printed
#' @method print bru_input
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname summary.bru_input
print.bru_input <- function(x, verbose = TRUE, ..., label.override = NULL) {
  cat(
    format(x,
      verbose = verbose,
      ...,
      label.override = label.override
    ),
    "\n",
    sep = ""
  )
  return(invisible(x))
}

#' @export
#' @rdname summary.bru_comp

print.summary_bru_input <- function(x, ...) {
  if (!is.null(x[["text"]])) {
    cat(x[["text"]], "\n", sep = "")
  }
  invisible(x)
}
