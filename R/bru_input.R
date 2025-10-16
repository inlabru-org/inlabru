#' @include deprecated.R

#' @title Store and evaluate component inputs
#' @description Store and evaluate component inputs, including support for
#'   non-standard-evaluation with `rlang`, automatic transfer to
#'   function calls, [eval_spatial()], and [MatrixModels::model.Matrix()].
#'
#'   These functions are normally only called internally, but users may find
#'   it useful to call them directly to experiment, and to make advanced
#'   model definitions using [ibm_input_set]/[ibm_input_new()].
#' @name bru_input
#' @rdname bru_input
NULL

#' @describeIn bru_input Create a `bru_input` object.
#' @param input An expression to be evaluated
#' @param label character; optional label used to identify the object in
#'   informational messages.
#' @param layer Optional expression that evaluates layer names or indices for
#'   use with [eval_spatial()]
#' @param selector character or integer; optional selector for use with
#'   use with [eval_spatial()]
#' @param \dots Passed on to sub-methods.
#' @export
#' @examples
#' (inp <- new_bru_input(x, "LABEL"))
#' bru_input(inp, data.frame(x = 1:3))
#'
new_bru_input <- function(input,
                          label = NULL,
                          layer = NULL,
                          selector = NULL,
                          ...) {
  inp <- structure(
    list(
      input = rlang::enquo(input),
      label = label,
      layer = rlang::enquo(layer),
      selector = selector
    ),
    class = "bru_input"
  )
  inp
}

#' @describeIn bru_input Evaluate inputs
#'
#' @export
#' @seealso [summary.bru_input()], [bru_comp()]
bru_input <- function(...) {
  UseMethod("bru_input")
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
#' Interactions can be handled by a formula input and `model = "fixed"`:
#' `components = y ~ 0 + name(~ 1 + x:z, model = "fixed")`
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
    glue("bru_input.bru_comp({component$label})"),
    verbosity = 4
  )

  if (is.null(component[["mapper"]])) {
    stop(glue("Component mapper for '{component$label}' is NULL"))
  }
  return(bru_input(component[["mapper"]], data = data, ...))

  if (is.null(component[["mapper"]])) {
    part_names <- intersect(c("main", "group", "replicate"), names(component))
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

#' @param model A [bru_model] object.
#' @export
#' @rdname bru_input
bru_input.bru_model <- function(model, lhoods, ...) {
  bru_input(lhoods, components = as_bru_comp_list(model))
}

#' @param components A [bru_comp_list] object.
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
#' @param lhoods A [bru_obs_list] object.
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
    rlang::eval_tidy(layer, data = envir, env = enclos),
    error = function(e) {
      e
    }
  )
  if (inherits(input_layer, "error")) {
    stop(paste0(
      "Failed to evaluate 'layer' input '",
      glue_collapse(rlang::as_label(layer), sep = "\n"),
      glue("' for '{label}:layer'.")
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
#' 1. If `rlang::eval_tidy()` failed, return NULL or map everything to 1
#'    (see the `null.on.fail` argument). This should normally not
#'    happen, unless the component use logic is incorrect,
#'    (via `used`)
#'    leading to missing columns for a certain likelihood in a
#'    multi-`bru_obs()` model.
#' 2. If we obtain a function, apply the function to the data object
#' 3. If we obtain an object supported by [eval_spatial()], extract the values
#'    of that data frame at the point locations
#' 4. If we obtain a formula, call `ModelMatrix::model.Matrix()`
#' 5. Else we obtain a vector and return as-is. This happens when input
#'    references a column of the data points, or some other complete expression
#'
#' @param env environment in which to evaluate the input expression. Default
#'   is `parent.frame()`
#' @param null.on.fail logical; if `TRUE`, return `NULL` if the input evaluation
#'   fails. If `FALSE` (default), return a vector of 1s and issue a warning
#'   (only for deprecated use of `coordinates` without having loaded the `sp`
#'   package), or stop with an error.
#' @export
bru_input.bru_input <- function(input, data, env = NULL,
                                null.on.fail = FALSE, ...) {
  bru_log_message(
    glue("bru_input(bru_input) for ({input$label})"),
    verbosity = 5
  )
  # Evaluate the map with the data in an environment
  enclos <-
    if (is.null(env)) {
      parent.frame()
    } else {
      env
    }
  if (is.list(data)) {
    envir <- rlang::as_data_mask(data)
  } else if (inherits(data, "Spatial")) {
    envir <- rlang::as_data_mask(as.data.frame(data))
  } else {
    envir <- rlang::as_data_mask(tibble::as_tibble(data))
  }
  assign(".data.", data, envir = envir)

  inp <- rlang::as_quosure(input$input, env = enclos)
  e_input <- tryCatch(
    rlang::eval_tidy(inp, data = envir, env = enclos),
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
      input_string <- paste0(rlang::as_label(input$input), collapse = "\n")
      if (identical(input_string, "coordinates")) {
        warning(
          glue(
            "The input evaluation '{input_string}' for '{input$label}'",
            " failed. Perhaps you need to load the 'sp' package ",
            "with 'library(sp)'?",
            " Attempting 'sp::coordinates'."
          ),
          immediate. = TRUE
        )
        input$input <- rlang::as_quosure(sp::coordinates, env = enclos)
        return(bru_input(
          input,
          data = data,
          env = enclos,
          null.on.fail = null.on.fail,
          ...
        ))
      } else {
        stop(
          glue(
            "The input evaluation '{input_string}' for '{input$label}' failed.",
            "\n  Perhaps the data object doesn't contain the needed variables?"
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
    msg <- glue(
      "Model input '",
      glue_collapse(rlang::as_label(input$input), sep = "\n"),
      "' for '{input$label}' returned some NA values.\n",
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
#' @describeIn inlabru-deprecated `r lifecycle::badge("deprecated")` from `2.12.0.9023`.
#' Use `bru_input()` instead.
input_eval <- function(...) {
  lifecycle::deprecate_warn(
    "2.12.0.9023",
    "input_eval()",
    "bru_input()"
  )
  bru_input(...)
}

#' @describeIn inlabru-deprecated `r lifecycle::badge("deprecated")` since version
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
#' @param object Object to be summarised.
#' @param \dots Passed on to other summary methods.
#' @param verbose logical; If `TRUE`, includes more details of the
#' component definitions. When `FALSE`, only show basic component
#' definition information.  Default `TRUE`.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [bru_input()], [bru_comp()]
#' @param label.override character; If not `NULL`, use this label instead of
#' the object's label.
#' @param type character; if non-NULL, added to the output'; `label =
#'   type(input)`.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname summary.bru_input
#' @export
#' @method format bru_input
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
  } else {
    the_input <- glue_collapse(rlang::as_label(inp), sep = "\n")
    if (is.null(lab) || identical(lab, "")) {
      the_lab <- ""
    } else {
      the_lab <- glue("{lab} = ")
    }
    if (is.null(type)) {
      text <- glue(the_lab, the_input)
    } else {
      text <- glue(the_lab, "{type}(", the_input, ")")
    }
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





#' @title Interface between `bru_input` and `bru_mapper`
#' @description Associate [bru_input] objects with [bru_mapper] objects.
#' @name ibm_input
#' @seealso [bru_input()]
NULL

#' @describeIn ibm_input Add an existing `bru_input` to a `bru_mapper`.
#' @param mapper A `bru_mapper` object.
#' @param input A `bru_input` object, or `NULL` to remove any existing input.
#'   Alternatively, an existing `bru_mapper` object with or without associated
#'   input can be provided, in which case the input from that mapper
#'   is copied.
#' @export
#' @examples
#' (m <- bm_autodetect())
#' ibm_input_available(m)
#' (m <- ibm_input_new(m, cos(x)))
#' ibm_input_available(m)
#' ibm_input_set(bm_linear(), m)
#'
ibm_input_set <- function(mapper, input) {
  stopifnot(inherits(mapper, "bru_mapper"))
  inp <- input
  if (inherits(input, "bru_mapper")) {
    if (ibm_input_available(input)) {
      inp <- ibm_input_get(input)
    } else {
      inp <- NULL
    }
  }
  stopifnot(is.null(inp) || inherits(inp, "bru_input"))
  # Use a special name that is unlikely to conflict with user-defined names.
  mapper[[".input"]] <- inp
  mapper
}
#' @describeIn ibm_input Create and add a `bru_input` to a `bru_mapper`.
#' @param \dots Passed on to [new_bru_input()].
#' @export
ibm_input_new <- function(mapper, ...) {
  stopifnot(inherits(mapper, "bru_mapper"))
  ibm_input_set(mapper, new_bru_input(...))
}
#' @describeIn ibm_input Check if a `bru_input` is associated with a
#' `bru_mapper`.
#' @export
ibm_input_available <- function(mapper) {
  !is.null(mapper[[".input"]])
}
#' @describeIn ibm_input Get the `bru_input` associated with a `bru_mapper`.
#' @export
ibm_input_get <- function(mapper) {
  stopifnot(ibm_input_available(mapper))
  mapper[[".input"]]
}
#' @describeIn bru_input Evaluate the input associated with a `bru_mapper`.
#' @seealso [ibm_input]
#' @export
bru_input.bru_mapper <- function(input, ..., label = "<unknown>") {
  bru_log_message(
    glue("bru_input.bru_mapper({label})"),
    verbosity = 5
  )
  if (!ibm_input_available(input)) {
    stop(glue(
      "No input defined for mapper '{label}'.",
      "\n Use `ibm_input_new()` or `ibm_input_set()` to assign one,",
      "\n like `bm_<type>(...) |> ibm_input_new(<expr>, ...)`."
    ))
  }
  inp <- ibm_input_get(input)
  bru_input(inp, ..., label = label)
}
#' @describeIn bru_input Evaluate the inputs for each sub-mapper in a `bm_pipe`
#' object.
#' @export
bru_input.bm_pipe <- function(input, ..., label = "<unknown>") {
  bru_log_message(
    glue("bru_input.bm_pipe({label})"),
    verbosity = 5
  )
  indexing <- names(input$mappers)
  if (is.null(indexing)) {
    indexing <- seq_along(input$mappers)
  } else {
    names(indexing) <- indexing
  }
  lapply(indexing, function(idx) {
    bru_input(input$mappers[[idx]], ..., label = glue("{label}:{idx}"))
  })
}
#' @describeIn bru_input Evaluate the inputs for each sub-mapper in a `bm_multi`
#' object.
#' @export
bru_input.bm_multi <- function(input, ..., label = "<unknown>") {
  bru_log_message(
    glue("bru_input.bm_multi({label})"),
    verbosity = 5
  )
  if (ibm_input_available(input)) {
    inp <- ibm_input_get(input)
    return(bru_input(inp, ..., label = label))
  }
  indexing <- names(input$mappers)
  if (is.null(indexing)) {
    indexing <- seq_along(input$mappers)
  } else {
    names(indexing) <- indexing
  }
  lapply(indexing, function(idx) {
    bru_input(input$mappers[[idx]], ..., label = glue("{label}:{idx}"))
  })
}
#' @describeIn bru_input Evaluate the inputs for each sub-mapper in a
#'   `bm_collect` object.
#' @export
bru_input.bm_collect <- function(input, ..., label = "<unknown>") {
  bru_log_message(
    glue("bru_input.bm_collect({label})"),
    verbosity = 5
  )
  if (ibm_input_available(input)) {
    inp <- ibm_input_get(input)
    return(bru_input(inp, ..., label = label))
  }
  indexing <- names(input$mappers)
  if (is.null(indexing)) {
    indexing <- seq_along(input$mappers)
  } else {
    names(indexing) <- indexing
  }
  lapply(indexing, function(idx) {
    bru_input(input$mappers[[idx]], ..., label = glue("{label}:{idx}"))
  })
}
#' @describeIn bru_input Evaluate the inputs for the sub-mapper in a `bm_repeat`
#' object.
#' @export
bru_input.bm_repeat <- function(input, ..., label = "<unknown>") {
  bru_log_message(
    glue("bru_input.bm_repeat({label})"),
    verbosity = 5
  )
  if (ibm_input_available(input)) {
    inp <- ibm_input_get(input)
    return(bru_input(inp, ..., label = label))
  }
  bru_input(input[["mapper"]], ..., label = label)
}
#' @describeIn bru_input Evaluate the inputs for each sub-mapper in a `bm_sum`
#' object.
#' @export
bru_input.bm_sum <- function(input, ..., label = "<unknown>") {
  bru_log_message(
    glue("bru_input.bm_sum({label})"),
    verbosity = 5
  )
  if (ibm_input_available(input)) {
    inp <- ibm_input_get(input)
    return(bru_input(inp, ..., label = label))
  }
  indexing <- names(input$mappers)
  if (is.null(indexing)) {
    indexing <- seq_along(input$mappers)
  } else {
    names(indexing) <- indexing
  }
  lapply(indexing, function(idx) {
    bru_input(input$mappers[[idx]], ..., label = glue("{label}:{idx}"))
  })
}



#' @describeIn ibm_input Create a `bru_mapper` placeholder object of class
#'   `bm_autodetect`. The main purpose of this class is to attach [bru_input]
#'   information to it, which is later used to determine a suitable 'real'
#'   mapper type.
#' @export
bm_autodetect <- function() {
  bru_mapper_define(
    list(n = 1L),
    new_class = "bm_autodetect"
  )
}
