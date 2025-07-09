#' Create an inlabru model object from model components
#'
#' The [inlabru] syntax for model formulae is different from what
#' `INLA::inla` considers a valid.
#' In inla most of the effects are defined by adding an `f(...)` expression to
#' the formula.
#' In [inlabru] the `f` is replaced by an arbitrary (exceptions: `const` and
#' `offset`) string that will determine the label of the effect. See Details for
#' further information.
#'
#' @details
#' For instance
#'
#' `y ~ f(myspde, ...)`
#'
#' in INLA is equivalent to
#'
#' `y ~ myspde(...)`
#'
#' in inlabru.
#'
#' A disadvantage of the inla way is that there is no clear separation between
#' the name of the covariate and the label of the effect. Furthermore, for some
#' models like SPDE it is much more natural to use spatial coordinates as
#' covariates rather than an index into the SPDE vertices. For this purpose
#' [inlabru] provides the `main` argument. For convenience, the `main` argument
#' can be used like the first argument of the f function, e.g., and is the first
#' argument of the component definition.
#' The `INLA` model formula
#'
#' `y ~ f(temperature, model = 'linear')`
#'
#' is equivalent to the `inlabru` component and formula definition
#'
#' `y ~ temperature(temperature, model = 'linear')`
#' and
#' `y ~ temperature(main = temperature, model = 'linear')`
#' as well as
#' `y ~ temperature(model = 'linear')`
#' which sets `main = temperature`.
#'
#' On the other hand, `main` can also be a function mapping, e.g the
#' `sp::coordinates()` function:
#'
#' `y ~ mySPDE(coordinates, ...)`
#'
#' This extracts the coordinates from the data object, and maps it to the latent
#' field via the information given in the `mapper`, which by default is
#' extracted from the `model` object, in the case of `spde` model objects.
#'
#' Morevover, `main` can be any expression that evaluates within your data as an
#' environment.
#' For instance, if your data has columns 'a' and 'b', you can create a fixed
#' effect of 'sin(a+b)' by setting `map` in the following way:
#'
#' `y ~ myEffect(sin(a+b))`
#'
#'
#' @export
#' @param components A [bru_comp_list] object
#' @param lhoods Either a [bru_obs()] object, a [bru_obs_list()] object, or
#'   a list of one or more [bru_obs()] or [bru_obs_list()] objects,
#'   or a list of arguments for a call to [bru_obs()].
#' @param options A [bru_options] options object or a list of options passed
#' on to [bru_options()]
#' @param .envir The environment in which the components are evaluated.
#'
#' @return A [bru_model] object
#' @keywords internal

bru_model <- function(components,
                      lhoods,
                      options = list(),
                      .envir = parent.frame()) {
  options <- bru_call_options(options)

  .response <- extract_response(components)

  # Turn input into a list of components (from existing list, or a special
  # formula)
  components <- bru_comp_list(components, .envir = .envir)

  lhoods <- bru_obs_list_construct(
    lhoods,
    options = options,
    .envir = .envir,
    .response = .response,
    .components = components
  )

  if (length(lhoods) == 0) {
    stop("No observation models provided.")
  }

  inputs <- bru_input(lhoods, components = components, null.on.fail = FALSE)

  # Back up environment
  env <- environment(components)

  # Create joint formula that will be used by inla
  formula <- BRU_response ~ -1
  included <- bru_used(lhoods)
  included <- union(included[["effect"]], included[["latent"]])

  # Complete the used component definitions based on data
  components <- bru_comp_list(
    components[included],
    lhoods = lhoods,
    inputs = inputs
  )

  for (cmp in included) {
    if (!(components[[cmp]][["main"]][["type"]] %in% c("offset", "const"))) {
      formula <- update.formula(formula, components[[cmp]]$inla.formula)
    }
  }

  # Restore environment
  environment(components) <- env
  environment(formula) <- env

  # Check linearity
  for (lh in seq_along(lhoods)) {
    if (lhoods[[lh]][["linear"]]) {
      used_lh <- bru_used(lhoods[[lh]])
      used_lh <- unique(c(used_lh$effect, used_lh$latent))
      lhoods[[lh]][["linear"]] <-
        all(vapply(components[used_lh], function(cmp) {
          ibm_is_linear(cmp$mapper)
        }, TRUE))
    }
  }

  # Make model
  mdl <- structure(
    list(
      effects = components,
      lhoods = lhoods,
      inputs = inputs,
      formula = formula
    ),
    class = "bru_model"
  )

  mdl
}


#' @export
#' @method summary bru_model
#' @param object Object to operate on
#' @param \dots Arguments passed on to other methods
#' @rdname bru_model
summary.bru_model <- function(object, ...) {
  result <- structure(
    list(
      components =
        summary(object[["effects"]], ...)
    ),
    class = "summary_bru_model"
  )
  result
}

#' @export
#' @param x An object to be printed
#' @rdname bru_model
print.summary_bru_model <- function(x, ...) {
  print(x[["components"]])
  invisible(x)
}

#' @export
#' @rdname bru_model
print.bru_model <- function(x, ...) {
  print(summary(x))
  invisible(x)
}




#' Evaluate or sample from a posterior result given a model and locations
#'
#' @export
#' @param model A [bru] model
#' @param state list of state lists, as generated by [evaluate_state()]
#' @param data A `list`, `data.frame`, or `Spatial*DataFrame`, with coordinates
#' and covariates needed to evaluate the predictor.
#' @param data_extra Additional data for the predictor evaluation
#' @param input Precomputed inputs list for the components
#' @param comp_simple Precomputed [bm_list] of simplified mappers for the
#' components
#' @param predictor A formula or an expression to be evaluated given the
#' posterior or for each sample thereof. The default (`NULL`) returns a
#' `data.frame` containing the sampled effects. In case of a formula the right
#' hand side is used for evaluation.
#' @param format character; determines the storage format of predictor output.
#' Available options:
#' * `"auto"` If the first evaluated result is a vector or single-column matrix,
#'   the "matrix" format is used, otherwise "list".
#' * `"matrix"` A matrix where each column contains the evaluated predictor
#' expression for a state.
#' * `"list"` A list where each element contains the evaluated predictor
#' expression for a state.
#' @param n Number of samples to draw.
#' @param seed If seed != 0L, the random seed
#' @param num.threads Specification of desired number of threads for parallel
#' computations. Default NULL, leaves it up to INLA.
#' When seed != 0, overridden to "1:1"
#' @param used A [bru_used()] object, or NULL (default)
#' @param n_pred integer. If provided, scalar predictor results are expanded to
#' vectors of length `n_pred`.
#' @param \dots Additional arguments passed on to `inla.posterior.sample`
#' @details * `evaluate_model` is a wrapper to evaluate model state, A-matrices,
#' effects, and predictor, all in one call.
#'
#' @keywords internal
#' @rdname evaluate_model
evaluate_model <- function(model,
                           state,
                           data = NULL,
                           data_extra = NULL,
                           input = NULL,
                           comp_simple = NULL,
                           predictor = NULL,
                           format = NULL,
                           used = NULL,
                           n_pred = NULL,
                           ...) {
  used <- bru_used(used, labels = names(model$effects))

  if (is.null(state)) {
    stop("Not enough information to evaluate model states.")
  }
  if (is.null(input)) {
    input <- bru_input(
      components = model$effects[used$effect],
      data = data
    )
  }
  if (is.null(comp_simple) && !is.null(input)) {
    comp_simple <- ibm_simplify(
      model$effects[used$effect],
      input = input,
      inla_f = TRUE
    )
  }
  if (is.null(comp_simple)) {
    effects <- NULL
  } else {
    effects <-
      lapply(
        state,
        function(x) {
          evaluate_effect_single_state(
            comp_simple,
            state = x,
            input = input
          )
        }
      )
  }

  if (is.null(predictor)) {
    return(effects)
  }

  values <- evaluate_predictor(
    model,
    state = state,
    data = data,
    data_extra = data_extra,
    effects = effects,
    predictor = predictor,
    format = format,
    used = used,
    n_pred = n_pred
  )

  values
}


#' @details * `evaluate_state` evaluates model state properties or samples
#' @param result A `bru` object from [bru()] or [lgcp()]
#' @param property Property of the model components to obtain value from.
#' Default: "mode". Other options are "mean", "0.025quant", "0.975quant",
#' "sd" and "sample". In case of "sample" you will obtain samples from the
#' posterior (see `n` parameter). If `result` is `NULL`, all-zero vectors are
#' returned for each component.
#' @param internal_hyperpar logical; If `TRUE`, return hyperparameter properties
#' on the internal scale. Currently ignored when `property="sample"`.
#' Default is `FALSE`.
#' @export
#' @rdname evaluate_model
#' @keywords internal
evaluate_state <- function(model,
                           result,
                           property = "mode",
                           n = 1,
                           seed = 0L,
                           num.threads = NULL,
                           internal_hyperpar = FALSE,
                           ...) {
  # Evaluate random states, or a single property
  if (property == "sample") {
    state <- post.sample.structured(result,
      n = n, seed = seed,
      num.threads = num.threads,
      ...
    )
  } else if (is.null(result)) {
    state <- list(lapply(
      model[["effects"]],
      function(x) {
        rep(0.0, ibm_n(x[["mapper"]]))
      }
    ))
  } else {
    state <- list(extract_property(
      result = result,
      property = property,
      internal_hyperpar = internal_hyperpar
    ))
  }

  state
}




#' @export
#' @rdname evaluate_effect
evaluate_effect_single_state <- function(...) {
  UseMethod("evaluate_effect_single_state")
}

#' Evaluate a component effect
#'
#' Calculate latent component effects given some data and the state of the
#' component's internal random variables.
#'
#' @export
#' @keywords internal
#' @param component A [bru_mapper], [bru_comp], or
#'   [bm_list].
#' @param input Pre-evaluated component input
#' @param state Specification of one latent variable state:
#' * `evaluate_effect_single_state.bru_mapper`:
#'   A vector of the latent component state.
#' * `evaluate_effect_single_state.*_list`: list of named state vectors.
#' @param ... Optional additional parameters, e.g. `inla_f`. Normally unused.
#' @param label Option label used for any warning messages, specifying the
#' affected component.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com} and
#' Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname evaluate_effect
#' @keywords internal

evaluate_effect_single_state.bru_mapper <- function(component, input, state,
                                                    ...,
                                                    label = NULL) {
  values <- ibm_eval(component, input = input, state = state, ...)

  not_ok <- ibm_invalid_output(
    component,
    input = input,
    state = state
  )
  if (any(not_ok)) {
    if (is.null(label)) {
      warning("Inputs for a mapper give some invalid outputs.",
        immediate. = TRUE
      )
    } else {
      warning("Inputs for '", label, "' give some invalid outputs.",
        immediate. = TRUE
      )
    }
  }

  as.vector(as.matrix(values))
}

#' @return * `evaluate_effect_single_state.bm_list`: A list of
#'   evaluated component effect values
#' @export
#' @rdname evaluate_effect
#' @keywords internal
evaluate_effect_single_state.bm_list <- function(components,
                                                 input,
                                                 state,
                                                 ...) {
  result <- list()
  for (label in names(components)) {
    result[[label]] <- evaluate_effect_single_state(
      components[[label]],
      input = input[[label]],
      state = state[[label]],
      ...,
      label = label
    )
  }
  result
}

#' @export
#' @rdname evaluate_effect
#' @keywords internal
evaluate_effect_single_state.bru_comp_list <- function(components,
                                                       input,
                                                       state,
                                                       ...) {
  comp_simple <- ibm_simplify(components, input = input, state = state, ...)
  evaluate_effect_single_state(comp_simple, input = input, state = state, ...)
}




#' Evaluate component effects or expressions
#'
#' Evaluate component effects or expressions, based on a bru model and one or
#' several states of the latent variables and hyperparameters.
#'
#' @param data A `list`, `data.frame`, or `Spatial*DataFrame`, with coordinates
#' and covariates needed to evaluate the model.
#' @param data_extra Additional data for the predictor evaluation. Variables
#' with the same name as in `data` will be ignored, unless accessed via
#' `.data_extra.[["name"]]` or `.data_extra.$name`.
#' @param state A list where each element is a list of named latent state
#' information, as produced by [evaluate_state()]
#' @param effects A list where each element is list of named evaluated effects,
#' each computed by [evaluate_effect_single_state.bru_comp_list()]
#' @param predictor Either a formula or expression
#' @param used A [bru_used()] object, or NULL (default)
#' @param format character; determines the storage format of the output.
#' Available options:
#' * `"auto"` If the first evaluated result is a vector or single-column matrix,
#'   the "matrix" format is used, otherwise "list".
#' * `"matrix"` A matrix where each column contains the evaluated predictor
#' expression for a state.
#' * `"list"` A list where each column contains the evaluated predictor
#' expression for a state.
#'
#' Default: "auto"
#' @param n_pred integer. If provided, scalar predictor results are expanded to
#' vectors of length `n_pred`.
#' @details For each component, e.g. "name", the state values are available as
#'   `name_latent`, and arbitrary evaluation can be done with `name_eval(...)`,
#'   see [bru_comp_eval()].
#' @return A list or matrix is returned, as specified by `format`
#' @keywords internal
#' @rdname evaluate_predictor
evaluate_predictor <- function(model,
                               state,
                               data,
                               data_extra,
                               effects,
                               predictor,
                               used = NULL,
                               format = "auto",
                               n_pred = NULL) {
  stopifnot(inherits(model, "bru_model"))
  format <- match.arg(format, c("auto", "matrix", "list"))
  pred.envir <- environment(predictor)
  if (inherits(predictor, "formula")) {
    predictor <- parse(
      text = as.character(predictor)[length(as.character(predictor))]
    )
  }
  formula.envir <- environment(model$formula)
  enclos <-
    if (!is.null(pred.envir)) {
      pred.envir
    } else if (!is.null(formula.envir)) {
      formula.envir
    } else {
      parent.frame()
    }

  used <- bru_used(used, labels = names(model$effects))

  # General evaluation environment
  envir <- new.env(parent = enclos)
  # Find .data. first,
  # then data variables,
  # then pred.envir variables (via enclos),
  # then formula.envir (via enclos if pred.envir is NULL):
  #  for (nm in names(pred.envir)) {
  #    assign(nm, pred.envir[[nm]], envir = envir)
  #  }
  #
  # Note: Since 2.7.0.9019, no longer converts Spatial*DataFrame to data frame
  # here; coordinates must be accessed via sp::coordinates() if needed.
  if (!is.null(data)) {
    if (inherits(data, "Spatial")) {
      for (nm in names(data)) {
        assign(nm, data[[nm]], envir = envir)
      }
    } else {
      list2env(data, envir = envir)
    }
  }
  assign(".data.", data, envir = envir)

  nms <- setdiff(names(data_extra), names(data))
  if (!is.null(data_extra[nms])) {
    if (inherits(data_extra, "Spatial")) {
      for (nm in nms) {
        assign(nm, data_extra[[nm]], envir = envir)
      }
    } else {
      list2env(data_extra[nms], envir = envir)
    }
  }
  assign(".data_extra.", data_extra, envir = envir)

  # Rename component states from label to label_latent
  state_names <- as.list(expand_labels(
    names(state[[1]]),
    names(model$effects),
    suffix = "_latent"
  ))
  names(state_names) <- names(state[[1]])

  # Construct _eval function names
  eval_names <- as.list(expand_labels(
    intersect(names(state[[1]]), names(model$effects)),
    intersect(names(state[[1]]), names(model$effects)),
    suffix = "_eval"
  ))
  names(eval_names) <- intersect(names(state[[1]]), names(model$effects))

  eval_fun_factory <-
    function(.comp, .envir, .enclos) {
      .is_offset <- .comp$main$type %in% c("offset", "const")
      .is_iid <- .comp$main$type %in% c("iid")
      .mapper <- .comp$mapper
      .label <- paste0(.comp$label, "_latent")
      .iid_precision <- paste0("Precision_for_", .comp$label)
      .iid_cache <- list()
      .iid_cache_index <- NULL
      eval_fun <- function(main, group = NULL,
                           replicate = NULL,
                           weights = NULL,
                           .state = NULL) {
        n_input <- ibm_n_output(
          .mapper[["mappers"]][["mapper"]][["mappers"]][["main"]],
          input = main
        )
        if (is.null(group)) {
          group <- rep(1, n_input)
        }
        if (is.null(replicate)) {
          replicate <- rep(1, n_input)
        }
        if (!.is_offset && is.null(.state)) {
          .state <- eval(
            parse(text = .label),
            envir = .envir,
            enclos = .enclos
          )
        }
        .input <- list(
          mapper = list(
            main = main,
            group = group,
            replicate = replicate
          ),
          scale = weights
        )

        .values <- ibm_eval(
          .mapper,
          input = .input,
          state = .state
        )
        if (!.is_iid) {
          not_ok <- ibm_invalid_output(
            .mapper[["mappers"]][[1]],
            input = .input[[1]],
            state = .state
          )
          if (any(not_ok)) {
            warning(
              "Inputs for `ibm_eval()` for '",
              .comp[["label"]],
              '" give some invalid outputs.'
            )
          }
        } else { # .is_iid, invalid indices give new samples
          # Check for known invalid output elements, based on the
          # initial mapper (subsequent mappers in the component pipe
          # are assumed to keep the same length and validity)
          not_ok <- ibm_invalid_output(
            .mapper[["mappers"]][[1]],
            input = .input[[1]],
            state = .state
          )
          if (any(not_ok)) {
            .cache_state_index <- eval(
              parse(text = ".cache_state_index"),
              envir = .envir,
              enclos = .enclos
            )
            if (!identical(.cache_state_index, .iid_cache_index)) {
              .iid_cache_index <<- .cache_state_index
              .iid_cache <<- list()
            }
            key <- as.character(main[not_ok])
            not_cached <- !(key %in% names(.iid_cache))
            if (any(not_cached)) {
              .prec <- eval(
                parse(text = .iid_precision),
                envir = .envir,
                enclos = .enclos
              )
              for (k in unique(key[not_cached])) {
                .iid_cache[k] <<- rnorm(1, mean = 0, sd = .prec^-0.5)
              }
            }
            .values[not_ok] <- vapply(
              key,
              function(k) .iid_cache[[k]],
              0.0
            )
          }
        }

        as.matrix(.values)
      }
      eval_fun
    }
  for (nm in names(eval_names)) {
    assign(
      eval_names[[nm]],
      eval_fun_factory(
        model$effects[[nm]],
        .envir = envir,
        .enclos = enclos
      ),
      envir = envir
    )
  }

  # Remove problematic objects:
  problems <- c(".Random.seed")
  remove(list = intersect(names(envir), problems), envir = envir)

  n <- length(state)
  for (k in seq_len(n)) {
    # Keep track of the iteration index so the iid cache can be invalidated
    assign(".cache_state_index", k, envir = envir)

    for (nm in names(state[[k]])) {
      assign(state_names[[nm]], state[[k]][[nm]], envir = envir)
    }
    for (nm in names(effects[[k]])) {
      assign(nm, effects[[k]][[nm]], envir = envir)
    }

    result_ <- eval(predictor, envir = envir, enclos = enclos)
    if (!is.null(n_pred) && is.numeric(result_) && length(result_) == 1) {
      result_ <- rep(result_, n_pred)
    }
    if (k == 1) {
      if (identical(format, "auto")) {
        if ((is.vector(result_) && !is.list(result_)) ||
          (is.matrix(result_) && (NCOL(result_) == 1))) {
          format <- "matrix"
        } else {
          format <- "list"
        }
      }
      if (identical(format, "matrix")) {
        result <- matrix(0.0, NROW(result_), n)
        rownames(result) <- row.names(as.matrix(result_))
      } else if (identical(format, "list")) {
        result <- vector("list", n)
      }
    }
    if (identical(format, "list")) {
      result[[k]] <- result_
    } else {
      result[, k] <- result_
    }
  }

  result
}



#' Evaluate component values in predictor expressions
#'
#' In predictor expressions, `name_eval(...)` can be used to evaluate
#' the effect of a component called "name".
#'
#' @aliases bru_component_eval
#' @param main,group,replicate,weights Specification of where to evaluate a
#'   component. The four inputs are passed on to the joint `bru_mapper` for the
#'   component, as
#'  ```
#'  list(mapper = list(
#'         main = main,
#'         group = group,
#'         replicate = replicate),
#'       scale = weights)
#'  ```
#'  NOTE: If you have model component with the same name as a data variable you
#'  want to supply as input to `name_eval()`, you need to use
#'  `.data.[["myvar"]]` to access it. Otherwise, it will try to use the other
#'  component effect as input, which is ill-defined.
#' @param .state The internal component state. Normally supplied automatically
#' by the internal methods for evaluating inlabru predictor expressions.
#' @return A vector of values for a component
#' @examples
#' if (bru_safe_inla() &&
#'   require("sf", quietly = TRUE) &&
#'   requireNamespace("sn", quietly = TRUE)) {
#'   mesh <- fmesher::fm_mesh_2d_inla(
#'     cbind(0, 0),
#'     offset = 2,
#'     max.edge = 2.5
#'   )
#'   spde <- INLA::inla.spde2.pcmatern(
#'     mesh,
#'     prior.range = c(1, NA),
#'     prior.sigma = c(0.2, NA)
#'   )
#'   set.seed(12345L)
#'   data <- sf::st_as_sf(
#'     data.frame(
#'       x = runif(50),
#'       y = runif(50),
#'       z = rnorm(50)
#'     ),
#'     coords = c("x", "y")
#'   )
#'   fit <- bru(
#'     z ~ -1 + field(geometry, model = spde),
#'     family = "gaussian", data = data,
#'     options = list(control.inla = list(int.strategy = "eb"))
#'   )
#'   pred <- generate(
#'     fit,
#'     newdata = data.frame(A = 0.5, B = 0.5),
#'     formula = ~ field_eval(cbind(A, B)),
#'     n.samples = 1L
#'   )
#' }
bru_comp_eval <- function(main,
                          group = NULL,
                          replicate = NULL,
                          weights = NULL,
                          .state = NULL) {
  stop(paste0(
    "In your predictor expression, use 'mylabel_eval(...)' instead of\n",
    "'bru_comp_eval(...)'.  See ?bru_comp_eval for more information."
  ))
}







#' @include mappers.R

#' @title Mapper methods for model objects
#' @description
#' Methods for the `ibm_linear()` and `ibm_simplify()` methods for
#' [bru] model objects and related classes.
#'
#' @inheritParams bru_mapper_generics
#'
#' @name bru_model_mapper_methods
#' @rdname bru_model_mapper_methods
NULL

#' @describeIn bru_model_mapper_methods Returns a list (one element per
#'   observation model) of [bm_list] objects, each with one [bm_taylor]
#'   entry for each included component.
#' @export
ibm_linear.bru_model <- function(mapper, input, state = NULL, ...) {
  model <- mapper
  stopifnot(inherits(model, "bru_model"))
  bru_log_message(
    paste0("Linearise components for each observation model"),
    verbosity = 3
  )
  mappers <-
    lapply(
      input,
      function(inp) {
        ibm_linear(
          model[["effects"]],
          input = inp,
          state = state,
          ...
        )
      }
    )

  mappers
}

#' @rdname bru_model_mapper_methods
#' @export
ibm_linear.bru_comp_list <- function(mapper, input, state = NULL, ...) {
  comp <- mapper
  included <- parse_inclusion(
    names(comp),
    names(input),
    NULL
  )

  mappers <- ibm_linear(
    as_bm_list(comp[included]),
    input = input[included],
    state = state[included],
    ...
  )

  mappers
}

#' @rdname bru_model_mapper_methods
#' @export
ibm_linear.bru_comp <- function(mapper,
                                input,
                                state = NULL,
                                ...) {
  bru_log_message(
    paste0("Linearise component '", mapper[["label"]], "'"),
    verbosity = 5
  )
  if (is.null(state)) {
    state <- rep(0, ibm_n(mapper[["mapper"]]))
  }
  ibm_linear(mapper[["mapper"]], input = input, state = state, ...)
}

#' @describeIn bru_model_mapper_methods Returns a list (one element per
#'   observation model) of [bm_list] objects, each with one [bru_mapper]
#'   entry for each included component.
#' @export
ibm_simplify.bru_model <- function(mapper, input = NULL, state = NULL, ...) {
  model <- mapper
  bru_log_message(
    paste0("Simplify component mappers for each observation model"),
    verbosity = 3
  )
  mappers <-
    lapply(
      input,
      function(inp) {
        ibm_simplify(
          model[["effects"]],
          input = inp,
          state = state,
          ...
        )
      }
    )

  mappers
}

#' @rdname bru_model_mapper_methods
#' @export
ibm_simplify.bru_comp <- function(mapper,
                                  input = NULL,
                                  state = NULL,
                                  ...) {
  bru_log_message(
    paste0("Linearise component '", mapper[["label"]], "'"),
    verbosity = 5
  )
  if (is.null(state)) {
    state <- rep(0, ibm_n(mapper[["mapper"]]))
  }
  ibm_simplify(mapper[["mapper"]], input = input, state = state, ...)
}


#' @rdname bru_model_mapper_methods
#' @export
ibm_simplify.bru_comp_list <- function(mapper,
                                       input = NULL,
                                       state = NULL,
                                       ...) {
  comp <- mapper
  included <- parse_inclusion(names(comp), names(input), NULL)

  mappers <- ibm_simplify(
    as_bm_list(comp[included]),
    input = input[included],
    state = state[included],
    ...
  )

  mappers
}


# @title Mapper methods for model objects
# @description
# Methods for the `ibm_linear()` and `ibm_simplify()` methods for
# [bru] model objects and related classes.
#
#' @inheritParams bru_mapper_generics
#'
#' @export
#' @rdname bm_list
#' @export
ibm_linear.bm_list <- function(mapper, input, state = NULL, ...) {
  label <- names(mapper)
  if (is.null(label)) {
    label <- as.character(seq_along(mapper))
  }
  mappers <- lapply(
    seq_along(mapper),
    function(k) {
      bru_log_message(
        paste0("Linearise component '", label[k], "'"),
        verbosity = 4
      )
      ibm_linear(
        mapper[[k]],
        input[[k]],
        state = NULL,
        ...
      )
    }
  )
  names(mappers) <- names(mapper)
  as_bm_list(mappers)
}

#' @rdname bm_list
#' @export
ibm_simplify.bm_list <- function(mapper, input = NULL, state = NULL, ...) {
  label <- names(mapper)
  if (is.null(label)) {
    label <- as.character(seq_along(mapper))
  }
  mappers <- lapply(
    seq_along(mapper),
    function(k) {
      bru_log_message(
        paste0("Simplify component '", label[k], "'"),
        verbosity = 4
      )
      ibm_simplify(
        mapper[[k]],
        input[[k]],
        state = NULL,
        ...
      )
    }
  )
  names(mappers) <- names(mapper)
  as_bm_list(mappers)
}
