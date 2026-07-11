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

bru_model <- function(
  components,
  lhoods,
  options = list(),
  .envir = parent.frame()
) {
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

  # Evaluate inputs for non-copy components
  is_copy <- vapply(components, function(x) !is.null(x[["copy"]]), TRUE)
  inputs <- bru_input(
    lhoods,
    components = components[!is_copy],
    null.on.fail = FALSE
  )

  # Back up environment
  env <- environment(components)

  # Complete the used component definitions based on data
  included <- bru_used(lhoods)
  included <- union(included[["effect"]], included[["latent"]])
  components <- bru_comp_list(
    components[included],
    lhoods = lhoods,
    inputs = inputs
  )

  # Evaluate inputs for copy components
  is_copy <- vapply(components, function(x) !is.null(x[["copy"]]), TRUE)
  if (any(is_copy)) {
    inputs_copy <- bru_input(
      lhoods,
      components = components[is_copy],
      null.on.fail = FALSE
    )
    inputs <- lapply(
      seq_along(inputs),
      function(k) {
        if (length(inputs_copy[[k]]) == 0) {
          return(inputs[[k]])
        }
        modifyList(inputs[[k]], inputs_copy[[k]], keep.null = TRUE)
      }
    )
  }

  # Create joint formula that will be used by inla
  formula <- bru_inla_formula(components)
  formula <- update.formula(formula, BRU_response ~ .)

  # Restore environment
  environment(components) <- env
  environment(formula) <- env

  # Check linearity
  for (lh in seq_along(lhoods)) {
    if (bru_is_linear(lhoods[[lh]])) {
      used_lh <- bru_used(lhoods[[lh]])
      lhoods[[lh]][["pred_expr"]][["is_linear"]] <-
        (length(used_lh$latent) == 0) &&
        all(vapply(
          components[used_lh$effect],
          function(cmp) {
            ibm_is_linear(cmp$mapper)
          },
          TRUE
        ))
    }
    lhoods[[lh]] <- bru_compat_pre_2_14_bru_obs(lhoods[[lh]])
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
      components = summary(as_bru_comp_list(object), ...),
      lhoods = {
        # Once lhoods is moved into bru_model, the NULL part will no longer be
        # needed, and as_bru_obs_list() can be used instead.
        bru_obs_lst <- object[["lhoods"]]
        if (is.null(bru_obs_lst)) {
          NULL
        } else {
          summary(bru_obs_lst, ...)
        }
      }
    ),
    class = "summary_bru_model"
  )
  result
}

#' @export
#' @param x An object to be printed
#' @rdname bru_model
print.summary_bru_model <- function(x, ...) {
  cat("Latent components:\n")
  print(x[["components"]])
  if (!is.null(x[["lhoods"]])) {
    cat("Observation models:\n")
    print(x[["lhoods"]])
  }
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
#' @param predictor A formula or a [bru_pred_expr] expression to be evaluated
#'   given the posterior or for each sample thereof. The default (`NULL`)
#'   returns a `data.frame` containing the sampled effects. In case of a formula
#'   the right hand side is used for evaluation.
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
#' When seed != 0, overridden to "1:1:1"
#' @param used A [bru_used()] object, or NULL (default)
#' @param n_pred integer. If provided, scalar predictor results are expanded to
#' vectors of length `n_pred`.
#' @param \dots Additional arguments passed on to `inla.posterior.sample`
#' @details * `evaluate_model` is a wrapper to evaluate model state, A-matrices,
#' effects, and predictor, all in one call.
#'
#' @keywords internal
#' @rdname evaluate_model
evaluate_model <- function(
  model,
  state,
  data = NULL,
  data_extra = NULL,
  input = NULL,
  comp_simple = NULL,
  predictor = NULL,
  format = NULL,
  used = NULL,
  n_pred = NULL,
  ...
) {
  comp_lst <- as_bru_comp_list(model)
  if (inherits(predictor, "bru_pred_expr")) {
    if (!is.null(used)) {
      warning(
        paste0(
          "Overriding `used` argument to `evaluate_model()` ",
          "in favour of `bru_used(predictor)`."
        )
      )
    }
    used <- bru_used(predictor, labels = names(comp_lst))
    predictor <- bru_pred_expr(predictor, format = "formula")
  } else {
    used <- bru_used(used, labels = names(comp_lst))
  }

  if (is.null(state)) {
    stop("Not enough information to evaluate model states.")
  }
  if (is.null(input)) {
    input <- bru_input(
      comp_lst[used$effect],
      data = data
    )
  }
  if (is.null(comp_simple) && !is.null(input)) {
    comp_simple <- ibm_simplify(
      comp_lst[used$effect],
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
evaluate_state <- function(
  model,
  result,
  property = "mode",
  n = 1,
  seed = 0L,
  num.threads = NULL,
  internal_hyperpar = FALSE,
  ...
) {
  # Evaluate random states, or a single property
  if (property == "sample") {
    state <- post.sample.structured(
      result,
      n = n,
      seed = seed,
      num.threads = num.threads,
      ...
    )
  } else if (is.null(result)) {
    state <- list(lapply(
      as_bru_comp_list(model),
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
#' @param \dots Optional additional parameters, e.g. `inla_f`. Normally unused.
#' @param label Option label used for any warning messages, specifying the
#' affected component.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com} and
#' Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname evaluate_effect
#' @keywords internal

evaluate_effect_single_state.bru_mapper <- function(
  component,
  input,
  state,
  ...,
  label = NULL
) {
  values <- ibm_eval(component, input = input, state = state, ...)

  not_ok <- ibm_invalid_output(
    component,
    input = input,
    state = state
  )
  if (any(not_ok)) {
    if (is.null(label)) {
      warning(
        "Inputs for a mapper give some invalid outputs.",
        immediate. = TRUE
      )
    } else {
      warning(
        "Inputs for '",
        label,
        "' give some invalid outputs.",
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
evaluate_effect_single_state.bm_list <- function(
  components,
  input,
  state,
  ...
) {
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
evaluate_effect_single_state.bru_comp_list <- function(
  components,
  input,
  state,
  ...
) {
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
#' `.data_extra.[["name"]]` or `.data_extra.$name`, or via pronouns;
#'  see Details.
#' @param state A list where each element is a list of named latent state
#' information, as produced by [evaluate_state()]
#' @param effects A list where each element is list of named evaluated effects,
#' each computed by [evaluate_effect_single_state.bru_comp_list()]
#' @param predictor Either a formula or [bru_pred_expr] expression
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
#' @details For each component, e.g. "name", the latent state values are
#'   available as `name_latent`, and arbitrary evaluation can be done with
#'   `name_eval(...)`, see [bru_comp_eval()].
#'
#'   The evaluation supports several [rlang::as_data_pronoun()] data masking
#'   pronouns, to access variables from different data sources, and some of
#'   these also have corresponding full objects, with an appended `.` in the
#'   name. The full objects can be passed as arguments to functions.
#'   \describe{
#'   \item{.effect/.effect.}{refers to the `effects` vectors}
#'   \item{.latent/.latent.}{refers to the latent state vectors}
#'   \item{.data/.data.}{refers to the main `data` argument}
#'   \item{.data_extra/.data_extra.}{refers to the `data_extra` argument}
#'   \item{.env}{refers to the evaluation environment of the predictor}
#'   }
#' @return A list or matrix is returned, as specified by `format`
#' @keywords internal
#' @rdname evaluate_predictor
evaluate_predictor <- function(
  model,
  state,
  data,
  data_extra,
  effects,
  predictor,
  used = NULL,
  format = "auto",
  n_pred = NULL
) {
  stopifnot(inherits(model, "bru_model"))
  format <- match.arg(format, c("auto", "matrix", "list"))
  comp_lst <- as_bru_comp_list(model)
  if (inherits(predictor, "bru_pred_expr")) {
    if (!is.null(used)) {
      warning(
        paste0(
          "Overriding `used` argument to `evaluate_predictor()` ",
          "in favour of `bru_used(predictor)`."
        )
      )
    }
    used <- bru_used(predictor, labels = names(comp_lst))
    predictor <- bru_pred_expr(predictor, format = "formula")
  }
  pred.envir <- environment(predictor)
  if (inherits(predictor, "formula")) {
    pred_text <- as.character(predictor)
    pred_text <- pred_text[length(pred_text)]
    predictor <- rlang::parse_expr(pred_text)
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

  used <- bru_used(used, labels = names(comp_lst))

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

  # Rename component states from label to label_latent
  state_names <- as.list(expand_labels(
    names(state[[1]]),
    names(comp_lst),
    suffix = "_latent"
  ))
  names(state_names) <- names(state[[1]])

  # Construct _eval function names
  eval_names <- as.list(expand_labels(
    intersect(names(state[[1]]), names(comp_lst)),
    intersect(names(state[[1]]), names(comp_lst)),
    suffix = "_eval"
  ))
  names(eval_names) <- intersect(names(state[[1]]), names(comp_lst))

  eval_fun_factory <-
    function(.comp, .envir, .enclos) {
      .is_offset <- .comp$main$type %in% c("offset", "const")
      .is_iid <- .comp$main$type %in% c("iid")
      .mapper <- .comp$mapper
      .label <- paste0(.comp$label, "_latent")
      .iid_precision <- paste0("Precision_for_", .comp$label)
      .iid_cache <- list()
      .iid_cache_index <- NULL
      eval_fun <- function(
        main,
        group = NULL,
        replicate = NULL,
        weights = NULL,
        .state = NULL
      ) {
        n_input <- ibm_n_output(
          .mapper[["mappers"]][["core"]][["mappers"]][["main"]],
          input = main
        )
        if (is.null(group)) {
          group <- rep(1, n_input)
        }
        if (is.null(replicate)) {
          replicate <- rep(1, n_input)
        }
        if (!.is_offset && is.null(.state)) {
          .state <- rlang::eval_tidy(
            rlang::parse_expr(.label),
            data = data_mask,
            env = .envir
          )
        }
        .input <- list(
          core = list(
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
        } else {
          # .is_iid, invalid indices give new samples
          # Check for known invalid output elements, based on the
          # initial mapper (subsequent mappers in the component pipe
          # are assumed to keep the same length and validity)
          not_ok <- ibm_invalid_output(
            .mapper[["mappers"]][[1]],
            input = .input[[1]],
            state = .state
          )
          if (any(not_ok)) {
            .cache_state_index <- rlang::eval_tidy(
              rlang::parse_expr(".cache_state_index"),
              data = data_mask,
              env = .envir
            )
            if (!identical(.cache_state_index, .iid_cache_index)) {
              .iid_cache_index <<- .cache_state_index
              .iid_cache <<- list()
            }
            key <- as.character(main[not_ok])
            not_cached <- !(key %in% names(.iid_cache))
            if (any(not_cached)) {
              .prec <- rlang::eval_tidy(
                rlang::parse_expr(.iid_precision),
                data = data_mask,
                env = .envir
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
  eval_list <- list()
  for (nm in names(eval_names)) {
    eval_list[[eval_names[[nm]]]] <-
      eval_fun_factory(
        comp_lst[[nm]],
        .envir = envir,
        .enclos = enclos
      )
  }

  # Remove problematic objects:
  problems <- c(".Random.seed")
  remove(list = intersect(names(envir), problems), envir = envir)

  n <- length(state)
  for (k in seq_len(n)) {
    state_df <- stats::setNames(state[[k]], state_names[names(state[[k]])])
    data_mask <- bru_data_mask(
      list(
        effect = effects[[k]],
        state_df,
        latent = state[[k]],
        data = data,
        data_extra = data_extra,
        eval_list,
        # Keep track of the iteration index so the iid cache can be
        # invalidated
        list(.cache_state_index = k)
      )
    )

    result_ <- rlang::eval_tidy(predictor, data = data_mask, env = envir)
    if (!is.null(n_pred) && is.numeric(result_) && length(result_) == 1) {
      result_ <- rep(result_, n_pred)
    }
    if (k == 1) {
      if (identical(format, "auto")) {
        if (
          (is.vector(result_) && !is.list(result_)) ||
            (is.matrix(result_) && (NCOL(result_) == 1))
        ) {
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
#'  list(core = list(
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
#' \donttest{
#' if (bru_safe_inla() &&
#'   requireNamespace("sn", quietly = TRUE)) {
#'   set.seed(12345L)
#'   data <-
#'     data.frame(
#'       x = runif(5),
#'       idx = seq_len(5)
#'     )
#'   data$y <- data$x + rnorm(5, sd = 0.1)
#'   fit <- bru(
#'     y ~ 0 + x + field(idx, model = "iid", mapper = bm_index(5L)),
#'     family = "gaussian", data = data,
#'     control.family = list(
#'       hyper = list(prec = list(initial = 9, fixed = TRUE))
#'     ),
#'     options = list(control.inla = list(int.strategy = "eb"))
#'   )
#'   pred <- generate(
#'     fit,
#'     newdata = list(),
#'     formula = ~ field_eval(c(seq_len(5), 6, 7, 6)),
#'     n.samples = 1L
#'   )
#' }
#' }
bru_comp_eval <- function(
  main,
  group = NULL,
  replicate = NULL,
  weights = NULL,
  .state = NULL
) {
  stop(paste0(
    "In your predictor expression, use 'mylabel_eval(...)' instead of\n",
    "'bru_comp_eval(...)'.  See ?bru_comp_eval for more information."
  ))
}


#' @include mappers.R

#' @title Mapper methods for model objects
#' @description
#' Methods for the `ibm_as_taylor()` and `ibm_simplify()` methods for
#' [bru] model objects and related classes.
#'
#' @inheritParams bru_mapper_generics
#' @inheritParams ibm_as_taylor
#' @inheritParams ibm_simplify
#'
#' @name bru_model_mapper_methods
#' @rdname bru_model_mapper_methods
NULL

#' @describeIn bru_model_mapper_methods Returns a list (one element per
#'   observation model) of [bm_list] objects, each with one [bm_taylor]
#'   entry for each included component. (autodiff == "pandemic")
#'
#'   If `autodiff` is not "pandemic", returns a list of [bm_taylor] objects, one
#'   for each observation model, with the offset and jacobians evaluated for the
#'   predictor of the observation model, and the component mappers passed on as
#'   `comp_mappers` for the evaluation of the jacobians.
#'   @param comp_mappers A [bm_list()] of mappers, either the original component
#'   mappers, or simplified mappers.
#'
#' @export
#' @param options A [bru_options] options object or a list of options passed
#' on to [bru_options()]
ibm_as_taylor.bru_model <- function(
  mapper,
  input,
  state = NULL,
  ...,
  options = NULL,
  comp_mappers = NULL,
  eval_fun = NULL
) {
  model <- mapper
  stopifnot(inherits(model, "bru_model"))
  options <- bru_call_options(options)
  if (identical(options[["bru_method"]][["autodiff"]], "pandemic")) {
    bru_log_message(
      paste0("Linearise components for each observation model"),
      verbosity = 3
    )
    mappers <-
      lapply(
        input,
        function(inp) {
          ibm_as_taylor(
            as_bru_comp_list(model),
            input = inp[["comp"]],
            state = state,
            inla_f = TRUE,
            ...
          )
        }
      )
  } else {
    bru_log_message(
      paste0("Linearise predictor for each observation model"),
      verbosity = 3
    )
    if (is.null(comp_mappers)) {
      comp_mappers <- ibm_simplify(
        model,
        input = input,
        inla_f = TRUE
      )
    }
    if (is.null(eval_fun)) {
      eval_fun <- bru_eval_fun(model)
    }
    mappers <- ibm_as_taylor(
      as_bru_obs_list(model),
      input = input,
      state = state,
      multi = TRUE,
      ...,
      comp_mappers = comp_mappers,
      eval_fun = eval_fun
    )
  }
  mappers
}

#' @rdname bru_model_mapper_methods
#' @export
ibm_as_taylor.bru_comp_list <- function(mapper, input, state = NULL, ...) {
  comp <- mapper
  included <- parse_inclusion(
    names(comp),
    names(input),
    NULL
  )

  mappers <- ibm_as_taylor(
    as_bm_list(comp[included]),
    input = input[included],
    state = state[included],
    ...
  )

  mappers
}

#' @rdname ibm_as_taylor
#' @export
#'
ibm_as_taylor.bru_comp <- function(mapper, input, state = NULL, ...) {
  bru_log_message(
    paste0("Linearise component '", mapper[["label"]], "'"),
    verbosity = 5
  )
  if (is.null(state)) {
    state <- rep(0, ibm_n(mapper[["mapper"]]))
  }
  ibm_as_taylor(mapper[["mapper"]], input = input, state = state, ...)
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
          as_bru_comp_list(model),
          input = inp[["comp"]],
          state = state,
          ...
        )
      }
    )

  mappers
}

#' @rdname bru_model_mapper_methods
#' @export
ibm_simplify.bru_comp <- function(mapper, input = NULL, state = NULL, ...) {
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
#'
ibm_simplify.bru_comp_list <- function(
  mapper,
  input = NULL,
  state = NULL,
  ...
) {
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
# Methods for the `ibm_as_taylor()` and `ibm_simplify()` methods for
# [bru] model objects and related classes.
#
#' @export
#' @rdname bru_model_mapper_methods
#' @export
#'
ibm_as_taylor.bm_list <- function(mapper, input, state = NULL, ...) {
  label <- names(mapper)
  if (is.null(label)) {
    idx <- seq_along(mapper)
    label <- as.character(seq_along(mapper))
  } else {
    idx <- setNames(names(mapper), names(mapper))
  }
  mappers <- lapply(
    idx,
    function(k) {
      bru_log_message(
        paste0("Linearise component '", label[k], "'"),
        verbosity = 4
      )
      ibm_as_taylor(
        mapper[[k]],
        input[[k]],
        state = state[[k]],
        ...
      )
    }
  )
  names(mappers) <- names(mapper)
  as_bm_list(mappers)
}

#' @rdname bru_model_mapper_methods
#' @export
ibm_simplify.bm_list <- function(mapper, input = NULL, state = NULL, ...) {
  label <- names(mapper)
  if (is.null(label)) {
    idx <- seq_along(mapper)
    label <- as.character(seq_along(mapper))
  } else {
    idx <- setNames(names(mapper), names(mapper))
  }
  mappers <- lapply(
    idx,
    function(k) {
      bru_log_message(
        paste0("Simplify component '", label[k], "'"),
        verbosity = 4
      )
      ibm_simplify(
        mapper[[k]],
        input[[k]],
        state = state[[k]],
        ...
      )
    }
  )
  names(mappers) <- names(mapper)
  as_bm_list(mappers)
}


#' @rdname bru_model_mapper_methods
#' @export
ibm_eval2.bru_model <- function(
  mapper,
  input,
  state = NULL,
  ...,
  options = NULL,
  comp_mappers = NULL,
  eval_fun = NULL
) {
  model <- mapper
  stopifnot(inherits(model, "bru_model"))
  options <- bru_call_options(options)
  if (identical(options[["bru_method"]][["autodiff"]], "pandemic")) {
    bru_log_abort(
      paste0(
        "ibm_eval2<bru_model> not supported for",
        " bru_method$autodiff == 'pandemic'"
      ),
      verbosity = 1
    )
  }
  if (is.null(comp_mappers)) {
    comp_mappers <- ibm_simplify(
      model,
      input = input,
      inla_f = TRUE
    )
  }
  if (is.null(eval_fun)) {
    eval_fun <- bru_eval_fun(model)
  }
  result <- ibm_eval2(
    as_bru_obs_list(model),
    input = input,
    state = state,
    multi = TRUE,
    ...,
    comp_mappers = comp_mappers,
    eval_fun = eval_fun
  )
  result
}

#' @rdname bru_model_mapper_methods
#' @export
ibm_eval.bru_model <- function(
  mapper,
  input,
  state = NULL,
  ...,
  options = NULL,
  comp_mappers = NULL,
  eval_fun = NULL
) {
  model <- mapper
  stopifnot(inherits(model, "bru_model"))
  options <- bru_call_options(options)
  if (identical(options[["bru_method"]][["autodiff"]], "pandemic")) {
    bru_log_abort(
      paste0(
        "ibm_eval2<bru_model> not supported for",
        " bru_method$autodiff == 'pandemic'"
      ),
      verbosity = 1
    )
  }
  if (is.null(comp_mappers)) {
    comp_mappers <- ibm_simplify(
      model,
      input = input,
      inla_f = TRUE
    )
  }
  if (is.null(eval_fun)) {
    eval_fun <- bru_eval_fun(model)
  }
  result <- ibm_eval(
    as_bru_obs_list(model),
    input = input,
    state = state,
    multi = TRUE,
    ...,
    comp_mappers = comp_mappers,
    eval_fun = eval_fun
  )
  result
}

#' @rdname bru_model_mapper_methods
#' @export
ibm_jacobian.bru_model <- function(
  mapper,
  input,
  state = NULL,
  ...,
  options = NULL,
  comp_mappers = NULL,
  eval_fun = NULL
) {
  model <- mapper
  stopifnot(inherits(model, "bru_model"))
  options <- bru_call_options(options)
  if (identical(options[["bru_method"]][["autodiff"]], "pandemic")) {
    bru_log_abort(
      paste0(
        "ibm_eval2<bru_model> not supported for",
        " bru_method$autodiff == 'pandemic'"
      ),
      verbosity = 1
    )
  }
  if (is.null(comp_mappers)) {
    comp_mappers <- ibm_simplify(
      model,
      input = input,
      inla_f = TRUE
    )
  }
  if (is.null(eval_fun)) {
    eval_fun <- bru_eval_fun(model)
  }
  result <- ibm_jacobian(
    as_bru_obs_list(model),
    input = input,
    state = state,
    multi = TRUE,
    ...,
    comp_mappers = comp_mappers,
    eval_fun = eval_fun
  )
  result
}


#' @title Evaluate or sample from a posterior result given a model and locations
#'
#' @description Evaluate
#'
#' @export
#' @param model A [bru] model
#' @param state list of state lists, as generated by [bru_state()]
#' @param data A `list`, `data.frame`, `af`, or `Spatial*DataFrame`, with
#'   coordinates and covariates needed to evaluate the predictor components.
#' @param data_extra Additional data for the predictor evaluation
#' @param input_comp Precomputed inputs list for the components
#' @param comp_mappers Precomputed [bm_list] of full
#'   (`as_bm_list(as_bru_comp_list(model))`) , simplified, or linearised mappers
#'   for the components
#' @param predictor A formula or a [bru_pred_expr] expression to be evaluated
#'   given the posterior or for each sample thereof. The default (`NULL`)
#'   returns a `data.frame` containing the sampled effects. In case of a formula
#'   the right hand side is used for evaluation.
#' @param format character; determines the storage format of predictor output.
#' Available options:
#' * `"auto"` If the first evaluated result is a vector or single-column matrix,
#'   the "matrix" format is used, otherwise "list".
#' * `"matrix"` A matrix where each column contains the evaluated predictor
#' expression for a state.
#' * `"list"` A list where each element contains the evaluated predictor
#' expression for a state.
#' @param used A [bru_used()] object, or NULL (default)
#' @param \dots Additional arguments, unused.
#'
# @keywords internal
#' @rdname bru_eval
#' @examples
#' if (bru_safe_inla()) {
#'   model <- bru_model(
#'     as_bru_comp_list(~ 1 + x + field(idx,
#'       model = "iid",
#'       mapper = bm_index(5L)
#'     )),
#'     bru_obs(
#'       y ~ Intercept + x + field,
#'       family = "poisson",
#'       data = data.frame(
#'         y = rpois(5, 10), x = rnorm(5),
#'         idx = c(1, 2, 3, 4, 5)
#'       )
#'     )
#'   )
#'   bru_eval(
#'     model,
#'     state = list(bru_state(model, property = "zeros")),
#'     data = data.frame(x = rnorm(4), idx = 4:1)
#'   )
#'   bru_eval(
#'     model,
#'     state = list(list(
#'       Intercept = 1, x = 2, field = c(0, 0, 0, 0, 0),
#'       Precision_for_field = 1
#'     )),
#'     data = data.frame(x = rnorm(4), idx = 4:1)
#'   )
#'   bru_eval(
#'     model,
#'     state = list(list(
#'       Intercept = 1, x = 2, field = c(0, 0, 0, 0, 0),
#'       Precision_for_field = 1
#'     )),
#'     data = data.frame(x = rnorm(4), idx = 4:1),
#'     predictor = new_bru_pred_expr(
#'       ~ list(a = Intercept + x + field, b = field_eval(c(5:7, 5:7)))
#'     )
#'   )
#' }
#'
bru_eval <- function(
  model,
  state,
  data = NULL,
  data_extra = NULL,
  input_comp = NULL,
  comp_mappers = NULL,
  predictor = NULL,
  format = NULL,
  used = NULL,
  ...
) {
  format <- match.arg(format, c("auto", "matrix", "list"))
  comp_lst <- as_bru_comp_list(model)
  if (inherits(predictor, "bru_pred_expr")) {
    if (!is.null(used)) {
      warning(
        paste0(
          "Overriding `used` argument to `bru_eval()` ",
          "in favour of `bru_used(predictor)`."
        )
      )
    }
    used <- bru_used(predictor, labels = names(comp_lst))
  } else {
    used <- bru_used(used, labels = names(comp_lst))
  }

  if (is.null(state)) {
    stop("Not enough information to evaluate model states.")
  }
  if (is.null(input_comp)) {
    input_comp <- bru_input(
      comp_lst[used$effect],
      data = data
    )
  }
  if (is.null(comp_mappers) && !is.null(input_comp)) {
    comp_mappers <- ibm_simplify(
      comp_lst[used$effect],
      input = input_comp,
      inla_f = TRUE
    )
  }

  if (!is.null(predictor)) {
    pred_quo <- bru_pred_expr(predictor, format = "quo")
    expr_mapper <- bm_expr(
      pred_quo,
      labels = list(
        root = "latent",
        derived = "effects",
        suffix = "_latent"
      ),
      assume = c(
        if (isTRUE(bru_is_rowwise(predictor))) "rowwise" else NULL,
        if (isTRUE(bru_is_linear(predictor))) "linear" else NULL,
        if (isTRUE(bru_is_additive(predictor))) "additive" else NULL,
        if (isTRUE(length(used[["latent"]]) == 0L)) "no_root" else NULL
      )
    )
  }

  nms <- unique(c(used$effect, used$latent))
  param_nms <- setdiff(names(state[[1]]), nms)
  state_nms <- setdiff(names(state[[1]]), param_nms)
  effect_nms <- intersect(nms, used$effect)

  eval_fun <- bru_eval_fun(model)
  derived <- list()

  values <- list()
  for (k in seq_along(state)) {
    state_k <- state[[k]]
    # Evaluate component effects
    for (nm in effect_nms) {
      res <- ibm_eval(
        comp_mappers[[nm]],
        input = input_comp[[nm]],
        state = state_k[[nm]],
        ...
      )
      derived[[nm]] <- res
    }
    if (is.null(predictor)) {
      values[[k]] <- derived
    } else {
      # Collect context data, including hyperparameters, if any
      data_list <- list(
        param = state_k[param_nms],
        data = data,
        data_extra = data_extra,
        eval_fun = eval_fun,
        # Keep track of the iteration index so the iid cache can be
        # invalidated
        list(.cache_state_index = k)
      )
      # Compute the expression
      values_ <-
        ibm_eval(
          expr_mapper,
          input = list(),
          state = state_k[state_nms],
          derived = derived,
          data = data_list
        )

      if (k == 1) {
        if (identical(format, "auto")) {
          if (
            (is.vector(values_) && !is.list(values_)) ||
              (is.matrix(values_) && (NCOL(values_) == 1))
          ) {
            format <- "matrix"
          } else {
            format <- "list"
          }
        }
        if (identical(format, "matrix")) {
          values <- matrix(0.0, NROW(values_), length(state))
          rownames(values) <- row.names(as.matrix(values_))
        } else if (identical(format, "list")) {
          values <- vector("list", length(state))
        }
      }
      if (identical(format, "list")) {
        values[[k]] <- values_
      } else {
        values[, k] <- values_
      }
    }
  }

  values
}

#' @title Extract model state properties or samples
#' @description Extracts model state summary values or generate samples
#' @param x Object
#' @param property Property of the model components to obtain value from.
#' Default: "mode". Other options are "mean", "0.025quant", "0.975quant",
#' "sd" and "sample", "zeros", "joint_mode", "predictor_sd".
#' In case of "sample" you will obtain samples from the
#' posterior. For "zeros", all-zero vectors are
#' returned for each component.
#' @param n The number of samples to generate. Default: 1L
#' @param seed Integer seed for random number generation. Default: 0L
#' @param num.threads `num.threads` setting for INLA. Default: NULL, which
#' leaves it up to INLA.
#' @param internal_hyperpar logical; If `TRUE`, return hyperparameter properties
#' on the internal scale. Currently ignored when `property="sample"`.
#' Default is `FALSE`.
#' @param \dots Further arguments, passed on to `post.sample.structured()` when
#'   `property="sample"`.
#' @returns A list of length `n`
#' @export
#' @rdname bru_state
# @keywords internal
#' @examples
#' if (bru_safe_inla()) {}
bru_state <- function(x, ...) {
  UseMethod("bru_state", x)
}

#' @rdname bru_state
#' @export
bru_state.bru_model <- function(x, ...) {
  state <- list(lapply(
    as_bru_comp_list(x),
    function(xx) {
      rep(0.0, ibm_n(xx[["mapper"]]))
    }
  ))

  state
}

#' @rdname bru_state
#' @export
bru_state.bru <- function(
  x,
  property = "mode",
  n = 1,
  seed = 0L,
  num.threads = NULL,
  internal_hyperpar = FALSE,
  ...
) {
  # Evaluate random states, or a single property
  property <- match.arg(
    property,
    c(
      "sample",
      "zeros",
      "mode",
      "mean",
      "0.025quant",
      "0.975quant",
      "sd",
      "joint_mode",
      "predictor_sd"
    )
  )
  if (property == "sample") {
    state <- post.sample.structured(
      x,
      n = n,
      seed = seed,
      num.threads = num.threads,
      ...
    )
    return(state)
  }
  if (property == "zeros") {
    return(bru_state(as_bru_model(x)))
  }
  state <- list(extract_property(
    result = x,
    property = property,
    internal_hyperpar = internal_hyperpar
  ))

  state
}
