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

  for (lh in seq_along(lhoods)) {
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
        bru_obs_lst <- as_bru_obs_list(object)
        summary(bru_obs_lst, ...)
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
#' @param comp_mappers A [bm_list()] of mappers, either the original component
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
ibm_eval2.bm_list <- function(
  mapper,
  input,
  state = NULL,
  ...
) {
  result <- list()
  for (nm in names(mapper)) {
    result[[nm]] <- ibm_eval2(
      mapper[[nm]],
      input = input[[nm]],
      state = state[[nm]],
      ...
    )
  }

  result
}

#' @rdname bru_model_mapper_methods
#' @export
ibm_eval.bm_list <- function(
  mapper,
  input,
  state = NULL,
  ...
) {
  result <- list()
  for (nm in names(mapper)) {
    result[[nm]] <- ibm_eval(
      mapper[[nm]],
      input = input[[nm]],
      state = state[[nm]],
      ...
    )
  }

  result
}

#' @rdname bru_model_mapper_methods
#' @export
ibm_jacobian.bm_list <- function(
  mapper,
  input,
  state = NULL,
  ...
) {
  result <- list()
  for (nm in names(mapper)) {
    result[[nm]] <- ibm_jacobian(
      mapper[[nm]],
      input = input[[nm]],
      state = state[[nm]],
      ...
    )
  }

  result
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
ibm_eval.bru_comp_list <- function(
  mapper,
  input,
  state = NULL,
  ...,
  comp_mappers = NULL
) {
  if (is.null(comp_mappers)) {
    comp_mappers <- ibm_simplify(
      mapper,
      input = input,
      inla_f = TRUE
    )
  }
  result <- ibm_eval(
    comp_mappers,
    input = input,
    state = state,
    ...
  )
  result
}

#' @rdname bru_model_mapper_methods
#' @export
ibm_eval.bru_comp <- function(
  mapper,
  input,
  state = NULL,
  ...
) {
  comp_mapper <- ibm_simplify(
    mapper,
    input = input,
    inla_f = TRUE
  )
  result <- ibm_eval(
    comp_mapper,
    input = input,
    state = state,
    ...
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
#' @returns A `matrix` or `list` (see `format`) of evaluated predictor
#'   expressions.
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
#' @param state A list of initial state vectors, one for each component. If
#' `NULL`, all-zero vectors are returned for each component.
#' @export
bru_state.bru_model <- function(x, ..., state = NULL) {
  comp_lst <- as_bru_comp_list(x)
  if (is.null(state)) {
    state <- list()
  }
  state <- state[intersect(names(comp_lst), names(state))]
  for (lab in names(comp_lst)) {
    if (is.null(state[[lab]])) {
      state[[lab]] <- rep(0, ibm_n(comp_lst[[lab]][["mapper"]]))
    } else if (length(state[[lab]]) == 1) {
      state[[lab]] <-
        rep(
          state[[lab]],
          ibm_n(comp_lst[[lab]][["mapper"]])
        )
    }
  }

  list(state)
}

#' @rdname bru_state
#' @export
bru_state.inla <- function(
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
  state <- list(extract_property(
    result = x,
    property = property,
    internal_hyperpar = internal_hyperpar
  ))

  state
}

#' @rdname bru_state
#' @export
bru_state.bru <- function(
  x,
  property = "mode",
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
  if (property == "zeros") {
    return(bru_state(as_bru_model(x)))
  }
  # Invoke method for "inla"
  NextMethod()
}
