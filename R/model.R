#' Internal [inlabru] model structure
#'
#' See [make.model].
#'
#' @name bru_model
#' @keywords internal
NULL


#' Create an inlabru model object from a component formula
#'
#' The [inlabru] syntax for model formulae is different from what
#' `INLA::inla` considers a valid.
#' In inla most of the effects are defined by adding an `f(...)` expression to the formula.
#' In [inlabru] the `f` is replaced by an arbitrary (exception: `offset`) string that will
#' determine the label of the effect. See Details for further information.
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
#' A disadvantage of the inla way is that there is no clear separation between the name of the covariate
#' and the label of the effect. Furthermore, for some models like SPDE it is much more natural to
#' use spatial coordinates as covariates rather than an index into the SPDE vertices. For this purpose
#' [inlabru] provides the new `main` agument. For convenience, the `main` argument can be used
#' like the first argument of the f function, e.g., and is the first argument of the component definition.
#'
#' `y ~ f(temperature, model = 'linear')`
#'
#' is equivalent to
#'
#' `y ~ temperature(temperature, model = 'linear')`
#' and
#' `y ~ temperature(main = temperature, model = 'linear')`
#' as well as
#' `y ~ temperature(model = 'linear')`
#' which sets `main = temperature`.
#'
#' On the other hand, map can also be a function mapping, e.g the [coordinates] function of the
#' [sp] package :
#'
#' `y ~ mySPDE(coordinates, ...)`
#'
#' This exctract the coordinates from the data object, and maps it to the latent
#' field via the information given in the `mapper`, which by default is
#' extracted from the `model` object, in the case of `spde` model
#' objects.
#'
#' Morevover, `main` can be any expression that evaluates within your data as an environment.
#' For instance, if your data has columns 'a' and 'b', you can create a fixed effect of 'sin(a+b)' by
#' setting `map` in the following way:
#'
#' `y ~ myEffect(sin(a+b))`
#'
#'
#' @export
#' @param components A component specification formula
#' @param lhoods A list of one or more `lhood` objects
#' @return A [bru_model] object
#' @keywords internal

make.model <- function(components, lhoods) {

  # Automatically add Intercept and -1 to components unless -1 is in components formula
  components <- auto.intercept(components)

  # Back up environment
  env <- environment(components)

  # Create effects
  effects <- component(components, lhoods)

  # Create joint formula that will be used by inla
  formula <- BRU_response ~ -1
  for (fm in lapply(effects, function(eff) {
    eff$inla.formula
  })) {
    formula <- update.formula(formula, fm)
  }

  # Restore environment
  environment(formula) <- env


  # Make model
  mdl <- list(effects = effects, formula = formula)
  class(mdl) <- c("bru_model", "list")
  return(mdl)
}





#' Evaluate or sample from a posterior result given a model and locations
#'
#' @export
#' @param model An [bru] model
#' @param result Posterior of an [bru] or [lgcp] run.
#' @param data A `list`, `data.frame`, or `Spatial*DataFrame`, with coordinates
#' and covariates needed to evaluate the model.
#' @param predictor A formula or an expression to be evaluated given the
#' posterior or for each sample thereof. The default (`NULL`) returns a
#' `data.frame` containing the sampled effects. In case of a formula the right
#' hand side is used for evaluation.
#' @param property Property of the model components to obtain value from.
#' Default: "mode". Other options are "mean", "0.025quant", "0.975quant",
#' "sd" and "sample". In case of "sample" you will obtain samples from the
#' posterior (see `n` parameter).
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
#' @param \dots Additional arguments passed on to `inla.posterior.sample`
#'
#' @keywords internal
#' @aliases  evaluate
#' @rdname evaluate
evaluate_model <- function(model,
                           result,
                           data,
                           state = NULL,
                           A = NULL,
                           predictor = NULL,
                           property = "mode",
                           format = NULL,
                           n = 1,
                           seed = 0L,
                           num.threads = NULL,
                           ...) {
  if (is.null(state) && !is.null(result)) {
    state <- evaluate_state(model,
      result = result,
      property = property,
      n = n, seed = seed, num.threads = num.threads,
      ...
    )
  }
  if (is.null(state)) {
    stop("Not enough information to evaluate model states.")
  }
  if (is.null(A) && !is.null(data)) {
    A <- amatrix_eval(model$effects, data = data)
  }
  if (is.null(A)) {
    effects <- NULL
  } else {
    effects <- evaluate_effect(model$effects,
      data = data,
      state = state, A = A
    )
  }

  if (is.null(predictor)) {
    return(effects)
  }

  values <- evaluate_predictor(model,
    data = data,
    state = state, effects = effects,
    predictor = predictor, format = format
  )

  values
}



#' @rdname evaluate
evaluate_state <- function(model,
                           result,
                           property = "mode",
                           n = 1,
                           seed = 0L,
                           num.threads = NULL,
                           ...) {
  # Evaluate random states, or a single property
  if (property == "sample") {
    state <- inla.posterior.sample.structured(result,
      n = n, seed = seed,
      num.threads = num.threads,
      ...
    )
  } else {
    state <- list(extract_property(result, property))
  }

  state
}




#' @export
#' @rdname evaluate
evaluate_effect <- function(...) {
  UseMethod("evaluate_effect")
}

#' Evaluate an component effect
#'
#' Calculates a latent component given some data and the state of the
#' component's internal random variables.
#'
#' TODO: Improve speed for iterated calls by making 'mapped' a parameter
#'
#' @export
#' @keywords internal
#' @param component An component.
#' @param data A `data.frame` or Spatial* object of covariates and/or point locations.
#' @param state A numeric vector of the latent variable state for a component
#' (for `evaluate_effect.component`) or
#' a list where each element is named list of such state vectors
#' (for `evaluate_effect.component_list`).
#' @param A A matrix overriding the default projection matrix or matrices
#' (named list of matrices for `evaluate_effect.component_list`)
#' @param ... Unused.
#' @return A numeric vector of the component effect values
#' state.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com} and
#' Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname evaluate

evaluate_effect.component <- function(component, data, state, A = NULL, ...) {
  # Make A-matrix (if not provided)
  if (is.null(A)) {
    A <- amatrix_eval(component, data)
  }

  # Determine component depending on the type of latent model
  if (component$main$type %in% c("offset")) {
    values <- input_eval(component, part = "main", data = data, env = component$env)
  } else {
    values <- A %*% state
  }

  as.matrix(values)
}


#' @export
#' @rdname evaluate
#' @keywords internal
evaluate_effect_single_state <- function(state, components, data,
                                         A = NULL, ...) {
  if (is.null(A)) {
    A <- amatrix_eval(components, data = data)
  }
  result <- list()
  for (label in names(components)) {
    result[[label]] <- evaluate_effect(
      components[[label]],
      data = data,
      state = state[[label]],
      A = A[[label]]
    )
  }
  as.data.frame(result)
}

#' @export
#' @rdname evaluate
#' @keywords internal
evaluate_effect.component_list <- function(components, data, state,
                                           A = NULL, ...) {
  if (is.null(A)) {
    A <- amatrix_eval(components, data = data)
  }
  lapply(
    state,
    evaluate_effect_single_state,
    components = components,
    data = data,
    A = A,
    ...
  )
}




#' Evaluate component effects or expressions
#'
#' Evaluate component effects or expressions, based on a bru model and one or
#' several states of the latent variables and hyperparameters.
#'
#' @param data A `list`, `data.frame`, or `Spatial*DataFrame`, with coordinates
#' and covariates needed to evaluate the model.
#' @param state A list where each element is a named list of latent state
#' information, as produced by [evaluate_state()]
#' @param effects A list where each element is named list of evaluated effects,
#' as computed by [evaluate_effect.component_list()]
#' @param predictor Either a formula or expression
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
#' @return A list is returned, as specified by `format`
#' @keywords internal
#' @rdname evaluate
evaluate_predictor <- function(model,
                               data,
                               state,
                               effects,
                               predictor,
                               format = "auto") {
  stopifnot(inherits(model, "bru_model"))
  format <- match.arg(format, c("auto", "data.frame", "matrix", "list"))
  pred.envir <- environment(predictor)
  if (inherits(predictor, "formula")) {
    predictor <- parse(text = as.character(predictor)[length(as.character(predictor))])
  }

  common_vars <- c(
    if (is.list(data)) {
      data
    } else {
      as.data.frame(data)
    },
    as.list(pred.envir),
    as.list(environment(model$formula))
  )

  n <- length(state)
  for (k in seq_len(n)) {
    state_k <- state[[k]]
    # Rename component states from label to label_latent
    names(state_k) <- expand_labels(
      names(state_k),
      names(model$effects),
      suffix = "_latent"
    )

    envir <- c(
      effects[[k]],
      state_k,
      common_vars
    )

    result_ <- eval(predictor, envir = envir)
    if (k == 1) {
      if (identical(format, "auto")) {
        if (is.vector(result_) ||
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
