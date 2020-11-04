#' Internal [inlabru] model structure
#'
#' See [make.model].
#'
#' @name bru_model
#' @keywords internal
NULL


# GENERICS ----

evaluate <- function(...) {
  UseMethod("evaluate")
}

# Constructor ----

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
#' @param data Locations and covariates needed to evaluate the model.
#' @param predictor A formula or an expression to be evaluated given the
#' posterior or for each sample thereof. The default (`NULL`) returns a
#' `data.frame` containing the sampled effects. In case of a formula the right
#' hand side is used for evaluation.
#' @param property Property of the model components to obtain value from.
#' Default: "mode". Other options are "mean", "0.025quant", "0.975quant",
#' "sd" and "sample". In case of "sample" you will obtain samples from the
#' posterior (see `n` parameter).
#' @param n Number of samples to draw.
#' @param seed If seed != 0L, the random seed
#' @param num.threads Specification of desired number of threads for parallel
#' computations. Default NULL, leaves it up to INLA.
#' When seed != 0, overridden to "1:1"
#' 
#' @keywords internal
evaluate.model <- function(model,
                           result,
                           data,
                           predictor = NULL,
                           property = "mode",
                           n = 1,
                           seed = 0L,
                           num.threads = NULL) {
  
  if (inherits(predictor, "formula")) {
    fml.envir <- as.list(environment(predictor))
    predictor <- parse(text = as.character(predictor)[length(as.character(predictor))])
  } else {
    fml.envir <- list()
  }

  # Do we obtain our values from sampling or from a property of a summary?
  if (property == "sample") {
    smp <- inla.posterior.sample.structured(result, n = n, seed = seed,
                                            num.threads = num.threads)
  } else {
    result$model <- model
    smp <- rep(list(extract.summary(result, property)), n)
  }

  # Which effects do we want? Remove effect from model that are not required for the evaluation
  if (is.null(predictor)) {
    vars <- setdiff(names(smp[[1]]), "Predictor")
  } else {
    effs <- model$effects
    model$effects <- effs[intersect(names(effs), all.vars(predictor))]
    vars <- intersect(names(smp[[1]]), all.vars(predictor))
  }

  if (!is.null(data)) {
    # Pre-calculate projection matrices
    As <- lapply(model$effects, amatrix_eval, data)
  }

  for (k in 1:n) {
    # Discard variables we do not need
    sm <- smp[[k]][vars]

    # Evaluate effects. Note that the expression may only contain hyper
    # parameters in which case there are no effects to evaluate.
    enm <- intersect(names(sm), names(model$effects))

    for (label in enm) {
      if (is.data.frame(sm[[label]])) {
        sm[[label]] <- sm[[label]]$value
      }
      if (!is.null(data)) {
        sm[[label]] <- value(model$effects[[label]], data = data,
                             state = sm[[label]], A = As[[label]])
      }
    }

    # If no predictor is provided simply return the samples.
    # Otherwise evaluate predictor with each sample as an environment
    if (is.null(predictor)) {
      if (is.null(data)) {
        stop(paste0("Both data and predictor are NULL.",
                    " Don't know what to predict."))
      } else {
        smp[[k]] <- data.frame(sm)
      }
    } else {
      # If data is NULL, the latent variables will be present in sm _instead_
      # of the effects, under the same name.
      # TODO: make latent variables consistently available, under special
      # names.
      envir <- c(sm, as.list(data.frame(data)), fml.envir,
                 as.list(environment(model$formula)))
      smp[[k]] <- eval(predictor, envir = envir)
    }
  }

  # Return
  if (property == "sample") {
    smp
  } else {
    smp[[1]]
  }
}
