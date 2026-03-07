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
#' @param \dots additional arguments affecting the samples produced.
#' @return The form of the value returned by `generate()` depends on the data
#' class and prediction formula. Normally, a data.frame is returned, or a list
#' of data.frames (if the prediction formula generates a list)
#' @example inst/examples/generate.bru.R

generate <- function(object, ...) {
  UseMethod("generate")
}

#' Methods for bru_info objects
#'
#' The `bru_info` class is used to store metadata about `bru` models.
#'
#' @export
#' @rdname bru_info
bru_info <- function(...) {
  UseMethod("bru_info")
}

#' @export
#' @param method character; The type of estimation method used
#' @param \dots Additional arguments to be stored in the `bru_info` object. For
#' `summary` and `print` methods, arguments passed on to submethods.
#' @param inlabru_version character; inlabru package version. Default: NULL, for
#' automatically detecting the version
#' @param INLA_version character; INLA package version. Default: NULL, for
#' automatically detecting the version
#' @describeIn bru_info Create a `bru_info` object
bru_info.character <- function(method,
                               ...,
                               inlabru_version = NULL,
                               INLA_version = NULL) {
  if (is.null(inlabru_version)) {
    inlabru_version <-
      tryCatch(
        expr = {
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
        expr = {
          getNamespaceVersion("INLA")
        },
        error = function(e) {
          "unknown"
        }
      )
  }
  object <- structure(
    list(
      method = method,
      ...,
      inlabru_version = inlabru_version,
      INLA_version = INLA_version
    ),
    class = "bru_info"
  )
  object
}

#' @export
#' @describeIn bru_info Extract the `bru_info` object from an estimated [bru()]
#'   result object. The default print method shows information about model
#'   components and observation models.
bru_info.bru <- function(object, ...) {
  object <- bru_check_object_bru(object)
  object[["bru_info"]]
}

#' @export
#' @describeIn bru_info Extract the `bru_info` object from an estimated [bru()]
#'   result object. The default print method shows information about model
#'   components and observation models.
as_bru_info <- function(object, ...) {
  UseMethod("as_bru_info")
}
#' @export
#' @describeIn bru_info Extract a `bru_info` object.
as_bru_info.bru_info <- function(object, ...) {
  object <- bru_info_upgrade(object)
  object
}
#' @export
#' @describeIn bru_info Extract the `bru_info` object from an estimated [bru()]
#'   result object.
as_bru_info.bru <- function(object, ...) {
  object <- bru_check_object_bru(object)
  object[["bru_info"]]
}


#' @title Extract INLA formula
#'
#' @description
#' Extracts the INLA formula used in a `bru` model
#'
#' @param x An object containing information about an INLA formula
#' @param \dots Additional arguments passed on to submethods
#' @export
#' @keywords internal
bru_inla_formula <- function(x, ...) {
  UseMethod("bru_inla_formula")
}
#' @rdname bru_inla_formula
#' @export
bru_inla_formula.bru <- function(x, ...) {
  bru_inla_formula(as_bru_info(x), ...)
}
#' @rdname bru_inla_formula
#' @export
bru_inla_formula.bru_info <- function(x, ...) {
  formula <- x[["model"]][["formula"]]
  if (is.null(formula)) {
    components <- x[["model"]][["effects"]]
    used <- bru_used(x)
    included <- union(used$effect, used$latent)
    formula <- bru_inla_formula(components[included])
    formula <- update.formula(formula, BRU_response ~ .)
  }
  formula
}

#' @rdname bru_inla_formula
#' @export
bru_inla_formula.bru_comp_list <- function(x, ...) {
  formula <- . ~ 0
  for (cmp in seq_along(x)) {
    formula <- update.formula(formula, x[[cmp]]$inla.formula)
  }
  environment(formula) <- environment(x)
  formula
}


#' @export
#' @method summary bru_info
#' @param object Object to operate on
#' @param verbose logical; If `TRUE`, include more details of the
#' component definitions. If `FALSE`, only show basic component
#' definition information. Default: `TRUE`
#' @param \dots Arguments passed on to other `summary` methods
#' @rdname bru_info
summary.bru_info <- function(object, verbose = TRUE, ...) {
  structure(
    list(
      inlabru_version = object[["inlabru_version"]],
      INLA_version = object[["INLA_version"]],
      components = summary(object[["model"]], verbose = verbose, ...),
      lhoods = summary(as_bru_obs_list(object), verbose = verbose, ...)
    ),
    class = "summary_bru_info"
  )
}

#' @export
#' @param x An  object to be printed
#' @rdname bru_info
print.summary_bru_info <- function(x, ...) {
  cat(glue("inlabru version: {x$inlabru_version}"), "\n")
  cat(glue("INLA version: {x$INLA_version}"), "\n")
  print(x$components)
  if (!is.null(x$lhoods)) {
    print(x$lhoods)
  }
  invisible(x)
}

#' @export
#' @rdname bru_info
print.bru_info <- function(x, ...) {
  print(summary(x))
  invisible(x)
}


#' @title Extract timing information from fitted [bru] object
#' @description
#' Extracts a data.frame or tibble with information about the `Time` (CPU),
#' `System`, and `Elapsed` time for each step of a `bru()` run.
#' @param object A fitted `bru` object
#' @param \dots unused
#' @export
bru_timings <- function(object, ...) {
  UseMethod("bru_timings")
}

#' @export
#' @rdname bru_timings
bru_timings.bru <- function(object, ...) {
  object <- bru_check_object_bru(object)
  object[["bru_timings"]]
}


# @param args Either a [bru_obs()] object, a [bru_obs_list()] object, or
#   a list of one or more [bru_obs()] or [bru_obs_list()] objects,
#   or a list of arguments for a call to [bru_obs()].
# Returns: A bru_obs_list object.
bru_obs_list_construct <- function(args, options, .envir = parent.frame(),
                                   .response = . ~ .,
                                   .components = NULL) {
  lhoods <- args
  if (inherits(lhoods, c("bru_obs", "bru_obs_list"))) {
    lhoods <- list(lhoods)
  }
  dot_is_lhood <- vapply(
    lhoods,
    function(lh) inherits(lh, "bru_obs"),
    TRUE
  )
  dot_is_lhood_list <- vapply(
    lhoods,
    function(lh) inherits(lh, "bru_obs_list"),
    TRUE
  )
  if (any(dot_is_lhood | dot_is_lhood_list)) {
    if (!all(dot_is_lhood | dot_is_lhood_list)) {
      stop(paste0(
        "Cannot mix `bru_obs()` parameters with `bru_obs` ",
        "and `bru_obs_list` objects.",
        "\n  Check if the argument(s) ",
        glue_collapse(
          glue("'{names(lhoods)[!(dot_is_lhood | dot_is_lhood_list)]}'"),
          sep = ", ", last = ", and "
        ),
        " were meant to be given to a call to 'bru_obs()',\n  or in the ",
        "'options' list argument instead."
      ))
    }
  } else {
    if (is.null(lhoods[["formula"]])) {
      lhoods[["formula"]] <- . ~ .
    }
    lhoods[["formula"]] <- auto_response(
      lhoods[["formula"]],
      .response
    )
    lhoods <- list(do.call(
      bru_obs,
      c(
        lhoods,
        list(
          options = options,
          .envir = .envir
        )
      )
    ))
    dot_is_lhood <- TRUE
    dot_is_lhood_list <- FALSE
  }
  lhoods <- do.call(c, lhoods[dot_is_lhood | dot_is_lhood_list])

  # Update include/exclude information to limit it to existing components
  # and check additivity
  lhoods <- bru_used_update(lhoods, labels = names(.components))

  lhoods
}


#' @title Convenient model fitting using (iterated) INLA
#'
#' @description This method is a wrapper for `INLA::inla` and provides
#'   multiple enhancements.
#'
#'   \itemize{
#'     \item
#'       Easy usage of spatial covariates and automatic construction of inla
#'       projection matrices for (spatial) SPDE models. This feature is
#'       accessible via the `components` parameter. Practical examples on how to
#'       use spatial data by means of the components parameter can also be found
#'       by looking at the [lgcp()] function's documentation.
#'     \item
#'       Constructing multiple observation models is straightforward. See
#'       [bru_obs()] for more information on how to provide additional
#'       models to `bru` using the `...` parameter list.
#'     \item
#'       Support for non-linear predictors. See example below.
#'     \item
#'       Log Gaussian Cox process (LGCP) inference is
#'       available by using the `"cp"` family or (even easier) by using the
#'       [lgcp()] function.
#'   }
#' @export
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'
#' @param components Latent component definitions, either as a [bru_comp_list()]
#'   object, or a `formula`-like specification.
#'   Also used to define a default linear additive predictor.  See
#'   [bru_comp()] for details.
#' @param \dots Observation models, each constructed by a calling [bru_obs()],
#'   or [bru_obs_list()].
#'
#'   Alternatively, for backwards compatibility, may be named parameters that
#'   can be passed to a single [bru_obs()] call. These arguments will be
#'   evaluated before calling [bru_obs()], in order to detect if they already
#'   are `bru_obs` objects. This means that special arguments that are only
#'   available in the context of `data` or `response_data` (such as `Ntrials`)
#'   will only work properly in direct calls to [bru_obs()].
#' @param .envir Environment for component evaluation (for when a non-formula
#' specification is used)
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
                options = list(),
                .envir = parent.frame()) {
  stopifnot(bru_safe_inla(multicore = TRUE))

  # Update default options
  options <- bru_call_options(options)
  bru_options_set_local(options, .reset = TRUE)

  bru_log_bookmark("bru")
  bru_log_message("bru: Preprocessing", verbosity = 1L)

  timing <- bru_timer_do(NULL, "Preprocess", 0L)

  # Turn model components and bru_obs objects into internal bru model
  bru.model <- bru_model(
    components,
    lhoods = list(...),
    options = options,
    .envir = .envir
  )

  # Set max iterations to 1 if all likelihood formulae are linear
  if (all(vapply(bru.model[["lhoods"]], bru_is_linear, TRUE))) {
    options$bru_max_iter <- 1
  }

  # Until the storage format is changed in a backwards incompatible way,
  # move lhoods and inputs out of bru.model:
  lhoods <- bru.model[["lhoods"]]
  inputs <- bru.model[["inputs"]]
  bru.model[["lhoods"]] <- NULL
  bru.model[["inputs"]] <- NULL

  info <- bru_info(
    method = "bru",
    model = bru.model,
    inputs = inputs,
    lhoods = lhoods,
    log = bru_log()["bru"],
    options = options
  )

  timing <- bru_timer_done(timing)

  # Run iterated INLA
  if (options$bru_run) {
    result <- iinla(
      model = info[["model"]],
      lhoods = as_bru_obs_list(info),
      inputs = info[["inputs"]],
      options = info[["options"]]
    )
  } else {
    result <- list()
  }

  result$bru_timings <-
    rbind(
      timing,
      result[["bru_iinla"]][["timings"]]
    )

  # Add bru information to the result
  result$bru_info <- info
  class(result) <- c("bru", class(result))
  result
}


#' @describeIn bru
#' Continue the optimisation from a previously computed estimate. The estimation
#' `options` list can be given new values to override the original settings.
#'
#' To rerun with a subset of the data (e.g. for cross validation or prior
#' sampling), use [bru_set_missing()] to set all or part of the response data
#' to `NA` before calling `bru_rerun()`.
#' @param result A previous estimation object of class `bru`
#'
#' @export
bru_rerun <- function(result, options = list()) {
  result <- bru_check_object_bru(result)

  info <- result[["bru_info"]]
  info[["options"]] <- bru_call_options(
    bru_options(
      info[["options"]],
      as.bru_options(options)
    )
  )
  bru_options_set_local(info[["options"]], .reset = TRUE)

  orig_timings <- result[["bru_timings"]]

  result <- iinla(
    model = info[["model"]],
    lhoods = as_bru_obs_list(info),
    inputs = info[["inputs"]],
    initial = result,
    options = info[["options"]]
  )

  new_timings <- result[["bru_iinla"]][["timings"]]$Iteration >
    max(orig_timings$Iteration)
  result$bru_timings <-
    rbind(
      orig_timings,
      result[["bru_iinla"]][["timings"]][new_timings, , drop = FALSE]
    )

  # Add bru information to the result
  result$bru_info <- info
  class(result) <- c("bru", class(result))
  result
}

#' @title Set missing values in observation models
#'
#' @description Set all or parts of the observation model response data
#' to `NA`, for example for use in cross validation (with [bru_rerun()])
#' or prior sampling (with [bru_rerun()] and [inlabru::generate()], but see
#' "Prior sampling caveats" below).
#'
#' @param object A `bru`, `bru_obs` or `bru_obs_list` object
#' @param keep For `bru_obs`, a single logical or an integer vector;
#'   If `TRUE`, keep all the response data, if `FALSE` (default),
#'   set all of it to `NA`. An integer vector determines which elements
#'   to keep (for positive values) or to set as missing (negative values).
#'
#'   For `bru` and `bru_obs_list`, a logical scalar or vector, or a list, see
#'   Details.
#' @param \dots Additional arguments passed on to the `bru_obs` and individual
#'   data class methods.  Currently unused.
#'
#' @details For `bru` and `bru_obs_list`,
#' \itemize{
#'   \item{`keep` must be either a single logical, which is expanded to a list,}
#'   \item{a logical vector, which is converted to a list,}
#'   \item{an unnamed list of the same length as the number of observation
#'     models, with elements compatible with the `bru_obs` method, or}
#'   \item{a named list with elements compatible with the `bru_obs` method,
#'     and only the named `bru_obs` models are acted upon, i.e. the elements
#'     not present in the list are treated as `keep = TRUE`.}
#' }
#'
#' E.g.: `keep = list(b = FALSE)` sets all observations in model `b` to missing,
#' and does not change model `a`.
#'
#' E.g.: `keep = list(a = 1:4, b = -(3:5))` keeps only observations `1:4` of
#' model `a`, marking the rest as missing, and sets observations `3:5` of model
#' `b` to missing.
#'
#' @section Prior sampling caveats:
#' Note that prior sampling requires special care for hyperparameters, as the
#' prior modes are not typically useful; in the future, we plan to have a
#' dedicated method that samples from the hyperparameters, and then uses
#' ```
#' bru_rerun(
#'   bru_set_missing(...),
#'   options = list(
#'     control.mode = list(
#'       theta = theta_sample,
#'       fixed = TRUE)
#'     )
#'   )
#' ```
#' for each sample.
#'
#' @export
#' @rdname bru_set_missing
#' @examples
#' obs <- c(
#'   A = bru_obs(y_A ~ ., data = data.frame(y_A = 1:6)),
#'   B = bru_obs(y_B ~ ., data = data.frame(y_B = 11:15))
#' )
#' bru_response_size(obs)
#' lapply(
#'   bru_set_missing(obs, keep = FALSE),
#'   function(x) {
#'     x[["response_data"]][[x[["response"]]]]
#'   }
#' )
#' lapply(
#'   bru_set_missing(obs, keep = list(B = FALSE)),
#'   function(x) {
#'     x[["response_data"]][[x[["response"]]]]
#'   }
#' )
#' lapply(
#'   bru_set_missing(obs, keep = list(1:4, -(3:5))),
#'   function(x) {
#'     x[["response_data"]][[x[["response"]]]]
#'   }
#' )
#'
bru_set_missing <- function(object, keep = FALSE, ...) {
  UseMethod("bru_set_missing")
}
#' @export
#' @param x Object on which to apply `bru_set_missing()`
#' @param value Value to be passed as `keep` to `bru_set_missing()`
#' @describeIn bru_set_missing Setter method for `bru_set_missing()`
`bru_set_missing<-` <- function(x, ..., value) {
  object <- bru_set_missing(x, keep = value, ...)
  object
}
#' @export
#' @rdname bru_set_missing
bru_set_missing.bru <- function(object, keep = FALSE, ...) {
  object <- bru_check_object_bru(object)
  bru_set_missing(object[["bru_info"]], ...) <- keep
  object
}
#' @export
#' @rdname bru_set_missing
bru_set_missing.bru_info <- function(object, keep = FALSE, ...) {
  bru_set_missing(object[["lhoods"]], ...) <- keep
  object
}
#' @export
#' @rdname bru_set_missing
bru_set_missing.bru_obs_list <- function(object, keep = FALSE, ...) {
  if (is.logical(keep)) {
    if (length(keep) == 1L) {
      keep <- rep(list(keep), length(object))
    } else {
      keep <- as.list(keep)
    }
  }
  if (is.null(names(keep))) {
    idxs <- seq_along(object)
    if (length(keep) != length(object)) {
      stop(paste0(
        "When `keep` is not named list or vector, it must have ",
        "length matching the number of observation models."
      ))
    }
  } else {
    idxs <- intersect(names(keep), names(object))
    if (length(setdiff(names(keep), names(object))) > 0) {
      stop(paste0(
        "Some observation models in `keep` were not found: ",
        glue_collapse(glue("'{setdiff(names(keep), names(object))}'"), ", ")
      ))
    }
  }
  for (idx in idxs) {
    bru_set_missing(object[[idx]], ...) <- keep[[idx]]
  }
  object
}
#' @export
#' @rdname bru_set_missing
bru_set_missing.bru_obs <- function(object, keep = FALSE, ...) {
  if (is.null(keep) || isTRUE(keep)) {
    return(object)
  }
  if (isFALSE(keep)) {
    keep <- -seq_len(bru_response_size(object))
  }
  bru_set_missing(
    object[["response_data"]][[object[["response"]]]],
    ...
  ) <- keep
  object
}

#' @export
#' @describeIn bru_set_missing From `> 2.13.0`, set missing values in any
#' object supporting [base::is.na<-()] for positive and negative indices.
bru_set_missing.default <- function(object, keep = FALSE, ...) {
  if (is.null(keep) || isTRUE(keep)) {
    return(object)
  }
  if (isFALSE(keep)) {
    keep <- -seq_len(bru_response_size(object))
  }
  is.na(object) <- -keep
  object
}

#' @export
#' @describeIn bru_set_missing From `> 2.13.0`, handles `data.frame`, tibbles,
#'   including `inla.mdata`.
#' @examplesIf bru_safe_inla()
#' (obs <- INLA::inla.mdata(y = 1:4, X = matrix(1:8, 4, 2)))
#' bru_set_missing(obs, keep = c(1, 4))
#' bru_set_missing(obs) <- -(1:2)
#' obs
#'
bru_set_missing.data.frame <- function(object, keep = FALSE, ...) {
  if (is.null(keep) || isTRUE(keep)) {
    return(object)
  }
  if (isFALSE(keep)) {
    keep <- -seq_len(bru_response_size(object))
  }
  object[-keep, ] <- NA
  object
}

#' @export
#' @describeIn bru_set_missing From `> 2.13.0`, handles `inla.surv`.
#' @examplesIf bru_safe_inla()
#' (obs <- INLA::inla.surv(time = 1:4, event = c(1, 0, 1, 0)))
#' bru_set_missing(obs, keep = c(1, 4))
#'
#' (obs <- INLA::inla.surv(
#'   time = 1:4,
#'   event = c(1, 0, 1, 0),
#'   cure = matrix(1:8, 4, 2)
#' ))
#' bru_set_missing(obs, keep = c(1, 4))
#'
#' (obs <- INLA::inla.surv(
#'   time = 1:4,
#'   event = c(1, 0, 1, 0),
#'   subject = c(1, 1, 2, 1)
#' ))
#' bru_set_missing(obs, keep = c(1, 4))
#'
bru_set_missing.inla.surv <- function(object, keep = FALSE, ...) {
  if (is.null(keep) || isTRUE(keep)) {
    return(object)
  }
  if (isFALSE(keep)) {
    keep <- -seq_len(bru_response_size(object))
  }
  cls <- class(object)
  if (!is.data.frame(object)) {
    # TODO: This block should be function for converting list-`inla.surv` to
    # `data.frame`-`inla.surv`, but inla.surv might standardise
    # to tibble or data.frame, so we should remove this when possible.
    dat <- object
    att <- attributes(dat)
    cure_null <- ("cure" %in% names(dat)) && is.null(dat[["cure"]])
    if (cure_null) {
      dat["cure"] <- NULL
    }
    dat <- tibble::as_tibble(unclass(dat))
    attr(dat, "names.ori") <- att[["names.ori"]]
    class(dat) <- c("inla.surv", class(dat))
    object <- dat
  }
  # Note: by only setting time,lower,upper to NA, printing of the object
  # still works without giving an NA indexing error.
  if ("subject" %in% names(object)) {
    object$time[-keep] <- NA
  } else {
    object$time[-keep] <- NA
    object$lower[-keep] <- NA
    object$upper[-keep] <- NA
  }
  object <- as.list(object)
  class(object) <- c("inla.surv", class(object))
  object
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

bru_data_mask <- function(data,
                          pronouns = NULL,
                          objects = NULL) {
  if (is.null(pronouns) || is.null(objects)) {
    nms <- setdiff(names(data), "")
    if (is.null(pronouns)) {
      pronouns <- nms
    }
    if (is.null(objects)) {
      objects <- nms
    }
  }
  data_orig <- data
  data <- lapply(data, function(x) {
    if (!is.null(x) && !is.list(x)) {
      x <- tibble::as_tibble(x)
    }
    x
  })
  top_envir <- rlang::new_environment()
  if (length(objects) > 0) {
    dot2_nms <- paste0(".", objects, ".")
    names(dot2_nms) <- objects
    for (nm in objects) {
      if (is.null(data_orig[[nm]])) {
        next
      }
      assign(dot2_nms[nm], data_orig[[nm]], envir = top_envir)
    }
  }
  bottom_envir <- top_envir
  for (k in rev(seq_along(data))) {
    if (is.null(data[[k]])) {
      next
    }
    bottom_envir <- rlang::new_environment(data[[k]], parent = bottom_envir)
  }
  mask <- rlang::new_data_mask(bottom_envir, top_envir)
  if (length(pronouns) > 0) {
    dot_nms <- paste0(".", pronouns)
    names(dot_nms) <- pronouns
    for (nm in pronouns) {
      if (is.null(data[[nm]])) {
        next
      }
      mask[[dot_nms[nm]]] <- rlang::as_data_pronoun(data[[nm]])
    }
  }
  mask[[".mask"]] <- rlang::as_data_pronoun(mask)
  class(mask) <- c("bru_data_mask", class(mask))
  mask
}

#' @title Evaluate expressions in data contexts
#' @description Evaluate an expression in a series of data contexts, also making
#'   the objects directly available as names surrounded by ".", stopping when
#'   the expression evaluation completes with no error, as well as `tidy`
#'   evaluation pronouns, e.g. `.data`.
#'
#'   This is an internal inlabru method, not intended for general use.
#' @param input An expression to be evaluated
#' @param data list of data objects in priority order. Named elements will be
#'   available as whole objects `.name.` as well as pronouns `.name` in the
#'   evaluation.
#' @param default Value used if the expression is evaluated as NULL. Default
#' NULL
#' @param .envir The evaluation environment
#' @return The result of expression evaluation
#' @keywords internal
#' @examples
#' # The A values come from the 'data' element, and the B values come from
#' # the 'response_data' element, as that is listed first.
#' bru_eval_in_data_context(
#'   list(A = .data$x, B = x),
#'   list(
#'     response_data = tibble::tibble(x = 1:5),
#'     data = tibble::tibble(x = 1:10)
#'   )
#' )
#' # Both A and B come from the 'data' element, as 'x' is found there,
#' # and the data objects are listed in order of precedence.
#' bru_eval_in_data_context(
#'   list(A = .data$x, B = x),
#'   list(
#'     data = tibble::tibble(x = 1:10),
#'     response_data = tibble::tibble(x = 1:5)
#'   )
#' )
#'
#' @export
bru_eval_in_data_context <- function(input,
                                     data = NULL,
                                     default = NULL,
                                     .envir = parent.frame()) {
  input <- rlang::enquo(input)
  deparse_input <- rlang::expr_text(input)
  if (inherits(data, "bru_data_mask")) {
    mask <- data
  } else {
    mask <- bru_data_mask(data)
  }
  success <- FALSE
  result <- NULL
  result <- try(
    rlang::eval_tidy(input, data = mask, env = .envir),
    silent = TRUE
  )
  success <- !inherits(result, "try-error")
  if (!success) {
    stop(glue::glue(
      "Input '{glue::glue_collapse(deparse_input, sep = '\n')}' could ",
      "not be evaluated."
    ))
  }
  if (is.null(result)) {
    result <- default
  }

  result
}


complete_coordnames <- function(data_coordnames, ips_coordnames) {
  new_coordnames <- character(
    max(
      length(ips_coordnames),
      length(data_coordnames)
    )
  )
  from_data <- which(!(data_coordnames %in% ""))
  new_coordnames[from_data] <- data_coordnames[from_data]
  from_ips <- setdiff(
    which(!(ips_coordnames %in% "")),
    from_data
  )
  new_coordnames[from_ips] <- ips_coordnames[from_ips]

  dummies <- seq_along(new_coordnames)
  dummies <- setdiff(dummies, c(from_data, from_ips))

  new_coordnames[dummies] <-
    glue("BRU_dummy_coordinate_{seq_along(new_coordnames)}")[dummies]

  list(
    data = new_coordnames[seq_along(data_coordnames)],
    ips = new_coordnames[seq_along(ips_coordnames)]
  )
}


# Extend bind_rows to handle XY/XYZ mismatches in one or more sfc columns
# TODO: detect non-point-data and either warn/stop or handle it properly.
# E.g. POLYGON + POINT data combinations currently fail with obscure error.
extended_bind_rows <- function(...) {
  dt <- list(...)
  names_ <- lapply(dt, names)
  is_sfc_ <- lapply(
    dt,
    function(data) {
      vapply(
        data,
        function(x) {
          inherits(x, "sfc")
        },
        TRUE
      )
    }
  )
  sfc_names_ <- unique(unlist(names_)[unlist(is_sfc_)])

  for (nm in sfc_names_) {
    # Which data objects have this column?
    sf_data_idx_ <-
      which(vapply(dt, function(data) !is.null(data[[nm]]), TRUE))

    # Unify CRS
    if (!inherits(dt[[sf_data_idx_[1]]][[nm]], "sfc")) {
      stop(glue(
        "Column '{nm}' is not a simple feature column in all data objects."
      ))
    }
    the_crs <- fm_crs(dt[[sf_data_idx_[1]]][[nm]])
    for (i in sf_data_idx_) {
      if (!inherits(dt[[i]][[nm]], "sfc")) {
        stop(glue(
          "Column '{nm}' is not a simple feature column in all data objects."
        ))
      }
      dt_crs <- fm_crs(dt[[i]][[nm]])
      if (!fm_crs_is_identical(dt_crs, the_crs)) {
        dt[[i]][[nm]] <- fm_transform(dt[[i]][[nm]], crs = the_crs)
      }
    }

    # Individual elements in sf columns may have different XY/XYZ properties,
    # so need to find out if any of them have Z or M, and then extend all others
    # to have Z too. M is handled by fmesher::fm_zm, but not sf::st_zm
    dt_nm <- lapply(dt[sf_data_idx_], function(x) x[[nm]])
    dt_nm <- fmesher::fm_zm(dt_nm)
    for (i in seq_along(sf_data_idx_)) {
      dt[sf_data_idx_[i]][[nm]] <- dt_nm[[i]][[nm]]
    }
  }

  # Technically, we only need this if at least one of the objects is a tibble
  dt <- lapply(dt, function(object) {
    if (!tibble::is_tibble(object)) {
      if (inherits(object, "sf")) {
        # Convert to tibble, preserving geometry column
        geometry_name <- attr(object, "sf_column", exact = TRUE)
        object <- tibble::as_tibble(object)
        object <- sf::st_as_sf(object, sf_column_name = geometry_name)
      } else {
        object <- tibble::as_tibble(object)
      }
    }
    object
  })
  result <- do.call(dplyr::bind_rows, dt)
  if (length(sfc_names_) > 0) {
    result <- sf::st_as_sf(result)
  }

  result
}


bru_get_parse_data <- function(x) {
  withr::with_options(
    list(keep.parse.data = TRUE),
    {
      utils::getParseData(parse(text = x, keep.source = TRUE))
    }
  )
}

# Returns TRUE if any argument is a Spatial* object, and gives deprecation
# messages.
check_sp_data_deprecation <- function(...) {
  .caller <- sys.call(-1)[[1]]
  obj <- list(...)
  nms <- names(obj)
  result <- vapply(
    nms,
    function(nm) {
      if (inherits(obj[[nm]], "Spatial")) {
        if (!is.null(.caller)) {
          lifecycle::deprecate_warn(
            "2.12.0.9023",
            as.character(glue::glue(
              "{deparse(.caller)}({nm} = ",
              "'has deprecated support for `Spatial` input')"
            )),
            I("`sf` input")
          )
        } else {
          lifecycle::deprecate_warn(
            "2.12.0.9023",
            I(as.character(glue::glue(
              "{nm} has deprecated support for `Spatial` input"
            ))),
            I("`sf` input")
          )
        }
        TRUE
      } else {
        FALSE
      }
    },
    logical(1)
  )
  invisible(any(result))
}


bru_agg_data <- function(
  data,
  ips = NULL,
  domain = NULL,
  samplers = NULL,
  options,
  .envir
) {
  if (is.null(ips)) {
    if (!is.null(domain)) {
      ips <- fm_int(
        domain = domain,
        samplers = samplers,
        int.args = options[["bru_int_args"]]
      )
    }
  } else {
    if (is.null(data)) {
      data <- ips
      ips <- NULL
    }
  }
  if (!is.null(ips)) {
    if (!is.null(data)) {
      ips$.block <- as.integer(ips$.block)
      if (".block" %in% names(data)) {
        ips <- dplyr::left_join(ips, data, by = ".block")
      } else {
        if (NROW(data) != max(ips$.block)) {
          stop(glue(
            "`data` has {NROW(data)} rows, but `ips` has ",
            "{max(ips$.block)} blocks. Cannot join."
          ))
        }
        ips <- dplyr::left_join(
          ips,
          dplyr::bind_cols(data, .block = seq_len(NROW(data))),
          by = ".block"
        )
      }
    }
    data <- ips
    ips <- NULL
  }

  data
}


bru_agg_input <- function(input, data_list, .envir) {
  aggregate_input <- bru_eval_in_data_context(
    {{ input }},
    data = data_list,
    default = NULL,
    .envir = .envir
  )
  if (is.null(aggregate_input)) {
    aggregate_input <- list()
  }
  if (is.null(aggregate_input[["block"]])) {
    aggregate_input[["block"]] <- bru_eval_in_data_context(
      .data[[".block"]],
      data = data_list,
      default = NULL,
      .envir = .envir
    )
  }
  if (is.null(aggregate_input[["weights"]])) {
    aggregate_input[["weights"]] <- bru_eval_in_data_context(
      .data[["weight"]],
      data = data_list,
      default = NULL,
      .envir = .envir
    )
  }
  if (is.null(aggregate_input[["n_block"]])) {
    agg_n_block_expr <- rlang::parse_expr(
      "bru_response_size(.response_data.)"
    )
    aggregate_input[["n_block"]] <- bru_eval_in_data_context(
      !!agg_n_block_expr,
      data = data_list,
      default = NULL,
      .envir = .envir
    )
  }
  lapply(c("block", "weights", "n_block"), function(nm) {
    if (is.null(aggregate_input[[nm]])) {
      bru_log_abort(
        glue::glue(
          "Aggregation requested, but `aggregate_input[['{nm}']]` ",
          "evaluates to NULL."
        )
      )
    }
  })

  if (!is.null(aggregate_input[["block_response"]])) {
    blk_resp <- aggregate_input[["block_response"]]
    if (is.null(data_list$response_data[[blk_resp]])) {
      msg <- glue::glue(
        '`block_response` variable "{blk_resp}" ',
        "not found in `response_data`."
      )

      bru_log_abort(msg)
    }
  } else {
    aggregate_input[["block_response"]] <- ".block"
  }
  blk_resp <- aggregate_input[["block_response"]]
  if (!is.null(data_list$response_data[[blk_resp]])) {
    agg_block_resp_expr <- rlang::parse_expr(
      glue::glue('.response_data[["{blk_resp}"]]')
    )
    block_response <- bru_eval_in_data_context(
      !!agg_block_resp_expr,
      data = data_list,
      default = NULL,
      .envir = .envir
    )

    aggregate_input[["block"]] <-
      match(aggregate_input[["block"]], block_response)
  }

  if (is.character(aggregate_input[["block"]])) {
    msg <- paste0(
      "'character' aggregation block information detected.\n",
      "Please use a numeric or integer vector for the block information ",
      "instead.\n",
      "If you want to use a character vector, supply a `block_response`\n",
      "argument with the name of a response variableto match the character ",
      "vector to."
    )

    bru_log_abort(msg)
  }

  if (anyNA(aggregate_input[["block"]])) {
    msg <- paste0(
      "'NA' aggregation block information detected.",
      "Either the `block` information was out of range, or a match was not ",
      "found when using `.response_data$block_response`."
    )

    bru_log_abort(msg)
  }

  aggregate_input
}


bru_obs_agg <- function(lh,
                        aggregate = NULL,
                        aggregate_input = NULL,
                        options = list(),
                        .envir = parent.frame()) {
  if (!is.null(aggregate) && is.character(aggregate)) {
    aggregate <- switch(aggregate,
      "none" = NULL,
      bm_aggregate(type = aggregate)
    )
  }
  if (!is.null(aggregate)) {
    lh$data <- bru_agg_data(
      data = lh$data,
      ips = lh$integration_info$ips,
      domain = lh$integration_info$domain,
      samplers = lh$integration_info$samplers,
      options = options,
      .envir = .envir
    )
    lh$integration_info$ips <- NULL

    aggregate_input <- bru_agg_input(
      {{ aggregate_input }},
      data_list = list(
        data = lh$data,
        response_data = lh$BRU_original_response_data,
        extra_data = lh$data_extra
      ),
      .envir = .envir
    )

    pred_text <- bru_pred_expr(lh, format = "text_raw")
    pred_text <- glue(
      "{{ibm_eval(BRU_aggregate_mapper, input = BRU_aggregate_input,",
      " state = {{{pred_text}}})}}"
    )
    lh$pred_expr <- new_bru_pred_expr(
      pred_text,
      used = lh$pred_expr$used,
      is_rowwise = FALSE,
      .envir = lh$pred_expr$.envir
    )
    lh$pred_expr$is_additive <- FALSE
    lh$pred_expr$is_linear <- FALSE
    lh$data_extra[["BRU_aggregate_mapper"]] <- aggregate
    lh$data_extra[["BRU_aggregate_input"]] <- aggregate_input

    lh <- bru_compat_pre_2_14_bru_obs(lh)
  }

  lh
}


bru_obs_family_cp_sp <- function(lh, options, .envir) {
  if (!is.null(lh[["aggregate"]])) {
    stop("The 'aggregate' feature cannot be used with family='cp'.")
  }

  response <- lh[["response_data"]][[lh[["response"]]]]
  orig_response_data <- lh[["BRU_original_response_data"]]
  response_data <- lh[["response_data"]]
  data <- lh[["data"]]
  formula <- lh[["formula"]]
  domain <- lh[["integration_info"]][["domain"]]
  samplers <- lh[["integration_info"]][["samplers"]]
  ips <- lh[["integration_info"]][["ips"]]

  response_expr_text <- as.character(formula)[2]

  # Catch and handle special cases:
  if (is.null(response) || !inherits(response, "list")) {
    domain_names <- trimws(strsplit(response_expr_text, split = "\\+")[[1]])
    if (!is.null(domain_names)) {
      # "a + b" conversion to list(a = a, b = b)
      domain_expr <- paste0(
        "list(",
        paste0(
          vapply(
            domain_names, function(x) {
              if (identical(x, "coordinates")) {
                glue("{x} = sp::{x}(.data.)")
              } else {
                glue("{x} = {x}")
              }
            },
            ""
          ),
          collapse = ", "
        ),
        ")"
      )
      response_expr <- rlang::parse_expr(domain_expr)
    }
    response <- tryCatch(
      expr = bru_eval_in_data_context(
        !!response_expr,
        data = list(response_data = orig_response_data, data = data),
        default = NULL,
        .envir = .envir
      ),
      error = function(e) {
        NULL
      }
    )
  }

  if (is.null(response)) {
    stop(paste0(
      "You called bru_obs() with family='cp' but the evaluated ",
      "response information is NULL"
    ))
  }

  if (is.null(ips)) {
    if (is.null(domain)) {
      stop(paste0(
        "The family='cp' model requires a 'domain' specification compatible ",
        "with 'fmesher::fm_int()'"
      ))
    }
    if (!setequal(names(response), names(domain))) {
      stop(glue("
          Mismatch between response and domain names:
            names(response) = ({glue_collapse(names(response), sep = ', ')})
            names(domain)   = ({glue_collapse(names(domain), sep = ', ')})"))
    }

    ips <- fm_int(
      domain = domain,
      samplers = samplers,
      int.args = options[["bru_int_args"]]
    )
    if ((inherits(samplers, "Spatial") ||
      inherits(data, "Spatial") ||
      inherits(response[["coordinates"]], "Spatial")) &&
      inherits(ips, "sf")) {
      ips <- sf::as_Spatial(ips)
    }
  }

  if (length(unique(response_data[["BRU_E"]])) > 1) {
    bru_log_warn(
      "Exposure/effort parameter E should be a scalar for likelihood 'cp'."
    )
  }

  ips_is_Spatial <- inherits(ips, "Spatial")
  if (ips_is_Spatial) {
    bru_safe_sp(force = TRUE)
    ips_coordnames <- sp::coordnames(ips)
    ips_crs <- fm_CRS(ips)
    # For backwards compatibility:
    data_crs <- fm_CRS(data)

    if ("coordinates" %in% names(response)) {
      data_coordnames <- colnames(response$coordinates)
      new_coordnames <- complete_coordnames(data_coordnames, ips_coordnames)
      colnames(response$coordinates) <- new_coordnames$data
      sp::coordnames(ips) <- new_coordnames$ips
    } else {
      if (inherits(response, c("sf", "sfc")) ||
        (is.list(response) &&
          any(vapply(
            response,
            function(x) inherits(x, c("sf", "sfc")), TRUE
          )))) {
        ips <- sf::st_as_sf(ips)
        ips_is_Spatial <- FALSE
      }
    }
  }
  # TODO: check that the crs info is the same

  # For non-Spatial models:
  # Use the response data list as the actual data object, since that's now the
  # canonical place where the point information is given.  This also allows
  # response_data to be used when constructing the response list. This makes
  # it a strict requirement that the predictor can be evaluated as a pure
  # function of the domain data.  When implementing sf support, might be able
  # to give explicit access to spatial coordinates, but otherwise the user can
  # extract it from the geometry with st_coordinates(geometry)[,1] or similar.
  # For Spatial models, keep the old behaviour for backwards compatibility for
  # now, but can likely realign that in the future after more testing.
  # Save general response data to add to response (precomputed covariates etc)
  if (!is.null(orig_response_data)) {
    data_ <- orig_response_data
  } else {
    data_ <- data
  }
  if (inherits(data_, "Spatial")) {
    data_ <- tibble::as_tibble(as.data.frame(data_))
  }
  if (ips_is_Spatial) {
    if ("coordinates" %in% names(response)) {
      idx <- names(response) %in% "coordinates"
      data <- as.data.frame(response$coordinates)
      if (!all(idx)) {
        data <- dplyr::bind_cols(data, as.data.frame(response[!idx]))
      }
    } else {
      data <- as.data.frame(response)
    }
  } else {
    data <- tibble::as_tibble(response)
    if (("geometry" %in% names(data)) &&
      inherits(data$geometry, "sfc")) {
      sf::st_geometry(data) <- "geometry"
    }
  }
  orig_response_data <- NULL
  n_block <- max(ips$.block)
  if (".block" %in% names(data_)) {
    n_block <- max(max(data_$.block), n_block)
    N_data <- base::tabulate(data_$.block, n_block)
  } else {
    N_data <- NROW(data)
  }

  # Add back additional data
  additional_data_names <- setdiff(names(data_), names(data))
  if ((length(additional_data_names) > 0) &&
    (NROW(data_) == sum(N_data))) {
    data <- cbind(data, tibble::as_tibble(data_)[additional_data_names])
  }
  if (is.null(data[[".block"]])) {
    data$.block <- 1L
  }

  if (ips_is_Spatial) {
    ips <- as.data.frame(ips)
  } else {
    if ("geometry" %in% names(ips)) {
      sf::st_geometry(ips) <- "geometry"
    }
  }

  # Use 'weights' for per-point weighting of eta
  if (length(response_data[["BRU_weights"]]) == 1L) {
    point_weights <- rep(response_data[["BRU_weights"]], sum(N_data))
  } else {
    stopifnot(length(response_data[["BRU_weights"]]) == sum(N_data))
    point_weights <- response_data[["BRU_weights"]]
  }
  response_data[["BRU_weights"]] <- 1L

  if (identical(options[["bru_compress_cp"]], TRUE)) {
    new_response_data <- tibble::tibble(
      BRU_E = c(
        rep(0, length(N_data[N_data > 0])),
        response_data[["BRU_E"]] * ips[["weight"]]
      ),
      BRU_response = c(
        N_data[N_data > 0],
        rep(0, NROW(ips))
      ),
      BRU_block = c(
        which(N_data > 0),
        ips[[".block"]]
      ),
      BRU_weights = response_data[["BRU_weights"]][1],
      BRU_Ntrials = response_data[["BRU_Ntrials"]][1],
      BRU_scale = response_data[["BRU_scale"]][1]
    )

    lh$data_extra[["BRU_cp_block_subset"]] <- which(N_data > 0)
    lh$data_extra[["BRU_cp_n_block"]] <- n_block
    lh$data_extra[["BRU_cp_agg"]] <-
      bm_aggregate(
        type = "average",
        n_block = n_block
      )

    pred_text <- bru_pred_expr(lh, format = "text_raw")
    pred_text <- glue("
        {{
          BRU_eta <- {{
            {pred_text}
          }}
          BRU_eta <- rep_len(BRU_eta, length(BRU_aggregate))
          c(ibm_eval(
              BRU_cp_agg,
              list(block = .block[BRU_aggregate]),
              state = BRU_eta[BRU_aggregate] *
                BRU_point_weights[BRU_aggregate])[BRU_cp_block_subset],
            BRU_eta[!BRU_aggregate])
        }}")
    old_pred_expr <- bru_pred_expr(lh)
    lh$pred_expr <- new_bru_pred_expr(
      pred_text,
      used = lh$pred_expr$used,
      is_rowwise = FALSE,
      .envir = lh$pred_expr$.envir
    )
    lh$pred_expr$is_additive <- FALSE
    lh$pred_expr$is_linear <- old_pred_expr$is_linear

    data <- extended_bind_rows(
      dplyr::bind_cols(data,
        BRU_aggregate = TRUE,
        BRU_point_weights = point_weights
      ),
      dplyr::bind_cols(ips,
        BRU_aggregate = FALSE,
        BRU_point_weights = 0.0
      )
    )

    if (!is.null(lh$control.gcpo)) {
      group_cv_block <- c(which(N_data > 0), ips$.block)
      # Would like:
      # group_cv_friends <- lapply(seq_len(max(ips$.block)), function(i) {
      #   which(group_cv_block == i)
      # })
      # Current inla.group.cv interface (2025-08-28)
      group_cv_friends <- lapply(seq_along(group_cv_block), function(i) {
        which(group_cv_block == group_cv_block[i])
      })

      lh$control.gcpo <- modifyList(
        lh$control.gcpo,
        list(friends = group_cv_friends)
      )
    }
  } else {
    new_response_data <- data.frame(
      BRU_E = c(
        rep(0, sum(N_data)),
        response_data[["BRU_E"]] * ips[["weight"]]
      ),
      BRU_response = c(
        point_weights,
        rep(0, NROW(ips))
      ),
      BRU_block = c(
        data[[".block"]],
        ips[[".block"]]
      ),
      BRU_weights = response_data[["BRU_weights"]][1],
      BRU_Ntrials = response_data[["BRU_Ntrials"]][1],
      BRU_scale = response_data[["BRU_scale"]][1]
    )
    data <- extended_bind_rows(data, ips)
  }
  if (ips_is_Spatial) {
    non_coordnames <- setdiff(names(data), data_coordnames)
    data <- sp::SpatialPointsDataFrame(
      coords = as.matrix(data[data_coordnames]),
      data = data[non_coordnames],
      proj4string = fm_CRS(data_crs),
      match.ID = FALSE
    )
  }

  lh$data <- data
  lh$response_data <- new_response_data
  lh$response <- "BRU_response"
  lh$inla.family <- "poisson"
  lh$integration_info$ips <- NULL

  lh <- bru_compat_pre_2_14_bru_obs(lh)

  lh
}

#
# > x<-sample(1:1000,size=1000000,replace=TRUE)
# > bench::mark(R = convert_group_cv_blocks_to_friends_list(x, method = "R"),
#               C = convert_group_cv_blocks_to_friends_list(x, method = "C"))
# # A tibble: 2 × 13
# expression      min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc
# total_time result memory time           gc
# <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>
# <bch:tm> <list> <list> <list>         <list>
# 1 R             4.23s    4.23s     0.236        NA    27.6      1   117
# 4.23s <list> <NULL> <bench_tm [1]> <tibble>
# 2 C            97.1ms  97.97ms    10.1          NA     1.68     6     1
# 595.38ms <list> <NULL> <bench_tm [6]> <tibble>
#
convert_group_cv_blocks_to_friends_list <- function(group_cv_block,
                                                    method = "C") {
  # Would like:
  # group_cv_friends <- lapply(seq_len(max(ips$.block)), function(i) {
  #   which(group_cv_block == i)
  # })
  # Current inla.group.cv interface (2025-08-28)
  # Inefficient implementation:
  # group_cv_friends <- lapply(seq_along(group_cv_block), function(i) {
  #   which(group_cv_block == group_cv_block[i])
  # })
  # Faster version:
  if (method == "R") {
    grp_index <- seq_len(max(group_cv_block))
    group_cv_friends_groups <- lapply(grp_index, function(i) {
      which(group_cv_block == i)
    })
    group_cv_friends <- lapply(group_cv_block, function(i) {
      group_cv_friends_groups[[i]]
    })
  } else if (method == "C") {
    group_cv_friends <-
      inlabru_group_cv_block_conversion(as.integer(group_cv_block),
        as.integer(max(group_cv_block)),
        per_node = as.logical(TRUE)
      )
  }
  group_cv_friends
}


bru_obs_family_cp <- function(lh, options, .envir) {
  if (inherits(lh[["data"]], "Spatial") ||
    inherits(lh[["BRU_original_response_data"]], "Spatial") ||
    inherits(lh[["integration_info"]][["samplers"]], "Spatial") ||
    inherits(lh[["integration_info"]][["ips"]], "Spatial")) {
    return(bru_obs_family_cp_sp(lh, options, .envir))
  }

  if (!is.null(lh[["aggregate"]])) {
    stop("The 'aggregate' feature cannot be used with family='cp'.")
  }

  response <- lh[["response_data"]][[lh[["response"]]]]
  orig_response_data <- lh[["BRU_original_response_data"]]
  response_data <- lh[["response_data"]]
  data <- lh[["data"]]
  domain <- lh[["integration_info"]][["domain"]]
  samplers <- lh[["integration_info"]][["samplers"]]
  ips <- lh[["integration_info"]][["ips"]]

  response_expr_text <- bru_pred_expr(lh, format = "resp_text")

  # Catch and handle special cases:
  if (is.null(response) || !inherits(response, "list")) {
    domain_names <- trimws(strsplit(response_expr_text, split = "\\+")[[1]])
    if (!is.null(domain_names)) {
      # "a + b" conversion to list(a = a, b = b)
      domain_expr <- paste0(
        "list(",
        paste0(
          vapply(
            domain_names, function(x) {
              if (identical(x, "coordinates")) {
                glue("{x} = sp::{x}(.data.)")
              } else {
                glue("{x} = {x}")
              }
            },
            ""
          ),
          collapse = ", "
        ),
        ")"
      )
      response_expr <- rlang::parse_expr(domain_expr)
    }
    response <- tryCatch(
      expr = bru_eval_in_data_context(
        !!response_expr,
        data = list(response_data = orig_response_data, data = data),
        default = NULL,
        .envir = .envir
      ),
      error = function(e) {
        NULL
      }
    )
  }

  if (is.null(response)) {
    stop(paste0(
      "You called bru_obs() with family='cp' but the evaluated ",
      "response information is NULL"
    ))
  }

  if (is.null(ips)) {
    if (is.null(domain)) {
      stop(paste0(
        "The family='cp' model requires a 'domain' specification compatible ",
        "with 'fmesher::fm_int()'"
      ))
    }
    if (!setequal(names(response), names(domain))) {
      stop(glue("
          Mismatch between response and domain names:
            names(response) = ({glue_collapse(names(response), sep = ', ')})
            names(domain)   = ({glue_collapse(names(domain), sep = ', ')})"))
    }

    ips <- fm_int(
      domain = domain,
      samplers = samplers,
      int.args = options[["bru_int_args"]]
    )
  }

  if (length(unique(response_data[["BRU_E"]])) > 1) {
    bru_log_warn(
      "Exposure/effort parameter E should be a scalar for likelihood 'cp'."
    )
  }

  # TODO: check that the crs info is the same

  # For non-Spatial models:
  # Use the response data list as the actual data object, since that's now the
  # canonical place where the point information is given.  This also allows
  # response_data to be used when constructing the response list. This makes
  # it a strict requirement that the predictor can be evaluated as a pure
  # function of the domain data.  When implementing sf support, might be able
  # to give explicit access to spatial coordinates, but otherwise the user can
  # extract it from the geometry with st_coordinates(geometry)[,1] or similar.
  # For Spatial models, keep the old behaviour for backwards compatibility for
  # now, but can likely realign that in the future after more testing.
  # Save general response data to add to response (precomputed covariates etc)
  if (!is.null(orig_response_data)) {
    data_ <- orig_response_data
  } else {
    data_ <- data
  }

  data <- tibble::as_tibble(response)
  if (("geometry" %in% names(data)) &&
    inherits(data$geometry, "sfc")) {
    sf::st_geometry(data) <- "geometry"
  }

  orig_response_data <- NULL
  n_block <- max(ips$.block)
  if (".block" %in% names(data_)) {
    n_block <- max(max(data_$.block), n_block)
    N_data <- base::tabulate(data_$.block, n_block)
  } else {
    N_data <- NROW(data)
  }

  # Add back additional data
  additional_data_names <- setdiff(names(data_), names(data))
  if ((length(additional_data_names) > 0) &&
    (NROW(data_) == sum(N_data))) {
    data <- cbind(data, tibble::as_tibble(data_)[additional_data_names])
  }
  if (is.null(data[[".block"]])) {
    data$.block <- 1L
  }

  if ("geometry" %in% names(ips)) {
    sf::st_geometry(ips) <- "geometry"
  }

  # Use 'weights' for per-point weighting of eta
  if (length(response_data[["BRU_weights"]]) == 1L) {
    point_weights <- rep(response_data[["BRU_weights"]], sum(N_data))
  } else {
    stopifnot(length(response_data[["BRU_weights"]]) == sum(N_data))
    point_weights <- response_data[["BRU_weights"]]
  }
  response_data[["BRU_weights"]] <- 1L

  if (identical(options[["bru_compress_cp"]], TRUE)) {
    new_response_data <- tibble::tibble(
      BRU_E = c(
        rep(0, length(N_data[N_data > 0])),
        response_data[["BRU_E"]] * ips[["weight"]]
      ),
      BRU_response = c(
        N_data[N_data > 0],
        rep(0, NROW(ips))
      ),
      BRU_block = c(
        which(N_data > 0),
        ips[[".block"]]
      ),
      BRU_weights = response_data[["BRU_weights"]][1],
      BRU_Ntrials = response_data[["BRU_Ntrials"]][1],
      BRU_scale = response_data[["BRU_scale"]][1]
    )

    lh$data_extra[["BRU_cp_block_subset"]] <- which(N_data > 0)
    lh$data_extra[["BRU_cp_n_block"]] <- n_block
    lh$data_extra[["BRU_cp_agg"]] <-
      bm_aggregate(
        type = "average",
        n_block = n_block
      )

    old_pred_expr <- bru_pred_expr(lh)
    pred_text <- bru_pred_expr(old_pred_expr, format = "text_raw")
    resp_text <- bru_pred_expr(old_pred_expr, format = "resp_text")
    lh$data_extra[["BRU_cp_predictor"]] <-
      function(eta, data, data_extra) {
        eta <- rep_len(eta, length(data$BRU_aggregate))
        c(
          ibm_eval(
            data_extra$BRU_cp_agg,
            list(block = data$.block[data$BRU_aggregate]),
            state = eta[data$BRU_aggregate] *
              data$BRU_point_weights[data$BRU_aggregate]
          )[data_extra$BRU_cp_block_subset],
          eta[!data$BRU_aggregate]
        )
      }
    pred_text <-
      glue("BRU_cp_predictor({{ {pred_text} }}, .data., .data_extra.)")
    lh$pred_expr <- new_bru_pred_expr(
      pred_text,
      used = lh$pred_expr$used,
      is_rowwise = FALSE,
      .envir = lh$pred_expr$.envir
    )
    lh$pred_expr$is_additive <- FALSE
    lh$pred_expr$is_linear <- old_pred_expr$is_linear
    lh$pred_expr$resp_text <- resp_text

    data <- extended_bind_rows(
      dplyr::bind_cols(data,
        BRU_aggregate = TRUE,
        BRU_point_weights = point_weights
      ),
      dplyr::bind_cols(ips,
        BRU_aggregate = FALSE,
        BRU_point_weights = 0.0
      )
    )
  } else {
    new_response_data <- data.frame(
      BRU_E = c(
        rep(0, length(point_weights)),
        response_data[["BRU_E"]] * ips[["weight"]]
      ),
      BRU_response = c(
        point_weights,
        rep(0, NROW(ips))
      ),
      BRU_block = c(
        data[[".block"]],
        ips[[".block"]]
      ),
      BRU_weights = response_data[["BRU_weights"]][1],
      BRU_Ntrials = response_data[["BRU_Ntrials"]][1],
      BRU_scale = response_data[["BRU_scale"]][1]
    )
    data <- extended_bind_rows(data, ips)
  }
  if (!is.null(lh$control.gcpo)) {
    group_cv_block <- new_response_data$BRU_block
    group_cv_friends <- convert_group_cv_blocks_to_friends_list(group_cv_block)
    lh$control.gcpo <- modifyList(
      lh$control.gcpo,
      list(friends = group_cv_friends)
    )
  }

  if (inherits(data, "sf")) {
    data <- fmesher::fm_zm(data)
  }
  lh$data <- data
  lh$response_data <- new_response_data
  lh$response <- "BRU_response"
  lh$inla.family <- "poisson"
  lh$integration_info$ips <- NULL

  lh <- bru_compat_pre_2_14_bru_obs(lh)

  lh
}


bru_obs_check_used_deprecation <- function(
  used,
  include,
  exclude,
  include_latent,
  expr_text,
  env = rlang::caller_env(),
  user_env = rlang::caller_env(2)
) {
  if (lifecycle::is_present(include) ||
    lifecycle::is_present(exclude) ||
    lifecycle::is_present(include_latent)) {
    if (lifecycle::is_present(include)) {
      bru_log_message(
        paste0(
          "The `include` argument of `bru_obs()` is deprecated ",
          "since inlabru 2.11.0 and will be ignored.\n\t",
          "If auto-detection doesn't work, use ",
          "`used = bru_used(effect = include)` instead."
        ),
        verbosity = 1L
      )
      lifecycle::deprecate_warn(
        "2.11.0",
        "bru_obs(include)",
        "bru_obs(used)",
        c(
          "The provided `include` value will be ignored.",
          "If auto-detection doesn't work, use `bru_used(effect = include)`"
        ),
        env = env,
        user_env = user_env
      )
    }
    if (lifecycle::is_present(exclude)) {
      bru_log_message(
        paste0(
          "The `exclude` argument of `bru_obs()` is deprecated ",
          "since inlabru 2.11.0 and may be ignored.\n\t",
          "If auto-detection doesn't work, use ",
          "`used = bru_used(effect_exclude = exclude)` instead."
        ),
        verbosity = 1L
      )
      lifecycle::deprecate_warn(
        "2.11.0",
        "bru_obs(exclude)",
        "bru_obs(used)",
        c(
          "The provided `exclude` value will be ignored.",
          paste0(
            "If auto-detection doesn't work, ",
            "use `bru_used(effect_exclude = exclude)`"
          )
        ),
        env = env,
        user_env = user_env
      )
    }
    if (lifecycle::is_present(include_latent)) {
      bru_log_message(
        paste0(
          "The `include_latent` argument of `bru_obs()` is deprecated ",
          "since inlabru 2.11.0 and will be ignored.\n\t",
          "If auto-detection doesn't work, use ",
          "`used = bru_used(latent = include_latent)` instead."
        ),
        verbosity = 1L
      )
      lifecycle::deprecate_warn(
        "2.11.0",
        "bru_obs(include_latent)",
        "bru_obs(used)",
        c(
          "The provided `include_latent` value will be ignored.",
          paste0(
            "If auto-detection doesn't work, ",
            "use `bru_used(latent = include_latent)`"
          )
        ),
        env = env,
        user_env = user_env
      )
    }
    if (is.null(used)) {
      used <- bru_used(expr_text)
    }
  }
  used
}

bru_obs_handle_is_rowwise <- function(pred_expr,
                                      is_rowwise = NULL,
                                      allow_combine = NULL,
                                      data,
                                      response_data,
                                      aggregate) {
  if (is.null(is_rowwise) && is.logical(allow_combine)) {
    is_rowwise <- !allow_combine
  }
  if (!is.logical(is_rowwise)) {
    if (!is.null(aggregate)) {
      is_rowwise <- FALSE
    } else if (!is.null(response_data)) {
      bru_log_warn(
        paste0(
          "Non-null response data supplied; ",
          "guessing is_rowwise=FALSE.",
          "\n  Specify is_rowwise explicitly to avoid this warning."
        )
      )
      is_rowwise <- FALSE
    } else if (is.list(data) && !is.data.frame(data)) {
      bru_log_warn(
        paste0(
          "Non data-frame list-like data supplied; ",
          "guessing is_rowwise=FALSE.\n",
          "  Specify is_rowwise explicitly to avoid this warning."
        )
      )
      is_rowwise <- FALSE
    } else {
      is_rowwise <- TRUE
    }
  }
  if (is.null(pred_expr$is_rowwise)) {
    pred_expr$is_rowwise <- is_rowwise
  } else {
    pred_expr$is_rowwise <- pred_expr$is_rowwise && is_rowwise
  }
  pred_expr
}


#' @title Observation model construction for usage with [bru()]
#'
#' @description Observation model construction for usage with [bru()].
#'
#' Note: Prior to version `2.12.0`, this function was called `like()`, and that
#' alias will remain for a while until examples etc have been updated and users
#' made aware of the change. The name change is to avoid issues with namespace
#' clashes, e.g. with `data.table::like()`, and also to signal that the function
#' defines observation models, not just likelihood functions.
#'
#' @rdname bru_obs
#' @export
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
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
#' @param data Predictor expression-specific data, as a `data.frame`, `tibble`,
#'  or `sf`.  Since `2.12.0.9023`, deprecated support for
#' `SpatialPoints[DataFrame]` objects.
#' @param response_data Observation/response-specific data for models that need
#'   different size/format for inputs and response variables, as a `data.frame`,
#'   `tibble`, or `sf`. Since `2.12.0.9023`, deprecated support for
#'   `SpatialPoints[DataFrame]` objects.
#' @param data_extra object convertible with `as.list()` with additional
#'   variables to be made available in predictor evaluations. Variables with
#'   the same names as the data object will be ignored, unless accessed via
#'   `.data_extra.[["name"]]` or `.data_extra.$name` in the formula.
#' @param E Exposure/effort parameter for family = 'poisson' passed on to
#'   `INLA::inla`. Special case if family is 'cp': rescale all integration
#'   weights by a scalar `E`. For sampler specific reweighting/effort, use a
#'   `weight` column in the `samplers` object instead, see [fmesher::fm_int()].
#'   Default taken from `options$E`, normally `1`.
#' @param Ntrials A vector containing the number of trials for the 'binomial'
#'   likelihood. Default taken from `options$Ntrials`, normally `1`.
#' @param weights Fixed (optional) weights parameters of the likelihood, so the
#'   log-likelihood`[i]` is changed into `weights[i] * log_likelihood[i]`.
#'   Default value is `1`. WARNING: The normalizing constant for the likelihood
#'   is NOT recomputed, so ALL marginals (and the marginal likelihood) must be
#'   interpreted with great care.
#'
#'   For `family = "cp"`, the weights are applied as `sum(weights * eta)` in
#'   the point location contribution part of the log-likelihood, where `eta` is
#'   the linear predictor, and do not affect the integration part of the
#'   likelihood. This can be used to implement approximative methods for point
#'   location uncertainty.
#' @param scale Fixed (optional) scale parameters of the precision for several
#'   models, such as Gaussian and student-t response models.
#' @param domain,samplers,ips Arguments used for `family="cp"` and `aggregate=`.
#' \describe{
#' \item{`domain`}{Named list of domain definitions, see [fmesher::fm_int()].}
#' \item{`samplers`}{Integration domain for `family="cp"` or subdomains for
#'   `aggregate=`, see [fmesher::fm_int()].}
#' \item{`ips`}{Integration points. Defaults
#'   to [fmesher::fm_int]`(domain, samplers)`. If explicitly given,
#'   overrides `domain` and `samplers`. `fmesher::new_fm_int()` (from fmesher
#'   `0.5.0.9013`) can be used for manually constructed integration schemes.}
#' }
#' @param used Either `NULL` (default) or a [bru_used()] object. When,
#'   `NULL`, the information about what effects and
#' latent vectors are made available to the predictor evaluation is defined by
#' `bru_used(formula)`, which will include all effects and latent vectors
#' used by the predictor expression.
#' @param is_rowwise,allow_combine logical; If `is_rowwise` is `FALSE`, the
#'   predictor expression may involve several rows of the input data to
#'   influence the same row. When `NULL`, it defaults to `TRUE`, unless
#'   `response_data` is non-`NULL`, or `data` is a `list`, or the likelihood
#'   construction requires it. The `allow_combine` argument is a deprecated and
#'   inverted version, such that `is_rowwise = !allow_combine`.
#' @param aggregate character ("none" or a valid name for the `type` argument of
#'   [bm_aggregate()]) or an aggregation `bru_mapper` object ([bm_aggregate()],
#'   [bm_logsumexp()], or [bm_logitaverage()]). Default `NULL`, interpreted as
#'   "none". `r lifecycle::badge("experimental")`, available from version
#'   `2.12.0.9013`.
#' @param aggregate_input `NULL` or an optional input list to the mapper
#'   defined by non-NULL `aggregate`, overriding the default,
#'   ```
#'   list(block = .data.[[".block"]],
#'        weights = .data.[["weight"]],
#'        n_block = bru_response_size(.response_data.),
#'        block_response = ".block")
#'   ```
#'   `r lifecycle::badge("experimental")`, available from version `2.12.0.9013`.
#'
#'   From `2.13.0.9016`, it will look for a `block_response` element in the
#'   list, which should be the name of a response variable in `response_data` to
#'   match the `block` information against, allowing character or factor
#'   aggregation block information to be used by replacing `block` with
#'   `match(block, block_response)`. If not supplied, the name
#'   `".block"` is tried. If that isn't available, the `block` information must
#'   be supplied directly as a numeric or integer vector, indexing into the rows
#'   of the response variable. Having no `block_response` variable is equivalent
#'   to having
#'   `response_data$.block = seq_len(bru_response_size(response_data))`.
#' @param control.family A optional `list` of `INLA::control.family` options
#' @param control.gcpo A optional `list` of `INLA::control.gcpo` options
#' @param tag character; Name that can be used to identify the relevant parts
#' of INLA predictor vector output, via [bru_index()].
#' @param options A [bru_options] options object or a list of options passed
#' on to [bru_options()]
#' @param .envir The evaluation environment to use for special arguments (`E`,
#'   `Ntrials`, `weights`, and `scale`) if not found in `response_data` or
#'   `data`. Defaults to the calling environment.
#' @param include,exclude,include_latent `r lifecycle::badge("deprecated")`, use
#'   `used` instead.
#'
#' @details
#' The `E`, `Ntrials`, `weights`, and `scale` arguments are evaluated in the
#' data context, with values from `response_data` taking precedence over `data`.
#'
#' @return A likelihood configuration which can be used to parameterise [bru()].
#' @seealso [bru_response_size()], [bru_used()], [bru_comp()],
#' [bru_comp_eval()]
#'
#' @example inst/examples/bru_obs.R
bru_obs <- function(formula = . ~ .,
                    family = "gaussian",
                    data = NULL,
                    response_data = NULL,
                    data_extra = NULL,
                    E = NULL,
                    Ntrials = NULL,
                    weights = NULL,
                    scale = NULL,
                    domain = NULL,
                    samplers = NULL,
                    ips = NULL,
                    used = NULL,
                    is_rowwise = NULL,
                    aggregate = NULL,
                    aggregate_input = NULL,
                    control.family = NULL,
                    control.gcpo = NULL,
                    tag = NULL,
                    options = list(),
                    .envir = parent.frame(),
                    include = deprecated(),
                    exclude = deprecated(),
                    include_latent = deprecated(),
                    allow_combine = deprecated()) {
  options <- bru_call_options(options)
  bru_options_set_local(options, .reset = TRUE)

  check_sp_data_deprecation(
    data = data,
    response_data = response_data,
    samplers = samplers,
    ips = ips
  )

  # Some defaults
  inla.family <- family

  # Only for deprecated include/exclude/include_latent argument handling:
  formula_char <- as.character(formula)
  used <- bru_obs_check_used_deprecation(
    used = used,
    include = rlang::maybe_missing(include),
    exclude = rlang::maybe_missing(exclude),
    include_latent = rlang::maybe_missing(include_latent),
    expr_text = formula_char[length(formula_char)]
  )

  # Build expressions from formula:
  pred_expr <- new_bru_pred_expr(formula,
    used = used,
    .envir = .envir
  )

  # Set the response name
  if (is.null(pred_expr$resp_text)) {
    stop("Missing response variable name or expression")
  }
  response_expr <- rlang::parse_expr(pred_expr$resp_text)
  data_list <- list(response_data = response_data, data = data)
  response <- tryCatch(
    expr = bru_eval_in_data_context(
      !!response_expr,
      data = data_list,
      default = NULL,
      .envir = .envir
    ),
    error = function(e) {
      NULL
    }
  )

  E <- bru_eval_in_data_context(
    {{ E }},
    data = data_list,
    default = options[["E"]],
    .envir = .envir
  )
  Ntrials <- bru_eval_in_data_context(
    {{ Ntrials }},
    data = data_list,
    default = options[["Ntrials"]],
    .envir = .envir
  )
  weights <- bru_eval_in_data_context(
    {{ weights }},
    data = data_list,
    default = 1,
    .envir = .envir
  )
  scale <- bru_eval_in_data_context(
    {{ scale }},
    data = data_list,
    default = 1,
    .envir = .envir
  )

  data_extra <- as.list(data_extra)

  # Decide on is_rowwise default
  pred_expr <- bru_obs_handle_is_rowwise(
    pred_expr = pred_expr,
    is_rowwise = is_rowwise,
    allow_combine = allow_combine,
    data = data,
    response_data = response_data,
    aggregate = aggregate
  )

  original_response_data <- response_data
  # Must be a list or tibble to allow inla.mdata responses
  # Until inla.surv is converted to tibble (or data.frame), must be a list
  response_data <- list(
    BRU_response = response,
    BRU_E = E,
    BRU_Ntrials = Ntrials,
    BRU_scale = scale,
    BRU_weights = weights
  )

  # Prototype object
  lh <- structure(
    list(
      family = family,
      pred_expr = pred_expr,
      response_data = response_data,
      data = data,
      data_extra = data_extra,
      integration_info = list(
        ips = ips,
        domain = domain,
        samplers = samplers
      ),
      response = "BRU_response",
      inla.family = inla.family,
      control.family = control.family,
      control.gcpo = control.gcpo,
      tag = tag,
      BRU_original_response_data = original_response_data
    ),
    class = "bru_obs"
  )

  lh$formula <- formula
  lh <- bru_compat_pre_2_14_bru_obs(lh)

  if (!is.null(aggregate)) {
    lh <- bru_obs_agg(lh,
      aggregate = aggregate,
      aggregate_input = {{ aggregate_input }},
      options = options,
      .envir = .envir
    )
  }

  # More on special bru likelihoods
  if (family == "cp") {
    lh <- bru_obs_family_cp(lh, options = options, .envir = .envir)
  }

  if (is.null(lh[["response_data"]][[lh[["response"]]]])) {
    stop("Response variable missing or could not be evaluated")
  }

  if (is.null(lh[["used"]])) {
    lh[["used"]] <- bru_used(formula)
  }

  # Remove temporary copy of input data
  lh$BRU_original_response_data <- NULL

  # Return likelihood
  lh
}

#' @describeIn bru_obs `r lifecycle::badge("deprecated")` Legacy `like()`
#' method for `inlabru` prior to version `2.12.0`. Use [bru_obs()] instead.
#' @param \dots Arguments passed from `like()` to `bru_obs()`
#' @export
like <- function(...,
                 E = NULL,
                 Ntrials = NULL,
                 weights = NULL,
                 scale = NULL,
                 options = list(),
                 .envir = parent.frame()) {
  options <- bru_call_options(options)
  bru_options_set_local(options, .reset = TRUE)
  bru_log_message(
    paste0(
      "The `like()` function has been deprecated in favour of `bru_obs()`, ",
      "since inlabru 2.12.0."
    ),
    verbosity = 1L
  )
  lifecycle::deprecate_warn(
    "2.12.0",
    "like()",
    "bru_obs()"
  )

  bru_obs(
    ...,
    E = {{ E }},
    Ntrials = {{ Ntrials }},
    weights = {{ weights }},
    scale = {{ scale }},
    options = options,
    .envir = .envir
  )
}


#' @title Response size queries
#'
#' @description
#' Extract the number of response values from `bru` and related objects.
#' @return An `integer` vector.
#' @param object An object from which to extract response size(s).
#' @seealso [like()]
#' @export
#' @examples
#' bru_response_size(
#'   bru_obs(y ~ 1, data = data.frame(y = rnorm(10)), family = "gaussian")
#' )
bru_response_size <- function(object) {
  UseMethod("bru_response_size")
}

#' @describeIn bru_response_size Extract the number of observations from an
#'   object supporting `NROW()`.
#' @export
bru_response_size.default <- function(object) {
  NROW(object)
}

#' @describeIn bru_response_size Extract the number of observations from an
#'   object supporting `NROW(object[[1]])`.
#' @export
bru_response_size.list <- function(object) {
  NROW(object[[1]])
}

#' @describeIn bru_response_size Extract the number of observations from an
#' `inla.surv` object.
#' @export
bru_response_size.inla.surv <- function(object) {
  # This special case handling is needed because inla.surv objects are
  # lists. When the INLA package changes representation of `inla.surv`
  # to a tibble or data.frame, this special case won't be needed anymore.
  NROW(object[[1]])
}

#' @describeIn bru_response_size Extract the number of observations from a
#' `bru_obs` object.
#' @export
bru_response_size.bru_obs <- function(object) {
  bru_response_size(object[["response_data"]][[object[["response"]]]])
}

#' @describeIn bru_response_size Extract the number of observations from a
#' `bru_obs_list` object, as a vector with one value per observation model.
#' @export
bru_response_size.bru_obs_list <- function(object) {
  vapply(object, bru_response_size, 1L)
}

#' @describeIn bru_response_size Extract the number of observations from a
#' `bru_info` object, as a vector with one value per observation model.
#' @export
bru_response_size.bru_info <- function(object) {
  bru_response_size(object[["lhoods"]])
}

#' @describeIn bru_response_size Extract the number of observations from a
#' `bru` object, as a vector with one value per observation model.
#' @export
bru_response_size.bru <- function(object) {
  bru_response_size(object[["bru_info"]])
}


set_list_names <- function(x, tag, priority = "immutable") {
  priority <- match.arg(priority, c("immutable", "tag", "name"))
  list_names <- names(x)
  if (is.null(list_names)) {
    list_names <- rep(NA_character_, length(x))
  } else {
    list_names <- vapply(
      list_names,
      function(xx) {
        if (identical(xx, "")) {
          NA_character_
        } else {
          xx
        }
      },
      ""
    )
  }
  tag_names <- vapply(
    x,
    function(xx) {
      if (is.null(xx[[tag]]) ||
        identical(xx[[tag]], "")) {
        NA_character_
      } else {
        xx[[tag]]
      }
    },
    ""
  )
  new_names <- vapply(
    seq_along(x),
    function(k) {
      if (priority == "immutable") {
        if (is.na(list_names[k])) {
          return(tag_names[k])
        }
        if (is.na(tag_names[k])) {
          return(list_names[k])
        }
        if (identical(list_names[k], tag_names[k])) {
          return(list_names[k])
        }
        stop(glue("
          When `priority = 'immutable'`, the tag and name must match.
          Cannot combine objects with mismatching tags and names:
            Tag: {tag_names[k]}
            Name: {list_names[k]}
          Use `NA` for either the tag or name, or use matching tags/names."))
      } else if (priority == "tag") {
        if (!is.na(tag_names[k])) {
          return(tag_names[k])
        }
        if (!is.na(list_names[k])) {
          return(list_names[k])
        }
        NA_character_
      } else if (priority == "name") {
        if (!is.na(list_names[k])) {
          return(list_names[k])
        }
        if (!is.na(tag_names[k])) {
          return(tag_names[k])
        }
        NA_character_
      }
    },
    ""
  )
  names(x) <- new_names
  # Set tags on missing tags
  for (k in which(is.na(tag_names))) {
    x[[k]][[tag]] <- new_names[k]
  }

  x
}

#' @describeIn bru_obs
#' Combine `bru_obs` observation model object into a `bru_obs_list` object
#' @param \dots For `bru_obs_list.bru_obs`, one or more `bru_obs` objects
#' @export
bru_obs_list <- function(..., .tag = NULL) {
  UseMethod("bru_obs_list")
}

#' @describeIn bru_obs
#' Combine one or more lists of `bru_obs` observation model objects
#' into a `bru_obs_list` object
#' @param object A list of `bru_obs` and/or `bru_obs_list` objects
#' @export
bru_obs_list.list <- function(object, ..., .tag = NULL) {
  if (length(list(...)) > 0) {
    return(bru_obs_list(list(object, ...)))
  }
  if (length(object) == 0) {
    return(
      structure(
        list(),
        class = c("bru_obs_list", "list")
      )
    )
  }
  is_obs <- vapply(object, function(x) inherits(x, "bru_obs"), TRUE)
  if (!all(is_obs)) {
    is_obs_list <- vapply(object, function(x) inherits(x, "bru_obs_list"), TRUE)
    if (!all(is_obs_list)) {
      object <- lapply(seq_along(object), function(k) {
        bru_obs_list(object[[k]], .tag = names(object)[k])
      })
    }
    object <- unlist(object, recursive = FALSE)
  }
  object <- structure(
    object,
    class = c("bru_obs_list", "list")
  )

  object <- set_list_names(object, tag = "tag", priority = "immutable")

  object
}

#' @describeIn bru_obs
#' Combine one or more lists of `bru_obs` observation model objects
#' into a `bru_obs_list` object
#' @param .tag Optional name to assign to a single `bru_obs` object. Reserved
#'   for internal use.
#' @export
bru_obs_list.bru_obs <- function(..., .tag = NULL) {
  if (length(list(...)) != 1L) {
    return(bru_obs_list(list(...)))
  }

  object <- structure(
    list(...),
    class = c("bru_obs_list", "list")
  )
  if (!is.null(.tag)) {
    names(object) <- .tag
  }
  object <- set_list_names(object, tag = "tag", priority = "immutable")
  object
}

#' @describeIn bru_obs
#' Combine one or more `bru_obs_list` objects into a `bru_obs_list` object
#' @export
bru_obs_list.bru_obs_list <- function(..., .tag = NULL) {
  if (length(list(...)) != 1L) {
    return(bru_obs_list(list(...)))
  }
  object <- list(...)[[1]]
  object <- set_list_names(object, tag = "tag", priority = "immutable")
  object
}


#' @describeIn bru_obs
#' Combine several `bru_obs` objects into a `bru_obs_list` object
#' @export
c.bru_obs <- function(...) {
  bru_obs_list(list(...))
}


#' @describeIn bru_obs
#' Combine several `bru_obs_list` objects into a `bru_obs_list` object
#' @export
c.bru_obs_list <- function(...) {
  bru_obs_list(list(...))
}


#' @export
#' @param x `bru_obs_list` object from which to extract element(s)
#' @param i indices specifying elements to extract
#' @rdname bru_obs
#' @seealso [summary.bru_obs()]
`[.bru_obs_list` <- function(x, i) {
  env <- environment(x)
  object <- NextMethod()
  class(object) <- c("bru_obs_list", "list")
  environment(object) <- env
  object
}


#' Summary and print methods for observation models
#'
#' @rdname bru_obs_print
#' @seealso [bru_obs()]
#' @param object Object to operate on
#' @param verbose logical; If `TRUE`, include more details of the
#' component definitions. If `FALSE`, only show basic component
#' definition information. Default: `TRUE`
#' @param \dots Arguments passed on to other `summary` methods
#' @param x Object to be printed
#' @method summary bru_obs
#' @export
#' @examples
#' obs <- bru_obs(y ~ ., data = data.frame(y = rnorm(10)))
#' summary(obs)
#' print(obs)
#'
summary.bru_obs <- function(object, verbose = TRUE, ...) {
  structure(
    list(
      family = object[["family"]],
      data_class = class(object[["data"]]),
      response_class = class(object[["response_data"]][[object[["response"]]]]),
      predictor =
        glue::glue("{bru_pred_expr(object, format = 'formula_text')}"),
      is_additive = bru_is_additive(object),
      is_linear = bru_is_linear(object),
      is_rowwise = bru_is_rowwise(object),
      used = bru_used(object),
      tag = object[["tag"]]
    ),
    class = "summary_bru_obs"
  )
}

#' @describeIn bru_obs `r lifecycle::badge("deprecated")`
#' Backwards compatibility for versions `<= 2.12.0`. For later versions, use
#' `as_bru_obs_list()`, `bru_obs_list()`, or `c()`.
#' @export
like_list <- function(...) {
  lifecycle::deprecate_warn(
    "2.12.0",
    "like_list()",
    "bru_obs_list()",
    details = paste0(
      "Use `as_bru_obs_list()`, `bru_obs_list(...)` or ",
      "`c(...)` to construct observation model lists."
    )
  )
  bru_obs_list(list(...))
}

#' @describeIn bru_obs `r lifecycle::badge("deprecated")`
#' Backwards compatibility for versions `<= 2.12.0.9017`. For later versions,
#' use `as_bru_obs_list()`, `bru_obs_list()` or `c()`.
#' @export
bru_like_list <- function(...) {
  lifecycle::deprecate_warn(
    "2.12.0.9017",
    "bru_like_list()",
    "bru_obs_list()",
    details = paste0(
      "Use `as_bru_obs_list()`, `bru_obs_list(...)` or ",
      "`c(...)` to construct observation model lists."
    )
  )
  bru_obs_list(list(...))
}

#' @rdname bru_obs_print
#' @method summary bru_obs_list
#' @export
summary.bru_obs_list <- function(object, verbose = TRUE, ...) {
  structure(
    lapply(
      object,
      function(x) {
        summary(x, verbose = verbose, ...)
      }
    ),
    class = "summary_bru_obs_list"
  )
}


#' @rdname bru_obs_print
#' @export
#' @importFrom glue glue
print.summary_bru_obs <- function(x, ...) {
  lh <- x
  cat(
    glue_data(
      lh,
      "  Model tag: {tag}\n",
      "    Family: '{family}'\n",
      "    Data class: {data_class}\n",
      "    Response class: {response_class}\n",
      "    Predictor: {predictor}\n",
      "    Additive/Linear/Rowwise: {is_additive}/{is_linear}/{is_rowwise}\n",
      "    Used components: {format(lh$used)}",
      .trim = FALSE,
      tag = if (is.null(lh$tag) || is.na(lh$tag)) {
        "<No tag>"
      } else {
        glue::glue_collapse(glue::glue_data(lh, "'{tag}'"), sep = ", ")
      },
      data_class = glue::glue_collapse(
        glue::glue_data(lh, "'{data_class}'"),
        sep = ", "
      ),
      response_class = glue::glue_collapse(
        glue::glue_data(lh, "'{response_class}'"),
        sep = ", "
      ),
      predictor =
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
    ),
    "\n"
  )
  invisible(x)
}

#' @rdname bru_obs_print
#' @export
print.summary_bru_obs_list <- function(x, ...) {
  cat("Observation models:\n")
  for (lh in x) {
    print(lh)
  }
  invisible(x)
}

#' @export
#' @rdname bru_obs_print
print.bru_obs <- function(x, ...) {
  print(summary(x))
  invisible(x)
}

#' @export
#' @rdname bru_obs_print
print.bru_obs_list <- function(x, ...) {
  print(summary(x))
  invisible(x)
}


#' Utility functions for bru observation model objects
#' @param x Object of `bru_obs` or `bru_obs_list` type
#' @param \dots Further arguments passed on to the submethods
#' @export
#' @keywords internal
#' @returns * `bru_obs_inla_family()` returns a string or vector of strings
#' @rdname bru_obs_methods
#' @name bru_obs_methods
#' @seealso [summary.bru_obs()]
bru_obs_inla_family <- function(x, ...) {
  UseMethod("bru_obs_inla_family")
}
#' @export
#' @rdname bru_obs_methods
bru_obs_inla_family.bru_obs <- function(x, ...) {
  x[["inla.family"]]
}
#' @export
#' @rdname bru_obs_methods
bru_obs_inla_family.bru_obs_list <- function(x, ...) {
  vapply(x, bru_obs_inla_family, "")
}

#' @param control.family list of INLA `control.family` options to override
#' @export
#' @keywords internal
#' @returns * `bru_obs_control_family()` returns a list with
#'   `INLA::control.family` options, or a list of such lists, with one element
#'   per observation model
#' @rdname bru_obs_methods
bru_obs_control_family <- function(x, control.family = NULL, ...) {
  UseMethod("bru_obs_control_family")
}

#' @export
#' @rdname bru_obs_methods
bru_obs_control_family.bru_obs <- function(x, control.family = NULL, ...) {
  if (!is.null(control.family)) {
    control.family
  } else if (is.null(x[["control.family"]])) {
    list()
  } else {
    x[["control.family"]]
  }
}

#' @export
#' @rdname bru_obs_methods
bru_obs_control_family.bru_obs_list <- function(x,
                                                control.family = NULL,
                                                ...) {
  # Extract the control.family information for each likelihood
  if (!is.null(control.family)) {
    if (length(control.family) != length(x)) {
      stop(paste0(
        "control.family supplied as option, but format doesn't match ",
        "the number of likelihoods"
      ))
    }
    like_has_cf <- vapply(
      seq_along(x),
      function(k) !is.null(x[[k]][["control.family"]]),
      TRUE
    )
    if (any(like_has_cf)) {
      bru_log_warn(
        glue(
          "Global control.family option overrides settings in likelihood(s) ",
          "{glue_collapse(which(like_has_cf), sep = ', ', last = ' and ')}",
        )
      )
    }
  } else {
    control.family <- lapply(x, bru_obs_control_family)
    # inla() requires a unnamed list of lists
    names(control.family) <- NULL
  }
  control.family
}

#' @export
#' @keywords internal
#' @returns * `bru_obs_control_gcpo()` returns a list with
#'   `INLA::control.gcpo` options, with predictor/response variable indices
#'   unified for multi-observation models. If a `bru_obs` model has a NULL
#'   `control.gcpo` argument, an empty list is returned for that model.
#' @rdname bru_obs_methods
bru_obs_control_gcpo <- function(x, ...) {
  UseMethod("bru_obs_control_gcpo")
}

#' @param index_offset integer; offset to add to indices in `control.gcpo`
#' @param index_length integer; length of the response vector for the
#'   observation model. Defaults to `bru_response_size(x)`.
#' @param force_weights logical; if `TRUE`, ensure that `control.gcpo$weights`
#'   is populated. This is needed if any of the observation models in a
#'   `bru_obs_list` has non-null `control.gcpo$weights`.
#' @export
#' @rdname bru_obs_methods
bru_obs_control_gcpo.bru_obs <- function(x,
                                         index_offset,
                                         index_length = bru_response_size(x),
                                         force_weights,
                                         ...) {
  c.gcpo <- x[["control.gcpo"]]
  if (is.null(c.gcpo)) {
    return(list())
  }
  for (nm in intersect(
    c("groups", "selection", "group.selection", "friends"),
    names(c.gcpo)
  )) {
    if (nm == "groups") {
      c.gcpo[[nm]] <- lapply(c.gcpo[[nm]], function(v) {
        v$idx <- v$idx + index_offset
        v
      })
    } else {
      c.gcpo[[nm]] <- lapply(c.gcpo[[nm]], function(v) v + index_offset)
    }
  }
  if (!("friends" %in% names(c.gcpo))) {
    c.gcpo[["friends"]] <- as.list(index_offset + seq_len(index_length))
  }
  if (force_weights && !("weights" %in% names(c.gcpo))) {
    c.gcpo[["weights"]] <- rep(1.0, index_length)
  } else if ("weights" %in% names(c.gcpo)) {
    if (length(c.gcpo[["weights"]]) != index_length) {
      stop(glue::glue(
        "Length of control.gcpo$weights ({length(c.gcpo[['weights']])}) ",
        "does not match response length ({index_length})"
      ))
    }
  }
  c.gcpo
}

#' @param control.gcpo list of INLA `control.gcpo` default options
#' @export
#' @rdname bru_obs_methods
bru_obs_control_gcpo.bru_obs_list <- function(x,
                                              control.gcpo = NULL,
                                              ...) {
  # Update the control.gcpo information for each likelihood
  response_sizes <- bru_response_size(x)
  any_element <- vapply(
    c("groups", "selection", "group.selection", "friends", "weights"),
    function(nm) {
      any(vapply(x, function(lh) {
        !is.null(lh[["control.gcpo"]][[nm]])
      }, TRUE))
    },
    TRUE
  )
  all_element <- vapply(
    c("groups", "selection", "group.selection", "friends", "weights"),
    function(nm) {
      all(vapply(x, function(lh) {
        !is.null(lh[["control.gcpo"]][[nm]])
      }, TRUE))
    },
    TRUE
  )
  c.gcpo <- lapply(
    seq_along(x),
    function(k) {
      bru_obs_control_gcpo(x[[k]],
        index_offset = sum(response_sizes[seq_len(k - 1)]),
        index_length = response_sizes[k],
        force_weights = any_element["weights"]
      )
    }
  )

  # If given in one model, must be given in all models
  c.gcpo.combined <- list()
  for (nm in c("groups", "selection", "group.selection")) {
    if (any_element[nm]) {
      if (!all_element[nm]) {
        stop(glue::glue(
          "control.gcpo${nm} given in some, but not all, observation models"
        ))
      }
      c.gcpo.combined[[nm]] <-
        do.call("c", lapply(c.gcpo, function(lh) lh[[nm]]))
    }
  }
  # Combine list/vector
  c.gcpo.combined[["friends"]] <-
    do.call("c", lapply(c.gcpo, function(lh) lh[["friends"]]))
  if (any_element[["weights"]]) {
    c.gcpo.combined[["weights"]] <-
      unlist(lapply(c.gcpo, function(lh) lh[["weights"]]))
  }

  merge_list_check <- function(a, b, exclude = character(0)) {
    for (nm in setdiff(names(b), exclude)) {
      if (nm %in% names(a)) {
        if (identical(a[[nm]], b[[nm]])) {
          next
        }
        stop(glue::glue(
          "Cannot merge control.gcpo lists: ",
          "both contain element '{nm}' with conflicting content."
        ))
      }
      a[[nm]] <- b[[nm]]
    }
    a
  }

  exc <- names(c.gcpo.combined)
  if (!is.null(control.gcpo)) {
    c.gcpo.combined <- merge_list_check(c.gcpo.combined,
      control.gcpo,
      exclude = exc
    )
  }
  for (idx in seq_along(c.gcpo)) {
    c.gcpo.combined <- merge_list_check(c.gcpo.combined,
      c.gcpo[[idx]],
      exclude = exc
    )
  }

  c.gcpo.combined
}


#' @title Log Gaussian Cox process (LGCP) inference using INLA
#'
#' @description
#' This function performs inference on a LGCP observed via points residing
#' possibly multiple dimensions. These dimensions are defined via the left hand
#' side of the formula provided via the model parameter. The left hand side
#' determines the intensity function that is assumed to drive the LGCP. This may
#' include effects that lead to a thinning (filtering) of the point process. By
#' default, the log intensity is assumed to be a linear combination of the
#' effects defined by the formula's RHS.
#'
#' More sophisticated models, e.g.
#' non-linear thinning, can be achieved by using the predictor argument. The
#' latter requires multiple runs of INLA for improving the required
#' approximation of the predictor. In many applications the LGCP is only
#' observed through subsets of the dimensions the process is living in. For
#' example, spatial point realizations may only be known in sub-areas of the
#' modelled space. These observed subsets of the LGCP domain are called samplers
#' and can be provided via the respective parameter. If samplers is NULL it is
#' assumed that all of the LGCP's dimensions have been observed completely.
#'
#' @export
#' @inheritParams bru_obs
#' @inheritParams bru
#' @param \dots Further arguments passed on to [bru_obs()].
#' @details
#' The `E` and `weights` arguments are evaluated in the data context, like for
#' [bru_obs()].
#' @return An [bru()] object
#' @examples
#' \donttest{
#' if (bru_safe_inla() &&
#'   require(ggplot2, quietly = TRUE) &&
#'   require(fmesher, quietly = TRUE) &&
#'   require(sn, quietly = TRUE)) {
#'   # Load the Gorilla data
#'   data <- gorillas_sf
#'
#'   # Plot the Gorilla nests, the mesh and the survey boundary
#'   ggplot() +
#'     geom_fm(data = data$mesh) +
#'     gg(data$boundary, fill = "blue", alpha = 0.2) +
#'     gg(data$nests, col = "red", alpha = 0.2)
#'
#'   # Define SPDE prior
#'   matern <- INLA::inla.spde2.pcmatern(
#'     data$mesh,
#'     prior.sigma = c(0.1, 0.01),
#'     prior.range = c(0.1, 0.01)
#'   )
#'
#'   # Define domain of the LGCP as well as the model components (spatial SPDE
#'   # effect and Intercept)
#'   cmp <- geometry ~ field(geometry, model = matern) + Intercept(1)
#'
#'   # Fit the model (with int.strategy="eb" to make the example take less time)
#'   fit <- lgcp(cmp, data$nests,
#'     samplers = data$boundary,
#'     domain = list(geometry = data$mesh),
#'     options = list(control.inla = list(int.strategy = "eb"))
#'   )
#'
#'   # Predict the spatial intensity surface
#'   lambda <- predict(
#'     fit,
#'     fm_pixels(data$mesh, mask = data$boundary),
#'     ~ exp(field + Intercept)
#'   )
#'
#'   # Plot the intensity
#'   ggplot() +
#'     gg(lambda, geom = "tile") +
#'     gg(data$nests, col = "red", alpha = 0.2)
#' }
#' }
#'
lgcp <- function(components,
                 data,
                 domain = NULL,
                 samplers = NULL,
                 ips = NULL,
                 formula = . ~ .,
                 E = NULL,
                 weights = NULL,
                 ...,
                 options = list(),
                 .envir = parent.frame()) {
  # If formula response missing, copy from components (for formula input)
  if (inherits(components, "formula")) {
    response <- extract_response(components)
    formula <- auto_response(formula, response)
  }
  lik <- bru_obs(
    family = "cp",
    formula = formula, data = data,
    domain = domain, samplers = samplers, ips = ips,
    ...,
    E = {{ E }},
    weights = {{ weights }},
    options = options,
    .envir = .envir
  )
  bru(components, lik, options = options, .envir = .envir)
}


expand_to_dataframe <- function(x, data = NULL) {
  if (is.null(data)) {
    data <- data.frame(matrix(nrow = NROW(x), ncol = 0))
  }
  only_x <- setdiff(names(x), names(data))
  if (length(only_x) < length(names(x))) {
    x <- x[!(names(x) %in% names(data))]
  }
  if (inherits(x, "Spatial")) {
    bru_safe_sp(force = TRUE)
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
    } else {
      result <- sp::cbind.Spatial(x, data)
    }
  } else if (is.data.frame(x)) {
    result <- cbind(x, data)
  } else {
    result <- c(x, as.list(data))
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
#' @export
#' @param object An object obtained by calling [bru()] or [lgcp()].
#' @param newdata A `data.frame` or `SpatialPointsDataFrame` of covariates
#'   needed for the prediction.
#' @param formula A formula where the right hand side defines an R expression
#' to evaluate for each generated sample. If `NULL`, the latent and
#' hyperparameter states are returned as named list elements.
#' See Details for more information.
#' @param n.samples Integer setting the number of samples to draw in order to
#' calculate the posterior statistics. The default is rather low but provides
#' a quick approximate result.
#' @param seed Random number generator seed passed on to `inla.posterior.sample`
#' @param probs A numeric vector of probabilities with values in `[0, 1]`,
#'   passed to `stats::quantile`
#' @param num.threads Specification of desired number of threads for parallel
#' computations. Default NULL, leaves it up to INLA.
#' When seed != 0, overridden to "1:1"
#' @param used Either `NULL` or a [bru_used()] object.
#'   Default, `NULL`, uses auto-detection of used variables in the formula.
#' @param drop logical; If `drop=FALSE`, and
#'   the prediction summary has the same number of rows as `newdata`, then the
#'   output is a joined object. Default `FALSE`.
#' @param \dots Additional arguments passed on to `inla.posterior.sample()`
#' @param include,exclude `r lifecycle::badge("deprecated")` If auto-detection
#' of used variables fails, use `used` instead.
#' @details
#' In addition to the component names (that give the effect of each component
#' evaluated for the input data), the suffix `_latent` variable name can be used
#' to directly access the latent state for a component, and the suffix function
#' `_eval` can be used to evaluate a component at other input values than the
#' expressions defined in the component definition itself, e.g.
#' `field_eval(cbind(x, y))` for a component that was defined with
#' `field(coordinates, ...)` (see also [bru_comp_eval()]).
#'
#' For "iid" models with `mapper = bm_index(n)`, `rnorm()` is used to
#' generate new realisations for indices greater than `n`, if accessed
#' via `<name>_eval(...)`.
#'
#' @return a `data.frame`, `sf`, or `Spatial*` object with predicted mean values
#'   and other summary statistics attached. Non-S4 object outputs have the class
#'   "bru_prediction" added at the front of the class list.
#' @example inst/examples/predict.bru.R

predict.bru <- function(object,
                        newdata = NULL,
                        formula = NULL,
                        n.samples = 100,
                        seed = 0L,
                        probs = c(0.025, 0.5, 0.975),
                        num.threads = NULL,
                        used = NULL,
                        drop = FALSE,
                        ...,
                        include = deprecated(),
                        exclude = deprecated()) {
  object <- bru_check_object_bru(object)

  # Convert data into list, data.frame or a Spatial object if not provided as
  # such
  if (is.character(newdata)) {
    newdata <- as.list(setNames(newdata, newdata))
  } else if (inherits(newdata, "formula")) {
    stop(paste0(
      "Formula supplied as data to predict.bru(). ",
      "Please check your argument order/names."
    ))
  }

  if (inherits(newdata, "sf")) {
    # Make sure coordinates have unified XY/XYZ/XYM/XYZM structure:
    newdata <- fmesher::fm_zm(newdata)
  }

  vals <- generate(
    object,
    newdata = newdata,
    formula = formula,
    n.samples = n.samples,
    seed = seed,
    num.threads = num.threads,
    include = include,
    exclude = exclude,
    used = used,
    ...
  )

  # Summarise

  newdata <- expand_to_dataframe(newdata)
  if (is.data.frame(vals[[1]])) {
    vals.names <- names(vals[[1]])
    covar <- intersect(vals.names, names(newdata))
    estim <- setdiff(vals.names, covar)
    smy <- list()

    for (nm in estim) {
      smy[[nm]] <-
        bru_summarise(
          lapply(
            vals,
            function(v) v[[nm]]
          ),
          x = vals[[1]][, covar, drop = FALSE],
          probs = probs
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
          if (NROW(newdata) == NROW(tmp)) {
            expand_to_dataframe(newdata, tmp)
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
      bru_log_warn(
        "Some generated list elements are unnamed"
      )
    }
    smy <- list()
    for (nm in vals.names) {
      tmp <-
        bru_summarise(
          data = lapply(
            vals,
            function(v) v[[nm]]
          ),
          probs = probs
        )
      if (!drop &&
        (NROW(newdata) == NROW(tmp))) {
        smy[[nm]] <- expand_to_dataframe(newdata, tmp)
      } else {
        smy[[nm]] <- tmp
      }
    }
  } else {
    tmp <- bru_summarise(data = vals, probs = probs)
    if (!drop &&
      (NROW(newdata) == NROW(tmp))) {
      smy <- expand_to_dataframe(newdata, tmp)
    } else {
      smy <- tmp
    }
  }

  if (!isS4(smy)) {
    class(smy) <- c("bru_prediction", class(smy))
  }
  smy
}


bru_generate_check_used_deprecation <- function(
  include,
  exclude,
  env = rlang::caller_env(),
  user_env = rlang::caller_env(2)
) {
  if (lifecycle::is_present(include)) {
    bru_log_message(
      paste0(
        "The `include` argument to `generate.bru()` is deprecated. ",
        "If auto-detection doesn't work, use ",
        "`used = bru_used(effect = include)`."
      ),
      verbosity = 1L
    )
    lifecycle::deprecate_warn(
      "2.12.0.9003",
      "generate(include)",
      "generate(used)",
      "If auto-detection doesn't work, use `bru_used(effect = include)`",
      env = env,
      user_env = user_env
    )
  } else {
    include <- NULL
  }
  if (lifecycle::is_present(exclude)) {
    bru_log_message(
      paste0(
        "The `exclude` argument to `generate.bru()` is deprecated. ",
        "If auto-detection doesn't work, use ",
        "`used = bru_used(effect_exclude = exclude)`."
      ),
      verbosity = 1L
    )
    lifecycle::deprecate_warn(
      "2.12.0.9003",
      "generate(exclude)",
      "generate(used)",
      paste0(
        "If auto-detection doesn't work, ",
        "use `bru_used(effect_exclude = exclude)`"
      ),
      env = env,
      user_env = user_env
    )
  } else {
    exclude <- NULL
  }
  list(include = include, exclude = exclude)
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
#' @export
#' @family sample generators
#' @param object A `bru` object obtained by calling [bru()].
#' @param newdata A `data.frame` or `SpatialPointsDataFrame` of covariates
#'   needed for sampling.
#' @param formula A formula where the right hand side defines an R expression
#' to evaluate for each generated sample. If `NULL`, the latent and
#' hyperparameter states are returned as named list elements.
#' See Details for more information.
#' @param n.samples Integer setting the number of samples to draw in order to
#' calculate the posterior statistics.
#' The default, 100, is rather low but provides a quick approximate result.
#' @param seed Random number generator seed passed on to
#'   `INLA::inla.posterior.sample`
#' @param num.threads Specification of desired number of threads for parallel
#' computations. Default NULL, leaves it up to INLA.
#' When seed != 0, overridden to "1:1"
#' @param used Either `NULL` or a [bru_used()] object.
#'   Default, `NULL`, uses auto-detection of used variables in the formula.
#' @param \dots additional, unused arguments.
#' @param include,exclude `r lifecycle::badge("deprecated")` If auto-detection
#' of used variables fails, use `used` instead.
#' @details
#' In addition to the component names (that give the effect of each component
#' evaluated for the input data), the suffix `_latent` variable name can be used
#' to directly access the latent state for a component, and the suffix function
#' `_eval` can be used to evaluate a component at other input values than the
#' expressions defined in the component definition itself, e.g.
#' `field_eval(cbind(x, y))` for a component that was defined with
#' `field(coordinates, ...)` (see also [bru_comp_eval()]).
#'
#' For "iid" models with `mapper = bm_index(n)`, `rnorm()` is used to
#' generate new realisations for indices greater than `n`, if accessed
#' via `<name>_eval(...)`.
#'
#' @return List of generated samples
#' @seealso [predict.bru]
#' @rdname generate

generate.bru <- function(object,
                         newdata = NULL,
                         formula = NULL,
                         n.samples = 100,
                         seed = 0L,
                         num.threads = NULL,
                         used = NULL,
                         ...,
                         include = deprecated(),
                         exclude = deprecated()) {
  object <- bru_check_object_bru(object)

  # Convert data into list, data.frame or a Spatial object if not provided as
  # such
  if (is.character(newdata)) {
    newdata <- as.list(setNames(newdata, newdata))
  } else if (inherits(newdata, "formula")) {
    stop(paste0(
      "Formula supplied as data to generate.bru(). ",
      "Please check your argument order/names."
    ))
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

    inc_exc <- bru_generate_check_used_deprecation(
      include,
      exclude
    )
    if (is.null(used)) {
      used <-
        bru_used(
          formula,
          effect = if (is.null(formula)) {
            character(0)
          } else {
            inc_exc$include
          },
          effect_exclude = inc_exc$exclude,
          latent = NULL
        )
    }
    used <- bru_used_update(used, labels = names(as_bru_comp_list(object)))

    pred <- new_bru_pred_expr(formula, used = used)

    if (inherits(newdata, "sf")) {
      # Make sure coordinates have unified XY/XYZ/XYM/XYZM structure:
      newdata <- fmesher::fm_zm(newdata)
    }

    vals <- evaluate_model(
      model = object$bru_info$model,
      state = state,
      data = newdata,
      predictor = pred
    )
    vals
  }
}


# Monte Carlo method for estimating posterior
#
# @export
# @param dfun A function returning a density for given x
# @param sfun A function providing samples from a posterior
# @param x Inital domain on which to perform the estimation. This will be
# adjusted as more samples are generated.
# @param samples An initial set of samples. Not required but will be used to
# estimate the inital domain \code{x} if \code{x} is \code{NULL}
# @param mcerr Monte Carlo error at which to stop the chain
# @param n Inital number of samples. This will be doubled for each iteration.
# @param discrete Set this to \code{TRUE} if the density is only defined for
#   integer \code{x}
# @param verbose Be verbose?

montecarlo.posterior <- function(dfun, sfun, x = NULL, samples = NULL,
                                 mcerr = 0.01, n = 100, discrete = FALSE,
                                 verbose = FALSE) {
  .Deprecated()
  xmaker <- function(hpd) {
    mid <- (hpd[2] + hpd[1]) / 2
    rg <- (hpd[2] - hpd[1]) / 2
    seq(mid - 1.2 * rg, mid + 1.2 * rg, length.out = 256)
  }
  xmaker2 <- function(hpd) {
    seq(hpd[1], hpd[2], length.out = 256)
  }

  initial.xmaker <- function(smp) {
    mid <- median(smp)
    rg <- (quantile(smp, 0.975) - quantile(smp, 0.25)) / 2
    seq(mid - 3 * rg, mid + 3 * rg, length.out = 256)
  }

  # Inital samples
  if (is.null(samples)) {
    samples <- sfun(n)
  }

  # Inital HPD
  if (is.null(x)) {
    x <- initial.xmaker(as.vector(unlist(samples)))
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
      cat(glue("hpd: ({min(x)}, {max(x)}), err = {err}, n = {n}"), "\n")
      # plot(x, lest, type = "l") ; lines(x, est, type = "l", col = "red")
    }

    # Convergence?
    if (err < mcerr) {
      converged <- TRUE
    } else {
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


#' Summarise and annotate data
#'
#' @export
#' @param data A list of samples, each either numeric or a `data.frame`
#' @param probs A numeric vector of probabilities with values in `[0, 1]`,
#'   passed to `stats::quantile`
#' @param x A `data.frame` of data columns that should be added to the summary
#'   data frame
#' @param cbind.only If TRUE, only `cbind` the samples and return a matrix where
#'   each column is a sample
#' @param max_moment integer, at least 2. Determines the largest moment
#'   order information to include in the output. If `max_moment > 2`,
#'   includes "skew" (skewness, `E[(x-m)^3/s^3]`), and
#'   if `max_moment > 3`, includes
#'   "ekurtosis" (excess kurtosis, `E[(x-m)^4/s^4] - 3`). Default 2.
#'   Note that the Monte Carlo variability of the `ekurtois` estimate may be
#'   large.
#' @return A `data.frame` or `Spatial[Points/Pixels]DataFrame` with summary
#'   statistics, "mean", "sd", `paste0("q", probs)`, "mean.mc_std_err",
#'   "sd.mc_std_err"
#'
#' @examples
#' bru_summarise(matrix(rexp(10000), 10, 1000), max_moment = 4, probs = NULL)
#'
#' @keywords internal
bru_summarise <- function(data, probs = c(0.025, 0.5, 0.975),
                          x = NULL, cbind.only = FALSE,
                          max_moment = 2) {
  if (is.list(data)) {
    data <- do.call(cbind, data)
  }
  N <- NCOL(data)
  if (cbind.only) {
    smy <- data.frame(data)
    colnames(smy) <- paste0("sample.", seq_len(N))
  } else {
    if (length(probs) == 0) {
      smy <- data.frame(
        rowMeans(data, na.rm = TRUE),
        apply(data, MARGIN = 1, sd, na.rm = TRUE)
      )
      qs_names <- NULL
    } else if (length(probs) == 1) {
      smy <- data.frame(
        rowMeans(data, na.rm = TRUE),
        apply(data, MARGIN = 1, sd, na.rm = TRUE),
        as.matrix(apply(data,
          MARGIN = 1, quantile, probs = probs, na.rm = TRUE,
          names = FALSE
        ))
      )
      qs_names <- paste0("q", probs)
    } else {
      smy <- data.frame(
        rowMeans(data, na.rm = TRUE),
        apply(data, MARGIN = 1, sd, na.rm = TRUE),
        t(apply(data,
          MARGIN = 1, quantile, probs = probs, na.rm = TRUE,
          names = FALSE
        ))
      )
      qs_names <- paste0("q", probs)
    }
    colnames(smy) <- c("mean", "sd", qs_names)
    # For backwards compatibility, add a median column:
    if (any(qs_names == "q0.5")) {
      smy[["median"]] <- smy[["q0.5"]]
    }

    # Get 3rd and 4rt central moments
    # Use 1/N normalisation of the sample sd
    skew <- rowMeans(
      ((data - smy$mean) / smy$sd)^3 * (N / (N - 1))^3,
      na.rm = TRUE
    )
    if (max_moment >= 3) {
      smy[["skew"]] <- skew
    }
    # eK + 3 >= skew^2 + 1
    # eK >= skew^2 - 2
    # Use 1/N normalisation of the sample sd
    ekurtosis <- pmax(
      skew^2 - 2,
      rowMeans(
        ((data - smy$mean) / smy$sd)^4 * (N / (N - 1))^4 - 3,
        na.rm = TRUE
      )
    )

    # Add Monte Carlo standard errors
    # Var(s) \approx (eK + 2) \sigma^2 / (4 n):
    # +2 replaced by 3-(n-3)/(n-1) = 2n/(n-1), from Rao 1973, p438
    sd.mc_std_err <-
      sqrt(pmax(0, ekurtosis + 2 * N / (N - 1))) *
        smy[["sd"]] / sqrt(4 * N)
    # Include sd MC error in estimate of mean MC error:
    smy[["mean.mc_std_err"]] <- (smy[["sd"]] + sd.mc_std_err * 2) / sqrt(N)
    smy[["sd.mc_std_err"]] <- sd.mc_std_err
    if (max_moment >= 3) {
      smy[["skew"]] <- skew
    }
    if (max_moment >= 4) {
      smy[["ekurtosis"]] <- ekurtosis
    }
  }
  if (!is.null(x)) {
    smy <- expand_to_dataframe(x, smy)
  }
  smy
}


lin_predictor <- function(lin, state) {
  do.call(
    c,
    lapply(
      lin,
      function(x) {
        as.vector(ibm_eval(x, state = state))
      }
    )
  )
}
nonlin_predictor <- function(param, state) {
  do.call(
    c,
    lapply(
      seq_along(param[["lhoods"]]),
      function(lh_idx) {
        as.vector(
          evaluate_model(
            model = param[["model"]],
            data = param[["lhoods"]][[lh_idx]][["data"]],
            data_extra = param[["lhoods"]][[lh_idx]][["data_extra"]],
            input = param[["input"]][[lh_idx]],
            state = list(state),
            comp_simple = param[["comp_simple"]][[lh_idx]],
            predictor = bru_pred_expr(param[["lhoods"]][[lh_idx]],
              format = "expr"
            ),
            format = "matrix",
            n_pred = bru_response_size(param[["lhoods"]][[lh_idx]])
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

line_search_opt_target <- function(x, param) {
  (x - 1)^2 * param[1] + 2 * (x - 1) * x^2 * param[2] + x^4 * param[3]
}

line_search_opt_target_exact <- function(x, param, nonlin_param) {
  state <- scale_state(param$state0, param$state1, x)
  nonlin <- nonlin_predictor(
    param = nonlin_param,
    state = state
  )
  sum((nonlin - param$lin)^2 * param$weights)
}


# @title FUNCTION_TITLE
# @description FUNCTION_DESCRIPTION
# @param model PARAM_DESCRIPTION
# @param lhoods PARAM_DESCRIPTION
# @param lin PARAM_DESCRIPTION
# @param state0 PARAM_DESCRIPTION
# @param state PARAM_DESCRIPTION
# @param A PARAM_DESCRIPTION
# @param weights PARAM_DESCRIPTION
# @param options PARAM_DESCRIPTION
# @return OUTPUT_DESCRIPTION
# @details DETAILS
# @examples
# \dontrun{
# if (interactive()) {
#  #EXAMPLE1
#  }
# }
# @export
# @rdname bru_line_search

bru_line_search <- function(model,
                            lhoods,
                            lin,
                            state0,
                            state,
                            input,
                            comp_lin,
                            comp_simple,
                            weights = 1,
                            options) {
  # Metrics ----
  pred_scalprod <- function(delta1, delta2) {
    if (!all(is.finite(delta1)) || !all(is.finite(delta2))) {
      Inf
    } else {
      sum(delta1 * delta2 * weights)
    }
  }
  pred_norm2 <- function(delta) {
    pred_scalprod(delta, delta)
  }
  pred_norm <- function(delta) {
    pred_scalprod(delta, delta)^0.5
  }

  # Is line search enabled? ----
  if (length(options$bru_method$search) == 0) {
    return(
      list(
        active = FALSE,
        step_scaling = 1,
        state = state
      )
    )
  }

  # Initialise ----
  if (is.null(weights)) {
    bru_log_warn(
      paste0(
        "NULL weights detected for line search. Using weights = 1 instead.",
        "\n\tThis is a bug in the inlabru package. Please notify the developer."
      )
    )
    weights <- 1
  }

  fact <- options$bru_method$factor

  nonlin_param <- list(
    model = model,
    lhoods = lhoods,
    input = input,
    comp_simple = comp_simple
  )

  state1 <- state
  lin_pred0 <- lin_predictor(lin, state0)
  lin_pred1 <- lin_predictor(lin, state1)
  nonlin_pred <- nonlin_predictor(
    param = nonlin_param,
    state = state1
  )

  if (length(lin_pred1) != length(nonlin_pred)) {
    bru_log_warn(
      paste0(
        "Please notify the inlabru package developer:",
        "\n\tThe line search linear and nonlinear predictors have ",
        "different lengths.",
        "\n\tThis should not happen!"
      )
    )
  }

  step_scaling <- 1

  if (isTRUE(options[["bru_debug"]])) {
    df_debug <- data.frame(
      idx = seq_along(lin_pred0),
      lin0 = lin_pred0,
      lin1 = lin_pred1,
      nonlin1 = nonlin_pred,
      weights = weights
    )
  }

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
  maximum_step <- options$bru_method$max_step

  finite_active <- 0
  contract_active <- 0
  expand_active <- 0

  # Contraction methods ----
  if (do_finite || do_contract) {
    nonfin <- !all(is.finite(nonlin_pred))
    norm0 <- pred_norm(nonlin_pred - lin_pred0)
    norm1 <- pred_norm(nonlin_pred - lin_pred1)

    while ((do_finite && nonfin) ||
      (do_contract && ((norm0 > norm01 * fact) || (norm1 > norm01 * fact)))) {
      if (do_finite && nonfin) {
        finite_active <- finite_active - 1
      } else {
        contract_active <- contract_active - 1
      }
      step_scaling <- step_scaling / fact
      state <- scale_state(state0, state1, step_scaling)
      nonlin_pred <- nonlin_predictor(
        param = nonlin_param,
        state = state
      )
      nonfin <- !all(is.finite(nonlin_pred))
      norm0 <- pred_norm(nonlin_pred - lin_pred0)
      norm1 <- pred_norm(nonlin_pred - lin_pred1)

      bru_log_message(
        glue(
          "iinla: Step rescaling: {signif(100 * step_scaling, 4)}%, Contract",
          " (",
          "norm0 = {signif(norm0, 4)}, ",
          "norm1 = {signif(norm1, 4)}, ",
          "norm01 = {signif(norm01, 4)}",
          ")"
        ),
        verbosity = 3
      )

      if (step_scaling < 2^-10) {
        break
      }
    }
  }

  # Expansion method ----
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

      nonlin_pred <- nonlin_predictor(
        param = nonlin_param,
        state = state
      )
      norm1_prev <- norm1
      norm0 <- pred_norm(nonlin_pred - lin_pred0)
      norm1 <- pred_norm(nonlin_pred - lin_pred1)
      overstep <-
        (norm1 > norm1_prev) || !is.finite(norm1)

      bru_log_message(
        glue(
          "iinla: Step rescaling: {signif(100 * step_scaling, 3)}%, Expand",
          " (",
          "norm0 = {signif(norm0, 4)}, ",
          "norm1 = {signif(norm1, 4)}, ",
          "norm01 = {signif(norm01, 4)}",
          ")"
        ),
        verbosity = 3
      )

      if (step_scaling > 2^10) {
        break
      }
    }
    # Overstep ----
    if (((step_scaling > 1) && overstep) || !is.finite(norm1)) {
      expand_active <- expand_active - 1
      step_scaling <- step_scaling / fact
      state <- scale_state(state0, state1, step_scaling)
      nonlin_pred <- nonlin_predictor(
        param = nonlin_param,
        state = state
      )
      norm0 <- pred_norm(nonlin_pred - lin_pred0)
      norm1 <- pred_norm(nonlin_pred - lin_pred1)

      bru_log_message(
        glue(
          "iinla: Step rescaling: {signif(100 * step_scaling, 3)}%, Overstep",
          " (",
          "norm0 = {signif(norm0, 4)}, ",
          "norm1 = {signif(norm1, 4)}, ",
          "norm01 = {signif(norm01, 4)}",
          ")"
        ),
        verbosity = 3
      )
    }
  }

  # Optimisation method ----
  if (do_optimise) {
    lin_pred <- lin_predictor(lin, state)
    delta_lin <- (lin_pred1 - lin_pred0)
    delta_nonlin <- (nonlin_pred - lin_pred) / step_scaling^2
    alpha <-
      optimise(
        line_search_opt_target,
        step_scaling * c(1 / fact^2, fact),
        param = c(
          pred_norm2(delta_lin),
          pred_scalprod(delta_lin, delta_nonlin),
          pred_norm2(delta_nonlin)
        )
      )
    step_scaling_opt_approx <- alpha$minimum

    state_opt <- scale_state(state0, state1, step_scaling_opt_approx)
    nonlin_pred_opt <- nonlin_predictor(
      param = nonlin_param,
      state = state_opt
    )
    norm0_opt <- pred_norm(nonlin_pred_opt - lin_pred0)
    norm1_opt <- pred_norm(nonlin_pred_opt - lin_pred1)

    bru_log_message(
      glue(
        "iinla: Step rescaling: {signif(100 * step_scaling_opt_approx, 4)}%,",
        " Approx Optimisation",
        " (",
        "norm0 = {signif(norm0_opt, 4)}, ",
        "norm1 = {signif(norm1_opt, 4)}, ",
        "norm01 = {signif(norm01, 4)}",
        ")"
      ),
      verbosity = 3
    )

    if (norm1_opt > norm01) {
      bru_log_message(
        glue(
          "iinla: norm1_opt > |delta|: ",
          "{signif(norm1_opt, 4)}",
          " > ",
          "{signif(norm01, 4)}"
        ),
        verbosity = 3
      )
    }

    if ((norm1_opt > norm01) ||
      identical(options$bru_method$line_opt_method, "full")) {
      alpha <-
        optimise(
          line_search_opt_target_exact,
          step_scaling * c(0, fact),
          param = list(
            lin = lin_pred1,
            state0 = state0,
            state1 = state1,
            weights = weights
          ),
          nonlin_param = nonlin_param
        )

      step_scaling_opt <- alpha$minimum
      state_opt <- scale_state(state0, state1, step_scaling_opt)
      nonlin_pred_opt <- nonlin_predictor(
        param = nonlin_param,
        state = state_opt
      )
      norm0_opt <- pred_norm(nonlin_pred_opt - lin_pred0)
      norm1_opt <- pred_norm(nonlin_pred_opt - lin_pred1)

      bru_log_message(
        glue(
          "iinla: Step rescaling: {signif(100 * step_scaling_opt, 4)}%,",
          " Optimisation",
          " (",
          "norm0 = {signif(norm0_opt, 4)}, ",
          "norm1 = {signif(norm1_opt, 4)}, ",
          "norm01 = {signif(norm01, 4)}",
          ")"
        ),
        verbosity = 3
      )

      if (norm1_opt > norm01) {
        bru_log_message(
          glue(
            "iinla: norm1_opt > |delta|: ",
            "{signif(norm1_opt, 4)}",
            " > ",
            "{signif(norm01, 4)}"
          ),
          verbosity = 3
        )
      }
    } else {
      step_scaling_opt <- step_scaling_opt_approx
    }

    if (norm1_opt < norm1) {
      step_scaling <- step_scaling_opt
      state <- state_opt
      nonlin_pred <- nonlin_pred_opt
      norm0 <- norm0_opt
      norm1 <- norm1_opt
    } else {
      bru_log_message(
        paste0("iinla: Optimisation did not improve on previous solution."),
        verbosity = 3
      )
    }
  }

  bru_log_message(
    glue("
      iinla: |lin1-lin0| = {signif(pred_norm(lin_pred1 - lin_pred0), 4)}
        <eta-lin1,delta>/|delta| = {signif(rel1, 4)}
        |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = {signif(rel2, 4)}",
      rel1 = pred_scalprod(nonlin_pred - lin_pred1, lin_pred1 - lin_pred0) /
        pred_norm(lin_pred1 - lin_pred0),
      rel2 = pred_norm(nonlin_pred - lin_pred0 - (lin_pred1 - lin_pred0) *
        pred_scalprod(lin_pred1 - lin_pred0, nonlin_pred - lin_pred0) /
        pred_norm2(lin_pred1 - lin_pred0))
    ),
    verbosity = 4
  )

  # Prevent long steps ----
  if (step_scaling > maximum_step) {
    step_scaling <- maximum_step
    state <- scale_state(state0, state1, step_scaling)
    nonlin_pred <- nonlin_predictor(
      param = nonlin_param,
      state = state
    )
    norm0 <- pred_norm(nonlin_pred - lin_pred0)
    norm1 <- pred_norm(nonlin_pred - lin_pred1)

    bru_log_message(
      glue(
        "iinla: Step rescaling: ",
        "{signif(100 * step_scaling, 3)}%, ",
        "Maximum step length",
        " (",
        "norm0 = {signif(norm0, 4)}, ",
        "norm1 = {signif(norm1, 4)}, ",
        "norm01 = {signif(norm01, 4)}",
        ")"
      ),
      verbosity = 3
    )
  }

  active <-
    (finite_active + contract_active + expand_active != 0) ||
      (abs(step_scaling - 1) > 0.001)

  if (active && (options$bru_verbose <= 2)) {
    bru_log_message(
      glue(
        "iinla: Step rescaling: ",
        "{signif(100 * step_scaling, 3)}%",
        " (",
        "norm0 = {signif(norm0, 4)}, ",
        "norm1 = {signif(norm1, 4)}, ",
        "norm01 = {signif(norm01, 4)}",
        ")"
      ),
      verbosity = 2
    )
  }

  if (isTRUE(options[["bru_debug"]])) {
    df_debug$nonlinopt <- nonlin_pred

    requireNamespace("ggplot2")
    requireNamespace("patchwork")
    pl1 <- ggplot2::ggplot(df_debug) +
      ggplot2::geom_point(ggplot2::aes(.data$idx,
        (.data$lin0 - .data$lin1),
        col = "start"
      )) +
      ggplot2::geom_point(ggplot2::aes(.data$idx,
        (.data$lin1 - .data$lin1),
        col = "inla"
      )) +
      ggplot2::geom_point(ggplot2::aes(.data$idx,
        (.data$nonlin1 - .data$lin1),
        col = "full"
      )) +
      ggplot2::geom_point(ggplot2::aes(.data$idx,
        (.data$nonlinopt - .data$lin1),
        col = "opt"
      )) +
      ggplot2::geom_abline(slope = 0, intercept = 0) +
      ggplot2::scale_color_discrete(
        breaks = c("start", "inla", "full", "opt")
      ) +
      ggplot2::ylab("eta - lin1")
    pl2 <- ggplot2::ggplot(df_debug) +
      ggplot2::geom_point(ggplot2::aes(
        .data$idx,
        (.data$lin0 - .data$lin1) * .data$weights^0.5,
        col = "start"
      )) +
      ggplot2::geom_point(ggplot2::aes(
        .data$idx,
        (.data$lin1 - .data$lin1) * .data$weights^0.5,
        col = "inla"
      )) +
      ggplot2::geom_point(ggplot2::aes(
        .data$idx,
        (.data$nonlin1 - .data$lin1) * .data$weights^0.5,
        col = "full"
      )) +
      ggplot2::geom_point(ggplot2::aes(
        .data$idx,
        (.data$nonlinopt - .data$lin1) * .data$weights^0.5,
        col = "opt"
      )) +
      ggplot2::geom_abline(slope = 0, intercept = 0) +
      ggplot2::scale_color_discrete(
        breaks = c("start", "inla", "full", "opt")
      ) +
      ggplot2::ylab("(eta - lin1) * sqrt(w)")
    pl3 <- ggplot2::ggplot(df_debug) +
      ggplot2::geom_point(ggplot2::aes(.data$idx,
        .data$lin0,
        col = "start"
      )) +
      ggplot2::geom_point(ggplot2::aes(.data$idx,
        .data$lin1,
        col = "inla"
      )) +
      ggplot2::geom_point(ggplot2::aes(.data$idx,
        .data$nonlin1,
        col = "full"
      )) +
      ggplot2::geom_point(ggplot2::aes(.data$idx,
        .data$nonlinopt,
        col = "opt"
      )) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          as.numeric(.data$idx),
          ymin = .data$lin1 - 2 * .data$weights^-0.5,
          ymax = .data$lin1 + 2 * .data$weights^-0.5
        ),
        alpha = 0.1
      ) +
      ggplot2::geom_abline(slope = 0, intercept = 0) +
      ggplot2::scale_color_discrete(
        breaks = c("start", "inla", "full", "opt")
      ) +
      ggplot2::ylab("eta")
    pl4 <- ggplot2::ggplot(data.frame(
      idx = seq_along(unlist(state0)),
      state0 = unlist(state0),
      state1 = unlist(state1)[names(unlist(state0))],
      state = unlist(state)[names(unlist(state0))]
    )) +
      ggplot2::geom_point(
        ggplot2::aes(.data$idx, .data$state0, col = "start")
      ) +
      ggplot2::geom_point(
        ggplot2::aes(.data$idx, .data$state1, col = "inla")
      ) +
      ggplot2::geom_point(
        ggplot2::aes(.data$idx, .data$state1, col = "full")
      ) +
      ggplot2::geom_point(
        ggplot2::aes(.data$idx, .data$state, col = "opt")
      ) +
      ggplot2::scale_color_discrete(
        breaks = c("start", "inla", "full", "opt")
      ) +
      ggplot2::ylab("state")
    pl1234 <-
      ((pl1 | pl2) / (pl3 | pl4)) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "right")

    # Compute deviations in the delta = lin1-lin0 direction and the orthogonal
    # direction
    # delta = lin1-lin0
    delta <- lin_pred1 - lin_pred0
    delta_norm <- norm01

    # <eta-lin0, delta>/|delta|
    along_opt <- pred_scalprod(nonlin_pred - lin_pred0, delta) / delta_norm
    # |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>|
    ortho_opt <- pred_norm(
      nonlin_pred - lin_pred0 - delta * along_opt / delta_norm
    )

    step_scaling_range <- seq(0, step_scaling * fact, length.out = 20)
    along <- ortho <- distance <- numeric(length(step_scaling_range))
    for (iii in seq_along(step_scaling_range)) {
      step_scaling_ <- step_scaling_range[iii]
      state_ <- scale_state(state0, state1, step_scaling_)
      nonlin_pred_ <- nonlin_predictor(
        param = nonlin_param,
        state = state_
      )
      # <eta-lin0, delta>/|delta|
      along[iii] <- pred_scalprod(nonlin_pred_ - lin_pred0, delta) / delta_norm
      # |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>|
      ortho[iii] <- pred_norm(
        nonlin_pred_ - lin_pred0 - delta * along[iii] / delta_norm
      )
      distance[iii] <-
        line_search_opt_target_exact(
          step_scaling_,
          param = list(
            lin = lin_pred1,
            state0 = state0,
            state1 = state1,
            weights = weights
          ),
          nonlin_param = nonlin_param
        )
    }


    pl0 <-
      ggplot2::ggplot(
        data =
          data.frame(
            along = along,
            ortho = ortho,
            distance = distance^0.5,
            what = "eta"
          ),
        mapping = ggplot2::aes(.data$along, .data$ortho, col = .data$what)
      ) +
      ggplot2::geom_path() +
      ggplot2::geom_point(ggplot2::aes(size = .data$distance)) +
      ggplot2::geom_point(
        data = data.frame(
          along = c(0, along_opt, delta_norm),
          ortho = c(0, ortho_opt, 0),
          what = c("eta0", "eta_opt", "eta1")
        )
      )

    browser()
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
        glue("{x}.{seq_along(state[[x]])}")
      }
    }
  )
  names(nm) <- names(state)
  nm
}
tidy_state <- function(state, value_name = "value") {
  df <- data.frame(
    effect = rep(names(state), lengths(state)),
    index = unlist(lapply(lengths(state), function(x) seq_len(x))),
    THEVALUE = unlist(state)
  )
  nm <- names(df)
  nm[nm == "THEVALUE"] <- value_name
  names(df) <- nm
  rownames(df) <- NULL
  df
}
tidy_states <- function(states, value_name = "value", id_name = "iteration") {
  df <- lapply(states, function(x) tidy_state(x, value_name = value_name))
  id <- rep(seq_along(states), each = nrow(df[[1]]))
  df <- do.call(rbind, df)
  df[[id_name]] <- id
  df
}


bru_timer_do <- function(timer, task = NULL, iteration = NA_integer_) {
  if (is.null(timer)) {
    timer <- data.frame(
      Task = character(0),
      Iteration = integer(0),
      Time = as.difftime(numeric(0), units = "secs"),
      System = as.difftime(numeric(0), units = "secs"),
      Elapsed = as.difftime(numeric(0), units = "secs")
    )
  }
  tag <- list(
    time = proc.time(),
    task = task,
    iter = iteration
  )
  old_tag <- attr(timer, "bru_timer_tag", exact = TRUE)
  if (!is.null(old_tag)) {
    ptm <- tag$time - old_tag$time
    if (!is.na(ptm[4])) {
      ptm[1] <- ptm[1] + ptm[4]
    }
    if (!is.na(ptm[5])) {
      ptm[2] <- ptm[2] + ptm[5]
    }

    if (is.null(old_tag$task)) {
      old_tag$task <- "Unknown"
    }
    if (is.null(old_tag$iter)) {
      old_tag$iter <- NA_integer_
    }
    timer <- rbind(
      timer,
      data.frame(
        Task = old_tag$task,
        Iteration = old_tag$iter,
        Time = as.difftime(unname(ptm[1]), units = "secs"),
        System = as.difftime(unname(ptm[2]), units = "secs"),
        Elapsed = as.difftime(unname(ptm[3]), units = "secs")
      )
    )
  }
  attr(timer, "bru_timer_tag") <- if (is.null(task)) {
    NULL
  } else {
    tag
  }
  timer
}
bru_timer_done <- function(timer) {
  bru_timer_do(timer, task = NULL)
}


#' Iterated INLA
#'
#' This is an internal wrapper for iterated runs of `INLA::inla`.
#' For nonlinear models, a linearisation is done with
#' `bru_compute_linearisation`, with a line search method between each
#' iteration. The `INLA::inla.stack` information is setup by [bru_make_stack()].
#'
#' @export
#' @param model A [bru_model] object
#' @param lhoods A list of likelihood objects from [bru_obs()]
#' @param inputs Optional pre-computed  list of per-likelihood component
#'   evaluations, from [bru_input.bru_obs_list()].
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

iinla <- function(model, lhoods, inputs = NULL, initial = NULL, options) {
  options <- bru_call_options(options)
  bru_options_set_local(options, .reset = TRUE)

  timings <- bru_timer_do(NULL, "Preprocess", 1L)

  inla.options <- bru_options_inla(options)

  bru_log_bookmark("iinla")
  orig_timings <- NULL
  original_log <- character(0) # Updated further below
  # Local utility method for collecting information object:
  collect_misc_info <- function(...) {
    timings <- bru_timer_done(timings)
    if (is.null(orig_track)) {
      track_df <- list()
      for (label in names(states[[1]])) {
        if (length(states[[1]][[label]]) > 0) {
          track_df[[label]] <-
            data.frame(
              effect = label,
              index = seq_along(states[[1]][[label]]),
              iteration = 0L,
              mode = NA_real_,
              sd = NA_real_,
              new_linearisation = states[[1]][[label]]
            )
        }
      }
      orig_track <- do.call(rbind, track_df)
    }

    if (is.null(orig_inla_track)) {
      orig_inla_track <- tibble::tibble(
        iteration = integer(0),
        iteration_part = numeric(0),
        f = numeric(0),
        nfunc = integer(0),
        nfunc_total = integer(0)
      )
      nfunc_offset <- 0L
    } else {
      nfunc_offset <- max(orig_inla_track$nfunc_total)
    }

    track_names <- if (length(track) > 0) {
      names(track[[1]])
    } else {
      c(
        "effect",
        "index",
        "iteration",
        "mode",
        "sd",
        "new_linearisation"
      )
    }

    track <-
      if (is.null(orig_track) ||
        setequal(names(orig_track), track_names)) {
        do.call(dplyr::bind_rows, c(list(orig_track), track))
      } else {
        track <- do.call(dplyr::bind_rows, track)
        original_names <- names(orig_track)
        new_names <- names(track)
        for (nn in setdiff(new_names, original_names)) {
          orig_track[[nn]] <- NA
        }
        for (nn in setdiff(original_names, new_names)) {
          track[[nn]] <- NA
        }
        dplyr::bind_rows(orig_track, track)
      }

    inla_track <- {
      offsets <- cumsum(vapply(
        seq_along(inla_track),
        function(k) {
          if (nrow(inla_track[[k]]) > 0) {
            max(inla_track[[k]]$nfunc)
          } else {
            0L
          }
        },
        0L
      ))
      offsets <- nfunc_offset + c(0, offsets)
      for (k in seq_along(inla_track)) {
        if (nrow(inla_track[[k]]) > 0) {
          inla_track[[k]]$nfunc_total <- inla_track[[k]]$nfunc + offsets[k]
        }
      }
      ## In case the number of theta values changes during the iterations,
      ## e.g.
      ## bru_rerun(..., options = list(control.mode = list(fixed = TRUE)))
      orig_n <- NCOL(orig_inla_track[["theta"]])
      if (length(inla_track) > 0) {
        new_n <- vapply(inla_track, function(x) NCOL(x[["theta"]]), 0L)
      } else {
        new_n <- orig_n
      }
      if (max(new_n) > orig_n) {
        if (is.null(orig_inla_track[["theta"]])) {
          orig_inla_track[["theta"]] <- matrix(NA_real_,
            nrow = NROW(orig_inla_track),
            ncol = max(new_n)
          )
        } else {
          orig_inla_track[["theta"]] <-
            cbind(
              orig_inla_track[["theta"]],
              matrix(NA_real_,
                nrow = NROW(orig_inla_track),
                ncol = max(new_n) - orig_n
              )
            )
        }
      }
      orig_n <- NCOL(orig_inla_track[["theta"]])
      if (any(new_n != orig_n)) {
        inla_track <- lapply(
          inla_track,
          function(x) {
            if (NCOL(x[["theta"]]) < orig_n) {
              if (is.null(x[["theta"]])) {
                x[["theta"]] <-
                  matrix(NA_real_,
                    nrow = NROW(x),
                    ncol = orig_n
                  )
              } else {
                x[["theta"]] <- cbind(
                  x[["theta"]],
                  matrix(NA_real_,
                    nrow = NROW(x),
                    ncol = orig_n - NCOL(x[["theta"]])
                  )
                )
              }
            }
            x
          }
        )
      }
      inla_track <- dplyr::bind_rows(inla_track)
      inla_track$iteration_part <-
        (seq_len(NROW(inla_track)) - 1L) / NROW(inla_track)
      dplyr::bind_rows(orig_inla_track, inla_track)
    }

    list(
      log = c(original_log, bru_log()["iinla"]),
      states = states,
      inla_stack = stk,
      track = track,
      inla_track = inla_track,
      timings = {
        iteration_offset <- if (is.null(orig_timings)) {
          0L
        } else {
          max(c(0L, orig_timings$Iteration), na.rm = TRUE)
        }
        rbind(
          orig_timings,
          data.frame(
            Task = timings$Task,
            Iteration = timings$Iteration + iteration_offset,
            Time = timings$Time,
            System = timings$System,
            Elapsed = timings$Elapsed
          )
        )
      },
      ...
    )
  }

  if (!is.null(inla.options[["offset"]])) {
    stop(paste0(
      "An offset option was specified which may interfere with the ",
      "inlabru model construction.\n",
      "Please use an explicit constant component instead; ",
      "e.g. ~ myoffset(value, model = 'const')"
    ))
  }

  # Initialise required local options
  inla.options <- modifyList(
    inla.options,
    list(
      control.mode = list(),
      control.predictor = list(compute = TRUE),
      control.compute = list()
    )
  )

  # Extract the family of each likelihood
  family <- bru_obs_inla_family(lhoods)

  # Extract the control.family information for each likelihood
  inla.options[["control.family"]] <-
    bru_obs_control_family(lhoods, inla.options[["control.family"]])

  # Combine the control.gcpo information from all the likelihoods
  control.gcpo.combined <-
    bru_obs_control_gcpo(
      lhoods,
      inla.options[["control.compute"]][["control.gcpo"]]
    )
  inla.options[["control.compute"]][["control.gcpo"]] <- control.gcpo.combined

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
    comp_lst <- as_bru_comp_list(model)
    state <- initial[intersect(names(comp_lst), names(initial))]
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
    states <- list(state)
    result <- NULL
  } else {
    stop("Unknown previous result information class")
  }
  old.result <- result

  # Preserve old log output
  original_log <- bru_log(
    if (is.null(old.result)) {
      character(0)
    } else if (is.null(old.result[["bru_iinla"]])) {
      character(0)
    } else if (is.null(old.result[["bru_iinla"]][["log"]])) {
      character(0)
    } else {
      old.result[["bru_iinla"]][["log"]]
    }
  )
  bru_log_message(
    "iinla: Start",
    verbosity = 3
  )

  # Track variables
  track <- list()
  if (is.null(old.result[["bru_iinla"]][["track"]])) {
    orig_track <- NULL
    track_size <- 0
  } else {
    orig_track <- old.result[["bru_iinla"]][["track"]]
    track_size <- max(orig_track[["iteration"]])
  }
  inla_track <- list()
  if (is.null(old.result[["bru_iinla"]][["inla_track"]]) ||
    (NROW(old.result[["bru_iinla"]][["inla_track"]]) == 0L)) {
    orig_inla_track <- NULL
    inla_track_size <- 0L
  } else {
    orig_inla_track <- old.result[["bru_iinla"]][["inla_track"]]
    inla_track_size <- max(orig_inla_track[["iteration"]])
  }

  # Preserve old timings
  if (!is.null(old.result[["bru_iinla"]][["timings"]])) {
    orig_timings <- old.result[["bru_iinla"]][["timings"]]
  }

  if (is.null(inputs)) {
    bru_log_message(
      "iinla: Evaluate component inputs",
      verbosity = 3
    )
    inputs <- bru_input(model, lhoods = lhoods)
  }
  bru_log_message(
    "iinla: Evaluate component linearisations",
    verbosity = 3
  )
  comp_lin <- ibm_linear(
    model,
    input = inputs,
    state = states[[length(states)]],
    inla_f = TRUE
  )
  bru_log_message(
    "iinla: Evaluate component simplifications",
    verbosity = 3
  )
  comp_simple <- ibm_simplify(
    model,
    input = inputs,
    inla_f = TRUE
  )
  bru_log_message(
    "iinla: Evaluate predictor linearisation",
    verbosity = 3
  )
  lin <- bru_compute_linearisation(
    model,
    lhoods = lhoods,
    input = inputs,
    state = states[[length(states)]],
    comp_simple = comp_simple,
    options = options
  )

  do_line_search <- (length(options[["bru_method"]][["search"]]) > 0)
  if (do_line_search) {
    # Always compute linearisation (above)
  }

  bru_log_message(
    "iinla: Construct inla stack",
    verbosity = 3
  )
  # Initial stack
  idx <- bru_index(model, used = bru_used(lhoods))
  stk <- bru_make_stack(lhoods, lin, idx)

  if (utils::packageVersion("INLA") <= "24.06.02") {
    stk.data <- INLA::inla.stack.data(stk)
  } else {
    stk.data <- INLA::inla.stack.data(stk, .response.name = "BRU.response")
  }
  inla.options$control.predictor$A <- INLA::inla.stack.A(stk)

  N_response <- sum(bru_response_size(lhoods))
  N_predictor <- NROW(inla.options$control.predictor$A)
  if (N_response != N_predictor) {
    msg <- glue(
      "The total number of response values (N={N_response})",
      " and predictor values (N={N_predictor})",
      " do not match.\n",
      "  This is likely due to a mistake in the component or predictor",
      " constructions.",
      .trim = FALSE
    )
    if ((N_response > 1L) && (N_predictor == 1L)) {
      msg <- paste0(
        msg,
        "\nPerhaps you only have components with scalar inputs?"
      )
    }
    bru_log_message(msg, verbosity = 1)
    stop(msg)
  }

  bru_log_message(
    "iinla: Model initialisation completed",
    verbosity = 3
  )

  k <- 1L
  interrupt <- FALSE
  line_search <- list(
    active = FALSE,
    step_scaling = 1,
    state = states[[length(states)]]
  )

  if ((!is.null(result) && !is.null(result$mode))) {
    previous_x <- result$mode$x
  } else {
    # TODO: construct from the initial linearisation state instead
    previous_x <- 0
  }

  track_df <- NULL
  do_final_integration <- (options$bru_max_iter == 1)
  do_final_theta_no_restart <- FALSE
  while (!interrupt) {
    if (k > 1L) {
      timings <- bru_timer_do(timings, "Preprocess", k)
    }

    if ((k >= options$bru_max_iter) && !do_final_integration) {
      do_final_integration <- TRUE
      bru_log_message(
        "iinla: Maximum iterations reached, running final INLA integration.",
        verbosity = 1
      )
    }

    # When running multiple times propagate theta
    if ((k > 1) || (!is.null(result) && !is.null(result$mode))) {
      inla.options[["control.mode"]]$restart <- TRUE
      inla.options[["control.mode"]]$theta <- result$mode$theta
      inla.options[["control.mode"]]$x <- (
        previous_x +
          (result$mode$x - previous_x) * line_search[["step_scaling"]]
      )
      result_indexing <- inla_result_latent_idx(result)
      # Reset the predictor values, to avoid spurious iteration effects.
      inla.options[["control.mode"]]$x[result_indexing$APredictor] <- 0
      inla.options[["control.mode"]]$x[result_indexing$Predictor] <- 0

      inla.options[["control.inla"]]$use.directions <-
        result$misc$opt.directions

      # Further settings, from inla.rerun
      inla.options[["control.inla"]]$optimise.strategy <- "plain"

      inla.options[["control.inla"]]$step.factor <- 1
      inla.options[["control.inla"]]$tolerance.step <- 1e-10

      if (identical(
        result$.args$control.inla[["stencil"]],
        INLA::inla.set.control.inla.default()$stencil
      )) {
        inla.options$control.inla$stencil <- 9
      }

      h.def <- INLA::inla.set.control.inla.default()$h
      if (identical(
        result$.args$control.inla[["h"]], h.def
      )) {
        inla.options$control.inla$h <-
          0.2 * result$.args$control.inla[["h"]]
      } else {
        inla.options$control.inla$h <-
          max(0.1 * h.def, 0.75 * result$.args$control.inla[["h"]])
      }

      tol.def <- INLA::inla.set.control.inla.default()$tolerance
      if (identical(result$.args$control.inla[["tolerance"]], tol.def)) {
        inla.options$control.inla$tolerance <-
          result$.args$control.inla[["tolerance"]] * 0.01
      } else {
        inla.options$control.inla$tolerance <-
          max(tol.def^2, result$.args$control.inla[["tolerance"]] * 0.1)
      }
    }
    if ((!is.null(result) && !is.null(result$mode))) {
      previous_x <- result$mode$x
    }

    bru_log_message(
      glue("iinla: Iteration {k} [max: {options$bru_max_iter}]"),
      verbosity = 1L
    )

    # Return previous result if inla crashes, e.g. when connection to server is
    # lost
    old.result <- result
    result <- NULL

    inla.formula <- update.formula(model$formula, BRU.response ~ .)
    inla.data <-
      c(
        stk.data,
        do.call(c, c(
          lapply(
            as_bru_comp_list(model),
            function(xx) as.list(bru_comp_env_extra(xx))
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
          weights = stk.data[["BRU.weights"]],
          scale = stk.data[["BRU.scale"]],
          offset = stk.data[["BRU.offset"]],
          control.predictor = list(link = stk.data[["BRU.link"]])
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
      inla.options.merged[["control.compute"]][["control.gcpo"]] <-
        control.gcpo.combined
    } else {
      # Compute the minimal amount required
      # Note: configs = TRUE only because latent indexing into mode$x is
      #   defined in configs$contents
      inla.options.merged <-
        modifyList(
          inla.options.merged,
          list(
            control.inla = list(
              control.vb = list(enable = FALSE),
              strategy = "gaussian",
              int.strategy = "eb"
            ),
            control.compute = list(
              config = TRUE,
              dic = FALSE,
              waic = FALSE
            ),
            # Required for line search weights:
            control.predictor = list(compute = TRUE)
          )
        )
      inla.options.merged[["control.compute"]][["control.gcpo"]] <- NULL
    }

    timings <- bru_timer_do(timings, "Run inla()", k)

    result <- fm_try_callstack(
      do.call(
        INLA::inla,
        inla.options.merged,
        envir = environment(as_bru_comp_list(model))
      )
    )

    timings <- bru_timer_do(timings, "Postprocess", k)

    if (inherits(result, "try-error")) {
      bru_log_warn(
        paste0("iinla: Problem in inla:\n", result)
      )
      bru_log_message(
        paste0(
          "iinla: Giving up and returning last successfully obtained result ",
          "for diagnostic purposes."
        ),
        verbosity = 1L,
        verbose = TRUE
      )
      if (is.null(old.result)) {
        old.result <- list()
      }
      if (!inherits(old.result, "iinla")) {
        class(old.result) <- c("iinla", class(old.result))
      }

      old.result[["bru_iinla"]] <- collect_misc_info()
      old.result$error <- result
      return(old.result)
    }

    # Extract values tracked for estimating convergence
    # Note: The number of fixed effects may be zero, and strong
    # non-linearities that don't necessarily affect the fixed
    # effects may appear in the random effects, so we need to
    # track all of them.
    result_mode <- evaluate_state(
      model,
      result,
      property = "joint_mode"
    )[[1]]
    result_sd <- evaluate_state(
      model,
      result,
      property = "sd",
      internal_hyperpar = TRUE
    )[[1]]
    track_df <- list()
    for (label in names(result_mode)) {
      track_df[[label]] <-
        data.frame(
          effect = label,
          index = seq_along(result_mode[[label]]),
          iteration = track_size + k,
          mode = result_mode[[label]],
          sd = result_sd[[label]]
        )
    }
    # label <- "bru_offset"
    # track_df[[label]] <-
    #   data.frame(
    #     effect = label,
    #     index = seq_along(stk.data[["BRU.offset"]]),
    #     iteration = track_size + k,
    #     mode = stk.data[["BRU.offset"]],
    #     sd = NA
    #   )
    label <- "log.posterior.mode"
    track_df[[label]] <-
      data.frame(
        effect = label,
        index = 1L,
        iteration = track_size + k,
        mode = result[["misc"]][["configs"]][["max.log.posterior"]],
        sd = Inf
      )
    track[[k]] <- do.call(rbind, track_df)

    if (is.null(result[["misc"]][["opt.trace"]])) {
      inla_track[[k]] <-
        tibble::tibble(
          iteration = integer(0),
          f = numeric(0),
          nfunc = integer(0)
        )
    } else {
      inla_track[[k]] <-
        tibble::tibble(
          iteration = inla_track_size + k,
          f = result[["misc"]][["opt.trace"]][["f"]],
          nfunc = result[["misc"]][["opt.trace"]][["nfunc"]],
          theta = result[["misc"]][["opt.trace"]][["theta"]]
        )
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
          timings <- bru_timer_do(timings, "Line search", k)
          line_weights <-
            extract_property(
              result = result,
              property = "predictor_sd",
              internal_hyperpar = FALSE
            )^{
              -2
            }
          line_search <- bru_line_search(
            model = model,
            lhoods = lhoods,
            lin = lin,
            state0 = state0,
            state = state,
            input = inputs,
            comp_lin = comp_lin,
            comp_simple = comp_simple,
            weights = line_weights,
            options = options
          )
          state <- line_search[["state"]]
        }
        timings <- bru_timer_do(timings, "Linearise", k)

        bru_log_message(
          "iinla: Evaluate component linearisations",
          verbosity = 3
        )
        comp_lin <- ibm_linear(
          model,
          input = inputs,
          state = state,
          inla_f = TRUE
        )
        bru_log_message(
          "iinla: Evaluate predictor linearisation",
          verbosity = 3
        )
        lin <- bru_compute_linearisation(
          model,
          lhoods = lhoods,
          input = inputs,
          state = state,
          comp_simple = comp_simple,
          options = options
        )

        stk <- bru_make_stack(lhoods, lin, idx)

        if (utils::packageVersion("INLA") <= "24.06.02") {
          stk.data <- INLA::inla.stack.data(stk)
        } else {
          stk.data <- INLA::inla.stack.data(stk,
            .response.name = "BRU.response"
          )
        }
        inla.options$control.predictor$A <- INLA::inla.stack.A(stk)

        N_response <- sum(bru_response_size(lhoods))
        N_predictor <- NROW(inla.options$control.predictor$A)
        if (N_response != N_predictor) {
          msg <- glue(
            "The total number of response values (N={N_response})",
            " and predictor values (N={N_predictor})",
            " do not match.\n",
            "  This is likely due to a mistake in the component or predictor",
            " constructions.",
            .trim = FALSE
          )
          if ((N_response > 1L) && (N_predictor == 1L)) {
            msg <- paste0(
              msg,
              "\nPerhaps you only have components with scalar inputs?"
            )
          }
          bru_log_message(msg, verbosity = 1L)
          stop(msg)
        }
      }
      # Store the state
      states <- c(states, list(state))
    }

    # Stopping criteria
    if (do_final_integration || (k >= options$bru_max_iter)) {
      interrupt <- TRUE
    } else if (!interrupt && (k > 1)) {
      max.dev <- options$bru_method$rel_tol
      dev <- abs(track[[k - 1]]$mode - track[[k]]$mode) / track[[k]]$sd
      ## do.call(c, lapply(by(do.call(rbind, track),
      ##   as.factor(do.call(rbind, track)$effect),
      ##   identity),
      ##   function(X) { abs(X$mean[k-1] - X$mean[k])/X$sd[k] }))
      bru_log_message(
        glue(
          "iinla: Max deviation from previous: ",
          "{signif(100 * max(dev), 3)}",
          "% of SD, and line search is ",
          if (line_search[["active"]]) "active" else "inactive",
          "\n       [stop if: < {100 * max.dev}% and line search inactive]"
        ),
        verbosity = 1L
      )
      do_final_integration <- all(dev < max.dev) && (!line_search[["active"]])
      if (do_final_integration) {
        do_final_theta_no_restart <- TRUE
        bru_log_message(
          "iinla: Convergence criterion met.",
          "\n       Running final INLA integration step with known theta mode.",
          verbosity = 1L
        )
      }
    }
    track[[k]][["new_linearisation"]] <- NA_real_
    for (label in names(states[[1]])) {
      track[[k]][["new_linearisation"]][track[[k]][["effect"]] == label] <-
        unlist(states[[length(states)]][[label]])
    }

    k <- k + 1L
  }

  bru_log_message("iinla: Computation completed", verbosity = 3)

  timings <- bru_timer_done(timings)

  result[["bru_iinla"]] <- collect_misc_info()
  class(result) <- c("iinla", class(result))
  result
}


auto_intercept <- function(components) {
  if (!inherits(components, "formula")) {
    if (inherits(components, "bru_comp_list")) {
      return(components)
    }
    stop("`components` must be a formula to auto-add an intercept component")
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
        as.formula(glue(
          ". ~ . - {var_names[inter_var]} -1"
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
  stopifnot(inherits(components, "bru_comp_list"))

  env <- environment(formula)
  formula <- update.formula(
    formula,
    glue(". ~ {glue_collapse(names(components), sep = ' + ')}")
  )
  environment(formula) <- env
  formula
}

# Extract the LHS of a formula, as response ~ .
extract_response <- function(formula) {
  if (inherits(formula, "formula") &&
    (length(as.character(formula)) == 3)) {
    as.formula(glue("{as.character(formula)[2]} ~ ."))
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


# Returns a formula's environment as a list. Removes all variable that are of
# type inla, function or formula. Also removes all variables that are not
# variables of the formula.
list.data <- function(formula) {
  # Formula environment as list
  elist <- as.list(environment(formula))

  # Remove previous inla results. For some reason these slow down the next INLA
  # call.
  elist <- elist[unlist(lapply(elist, function(x) !inherits(x, "inla")))]

  # Remove functions. This can cause problems as well.
  # elist <- elist[unlist(lapply(elist, function(x) !is.function(x)))]

  # Remove formulae. This can cause problems as well.
  # elist <- elist[unlist(lapply(elist, function(x) !inherits(x, "formula")))]

  # The formula expression is too general for this to be reliable:
  # Keep only purse formula variables
  # elist <- elist[names(elist) %in% all.vars(formula)]

  elist
}


#' Summary for an inlabru fit
#'
#' Takes a fitted `bru` object produced by [bru()] or [lgcp()] and creates
#' various summaries from it.
#'
#' @export
#' @method summary bru
#' @param object An object obtained from a [bru()] or [lgcp()] call
#' @param verbose logical; If `TRUE`, include more details of the
#' component definitions. If `FALSE`, only show basic component
#' definition information. Default: `FALSE`
#' @param \dots arguments passed on to component summary functions, see
#' [summary.bru_comp()].
#' @example inst/examples/bru.R
#'

summary.bru <- function(object, verbose = FALSE, ...) {
  object <- bru_check_object_bru(object)

  result <- list(
    bru_info = summary(object[["bru_info"]], verbose = verbose, ...)
  )

  if (inherits(object, "inla")) {
    result$inla <- NextMethod("summary", object)
    result[["inla"]][["call"]] <- NULL
  }

  class(result) <- "summary_bru"
  result
}


#' @export
#' @param x An object to be printed
#' @rdname summary.bru

print.summary_bru <- function(x, ...) {
  print(x$bru_info, ...)
  print(x$inla, ...)
  invisible(x)
}

#' @export
#' @param x A `bru` object to be printed
#' @describeIn bru Print a summary of a `bru` object.

print.bru <- function(x, ...) {
  print(summary(x, ...), ...)
  invisible(x)
}
