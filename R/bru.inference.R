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
#' @return The form of the value returned by `generate()` depends on the data
#' class and prediction formula. Normally, a data.frame is returned, or a list
#' of data.frames (if the prediction formula generates a list)
#' @example inst/examples/generate.bru.R

generate <- function(object, ...) {
  UseMethod("generate")
}

bru_check_object_bru <- function(object,
                                 new_version = getNamespaceVersion("inlabru")) {
  object <-
    bru_info_upgrade(
      object,
      new_version = new_version
    )
  object
}

bru_info_upgrade <- function(object,
                             new_version = getNamespaceVersion("inlabru")) {
  object_full <- object
  object <- object[["bru_info"]]
  msg <- NULL
  if (!is.list(object)) {
    msg <- "Not a list"
  } else if (!inherits(object, "bru_info")) {
    msg <- "Not a bru_info object"
  } else if (is.null(object[["inlabru_version"]])) {
    msg <- "`inlabru_version` is missing"
  }
  if (!is.null(msg)) {
    stop(paste0(
      "bru_info part of the object can't be converted to `bru_info`; ",
      msg
    ))
  }

  old_ver <- object[["inlabru_version"]]
  if (utils::compareVersion(new_version, old_ver) > 0) {
    warning(
      "Old bru_info object version ",
      old_ver,
      " detected. Attempting upgrade to version ",
      new_version
    )

    message("Detected bru_info version ", old_ver)
    if (utils::compareVersion("2.1.14.901", old_ver) > 0) {
      message("Upgrading bru_info to 2.1.14.901")
      # Ensure INLA_version stored
      if (is.null(object[["INLA_version"]])) {
        object[["INLA_version"]] <- "0.0.0"
      }
      # Check for component$mapper as a bm_multi
      for (k in seq_along(object[["model"]][["effects"]])) {
        cmp <- object[["model"]][["effects"]][[k]]
        cmp[["mapper"]] <-
          bm_multi(list(
            main = cmp$main$mapper,
            group = cmp$group$mapper,
            replicate = cmp$replicate$mapper
          ))
        object[["model"]][["effects"]][[k]] <- cmp
      }
      object[["inlabru_version"]] <- "2.1.14.901"
    }
    if (utils::compareVersion("2.5.3.9003", old_ver) > 0) {
      message("Upgrading bru_info to 2.5.3.9003")
      # Check that likelihoods store 'weights'
      for (k in seq_along(object[["lhoods"]])) {
        lhood <- object[["lhoods"]][[k]]
        if (is.null(lhood[["weights"]])) {
          lhood[["weights"]] <- 1
          object[["lhoods"]][[k]] <- lhood
        }
      }
      object[["inlabru_version"]] <- "2.5.3.9003"
    }

    if (utils::compareVersion("2.5.3.9005", old_ver) > 0) {
      message("Upgrading bru_info to 2.5.3.9005")
      # Make sure component$mapper is a bm_pipe
      for (k in seq_along(object[["model"]][["effects"]])) {
        cmp <- object[["model"]][["effects"]][[k]]
        # Convert offset mappers to const mappers
        if (inherits(cmp$main$mapper, c("bm_offset", "bru_mapper_offset"))) {
          cmp$main$mapper <- bm_const()
        }
        cmp[["mapper"]] <-
          bm_pipe(
            list(
              mapper = bm_multi(list(
                main = cmp$main$mapper,
                group = cmp$group$mapper,
                replicate = cmp$replicate$mapper
              )),
              scale = bm_scale()
            )
          )
        object[["model"]][["effects"]][[k]] <- cmp
      }
      object[["inlabru_version"]] <- "2.5.3.9005"
    }

    if (utils::compareVersion("2.6.0.9000", old_ver) > 0) {
      message("Upgrading bru_info to 2.6.0.9000")
      # Make sure component$mapper is a bm_pipe
      for (k in seq_along(object[["model"]][["effects"]])) {
        cmp <- object[["model"]][["effects"]][[k]]
        cmp[["mapper"]] <-
          bm_pipe(
            list(
              mapper = bm_multi(list(
                main = cmp$main$mapper,
                group = cmp$group$mapper,
                replicate = cmp$replicate$mapper
              )),
              scale = bm_scale()
            )
          )
        object[["model"]][["effects"]][[k]] <- cmp
      }
      object[["inlabru_version"]] <- "2.6.0.9000"
    }

    if (utils::compareVersion("2.7.0.9010", old_ver) > 0) {
      message("Upgrading bru_info to 2.7.0.9010")
      # Make sure component$group/replicate$input isn't NULL
      for (k in seq_along(object[["model"]][["effects"]])) {
        cmp <- object[["model"]][["effects"]][[k]]
        if (identical(deparse(cmp[["group"]][["input"]][["input"]]), "NULL")) {
          cmp[["group"]][["input"]][["input"]] <- expression(1L)
        }
        if (identical(
          deparse(cmp[["replicate"]][["input"]][["input"]]),
          "NULL"
        )) {
          cmp[["replicate"]][["input"]][["input"]] <- expression(1L)
        }
        object[["model"]][["effects"]][[k]] <- cmp
      }
      object[["inlabru_version"]] <- "2.7.0.9010"
    }

    if (utils::compareVersion("2.7.0.9016", old_ver) > 0) {
      message("Upgrading bru_info to 2.7.0.9016")
      # Convert old style include_component/allow_latent to intermediate
      # include_component/include_latent format

      # Update include/exclude information to limit it to existing components
      for (k in seq_along(object[["lhoods"]])) {
        if (isTRUE(object[["lhoods"]][[k]][["allow_latent"]])) {
          # Force inclusion of all components
          object[["lhoods"]][[k]][["include_latent"]] <- NULL
        } else {
          if (is.null(object[["lhoods"]][[k]][["include_latent"]])) {
            object[["lhoods"]][[k]][["include_latent"]] <- character(0)
          }
        }
      }
      object[["lhoods"]] <-
        bru_used_upgrade(
          object[["lhoods"]],
          labels = names(object[["model"]][["effects"]])
        )
      object[["inlabru_version"]] <- "2.7.0.9016"
    }

    if (utils::compareVersion("2.7.0.9017", old_ver) > 0) {
      message("Upgrading bru_info to 2.7.0.9017")
      # Convert old style include_component/allow_latent to new
      # bru_used format

      # Update include/exclude information to limit it to existing components
      for (k in seq_along(object[["lhoods"]])) {
        if (is.null(object[["lhoods"]][[k]][["used_components"]])) {
          if (isTRUE(object[["lhoods"]][[k]][["allow_latent"]])) {
            # Force inclusion of all components
            object[["lhoods"]][[k]][["include_latent"]] <- NULL
          } else {
            if (is.null(object[["lhoods"]][[k]][["include_latent"]])) {
              object[["lhoods"]][[k]][["include_latent"]] <- character(0)
            }
          }
        }
      }
      object[["lhoods"]] <-
        bru_used_upgrade(
          object[["lhoods"]],
          labels = names(object[["model"]][["effects"]])
        )
      object[["inlabru_version"]] <- "2.7.0.9017"
    }

    if (utils::compareVersion("2.7.0.9021", old_ver) > 0) {
      message("Upgrading bru_info to 2.7.0.9021")
      # Make sure 'used' components format is properly stored

      if (is.null(object[["lhoods"]][["is_additive"]])) {
        object[["lhoods"]][["is_additive"]] <-
          object[["lhoods"]][["linear"]]
      }

      object[["lhoods"]] <-
        bru_used_upgrade(
          object[["lhoods"]],
          labels = names(object[["model"]][["effects"]])
        )
      object[["inlabru_version"]] <- "2.7.0.9021"
    }

    if (utils::compareVersion("2.10.1.9007", old_ver) > 0) {
      message("Upgrading bru_info to 2.10.1.9007")
      # Update timings info to difftime format

      warning(
        paste0(
          "From 2.10.1.9007, elapsed time is in Elapsed, Time is CPU time.",
          "\n  Copying old elapsed time to both Elapsed and Time, ",
          "setting System to zero;",
          "  Do not over-interpret."
        ),
        immediate. = TRUE
      )
      conversion <- as.difftime(
        object_full[["bru_timings"]][["Time"]],
        units = "secs"
      )
      object_full[["bru_timings"]][["Time"]] <- conversion
      object_full[["bru_timings"]][["System"]] <-
        as.difftime(rep(0.0, length(conversion)), units = "secs")
      object_full[["bru_timings"]][["Elapsed"]] <- conversion

      object[["inlabru_version"]] <- "2.10.1.9007"
    }

    if (utils::compareVersion("2.10.1.9012", old_ver) > 0) {
      message("Upgrading bru_info to 2.10.1.9012")
      # Update log format to new format

      object_full[["bru_iinla"]][["log"]] <-
        bru_log_new(
          object_full[["bru_iinla"]][["log"]][["log"]],
          bookmarks = object_full[["bru_iinla"]][["log"]][["bookmarks"]]
        )

      object[["inlabru_version"]] <- "2.10.1.9012"
    }

    if (utils::compareVersion("2.12.0.9014", old_ver) > 0) {
      message("Upgrading bru_info to 2.12.0.9014")
      # Update is_additive/linear storage

      object[["lhoods"]][["is_additive"]] <-
        object[["lhoods"]][["linear"]]

      # Update predictor expression
      object[["lhoods"]] <-
        bru_used_upgrade(
          object[["lhoods"]],
          labels = names(object[["model"]][["effects"]])
        )

      object[["inlabru_version"]] <- "2.12.0.9014"
    }

    if (utils::compareVersion("2.12.0.9017", old_ver) > 0) {
      message("Upgrading bru_info to 2.12.0.9017")

      # Update bru_like class names to bru_obs
      if (!is.null(object[["lhoods"]])) {
        object[["lhoods"]] <-
          lapply(object[["lhoods"]], function(x) {
            class(x) <- "bru_obs"
            x
          })
        class(object[["lhoods"]]) <- c("bru_obs_list", "list")
      }

      eff <- object[["model"]][["effects"]]
      for (k in seq_along(object[["model"]][["effects"]])) {
        class(eff[[k]][["main"]]) <- "bru_subcomp"
        class(eff[[k]][["group"]]) <- "bru_subcomp"
        class(eff[[k]][["replicate"]]) <- "bru_subcomp"
        class(eff[[k]]) <- "bru_comp"
      }
      class(eff) <- c("bru_comp_list", "list")
      eff <- object[["model"]][["effects"]] <- eff

      object[["inlabru_version"]] <- "2.12.0.9017"
    }

    object[["inlabru_version"]] <- new_version
    message(paste0("Upgraded bru_info to ", new_version))

    object_full[["bru_info"]] <- object
  }
  object_full
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
#'   result object. The default print method show information about model
#'   components and observation models.
bru_info.bru <- function(object, ...) {
  object <- bru_check_object_bru(object)
  object[["bru_info"]]
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
      lhoods = summary(object[["lhoods"]], verbose = verbose, ...)
    ),
    class = "summary_bru_info"
  )
}

#' @export
#' @param x An  object to be printed
#' @rdname bru_info
print.summary_bru_info <- function(x, ...) {
  cat(paste0("inlabru version: ", x$inlabru_version, "\n"))
  cat(paste0("INLA version: ", x$INLA_version, "\n"))
  cat(paste0("Components:\n"))
  print(x$components)
  cat(paste0("Observation models:\n"))
  print(x$lhoods)
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
#' @param ... unused
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
        paste0("'", names(lhoods)[!(dot_is_lhood | dot_is_lhood_list)], "'",
          collapse = ", "
        ),
        " were meant to be given to a call to 'bru_obs()',\n  or in the ",
        "'options' list argument instead."
      ))
    }
    if (!all(dot_is_lhood) && !all(dot_is_lhood_list)) {
      stop(paste0(
        "Cannot mix `bru_obs` and `bru_obs_list` objects in the `...` ",
        "argument of `bru()`."
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
#'       by looking at the [lgcp] function's documentation.
#'     \item
#'       Constructing multiple likelihoods is straight forward. See [like] for
#'       more information on how to provide additional likelihoods to `bru`
#'       using the `...` parameter list.
#'     \item
#'       Support for non-linear predictors. See example below.
#'     \item
#'       Log Gaussian Cox process (LGCP) inference is
#'       available by using the `cp` family or (even easier) by using the
#'       [lgcp] function.
#'   }
#' @export
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'
#' @param components A `formula`-like specification of latent components.
#'   Also used to define a default linear additive predictor.  See
#'   [bru_comp()] for details.
#' @param ... Obervation models, each constructed by a calling [bru_obs()], or
#'   named parameters that can be passed to a single [bru_obs()] call. Note that
#'   all the arguments will be evaluated before calling [bru_obs()] in order to
#'   detect if they are `like` objects. This means that special arguments that
#'   need to be evaluated in the context of `response_data` or `data` (such as
#'   `Ntrials`) may will only work that way in direct calls to [bru_obs()].
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

  timings_convert <- function(x) {
    if (!is.na(x[4])) {
      x[1] <- x[1] + x[4]
    }
    if (!is.na(x[5])) {
      x[2] <- x[2] + x[5]
    }
    x[1:3]
  }

  timings_collect <- function(Task, Iteration, time_diff) {
    data.frame(
      Task = Task,
      Iteration = Iteration,
      Time = as.difftime(unname(time_diff[1]), units = "secs"),
      System = as.difftime(unname(time_diff[2]), units = "secs"),
      Elapsed = as.difftime(unname(time_diff[3]), units = "secs")
    )
  }

  timing <- list(start = timings_convert(proc.time()))

  # Turn model components and bru_obs objects into internal bru model
  bru.model <- bru_model(
    components,
    lhoods = list(...),
    options = options,
    .envir = .envir
  )

  # Set max iterations to 1 if all likelihood formulae are linear
  if (all(vapply(
    bru.model[["lhoods"]],
    function(lh) isTRUE(lh[["linear"]]), TRUE
  ))) {
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

  timing$setup <- timings_convert(proc.time())

  # Run iterated INLA
  if (options$bru_run) {
    result <- iinla(
      model = info[["model"]],
      lhoods = info[["lhoods"]],
      inputs = info[["inputs"]],
      options = info[["options"]]
    )
  } else {
    result <- list()
  }

  timing$end <- timings_convert(proc.time())
  result$bru_timings <-
    rbind(
      timings_collect("Preprocess", 0L, timing$setup - timing$start),
      result[["bru_iinla"]][["timings"]]
    )

  # Add bru information to the result
  result$bru_info <- info
  class(result) <- c("bru", class(result))
  return(result)
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
    lhoods = info[["lhoods"]],
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
  return(result)
}

#' @title Set missing values in observation models
#'
#' @description Set all or parts of the observation model response data
#' to `NA`, for example for use in cross validation (with [bru_rerun()])
#' or prior sampling (with [bru_rerun()] and [inlabru::generate()]).
#'
#' @param object A `bru`, `bru_obs` or `bru_obs_list` object
#' @param keep For `bru_obs`, a single logical or an integer vector;
#'   If `TRUE`, keep all the response data, if `FALSE` (default),
#'   set all of it to `NA`. An integer vector determines which elements
#'   to keep (for positive values) or to set as missing (negative values).
#'
#'   For `bru` and `bru_obs_list`, a logical scalar or vector, or a list, see
#'   Details.
#' @param \dots Additional arguments passed on to the `bru_obs` method.
#'   Currently unused.
#'
#' @details For `bru` and `bru_obs_list`,
#' \itemize{
#'   \item{`keep` must be either a single logical, which is expanded to a list,}
#'   \item{a logical vector, which is converted to a list,}
#'   \item{an unnamed list of the same length as the number of observation
#'     models, with elements compatible with the `bru_obs` method, or}
#'   \item{a named list with elements compatible with the `bru_obs` method,
#'     and only the named `bro_obs` models are acted upon, i.e. the elements
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
#'     x[["response_data"]][["BRU_response"]]
#'   }
#' )
#' lapply(
#'   bru_set_missing(obs, keep = list(B = FALSE)),
#'   function(x) {
#'     x[["response_data"]][["BRU_response"]]
#'   }
#' )
#' lapply(
#'   bru_set_missing(obs, keep = list(1:4, -(3:5))),
#'   function(x) {
#'     x[["response_data"]][["BRU_response"]]
#'   }
#' )
bru_set_missing <- function(object, keep = FALSE, ...) {
  UseMethod("bru_set_missing")
}
#' @export
#' @rdname bru_set_missing
bru_set_missing.bru <- function(object, keep = FALSE, ...) {
  object <- bru_check_object_bru(object)
  object[["bru_info"]][["lhoods"]] <-
    bru_set_missing(object[["bru_info"]][["lhoods"]], keep = keep, ...)
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
        paste0("'", setdiff(names(keep), names(object)), "'", collapse = ", ")
      ))
    }
  }
  for (idx in idxs) {
    object[[idx]] <- bru_set_missing(object[[idx]],
      keep = keep[[idx]],
      ...
    )
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
  if (inherits(object[["response_data"]][["BRU_response"]], "inla.surv")) {
    if (!is.data.frame(object[["response_data"]][["BRU_response"]])) {
      # TODO: This block should be function, but inla.surv might standardise
      # to tibble or data.frame, so we should remove this when possible.
      dat <- object[["response_data"]][["BRU_response"]]
      cls <- class(dat)
      att <- attributes(dat)
      cure_null <- ("cure" %in% names(dat)) && is.null(dat[["cure"]])
      if (cure_null) {
        dat["cure"] <- NULL
      }
      dat <- as.data.frame(unclass(dat))
      attr(dat, "names.ori") <- att[["names.ori"]]
      class(dat) <- c("inla.surv", "data.frame")
      object[["response_data"]][["BRU_response"]] <- dat
    }
    # Note: by only setting time,lower,upper to NA, printing of the object
    # still works without giving an NA indexing error.
    object[["response_data"]][["BRU_response"]]$time[-keep] <- NA
    object[["response_data"]][["BRU_response"]]$lower[-keep] <- NA
    object[["response_data"]][["BRU_response"]]$upper[-keep] <- NA
  } else if (is.data.frame(object[["response_data"]][["BRU_response"]])) {
    # Handles data.frame and tibble, including inla.mdata
    object[["response_data"]][["BRU_response"]][-keep, ] <- NA
  } else {
    object[["response_data"]][["BRU_response"]][-keep] <- NA
  }
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

#' @title Evaluate expressions in data contexts
#' @description Evaluate an expression in a series of data contexts, also making
#'   the objects directly available as names surrounded by ".", stopping when
#'   the expression evaluation completes with no error.
#'
#'   This is an internal inlabru method, not intended for general use.
#' @param input An expression to be evaluated
#' @param data list of data objects in priority order. Named elements will
#' be available as `.name.` in the evaluation. The `input` expression is
#' evaluated with each non-NULL `data` object as `envir`, in order,
#' until success. If there are no non-NULL data objects, the expression is
#' evaluated in an empty environment, potentially falling back to enclosing
#' environment variables.
#' @param default Value used if the expression is evaluated as NULL. Default
#' NULL
#' @param .envir The evaluation environment
#' @return The result of expression evaluation
#' @keywords internal
#' @examples
#' # The A values come from the 'data' element, and the B values come from
#' # the 'response_data' element, as that is listed first.
#' bru_eval_in_data_context(
#'   quote(
#'     list(A = .data.$x, B = x)
#'   ),
#'   list(
#'     response_data = tibble::tibble(x = 1:5),
#'     data = tibble::tibble(x = 1:10)
#'   )
#' )
#' # Both A and B come from the 'data' element, as 'x' is found there,
#' # terminating the evaluation attempt.
#' bru_eval_in_data_context(
#'   quote(
#'     list(A = .data.$x, B = x)
#'   ),
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
  data_orig <- data
  data <- lapply(data, function(x) {
    if (!is.null(x) && !is.list(x)) {
      x <- tibble::as_tibble(x)
    }
    x
  })
  enclos_envir <- new.env(parent = .envir)
  nms <- names(data)
  for (nm in setdiff(nms, "")) {
    assign(paste0(".", nm, "."), data_orig[[nm]], envir = enclos_envir)
  }
  success <- FALSE
  result <- NULL
  for (k in seq_along(data)) {
    if (is.null(data[[k]])) {
      next
    }
    result <- try(
      eval(input, envir = data[[k]], enclos = enclos_envir),
      silent = TRUE
    )
    if (!inherits(result, "try-error")) {
      success <- TRUE
      break
    }
  }
  if (all(vapply(data, is.null, TRUE))) {
    result <- try(
      eval(input, envir = NULL, enclos = enclos_envir),
      silent = TRUE
    )
    if (!inherits(result, "try-error")) {
      success <- TRUE
    }
  }
  if (!success) {
    stop(paste0(
      "Input '",
      paste0(deparse(input), collapse = "\n"),
      "' could not be evaluated."
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

  dummies <- seq_len(length(new_coordnames))
  dummies <- setdiff(dummies, c(from_data, from_ips))

  new_coordnames[dummies] <-
    paste0(
      "BRU_dummy_coordinate_",
      seq_along(new_coordnames)
    )[dummies]

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
      stop(paste0(
        "Column '", nm, "' is not a simple feature column in all data objects."
      ))
    }
    the_crs <- fm_crs(dt[[sf_data_idx_[1]]][[nm]])
    for (i in sf_data_idx_) {
      if (!inherits(dt[[i]][[nm]], "sfc")) {
        stop(paste0(
          "Column '",
          nm,
          "' is not a simple feature column in all data objects."
        ))
      }
      dt_crs <- fm_crs(dt[[i]][[nm]])
      if (!fm_crs_is_identical(dt_crs, the_crs)) {
        dt[[i]][[nm]] <- fm_transform(dt[[i]][[nm]], crs = the_crs)
      }
    }

    ncol_ <- vapply(
      sf_data_idx_,
      function(i) {
        crds <- sf::st_coordinates(dt[[i]][[nm]])
        length(intersect(colnames(crds), c("X", "Y", "Z")))
      },
      0L
    )
    if (length(unique(ncol_)) > 0) {
      # Some dimension mismatch
      ncol_max <- max(ncol_)
      for (ii in which(ncol_ < ncol_max)) {
        i <- sf_data_idx_[ii]
        # Extend columns
        if (nrow(dt[[i]]) > 0) {
          dt[[i]][[nm]] <-
            sf::st_zm(dt[[i]][[nm]], drop = FALSE, what = "Z")
        } else {
          dt[[i]][[nm]] <-
            sf::st_as_sf(
              as.data.frame(
                matrix(0.0, 0L, ncol_max)
              ),
              coords = seq_len(ncol_max),
              crs = fm_crs(dt[[i]][[nm]])
            )$geometry
        }
      }
    }
  }
  # Technically, we only need this if at lest one of the objects is a tibble
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

#' @title Check for predictor expression additivity
#' @description Checks if a predictor expression is additive or not
#' @param x A predictor `expression`, `formula`, or parse information
#'   `data.frame`.
#' @param \dots Arguments passed on recursively.
#' @param verbose logical; if `TRUE`, print diagnostic parsing information.
#' @return `TRUE` if the expression is detected to be additive, `FALSE`
#'   otherwise.
#' @keywords internal
#' @export
bru_is_additive <- function(x, ...) {
  UseMethod("bru_is_additive")
}

#' @rdname bru_is_additive
#' @export
bru_is_additive.data.frame <- function(x, root_id = 0, ..., verbose = FALSE) {
  if (root_id == 0) {
    root_id <- x$id[x$parent == 0]
    if (length(root_id) != 1L) {
      stop("Cannot determine parser root id")
    }
  }

  x_root <- x[x$id == root_id, , drop = FALSE]
  if (x_root$token == "SYMBOL") {
    if (verbose) {
      message("SYMBOL found")
    }
    return(TRUE)
  }
  if (x_root$token == "expr") {
    if (verbose) {
      message("expr found, may be additive")
    }
    x_terms <- x[x$parent == root_id, , drop = FALSE]
    if (nrow(x_terms) < 1L) {
      if (verbose) {
        message("No terms found, assuming non-additive")
      }
      return(FALSE)
    }
    if (nrow(x_terms) == 1L) {
      if (verbose) {
        message("Only one term, checking for symbol")
      }
      add <- x_terms$token[1] == "SYMBOL"
      if (add && verbose) {
        message("SYMBOL found")
      }
      return(add)
    }
    if (nrow(x_terms) == 2L) {
      if (x_terms$token[1] != "'+'") {
        if (verbose) {
          message("Only two terms and no '+' found, assuming non-additive")
        }
        return(FALSE)
      }
      if (verbose) {
        message("'+' found, may be additive")
      }
      add <- bru_is_additive(
        x = x,
        root_id = x_terms$id[2],
        ...,
        verbose = verbose
      )
      return(all(add))
    }
    if (nrow(x_terms) >= 4L) {
      if (verbose) {
        message("More than 3 terms, assuming non-additive")
      }
      return(FALSE)
    }
    if (x_terms$token[2] == "'+'") {
      if (verbose) {
        message("'+' found, may be additive")
      }
      add <- c(
        bru_is_additive(x = x, root_id = x_terms$id[1], verbose = verbose),
        bru_is_additive(x = x, root_id = x_terms$id[3], verbose = verbose)
      )
      return(all(add))
    }
    if ((x_terms$token[1] == "'('") && (x_terms$token[3] == "')'")) {
      if (x_terms$token[2] == "expr") {
        if (verbose) {
          message("(expr) found, may be additive")
        }
        add <- bru_is_additive(
          x = x,
          root_id = x_terms$id[2],
          verbose = verbose
        )
        return(add)
      }
    }
  }

  if (verbose) {
    message("No known additive structure found, assuming non-additive")
  }

  return(FALSE)
}

#' @rdname bru_is_additive
#' @export
bru_is_additive.character <- function(x, ...) {
  bru_is_additive(
    bru_get_parse_data(x),
    root_id = 0,
    ...
  )
}

#' @rdname bru_is_additive
#' @export
bru_is_additive.expression <- function(x, ...) {
  bru_is_additive(
    bru_get_parse_data(as.character(x)),
    root_id = 0,
    ...
  )
}

#' @rdname bru_is_additive
#' @export
bru_is_additive.formula <- function(x, ...) {
  bru_is_additive(
    bru_get_parse_data(as.character(x)[length(x)]),
    root_id = 0,
    ...
  )
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
#'   to `fmesher::fm_int(domain, samplers)`. If explicitly given,
#'   overrides `domain` and `samplers`.}
#' }
#' @param used Either `NULL` (default) or a [bru_used()] object. When,
#'   `NULL`, the information about what effects and
#' latent vectors are made available to the predictor evaluation is defined by
#' `bru_used(formula)`, which will include all effects and latent vectors
#' used by the predictor expression.
#' @param allow_combine logical; If `TRUE`, the predictor expression may involve
#'   several rows of the input data to influence the same row. When `NULL`,
#'   defaults to `FALSE`, unless `response_data` is non-`NULL`, or `data` is a
#'   `list`, or the likelihood construction requires it.
#' @param aggregate character ("none", "sum", "average", "logsumexp", or
#'   "logaverageexp", as defined by `bm_aggregate(type = aggregate)`)
#'   or an aggregation `bru_mapper` object
#'   ([bm_aggregate()] or [bm_logsumexp()]). Default `NULL`,
#'   interpreted as "none". `r lifecycle::badge("experimental")`, available
#'   from version `2.12.0.9013`.
#' @param aggregate_input `NULL` or an optional input list to the mapper
#'   defined by non-NULL `aggregate`, overriding the default,
#'   ```
#'   list(block = .data.[[".block"]],
#'        weights = .data.[["weight"]],
#'        n_block = bru_response_size(.response_data.))
#'   ```
#'   `r lifecycle::badge("experimental")`, available from version `2.12.0.9013`.
#' @param control.family A optional `list` of `INLA::control.family` options
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
                    allow_combine = NULL,
                    aggregate = NULL,
                    aggregate_input = NULL,
                    control.family = NULL,
                    tag = NULL,
                    options = list(),
                    .envir = parent.frame(),
                    include = deprecated(),
                    exclude = deprecated(),
                    include_latent = deprecated()) {
  options <- bru_call_options(options)
  bru_options_set_local(options, .reset = TRUE)

  # Some defaults
  inla.family <- family

  formula_char <- as.character(formula)

  # Does the likelihood formula imply an additive predictor?
  is_additive <- bru_is_additive(formula_char[length(formula_char)])
  is_additive_dot <- formula_char[length(formula_char)] == "."

  # If not additive, set predictor expression according to the formula's RHS
  if (!is_additive) {
    expr <- parse(text = formula_char[length(formula_char)])
  } else {
    expr <- NULL
  }

  # Set the response name
  if (length(formula_char) < 3) {
    stop("Missing response variable names")
  }
  response_expr <- parse(text = formula_char[2])
  response <- tryCatch(
    expr = bru_eval_in_data_context(
      substitute(response_expr),
      data = list(response_data = response_data, data = data),
      default = NULL,
      .envir = .envir
    ),
    error = function(e) {
      NULL
    }
  )

  # Catch and handle special cases:
  if ((family == "cp") && (is.null(response) || !inherits(response, "list"))) {
    domain_names <- trimws(strsplit(formula_char[2], split = "\\+")[[1]])
    if (!is.null(domain_names)) {
      # "a + b" conversion to list(a = a, b = b)
      domain_expr <- paste0(
        "list(",
        paste0(
          vapply(
            domain_names, function(x) {
              if (identical(x, "coordinates")) {
                paste0(x, " = sp::", x, "(.data.)")
              } else {
                paste0(x, " = ", x)
              }
            },
            ""
          ),
          collapse = ", "
        ),
        ")"
      )
      response_expr <- parse(text = domain_expr)
    }
    response <- tryCatch(
      expr = bru_eval_in_data_context(
        substitute(response_expr),
        data = list(response_data = response_data, data = data),
        default = NULL,
        .envir = .envir
      ),
      error = function(e) {
        NULL
      }
    )
  }

  if (is.null(response)) {
    stop("Response variable missing or could not be evaluated")
  }

  E <- bru_eval_in_data_context(
    substitute(E),
    data = list(response_data = response_data, data = data),
    default = options[["E"]],
    .envir = .envir
  )
  Ntrials <- bru_eval_in_data_context(
    substitute(Ntrials),
    data = list(response_data = response_data, data = data),
    default = options[["Ntrials"]],
    .envir = .envir
  )
  weights <- bru_eval_in_data_context(
    substitute(weights),
    data = list(response_data = response_data, data = data),
    default = 1,
    .envir = .envir
  )
  scale <- bru_eval_in_data_context(
    substitute(scale),
    data = list(response_data = response_data, data = data),
    default = 1,
    .envir = .envir
  )
  if (!is.null(aggregate) && is.character(aggregate)) {
    aggregate <- match.arg(
      aggregate,
      c("none", "sum", "average", "logsumexp", "logaverageexp")
    )
    aggregate <- switch(aggregate,
      "none" = NULL,
      bm_aggregate(type = aggregate)
    )
  }
  if (!is.null(aggregate)) {
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
            stop(paste0(
              "`data` has ", NROW(data), " rows, but `ips` has ",
              max(ips$.block), " blocks. Cannot join."
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

    aggregate_input <- bru_eval_in_data_context(
      substitute(aggregate_input),
      data = list(data = data, response_data = response_data),
      default = NULL,
      .envir = .envir
    )
    if (is.null(aggregate_input)) {
      aggregate_input <- list()
    }
    if (is.null(aggregate_input[["block"]])) {
      aggregate_input[["block"]] <- bru_eval_in_data_context(
        quote(.data.[[".block"]]),
        data = list(data = data, response_data = response_data),
        default = NULL,
        .envir = .envir
      )
    }
    if (is.character(aggregate_input[["block"]])) {
      msg <- paste0(
        "'character' aggregation block information detected.\n",
        "Please use a numeric or integer vector for the block information ",
        "instead.\n",
        "If you want to use a character vector, ",
        "please convert it to a factor first,\n",
        "making sure the factor level order matches your intended order,\n",
        "and use `as.integer()`."
      )

      if (utils::packageVersion("fmesher") < "0.5.0") {
        msg <- c(msg, paste0(
          "You have fmesher < 0.5.0. From version 0.5.0,",
          "`fm_int()`/`fm_cprod()`\ncreates integer block information ",
          "automatically."
        ))
      }
      bru_log_abort(msg)
    }
    if (is.null(aggregate_input[["weights"]])) {
      aggregate_input[["weights"]] <- bru_eval_in_data_context(
        quote(.data.[["weight"]]),
        data = list(data = data, response_data = response_data),
        default = NULL,
        .envir = .envir
      )
    }
    if (is.null(aggregate_input[["n_block"]])) {
      aggregate_input[["n_block"]] <- bru_eval_in_data_context(
        quote(bru_response_size(.response_data.)),
        data = list(data = data, response_data = response_data),
        default = NULL,
        .envir = .envir
      )
    }
    if (is.null(aggregate_input[["block"]])) {
      stop(paste0(
        "Aggregation requested, but `aggregate_input[['block']]` ",
        "evaluates to NULL."
      ))
    }
    if (is.null(aggregate_input[["weights"]])) {
      stop(paste0(
        "Aggregation requested, but `aggregate_input[['weights']]` ",
        "evaluates to NULL."
      ))
    }
    if (is.null(aggregate_input[["n_block"]])) {
      stop(paste0(
        "Aggregation requested, but `aggregate_input[['n_block']]` ",
        "evaluates to NULL."
      ))
    }
  }

  data_extra <- as.list(data_extra)
  if (!is.null(aggregate)) {
    if (!is_additive) {
      expr_text <- formula_char[length(formula_char)]
    } else {
      expr_text <- "BRU_EXPRESSION"
    }
    expr_text <- paste0(
      "{ibm_eval(BRU_aggregate_mapper, input = BRU_aggregate_input,",
      " state = {", expr_text, "})}"
    )
    expr <- parse(text = expr_text)
    data_extra[["BRU_aggregate_mapper"]] <- aggregate
    data_extra[["BRU_aggregate_input"]] <- aggregate_input
    allow_combine <- TRUE
    is_additive <- FALSE
  }

  if (inherits(data, "Spatial")) {
    lifecycle::deprecate_warn(
      "2.12.0.9023",
      "bru_obs(data = 'has deprecated support for `Spatial` input')",
      I("`sf` input")
    )
  }
  if (inherits(response, "Spatial")) {
    lifecycle::deprecate_warn(
      "2.12.0.9023",
      I("`bru_obs() response objects of `Spatial` type"),
      I("`sf` input")
    )
  } else if (is.list(response) &&
    inherits(response[["coordinates"]], "Spatial")) {
    lifecycle::deprecate_warn(
      "2.12.0.9023",
      I("`bru_obs() response objects of `Spatial` type"),
      I("`sf` input")
    )
  }
  if (inherits(samplers, "Spatial")) {
    lifecycle::deprecate_warn(
      "2.12.0.9023",
      "bru_obs(samplers = 'has deprecated support for `Spatial` input')",
      I("`sf` input")
    )
  }
  if (inherits(ips, "Spatial")) {
    lifecycle::deprecate_warn(
      "2.12.0.9023",
      "bru_obs(ips = 'has deprecated support for `Spatial` input')",
      I("`sf` input")
    )
  }

  # More on special bru likelihoods
  if (family == "cp") {
    if (!is.null(aggregate)) {
      stop("The 'aggregate' feature cannot be used with family='cp'.")
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
        stop(paste0(
          "Mismatch between names of observed dimensions and the given ",
          "domain names:\n",
          "    names(response) = (",
          paste0(names(response), collapse = ", "),
          ")\n",
          "    names(domain)   = (",
          paste0(names(domain), collapse = ", "),
          ")"
        ))
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

    if (length(E) > 1) {
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
    if (!is.null(response_data)) {
      data_ <- response_data
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
        if (any(!idx)) {
          data <- dplyr::bind_cols(data, as.data.frame(response[!idx]))
        }
      } else {
        data <- as.data.frame(response)
      }
      response_data <- NULL
      N_data <- NROW(data)
    } else {
      data <- tibble::as_tibble(response)
      if (("geometry" %in% names(data)) &&
        inherits(data$geometry, "sfc")) {
        sf::st_geometry(data) <- "geometry"
      }
      response_data <- NULL
      N_data <- NROW(data)
    }

    # Add back additional data
    additional_data_names <- setdiff(names(data_), names(data))
    if ((length(additional_data_names) > 0) &&
      (NROW(data_) == N_data)) {
      data <- cbind(data, tibble::as_tibble(data_)[additional_data_names])
    }

    if (ips_is_Spatial) {
      ips <- as.data.frame(ips)
    } else {
      if ("geometry" %in% names(ips)) {
        sf::st_geometry(ips) <- "geometry"
      }
    }

    # Use 'weights' for per-point weighting of eta
    if (length(weights) == 1L) {
      point_weights <- rep(weights, N_data)
    } else {
      stopifnot(length(weights) == N_data)
      point_weights <- weights
    }
    weights <- 1L

    if (identical(options[["bru_compress_cp"]], TRUE)) {
      allow_combine <- TRUE
      response_data <- tibble::tibble(
        BRU_E = c(
          0,
          E * ips[["weight"]]
        ),
        BRU_response_cp = c(
          N_data,
          rep(0, NROW(ips))
        )
      )
      if (!is_additive) {
        expr_text <- formula_char[length(formula_char)]
      } else {
        expr_text <- "BRU_EXPRESSION"
      }
      expr_text <- paste0(
        "{\n",
        "  BRU_eta <- {", expr_text, "}\n",
        "  if (length(BRU_eta) == 1L) {\n",
        "    BRU_eta <- rep(BRU_eta, length(BRU_aggregate))\n",
        "  }\n",
        "  c(mean(BRU_point_weights[BRU_aggregate] *\n",
        "         BRU_eta[BRU_aggregate]),\n",
        "    BRU_eta[!BRU_aggregate])\n",
        "}"
      )
      expr <- parse(text = expr_text)

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
      if (!all(point_weights == 1)) {
        stop(
          "Point 'weights' are not supported for non-compressed Cox processes."
        )
      }
      response_data <- data.frame(
        BRU_E = c(
          rep(0, N_data),
          E * ips[["weight"]]
        ),
        BRU_response_cp = c(
          rep(1, N_data),
          rep(0, NROW(ips))
        )
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

    response <- "BRU_response_cp"
    inla.family <- "poisson"
    E <- response_data[["BRU_E"]]

    allow_combine <- TRUE
  } else {
    if (!is.logical(allow_combine)) {
      if (!is.null(response_data)) {
        bru_log_warn(
          paste0(
            "Non-null response data supplied; ",
            "guessing allow_combine=TRUE.",
            "\n  Specify allow_combine explicitly to avoid this warning."
          )
        )
        allow_combine <- TRUE
      } else if (is.list(data) && !is.data.frame(data)) {
        bru_log_warn(
          paste0(
            "Non data-frame list-like data supplied; ",
            "guessing allow_combine=TRUE.",
            "\n  Specify allow_combine explicitly to avoid this warning."
          )
        )
        allow_combine <- TRUE
      } else {
        allow_combine <- FALSE
      }
    }

    # Need to make a list instead of data.frame, to allow inla.mdata responses
    response_data <- list(
      BRU_response = response,
      BRU_E = E,
      BRU_Ntrials = Ntrials,
      BRU_scale = scale
    )
    response <- "BRU_response"
  }

  if (lifecycle::is_present(include) ||
    lifecycle::is_present(exclude) ||
    lifecycle::is_present(include_latent)) {
    if (!lifecycle::is_present(include)) {
      include <- NULL
    } else {
      bru_log_message(
        paste0(
          "The `include` argument of `bru_obs()` is deprecated ",
          "since inlabru 2.11.0 and may be ignored.\n\t",
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
          "The provided `include` value may be ignored.",
          "If auto-detection doesn't work, use `bru_used(effect = include)`"
        )
      )
    }
    if (!lifecycle::is_present(exclude)) {
      exclude <- NULL
    } else {
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
          "The provided `exclude` value may be ignored.",
          paste0(
            "If auto-detection doesn't work, ",
            "use `bru_used(effect_exclude = exclude)`"
          )
        )
      )
    }
    if (!lifecycle::is_present(include_latent)) {
      include_latent <- NULL
    } else {
      bru_log_message(
        paste0(
          "The `include_latent` argument of `bru_obs()` is deprecated ",
          "since inlabru 2.11.0 and may be ignored.\n\t",
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
          "The provided `include_latent` value may be ignored.",
          paste0(
            "If auto-detection doesn't work, ",
            "use `bru_used(latent = include_latent)`"
          )
        )
      )
    }
    if (is.null(used)) {
      used <- bru_used(formula,
        effect = include,
        effect_exclude = exclude,
        latent = include_latent
      )
    }
  }
  if (is.null(used)) {
    used <- bru_used(formula)
  }

  # The likelihood object that will be returned

  lh <- structure(
    list(
      family = family,
      formula = formula,
      response_data = response_data, # agg
      data = data,
      data_extra = data_extra,
      E = E,
      Ntrials = Ntrials,
      weights = weights,
      scale = scale,
      samplers = samplers,
      is_additive = is_additive,
      linear = is_additive, # Not quite correct, as depends on component defs
      expr = expr,
      response = response,
      inla.family = inla.family,
      domain = domain,
      used = used,
      allow_combine = allow_combine,
      control.family = control.family,
      tag = tag
    ),
    class = "bru_obs"
  )

  # Return likelihood
  lh
}

#' @describeIn bru_obs `r lifecycle::badge("deprecated")` Legacy `like()`
#' method for `inlabru` prior to version `2.12.0`. Use [bru_obs()] instead.
#' @param mesh `r lifecycle::badge("deprecated")` Ignored.
#' @export
like <- function(formula = . ~ .,
                 family = "gaussian",
                 data = NULL,
                 response_data = NULL,
                 E = NULL,
                 Ntrials = NULL,
                 weights = NULL,
                 scale = NULL,
                 domain = NULL,
                 samplers = NULL,
                 ips = NULL,
                 used = NULL,
                 allow_combine = NULL,
                 control.family = NULL,
                 tag = NULL,
                 options = list(),
                 .envir = parent.frame(),
                 mesh = deprecated(),
                 include = deprecated(),
                 exclude = deprecated(),
                 include_latent = deprecated()) {
  options <- bru_call_options(options)
  bru_options_set_local(options, .reset = TRUE)
  bru_log_message(
    paste0(
      "The `like()` function has been deprecated in favour of `bru_obs()`, ",
      "since inlabru 2.12.0."
    ),
    verbosity = 1L
  )
  lifecycle::deprecate_soft(
    "2.12.0",
    "like()",
    "bru_obs()"
  )

  E <- bru_eval_in_data_context(
    substitute(E),
    data = list(response_data = response_data, data = data),
    default = options[["E"]],
    .envir = .envir
  )
  Ntrials <- bru_eval_in_data_context(
    substitute(Ntrials),
    data = list(response_data = response_data, data = data),
    default = options[["Ntrials"]],
    .envir = .envir
  )
  weights <- bru_eval_in_data_context(
    substitute(weights),
    data = list(response_data = response_data, data = data),
    default = 1,
    .envir = .envir
  )
  scale <- bru_eval_in_data_context(
    substitute(scale),
    data = list(response_data = response_data, data = data),
    default = 1,
    .envir = .envir
  )

  special_env <- new.env(parent = .envir)
  assign("E", E, envir = special_env)
  assign("Ntrials", Ntrials, envir = special_env)
  assign("weights", weights, envir = special_env)
  assign("scale", scale, envir = special_env)

  bru_obs(
    formula = formula,
    family = family,
    data = data,
    response_data = response_data,
    E = E,
    Ntrials = Ntrials,
    weights = weights,
    scale = scale,
    domain = domain,
    samplers = samplers,
    ips = ips,
    include = include,
    exclude = exclude,
    include_latent = include_latent,
    used = used,
    allow_combine = allow_combine,
    control.family = control.family,
    tag = if (is.null(tag) || identical(tag, "")) {
      NA_character_
    } else {
      tag
    },
    options = options,
    .envir = special_env
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


#' @describeIn bru_obs
#' Combine `bru_obs` observation model object into a `bru_obs_list` object
#' @param \dots For `bru_obs_list.bru_obs`, one or more `bru_obs` objects
#' @export
bru_obs_list <- function(...) {
  UseMethod("bru_obs_list")
}

#' @describeIn bru_obs
#' Combine one or more lists of `bru_obs` observation model objects
#' into a `bru_obs_list` object
#' @param object A list of `bru_obs` objects
#' @export
bru_obs_list.list <- function(object, ..., .envir = NULL) {
  object <- lapply(object, as_bru_obs)
  class(object) <- c("bru_obs_list", "list")
  bru_obs_list(object, .envir = .envir)
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
        stop(paste0(
          "Cannot combine objects with mismatching tags and names:\n",
          "  Tag: ", tag_names[k], "\n",
          "  Name: ", list_names[k], "\n",
          "  Use `NA` for either the tag or name,",
          " or use matching tags and names."
        ))
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
#' Combine a list of `bru_obs` observation model objects
#' into a `bru_obs_list` object
#' @param object A list of `bru_obs` objects
#' @export
bru_obs_list.bru_obs_list <- function(..., .envir = NULL) {
  if (length(list(...)) > 1) {
    # If multiple objects are given, combine them into a single list
    # and then recall the method
    object <- lapply(list(...), as_bru_obs_list)
    object <- structure(
      unlist(object, recursive = FALSE),
      class = c("bru_obs_list", "list")
    )
  } else {
    object <- list(...)[[1]]
  }

  if (is.null(.envir)) {
    .envir <- environment(object)
  }

  environment(object) <- .envir
  object <- set_list_names(object, tag = "tag", priority = "immutable")
  object
}



#' @describeIn bru_obs
#' Combine several `bru_obs` objects into a `bru_obs_list` object
#' @export
c.bru_obs <- function(..., .envir = NULL) {
  bru_obs_list(list(...), .envir = .envir)
}


#' @describeIn bru_obs
#' Combine several `bru_obs_list` objects into a `bru_obs_list` object
#' @export
c.bru_obs_list <- function(..., .envir = NULL) {
  bru_obs_list(..., .envir = .envir)
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
      predictor = deparse(object[["formula"]]),
      is_additive = object[["is_additive"]],
      is_linear = object[["linear"]],
      used = object[["used"]],
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
  lifecycle::deprecate_soft(
    "2.12.0",
    "like_list()",
    "bru_obs_list()",
    details = paste0(
      "Use `as_bru_obs_list()`, `bru_obs_list(...)` or ",
      "`c(...)` to construct observation model lists."
    )
  )
  as_bru_obs_list(list(...))
}

#' @describeIn bru_obs `r lifecycle::badge("deprecated")`
#' Backwards compatibility for versions `<= 2.12.0.9017`. For later versions,
#' use `as_bru_obs_list()`, `bru_obs_list()` or `c()`.
#' @export
bru_like_list <- function(...) {
  lifecycle::deprecate_soft(
    "2.12.0.9017",
    "bru_like_list()",
    "bru_obs_list()",
    details = paste0(
      "Use `as_bru_obs_list()`, `bru_obs_list(...)` or ",
      "`c(...)` to construct observation model lists."
    )
  )
  as_bru_obs_list(list(...))
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
print.summary_bru_obs <- function(x, ...) {
  lh <- x
  cat(sprintf(
    paste0(
      "  Family: '%s'\n",
      "    Tag: %s\n",
      "    Data class: %s\n",
      "    Response class: %s\n",
      "    Predictor: %s\n",
      "    Additive/Linear: %s/%s\n",
      "    Used components: %s\n"
    ),
    lh$family,
    if (is.null(lh$tag) || is.na(lh$tag)) {
      "<No tag>"
    } else {
      paste0("'", lh$tag, "'", collapse = ", ")
    },
    paste0("'", lh$data_class, "'", collapse = ", "),
    paste0("'", lh$response_class, "'", collapse = ", "),
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
    },
    as.character(lh$is_additive),
    as.character(lh$is_linear),
    format(lh$used)
  ))
  invisible(x)
}

#' @rdname bru_obs_print
#' @export
print.summary_bru_obs_list <- function(x, ...) {
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
#' @export
#' @keywords internal
#' @returns * `bru_obs_inla_family()` returns a string or vector of strings
#' @rdname bru_obs_methods
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
        paste0(
          "Global control.family option overrides settings in likelihood(s) ",
          paste0(which(like_has_cf), collapse = ", "),
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

bru_obs_expr <- function(lhood, components) {
  if (is.null(lhood[["expr"]])) {
    # Only needed pre-2.12.0.9014
    # Later versions construct expressions for all models,
    # when calling bru_used_update.bru_obs()
    expr_text <- "BRU_EXPRESSION"
    if (utils::packageVersion("inlabru") >= "2.12.0.9014") {
      warning(
        paste0(
          "Code comment in `bru_obs_expr` claims lhood[['expr']] cannot ",
          "be null, but it is null."
        ),
        immediate. = TRUE
      )
    }
  } else {
    expr_text <- as.character(lhood[["expr"]])
  }
  if (grepl(
    pattern = "BRU_EXPRESSION",
    x = expr_text
  )) {
    included <- bru_used(lhood)[["effect"]]
    expr_text <-
      gsub(
        pattern = "BRU_EXPRESSION",
        replacement = paste0(included, collapse = " + "),
        x = expr_text
      )
  }
  parse(text = expr_text)
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
#' @param \dots Further arguments passed on to [bru_obs()]. In particular,
#'   optional `E`, a single numeric used rescale all integration weights by a
#'   fixed factor.
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
#'     geom_fm(data = data$mesh, alpha = 0, linewidth = 0.05) +
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
#' @param drop logical; If `keep=FALSE`, `newdata` is a `Spatial*DataFrame`, and
#'   the prediciton summary has the same number of rows as `newdata`, then the
#'   output is a `Spatial*DataFrame` object. Default `FALSE`.
#' @param \dots Additional arguments passed on to `inla.posterior.sample()`
#' @param data `r lifecycle::badge("deprecated")` Use `newdata` instead.
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
#' generate new realisations for indices greater than `n`.
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
                        data = deprecated(),
                        include = deprecated(),
                        exclude = deprecated()) {
  object <- bru_check_object_bru(object)
  if (lifecycle::is_present(data)) {
    if (is.null(newdata)) {
      lifecycle::deprecate_stop(
        "2.8.0",
        "predict(data)",
        "predict(newdata)",
        details =
          "Only `data` provided and not `newdata`. Use `newdata` only."
      )
    } else {
      lifecycle::deprecate_stop(
        "2.8.0",
        "predict(data)",
        "predict(newdata)",
        details = "Both `newdata` and `data` provided. Use `newdata` only."
      )
    }
  }

  # Convert data into list, data.frame or a Spatial object if not provided as
  # such
  if (is.character(newdata)) {
    newdata <- as.list(setNames(newdata, newdata))
  } else if (inherits(newdata, c("fm_mesh_2d", "inla.mesh"))) {
    lifecycle::deprecate_stop(
      "2.8.0",
      "predict(newdata = 'should not be an `fm_mesh_2d`/`inla.mesh` object')",
      details = paste0(
        "Use 'newdata = fm_vertices(mesh, format = ...)' ",
        "instead of 'newdata = mesh'"
      )
    )
  } else if (inherits(newdata, "formula")) {
    stop(paste0(
      "Formula supplied as data to predict.bru(). ",
      "Please check your argument order/names."
    ))
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
#' @param ... additional, unused arguments.
#' @param data `r lifecycle::badge("deprecated")` Use `newdata` instead.
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
#' generate new realisations for indices greater than `n`.
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
                         data = deprecated(),
                         include = deprecated(),
                         exclude = deprecated()) {
  object <- bru_check_object_bru(object)
  if (lifecycle::is_present(data)) {
    if (is.null(newdata)) {
      lifecycle::deprecate_stop(
        "2.8.0",
        "generate(data)",
        "generate(newdata)",
        details =
          "Only `data` provided and not `newdata`. Use `newdata` only."
      )
    } else {
      lifecycle::deprecate_stop(
        "2.8.0",
        "generate(data)",
        "generate(newdata)",
        details = "Both `newdata` and `data` provided. Use `newdata` only."
      )
    }
  }

  # Convert data into list, data.frame or a Spatial object if not provided as
  # such
  if (is.character(newdata)) {
    newdata <- as.list(setNames(newdata, newdata))
  } else if (inherits(newdata, c("fm_mesh_2d", "inla.mesh"))) {
    lifecycle::deprecate_stop(
      "2.8.0",
      "predict(newdata = 'should not be an `fm_mesh_2d`/`inla.mesh` object')",
      details = paste0(
        "Use 'newdata = fm_vertices(mesh, format = ...)' ",
        "instead of 'newdata = mesh'"
      )
    )
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

    if (is.null(used)) {
      if (lifecycle::is_present(include)) {
        bru_log_message(
          paste0(
            "The `include` argument to `generate.bru()` is deprecated. ",
            "If auto-detection doesn't work, use ",
            "`used = bru_used(effect = include)`."
          ),
          verbosity = 1L
        )
        lifecycle::deprecate_soft(
          "2.12.0.9003",
          "generate(include)",
          "generate(used)",
          "If auto-detection doesn't work, use `bru_used(effect = include)`"
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
        lifecycle::deprecate_soft(
          "2.12.0.9003",
          "generate(exclude)",
          "generate(used)",
          paste0(
            "If auto-detection doesn't work, ",
            "use `bru_used(effect_exclude = exclude)`"
          )
        )
      } else {
        exclude <- NULL
      }
      used <-
        bru_used(
          formula,
          effect = if (is.null(formula)) {
            character(0)
          } else {
            include
          },
          effect_exclude = exclude,
          latent = NULL
        )
    }
    used <- bru_used_update(used, labels = names(object$bru_info$model$effects))

    vals <- evaluate_model(
      model = object$bru_info$model,
      state = state,
      data = newdata,
      predictor = formula,
      used = used
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
      cat(paste0(
        "hpd:",
        min(x),
        " ",
        max(x),
        ", err = ",
        err,
        ", n = ",
        n,
        "\n"
      ))
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
        apply(data, MARGIN = 1, mean, na.rm = TRUE),
        apply(data, MARGIN = 1, sd, na.rm = TRUE)
      )
      qs_names <- NULL
    } else if (length(probs) == 1) {
      smy <- data.frame(
        apply(data, MARGIN = 1, mean, na.rm = TRUE),
        apply(data, MARGIN = 1, sd, na.rm = TRUE),
        as.matrix(apply(data,
          MARGIN = 1, quantile, probs = probs, na.rm = TRUE,
          names = FALSE
        ))
      )
      qs_names <- paste0("q", probs)
    } else {
      smy <- data.frame(
        apply(data, MARGIN = 1, mean, na.rm = TRUE),
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
    skew <- apply(((data - smy$mean) / smy$sd)^3 * (N / (N - 1))^3,
      MARGIN = 1, mean, na.rm = TRUE
    )
    if (max_moment >= 3) {
      smy[["skew"]] <- skew
    }
    # eK + 3 >= skew^2 + 1
    # eK >= skew^2 - 2
    # Use 1/N normalisation of the sample sd
    ekurtosis <- pmax(
      skew^2 - 2,
      apply(((data - smy$mean) / smy$sd)^4 * (N / (N - 1))^4 - 3,
        MARGIN = 1,
        mean,
        na.rm = TRUE
      )
    )
    if (max_moment >= 4) {
      smy[["ekurtosis"]] <- ekurtosis
    }

    # Add Monte Carlo standard errors
    # Var(s) \approx (eK + 2) \sigma^2 / (4 n):
    # +2 replaced by 3-(n-3)/(n-1) = 2n/(n-1), from Rao 1973, p438
    smy[["sd.mc_std_err"]] <-
      sqrt(pmax(0, ekurtosis + 2 * N / (N - 1))) *
        smy[["sd"]] / sqrt(4 * N)
    # Include sd MC error in estimate of mean MC error:
    smy[["mean.mc_std_err"]] <- (
      smy[["sd"]] +
        smy[["sd.mc_std_err"]] * 2
    ) / sqrt(N)
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
            predictor = bru_obs_expr(
              param[["lhoods"]][[lh_idx]],
              param[["model"]][["effects"]]
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
    if (any(!is.finite(delta1)) || any(!is.finite(delta2))) {
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
    nonfin <- any(!is.finite(nonlin_pred))
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
      nonfin <- any(!is.finite(nonlin_pred))
      norm0 <- pred_norm(nonlin_pred - lin_pred0)
      norm1 <- pred_norm(nonlin_pred - lin_pred1)

      bru_log_message(
        paste0(
          "iinla: Step rescaling: ",
          signif(100 * step_scaling, 4),
          "%, Contract",
          " (",
          "norm0 = ", signif(norm0, 4),
          ", norm1 = ", signif(norm1, 4),
          ", norm01 = ", signif(norm01, 4),
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
        paste0(
          "iinla: Step rescaling: ",
          signif(100 * step_scaling, 3),
          "%, Expand",
          " (",
          "norm0 = ", signif(norm0, 4),
          ", norm1 = ", signif(norm1, 4),
          ", norm01 = ", signif(norm01, 4),
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
        paste0(
          "iinla: Step rescaling: ",
          signif(100 * step_scaling, 3),
          "%, Overstep",
          " (",
          "norm0 = ", signif(norm0, 4),
          ", norm1 = ", signif(norm1, 4),
          ", norm01 = ", signif(norm01, 4),
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
      paste0(
        "iinla: Step rescaling: ",
        signif(100 * step_scaling_opt_approx, 4),
        "%, Approx Optimisation",
        " (",
        "norm0 = ", signif(norm0_opt, 4),
        ", norm1 = ", signif(norm1_opt, 4),
        ", norm01 = ", signif(norm01, 4),
        ")"
      ),
      verbosity = 3
    )

    if (norm1_opt > norm01) {
      bru_log_message(
        paste0(
          "iinla: norm1_opt > |delta|: ",
          signif(norm1_opt, 4),
          " > ",
          signif(norm01, 4)
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
        paste0(
          "iinla: Step rescaling: ",
          signif(100 * step_scaling_opt, 4),
          "%, Optimisation",
          " (",
          "norm0 = ", signif(norm0_opt, 4),
          ", norm1 = ", signif(norm1_opt, 4),
          ", norm01 = ", signif(norm01, 4),
          ")"
        ),
        verbosity = 3
      )

      if (norm1_opt > norm01) {
        bru_log_message(
          paste0(
            "iinla: norm1_opt > |delta|: ",
            signif(norm1_opt, 4),
            " > ",
            signif(norm01, 4)
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
    paste0(
      "iinla: |lin1-lin0| = ",
      signif(pred_norm(lin_pred1 - lin_pred0), 4),
      "\n       <eta-lin1,delta>/|delta| = ",
      signif(pred_scalprod(nonlin_pred - lin_pred1, lin_pred1 - lin_pred0) /
        pred_norm(lin_pred1 - lin_pred0), 4),
      "\n       |eta-lin0 - delta <delta,eta-lin0>/<delta,delta>| = ",
      signif(pred_norm(nonlin_pred - lin_pred0 - (lin_pred1 - lin_pred0) *
        pred_scalprod(lin_pred1 - lin_pred0, nonlin_pred - lin_pred0) /
        pred_norm2(lin_pred1 - lin_pred0)), 4)
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
      paste0(
        "iinla: Step rescaling: ",
        signif(100 * step_scaling, 3),
        "%, Maximum step length",
        " (",
        "norm0 = ", signif(norm0, 4),
        ", norm1 = ", signif(norm1, 4),
        ", norm01 = ", signif(norm01, 4),
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
      paste0(
        "iinla: Step rescaling: ",
        signif(100 * step_scaling, 3),
        "%",
        " (",
        "norm0 = ", signif(norm0, 4),
        ", norm1 = ", signif(norm1, 4),
        ", norm01 = ", signif(norm01, 4),
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
        paste0(x, ".", seq_len(length(state[[x]])))
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
  id <- rep(seq_len(length(states)), each = nrow(df[[1]]))
  df <- do.call(rbind, df)
  df[[id_name]] <- id
  df
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
  add_timing <- function(timings, task, iteration = NA_integer_) {
    AbsTime <- proc.time()
    if (!is.na(AbsTime[4])) {
      AbsTime[1] <- AbsTime[1] + AbsTime[4]
    }
    if (!is.na(AbsTime[5])) {
      AbsTime[2] <- AbsTime[2] + AbsTime[2]
    }
    return(rbind(
      timings,
      data.frame(
        Task = task,
        Iteration = iteration,
        Time = AbsTime[1],
        System = AbsTime[2],
        Elapsed = AbsTime[3]
      )
    ))
  }

  options <- bru_call_options(options)
  bru_options_set_local(options, .reset = TRUE)

  timings <- add_timing(NULL, "Start")

  inla.options <- bru_options_inla(options)

  bru_log_bookmark("iinla")
  orig_timings <- NULL
  original_log <- character(0) # Updated further below
  # Local utility method for collecting information object:
  collect_misc_info <- function(...) {
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

    list(
      log = c(original_log, bru_log()["iinla"]),
      states = states,
      inla_stack = stk,
      track = if (is.null(orig_track) ||
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
      },
      inla_track = {
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
        dplyr::bind_rows(orig_inla_track, inla_track)
      },
      timings = {
        iteration_offset <- if (is.null(orig_timings)) {
          0L
        } else {
          max(c(0L, orig_timings$Iteration), na.rm = TRUE)
        }
        rbind(
          orig_timings,
          data.frame(
            Task = timings$Task[-1],
            Iteration = timings$Iteration[-1] + iteration_offset,
            Time = as.difftime(diff(timings$Time), units = "secs"),
            System = as.difftime(diff(timings$System), units = "secs"),
            Elapsed = as.difftime(diff(timings$Elapsed), units = "secs")
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
      control.predictor = list(compute = TRUE)
    )
  )

  # Extract the family of each likelihood
  family <- bru_obs_inla_family(lhoods)

  # Extract the control.family information for each likelihood
  inla.options[["control.family"]] <-
    bru_obs_control_family(lhoods, inla.options[["control.family"]])

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
    msg <- paste0(
      "The total number of response values (N=",
      N_response,
      ") and predictor values (N=",
      N_predictor,
      ") do not match.\n",
      "  This is likely due to a mistake in the component or predictor ",
      "constructions."
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
      paste0("iinla: Iteration ", k, " [max:", options$bru_max_iter, "]"),
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
    }

    timings <- add_timing(timings, "Preprocess", k)

    result <- fm_try_callstack(
      do.call(
        INLA::inla,
        inla.options.merged,
        envir = environment(model$effects)
      )
    )

    timings <- add_timing(timings, "Run inla()", k)

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
          timings <- add_timing(timings, "Line search", k)
        }

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
          msg <- paste0(
            "The total number of response values (N=",
            N_response,
            ") and predictor values (N=",
            N_predictor,
            ") do not match.\n",
            "  This is likely due to a mistake in the component or predictor ",
            "constructions."
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
        paste0(
          "iinla: Max deviation from previous: ",
          signif(100 * max(dev), 3),
          "% of SD, and line search is ",
          if (line_search[["active"]]) "active" else "inactive",
          "\n       [stop if: <", 100 * max.dev, "% and line search inactive]"
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

  result[["bru_iinla"]] <- collect_misc_info()
  class(result) <- c("iinla", class(result))
  return(result)
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
  stopifnot(inherits(components, "bru_comp_list"))

  env <- environment(formula)
  formula <- update.formula(
    formula,
    paste0(
      ". ~ ",
      paste0(names(components), collapse = " + ")
    )
  )
  environment(formula) <- env
  formula
}

# Extract the LHS of a formula, as response ~ .
extract_response <- function(formula) {
  if (inherits(formula, "formula") &&
    (length(as.character(formula)) == 3)) {
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
  return(result)
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
