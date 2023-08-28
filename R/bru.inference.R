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
  object[["bru_info"]] <-
    bru_info_upgrade(
      object[["bru_info"]],
      old = old,
      new_version = new_version
    )
  object
}

bru_info_upgrade <- function(object,
                             old = FALSE,
                             new_version = getNamespaceVersion("inlabru")) {
  if (!is.list(object)) {
    stop("bru_info part of the object can't be converted to `bru_info`; not a list")
  }
  if (!inherits(object, "bru_info")) {
    old <- TRUE
    class(object) <- c("bru_info", "list")
  }
  if (is.null(object[["inlabru_version"]])) {
    object[["inlabru_version"]] <- "0.0.0"
    old <- TRUE
  }
  old_ver <- object[["inlabru_version"]]
  if (old || (utils::compareVersion(new_version, old_ver) > 0)) {
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
      # Make sure component$mapper is a bru_mapper_scale
      for (k in seq_along(object[["model"]][["effects"]])) {
        cmp <- object[["model"]][["effects"]][[k]]
        # Convert offset mappers to const mappers
        if (inherits(cmp$main$mapper, "bru_mapper_offset")) {
          cmp$main$mapper <- bru_mapper_const()
        }
        cmp[["mapper"]] <-
          bru_mapper_scale(
            bru_mapper_multi(list(
              main = cmp$main$mapper,
              group = cmp$group$mapper,
              replicate = cmp$replicate$mapper
            ))
          )
        object[["model"]][["effects"]][[k]] <- cmp
      }
      object[["inlabru_version"]] <- "2.5.3.9005"
    }

    if (utils::compareVersion("2.6.0.9000", old_ver) > 0) {
      message("Upgrading bru_info to 2.6.0.9000")
      # Make sure component$mapper is a bru_mapper_pipe
      for (k in seq_along(object[["model"]][["effects"]])) {
        cmp <- object[["model"]][["effects"]][[k]]
        cmp[["mapper"]] <-
          bru_mapper_pipe(
            list(
              mapper = bru_mapper_multi(list(
                main = cmp$main$mapper,
                group = cmp$group$mapper,
                replicate = cmp$replicate$mapper
              )),
              scale = bru_mapper_scale()
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
        if (identical(deparse(cmp[["replicate"]][["input"]][["input"]]), "NULL")) {
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

      object[["lhoods"]] <-
        bru_used_upgrade(
          object[["lhoods"]],
          labels = names(object[["model"]][["effects"]])
        )
      object[["inlabru_version"]] <- "2.7.0.9021"
    }

    object[["inlabru_version"]] <- new_version
    message(paste0("Upgraded bru_info to ", new_version))
  }
  object
}


#' Methods for bru_info objects
#' @export
#' @method summary bru_info
#' @param object Object to operate on
#' @param verbose logical; If `TRUE`, include more details of the
#' component definitions. If `FALSE`, only show basic component
#' definition information. Default: `TRUE`
#' @param \dots Arguments passed on to other methods
#' @rdname bru_info
summary.bru_info <- function(object, verbose = TRUE, ...) {
  result <- list(
    inlabru_version = object[["inlabru_version"]],
    INLA_version = object[["INLA_version"]],
    components = summary(object[["model"]], verbose = verbose, ...),
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

#' @export
#' @param x A `summary_bru_info` object to be printed
#' @rdname bru_info
print.summary_bru_info <- function(x, ...) {
  cat(paste0("inlabru version: ", x$inlabru_version, "\n"))
  cat(paste0("INLA version: ", x$INLA_version, "\n"))
  cat(paste0("Components:\n"))
  print(x$components)
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


# Used for upgrading from versions <= 2.7.0.9017 to >= 2.7.0.2021
bru_used_upgrade <- function(lhoods, labels) {
  for (k in seq_along(lhoods)) {
    if (is.null(lhoods[[k]][["used"]]) &&
      is.null(lhoods[[k]][["used_components"]])) {
      used <- bru_used(
        NULL,
        effect = lhoods[[k]][["include_components"]],
        latent = lhoods[[k]][["include_latent"]],
        effect_exclude = lhoods[[k]][["exclude_components"]],
      )

      lhoods[[k]][["include_components"]] <- NULL
      lhoods[[k]][["exclude_components"]] <- NULL
      lhoods[[k]][["include_latent"]] <- NULL
    } else if (!is.null(lhoods[[k]][["used_components"]])) {
      used <- bru_used(
        NULL,
        effect = lhoods[[k]][["effect"]],
        latent = lhoods[[k]][["latent"]]
      )
      lhoods[[k]][["used_components"]] <- NULL
    } else {
      used <- bru_used(lhoods[[k]])
    }

    used <- bru_used_update(used, labels = labels)
    lhoods[[k]][["used"]] <- used
  }
  lhoods
}


#' Update used_component information objects
#'
#' Merge available component labels information with used components
#' information.
#'
#' @param x Object to be updated
#' @param labels character vector of component labels
#' @param ... Unused
#'
#' @keywords internal
#' @export
#' @family bru_used
bru_used_update <- function(x, labels, ...) {
  UseMethod("bru_used_update")
}

#' @rdname bru_used_update
#' @export
bru_used_update.bru_like_list <- function(x, labels, ...) {
  for (k in seq_along(x)) {
    x[[k]] <- bru_used_update(x[[k]], labels = labels, ...)
  }
  x
}

#' @rdname bru_used_update
#' @export
bru_used_update.bru_like <- function(x, labels, ...) {
  x[["used"]] <-
    bru_used_update(bru_used(x), labels = labels, ...)
  x
}

#' @rdname bru_used_update
#' @export
bru_used_update.bru_used <- function(x, labels, ...) {
  used <- x
  used$effect <-
    parse_inclusion(
      labels,
      include = used[["effect"]],
      exclude = used[["effect_exclude"]]
    )
  used[["effect_exclude"]] <- NULL
  used$latent <-
    parse_inclusion(
      labels,
      include = used[["latent"]],
      exclude = NULL
    )
  used
}


#' List components used in a model
#'
#' Create or extract information about which components are used by a model, or
#' its individual observation models. If a non-NULL `labels` argument is
#' supplied, also calls [bru_used_update()] on the `bru_used` objects.
#'
#' @param x An object that contains information about used components
#' @param effect character; components used as effects. When `NULL`, auto-detect
#' components to include or include all components.
#' @param effect_exclude character; components to specifically exclude from
#' effect evaluation. When `NULL`, do not specifically exclude any components.
#' @param latent character; components used as `_latent` or `_eval()`. When
#' `NULL`, auto-detect components.
#' @param labels character; component labels passed on to
#' [bru_used_update()]

#' @param join Whether to join list output into a single object; Default
#' may depend on the input object class
#' @param ... Parameters passed on to the other methods
#' @returns A `bru_used` object (a list with elements `effect`
#' and `latent`), or a list of such objects
#' (for methods with `join = FALSE`)
#'
#' @examples
#' (used <- bru_used(~.))
#' bru_used(used, labels = c("a", "c"))
#' (used <- bru_used(~ a + b + c_latent + d_latent))
#' bru_used(used, labels = c("a", "c"))
#' (used <- bru_used(expression(a + b + c_latent + d_latent)))
#' bru_used(used, labels = c("a", "c"))
#'
#' @export
#' @keywords internal
#' @family bru_used
bru_used <- function(x = NULL, ...) {
  # Need to specify the dispatch object explicitly to handle the NULL case:
  UseMethod("bru_used", x)
}

#' @describeIn bru_used Create a `bru_used` object.
#' @export
bru_used.default <- function(x = NULL, ...,
                             effect = NULL,
                             effect_exclude = NULL,
                             latent = NULL,
                             labels = NULL) {
  if (!is.null(x)) {
    stop(paste0(
      "bru_used input class not supported: ",
      "'", paste0(class(x), collapse = "', '"), "'"
    ))
  }

  used <- list(
    effect = effect,
    latent = latent
  )
  used[["effect_exclude"]] <- effect_exclude
  class(used) <- c("bru_used", class(used))

  if (!is.null(labels)) {
    used <- bru_used_update(used, labels = labels)
  }

  used
}


# Function from
# https://stackoverflow.com/questions/63580260/is-there-a-way-to-stop-all-vars-returning-names-from-the-right-hand-side-of
# corrected to handle multiple $ correctly
replace_dollar <- function(expr) {
  if (!is.language(expr) || length(expr) == 1L) {
    return(expr)
  }
  if (expr[[1]] == quote(`$`)) {
    expr[[1]] <- quote(`[[`)
    expr[[3]] <- as.character(expr[[3]])
    expr[[2]] <- replace_dollar(expr[[2]])
    expr[[3]] <- replace_dollar(expr[[3]])
  } else {
    for (i in seq_along(expr)[-1]) {
      expr[[i]] <- replace_dollar(expr[[i]])
    }
  }
  expr
}

#' @title Extract basic variable names from expression
#'
#' @description
#' Extracts the variable names from an R expression by pre- and post-processing
#' around [all.vars()].
#' First replaces `$` with `[[` indexing, so that internal column/variable names
#' are ignored, then calls `all.vars()`.
#'
#' @param x A `formula`, `expression`, or `character`
#' @param functions logical; if TRUE, include function names
#'
#' @returns If successful, a character vector, otherwise `NULL`
#'
#' @examples
#' bru_used_vars(~.)
#' bru_used_vars(~ a + b + c_latent + d_eval())
#' bru_used_vars(expression(a + b + c_latent + d_eval()))
#'
#' bru_used_vars(~., functions = TRUE)
#' bru_used_vars(~ a + b + c_latent + d_eval(), functions = TRUE)
#' bru_used_vars(expression(a + b + c_latent + d_eval()), functions = TRUE)
#'
#' bru_used_vars(a ~ b)
#' bru_used_vars(expression(a ~ b))
#'
#' @keywords internal
#' @export
#' @family bru_used
bru_used_vars <- function(x, functions = FALSE) {
  UseMethod("bru_used_vars")
}

#' @rdname bru_used_vars
#' @export
bru_used_vars.character <- function(x, functions = FALSE) {
  ex <- paste0(x, collapse = "\n")
  ex <- str2lang(ex)
  ex <- replace_dollar(ex)
  vars <- all.vars(ex, functions = functions, unique = TRUE)
  if (identical(vars, ".") || identical(vars, character(0))) {
    vars <- NULL
  }
  vars
}

#' @rdname bru_used_vars
#' @export
bru_used_vars.expression <- function(x, functions = FALSE) {
  attributes(x) <- NULL
  ex <- deparse1(x, collapse = "\n")
  bru_used_vars(ex, functions = functions)
}

#' @describeIn bru_used_vars Only the right-hand side is used.
#' @export
bru_used_vars.formula <- function(x, functions = FALSE) {
  form <- x[[length(x)]]
  ex <- deparse1(form, collapse = "\n")
  bru_used_vars(ex, functions = functions)
}



#' @describeIn bru_used Create a `bru_used` object from a `character`
#' representation of an expression.
#' @export
bru_used.character <- function(x, ...,
                               effect = NULL,
                               effect_exclude = NULL,
                               latent = NULL,
                               labels = NULL) {
  form <- x
  if (is.null(effect)) {
    effect <- bru_used_vars(form, functions = FALSE)
    effect <- effect[
      !grepl("^.*_latent$", effect) &
        !grepl("^.*_eval$", effect)
    ]
  }
  if (is.null(latent)) {
    latent <- bru_used_vars(form, functions = TRUE)

    include_latent <- latent[grepl("^.*_latent$", latent)]
    include_latent <- gsub("_latent$", "", include_latent)
    include_eval <- latent[grepl("^.*_eval$", latent)]
    include_eval <- gsub("_eval$", "", include_eval)
    latent <- union(include_latent, include_eval)
    if (length(latent) == 0) {
      include_latent <- character(0)
    }
  }

  bru_used(
    x = NULL,
    ...,
    effect = effect,
    effect_exclude = effect_exclude,
    latent = latent,
    labels = labels
  )
}


#' @describeIn bru_used Create a `bru_used` object from an expression object.
#' @export
bru_used.expression <- function(x, ...,
                                effect = NULL,
                                effect_exclude = NULL,
                                latent = NULL,
                                labels = NULL) {
  attributes(x) <- NULL
  y <- deparse1(x, collapse = "\n")
  bru_used(
    x = y,
    ...,
    effect = effect,
    effect_exclude = effect_exclude,
    latent = latent,
    labels = labels
  )
}


#' @describeIn bru_used Create a `bru_used` object from a formula (only the
#' right-hand side is used).
#' @export
bru_used.formula <- function(x, ...,
                             effect = NULL,
                             effect_exclude = NULL,
                             latent = NULL,
                             labels = NULL) {
  form <- x[[length(x)]]
  ex <- deparse1(form, collapse = "\n")
  bru_used(
    x = ex,
    ...,
    effect = effect,
    effect_exclude = effect_exclude,
    latent = latent,
    labels = labels
  )
}

#' @rdname bru_used
#' @export
bru_used.bru <- function(x, ..., join = TRUE) {
  bru_used(x[["bru_info"]][["lhoods"]], ..., join = join)
}

#' @rdname bru_used
#' @export
bru_used.list <- function(x, ..., join = TRUE) {
  used <- lapply(x, function(y) bru_used(y, ...))
  if (join) {
    effect_exclude <- unique(unlist(lapply(used, function(y) y[["effect_exclude"]])))
    if (length(effect_exclude) > 0) {
      stop("Cannot join 'bru_used' objects with non-null 'effect_exclude' information.")
    }
    used <-
      bru_used(
        effect = unique(unlist(lapply(used, function(y) y[["effect"]]))),
        latent = unique(unlist(lapply(used, function(y) y[["latent"]])))
      )
  }
  used
}

#' @rdname bru_used
#' @export
bru_used.bru_like <- function(x, ...) {
  bru_used(x[["used"]], ...)
}

#' @describeIn bru_used Convenience method that takes
#' an existing `bru_used` object and calls [bru_used_update()]
#' if `labels` is non-NULL.
#' @export
bru_used.bru_used <- function(x, labels = NULL, ...) {
  if (!is.null(labels)) {
    x <- bru_used_update(x, labels = labels)
  }
  x
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
#' @aliases bru
#' @export
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'
#' @param components A `formula`-like specification of latent components.
#'   Also used to define a default linear additive predictor.  See
#'   [component()] for details.
#' @param ... Likelihoods, each constructed by a calling [like()], or named
#'   parameters that can be passed to a single [like()] call. Note that
#'   all the arguments will be evaluated before calling [like()] in order
#'   to detect if they are `like` objects. This means that
#'   special arguments that need to be evaluated in the context of
#'   `response_data` or `data` (such as Ntrials) may will only work that
#'   way in direct calls to [like()].
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
  stopifnot(bru_safe_inla())

  timing_start <- Sys.time()

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
    lhoods <- list(do.call(
      like,
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

  if (length(lhoods) == 0) {
    stop("No response likelihood models provided.")
  }

  # Turn input into a list of components (from existing list, or a special
  # formula)
  components <- component_list(components, .envir = .envir)

  # Update include/exclude information to limit it to existing components
  lhoods <- bru_used_update(lhoods, labels = names(components))

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

  timing_setup <- Sys.time()

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

  timing_end <- Sys.time()
  result$bru_timings <-
    rbind(
      data.frame(
        Task = c("Preparation"),
        Iteration = rep(NA_integer_, 1),
        Time = c(timing_setup - timing_start)
      ),
      result[["bru_iinla"]][["timings"]]
    )

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
  result <- bru_check_object_bru(result)

  info <- result[["bru_info"]]
  info[["options"]] <- bru_call_options(
    bru_options(
      info[["options"]],
      as.bru_options(options)
    )
  )

  original_timings <- result[["bru_timings"]]

  result <- iinla(
    model = info[["model"]],
    lhoods = info[["lhoods"]],
    initial = result,
    options = info[["options"]]
  )

  timing_end <- Sys.time()
  result$bru_timings <-
    rbind(
      original_timings[1, , drop = FALSE],
      result[["bru_iinla"]][["timings"]]
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

#' Evaluate expressions in the data context
#' @param input An expression to be evaluated
#' @param data Likelihood-specific data, as a `data.frame` or
#' `SpatialPoints[DataFrame]`
#'   object.
#' @param response_data Likelihood-specific data for models that need different
#'  size/format for inputs and response variables, as a `data.frame` or
#' `SpatialPoints[DataFrame]`
#'   object.
#' @param default Value used if the expression is evaluated as NULL. Default
#' NULL
#' @param .envir The evaluation environment
#' @return The result of expression evaluation
#' @keywords internal

eval_in_data_context <- function(input,
                                 data = NULL,
                                 response_data = NULL,
                                 default = NULL,
                                 .envir = parent.frame()) {
  response_data_orig <- response_data
  if (!is.null(response_data)) {
    if (is.list(response_data) && !is.data.frame(response_data)) {
    } else {
      response_data <- as.data.frame(response_data)
    }
  }
  if (!is.null(response_data)) {
    enclos_envir <- new.env(parent = .envir)
    assign(".data.", response_data_orig, envir = enclos_envir)
    result <- try(
      eval(input, envir = response_data, enclos = enclos_envir),
      silent = TRUE
    )
  }
  if (is.null(response_data) || inherits(result, "try-error")) {
    data_orig <- data
    if (!is.null(data)) {
      if (is.list(data) && !is.data.frame(data)) {
      } else {
        data <- as.data.frame(data)
      }
    }
    enclos_envir <- new.env(parent = .envir)
    assign(".data.", data_orig, envir = enclos_envir)
    result <- try(
      eval(input, envir = data, enclos = enclos_envir),
      silent = TRUE
    )
  }
  if (inherits(result, "try-error")) {
    stop(paste0(
      "Input '",
      deparse(input),
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
    the_crs <- fm_crs(dt[[sf_data_idx_[1]]][[nm]])
    for (i in sf_data_idx_) {
      dt_crs <- fm_crs(dt[[i]][[nm]])
      if (!fm_crs_is_identical(dt_crs, the_crs)) {
        dt[[i]][[nm]] <- fm_transform(dt[[i]][[nm]], crs = the_crs)
      }
    }

    ncol_ <- vapply(
      sf_data_idx_,
      function(i) {
        ncol(sf::st_coordinates(dt[[i]][[nm]]))
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
            sf::st_as_sf(
              as.data.frame(cbind(
                sf::st_coordinates(dt[[i]][[nm]]),
                matrix(0.0, nrow(dt[[i]]), ncol_max - ncol_[[i]])
              )),
              coords = seq_len(ncol_max),
              crs = fm_crs(dt[[i]][[nm]])
            )$geometry
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
  result <- do.call(dplyr::bind_rows, dt)
  if (length(sfc_names_) > 0) {
    result <- sf::st_as_sf(result)
  }

  result
}


#' Observation model construction for usage with [bru()]
#'
#' @aliases like
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
#' @param data Likelihood-specific data, as a `data.frame` or
#' `SpatialPoints[DataFrame]`
#'   object.
#' @param response_data Likelihood-specific data for models that need different
#'  size/format for inputs and response variables, as a `data.frame` or
#' `SpatialPoints[DataFrame]`
#'   object.
#' @param mesh Deprecated.
#' @param E Exposure parameter for family = 'poisson' passed on to
#'   `INLA::inla`. Special case if family is 'cp': rescale all integration
#'   weights by E. Default taken from `options$E`, normally `1`.
#' @param Ntrials A vector containing the number of trials for the 'binomial'
#'  likelihood. Default taken from `options$Ntrials`, normally `1`.
#' @param weights Fixed (optional) weights parameters of the likelihood,
#' so the log-likelihood`[i]` is changed into `weights[i] * log_likelihood[i]`.
#' Default value is `1`. WARNING: The normalizing constant for the likelihood
#' is NOT recomputed, so ALL marginals (and the marginal likelihood) must be
#' interpreted with great care.
#' @param samplers Integration domain for 'cp' family.
#' @param ips Integration points for 'cp' family. Overrides `samplers`.
#' @param domain Named list of domain definitions.
#' @param include Character vector of component labels that are used as effects
#'   by the
#'   predictor expression; Default: the result of `[all.vars()]` on the
#'   predictor expression, unless the expression is not ".", in which case
#'   `include=NULL`, to include all components that are not
#'   explicitly excluded. The [bru_used()] methods are used
#'   to extract the variable names, followed by removal of non-component names
#'   when the components are available.
#' @param exclude Character vector of component labels that are not used by the
#'   predictor expression. The exclusion list is applied to the list
#'   as determined by the `include` parameter; Default: NULL (do not remove
#'   any components from the inclusion list)
#' @param include_latent character vector.
#' Specifies which the latent state variables are
#' directly available to the predictor expression, with a `_latent` suffix.
#' This also makes evaluator functions with suffix `_eval` available, taking
#' parameters `main`, `group`, and `replicate`, taking values for where to
#' evaluate the component effect that are different than those defined in the
#' component definition itself (see [component_eval()]). Default `NULL`
#' auto-detects use of `_latent` and `_eval` in the predictor expression.
#' @param used Either `NULL` or a [bru_used()] object, overriding `include`,
#'   `exclude`, and `include_latent`.
#' @param allow_latent `r lifecycle::badge("deprecated")` logical, deprecated.
#'   Use `include_latent` instead.
#' @param allow_combine logical; If `TRUE`, the predictor expression may involve
#'   several rows of the input data to influence the same row. Default `FALSE`,
#'   but forced to `TRUE` if `response_data` is `NULL` or `data` is a `list`
#' @param control.family A optional `list` of `INLA::control.family` options
#' @param options A [bru_options] options object or a list of options passed
#' on to [bru_options()]
#' @param .envir The evaluation environment to use for special arguments (`E`,
#'   `Ntrials`, and `weights`) if not found in `response_data` or `data`.
#'   Defaults to the calling environment.
#'
#' @return A likelihood configuration which can be used to parameterise [bru()].
#'
#' @example inst/examples/like.R

like <- function(formula = . ~ ., family = "gaussian", data = NULL,
                 response_data = NULL, # agg
                 mesh = deprecated(), E = NULL, Ntrials = NULL, weights = NULL,
                 samplers = NULL, ips = NULL, domain = NULL,
                 include = NULL,
                 exclude = NULL,
                 include_latent = NULL,
                 used = NULL,
                 allow_latent = deprecated(),
                 allow_combine = NULL,
                 control.family = NULL,
                 options = list(),
                 .envir = parent.frame()) {
  options <- bru_call_options(options)

  # Some defaults
  inla.family <- family

  formula_char <- as.character(formula)

  # Does the likelihood formula imply a linear predictor?
  linear <- formula_char[length(formula_char)] == "."

  # If not linear, set predictor expression according to the formula's RHS
  if (!linear) {
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
    expr = eval_in_data_context(
      substitute(response_expr),
      data = data,
      response_data = response_data,
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
      expr = eval_in_data_context(
        substitute(response_expr),
        data = data,
        response_data = response_data,
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

  E <- eval_in_data_context(
    substitute(E),
    data = data,
    response_data = response_data,
    default = options[["E"]],
    .envir = .envir
  )
  Ntrials <- eval_in_data_context(
    substitute(Ntrials),
    data = data,
    response_data = response_data,
    default = options[["Ntrials"]],
    .envir = .envir
  )
  weights <- eval_in_data_context(
    substitute(weights),
    data = data,
    response_data = response_data,
    default = 1,
    .envir = .envir
  )

  # More on special bru likelihoods
  if (family == "cp") {
    if (is.null(response)) {
      stop("You called like() with family='cp' but the evaluated response information is NULL")
    }

    if (is.null(ips)) {
      ips <- fm_int(
        domain = domain,
        samplers = samplers,
        int.args = options[["bru_int_args"]]
      )
    }

    if (length(E) > 1) {
      warning("Exposure/effort parameter E should be a scalar for likelihood 'cp'.")
    }

    # TODO!!! ####
    ips_is_Spatial <- inherits(ips, "Spatial")
    if (ips_is_Spatial) {
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
            any(vapply(response, function(x) inherits(x, c("sf", "sfc")), TRUE)))) {
          ips <- sf::st_as_sf(ips)
          ips_is_Spatial <- FALSE
        }
      }
    }
    # TODO: check that the crs info is the same

    # For non-Spatial models:
    # Use the response data list as the actual data object, since that's now the
    # canonical place where the point information is given.  This also allows
    # response_data to be used when constructing the response list.
    # This makes it a strict requirement that the predictor can be evaluated as
    # a pure function of the domain data.  When implementing sf support, might
    # be able to give explicit access to spatial coordinates, but otherwise the
    # user can extract it from the geometry with st_coordinates(geometry)[,1] or
    # similar.
    # For Spatial models, keep the old behaviour for backwards compatibility for
    # now, but can likely realign that in the future after more testing.
    # Save general response data to add to response (precomputed covariates etc)
    if (!is.null(response_data)) {
      data_ <- response_data
    } else {
      data_ <- data
    }
    if (inherits(data_, "Spatial")) {
      data_ <- as.data.frame(data_)
    }
    if (ips_is_Spatial) {
      if ("coordinates" %in% names(response)) {
        idx <- names(response) %in% "coordinates"
        data <- as.data.frame(response$coordinates)
        if (any(!idx)) {
          data <- cbind(data, as.data.frame(response[!idx]))
        }
      } else {
        data <- as.data.frame(response)
      }
      response_data <- NULL
      N_data <- NROW(data)
    } else {
      data <- as.data.frame(response)
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

    if (identical(options[["bru_compress_cp"]], TRUE)) {
      allow_combine <- TRUE
      response_data <- data.frame(
        BRU_E = c(
          0,
          E * ips[["weight"]]
        ),
        BRU_response_cp = c(
          N_data,
          rep(0, NROW(ips))
        )
      )
      if (!linear) {
        expr_text <- formula_char[length(formula_char)]
        expr_text <- paste0(
          "{BRU_eta <- ", expr_text, "\n",
          " c(mean(BRU_eta[BRU_aggregate]), BRU_eta[!BRU_aggregate])}"
        )
      } else {
        expr_text <- paste0(
          "{BRU_eta <- BRU_EXPRESSION\n",
          " c(mean(BRU_eta[BRU_aggregate]), BRU_eta[!BRU_aggregate])}"
        )
      }
      expr <- parse(text = expr_text)

      data <- extended_bind_rows(
        dplyr::bind_cols(data, BRU_aggregate = TRUE),
        dplyr::bind_cols(ips, BRU_aggregate = FALSE)
      )
    } else {
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
    allow_combine <-
      if (!is.null(response_data) ||
        (is.list(data) && !is.data.frame(data))) {
        TRUE
      } else if (is.null(allow_combine)) {
        FALSE
      } else {
        allow_combine
      }

    # Need to make a list instead of data.frame, to allow inla.mdata responses
    response_data <- list(
      BRU_response = response,
      BRU_E = E,
      BRU_Ntrials = Ntrials
    )
    response <- "BRU_response"
  }

  if (is.null(used)) {
    used <- bru_used(
      formula,
      effect = include,
      effect_exclude = exclude,
      latent = include_latent
    )
  }
  if (lifecycle::is_present(allow_latent)) {
    lifecycle::deprecate_warn(
      "2.8.0",
      "like(allow_latent = 'is deprecated')",
      "like(include_latent)",
      details = "The default `like(..., include_latent = NULL)` auto-detects use of `_latent` and `_eval`."
    )
    if (allow_latent && is.null(include_latent)) {
      # Set to NULL so that later bru_used_update adds all components.
      used[["latent"]] <- NULL
    }
  }


  # The likelihood object that will be returned

  lh <- list(
    family = family,
    formula = formula,
    response_data = response_data, # agg
    data = data,
    E = E,
    Ntrials = Ntrials,
    weights = weights,
    samplers = samplers,
    linear = linear,
    expr = expr,
    response = response,
    inla.family = inla.family,
    domain = domain,
    used = used,
    allow_combine = allow_combine,
    control.family = control.family
  )

  class(lh) <- c("bru_like", "list")

  # Return likelihood
  lh
}

# Placeholder for future modularised version of like()
## @describeIn like
## Alias for `like()`, with `obs` standing for "observation model"
## @export
# bru_obs <- like


#' @describeIn like
#' Combine `bru_like` likelihoods into a `bru_like_list` object
#' @param \dots For `like_list.bru_like`, one or more `bru_like` objects
#' @export
like_list <- function(...) {
  UseMethod("like_list")
}

#' @describeIn like
#' Combine a list of `bru_like` likelihoods
#' into a `bru_like_list` object
#' @param object A list of `bru_like` objects
#' @param envir An optional environment for the new `bru_like_list` object
#' @export
like_list.list <- function(object, envir = NULL, ...) {
  if (is.null(envir)) {
    envir <- environment(object)
  }
  if (any(vapply(object, function(x) !inherits(x, "bru_like"), TRUE))) {
    if (any(vapply(object, function(x) inherits(x, "bru_like_list"), TRUE))) {
      stop(paste0(
        "All list elements must be of class 'bru_like'.\n",
        "To combine with 'bru_like_list' objects, use c(...)."
      ))
    }
    stop("All list elements must be of class 'bru_like'.")
  }

  class(object) <- c("bru_like_list", "list")
  environment(object) <- envir
  object
}

#' @describeIn like
#' Combine several `bru_like` likelihoods
#' into a `bru_like_list` object
#' @export
like_list.bru_like <- function(..., envir = NULL) {
  do.call(c, list(..., envir = envir))
}

#' @describeIn like
#' Combine several `bru_like` likelihoods and/or `bru_like_list`
#' objects into a `bru_like_list` object
#' @export
c.bru_like <- function(..., envir = NULL) {
  lst <- lapply(list(...), function(x) {
    if (inherits(x, "bru_like")) {
      list(x)
    } else if (inherits(x, "bru_like_list")) {
      x
    } else {
      stop("Can only combine 'bru_like' and 'bru_like_list' objects.")
    }
  })
  lst <- do.call(c, lst)
  like_list(lst, envir = envir)
}

#' @describeIn like
#' Combine several `bru_like` likelihoods and/or `bru_like_list`
#' objects into a `bru_like_list` object
#' @export
c.bru_like_list <- function(..., envir = NULL) {
  if (!all(vapply(
    list(...),
    function(xx) is.null(xx) || inherits(xx, "bru_like_list"),
    TRUE
  ))) {
    lst <- lapply(list(...), function(x) {
      if (inherits(x, "bru_like")) {
        structure(
          list(x),
          class = "bru_like_list"
        )
      } else if (inherits(x, "bru_like_list")) {
        x
      } else {
        stop("Can only combine 'bru_like' and 'bru_like_list' objects.")
      }
    })
    return(do.call("c", lst))
  }
  object <- NextMethod()
  like_list(object, envir = envir)
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


#' Utility functions for bru likelihood objects
#' @param x Object of `bru_like` or `bru_like_list` type
#' @export
#' @keywords internal
#' @returns * `bru_like_inla_family()` returns a string or vector of strings
#' @rdname bru_like_methods
bru_like_inla_family <- function(x, ...) {
  UseMethod("bru_like_inla_family")
}
#' @export
#' @rdname bru_like_methods
bru_like_inla_family.bru_like <- function(x, ...) {
  x[["inla.family"]]
}
#' @export
#' @rdname bru_like_methods
bru_like_inla_family.bru_like_list <- function(x, ...) {
  vapply(x, bru_like_inla_family, "")
}

#' @param control.family list of INLA `control.family` options to override
#' @export
#' @keywords internal
#' @returns * `bru_like_control_family()` returns a list with `INLA::control.family` options,
#' or a list of such lists, with one element per observation model
#' @rdname bru_like_methods
bru_like_control_family <- function(x, control.family = NULL, ...) {
  UseMethod("bru_like_control_family")
}
#' @export
#' @rdname bru_like_methods
bru_like_control_family.bru_like <- function(x, control.family = NULL, ...) {
  if (!is.null(control.family)) {
    control.family
  } else if (is.null(x[["control.family"]])) {
    list()
  } else {
    x[["control.family"]]
  }
}
#' @export
#' @rdname bru_like_methods
bru_like_control_family.bru_like_list <- function(x, control.family = NULL, ...) {
  # Extract the control.family information for each likelihood
  if (!is.null(control.family)) {
    if (length(control.family) != length(x)) {
      stop("control.family supplied as option, but format doesn't match the number of likelihoods")
    }
    like_has_cf <- vapply(
      seq_along(x),
      function(k) !is.null(x[[k]][["control.family"]]),
      TRUE
    )
    if (any(like_has_cf)) {
      warning("Global control.family option overrides settings in likelihood(s) ",
        paste0(which(like_has_cf)),
        collapse = ", "
      )
    }
  } else {
    control.family <- lapply(x, bru_like_control_family)
  }
  control.family
}

bru_like_expr <- function(lhood, components) {
  if (is.null(lhood[["expr"]])) {
    expr_text <- "BRU_EXPRESSION"
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
#' @aliases lgcp
#' @export
#' @param components A formula describing the latent components
#' @param data A data frame or `SpatialPoints(DataFrame)` object
#' @param samplers A data frame or `Spatial[Points/Lines/Polygons]DataFrame`
#'   objects
#' @param domain Named list of domain definitions
#' @param ips Integration points (overrides `samplers`)
#' @param formula If NULL, the linear combination implied by the `components` is
#'   used as a predictor for the point location intensity. If a (possibly
#'   non-linear) expression is provided the respective Taylor approximation is
#'   used as a predictor. Multiple runs of INLA are then required for a better
#'   approximation of the posterior.
#' @param \dots Further arguments passed on to [like()]. In particular,
#' optional `E`, a single numeric used rescale all integration weights by a fixed
#'   factor.
#' @param options See [bru_options_set()]
#' @param .envir The evaluation environment to use for special arguments
#' (`E`, `Ntrials`, and `weights`) if not found in response_data or data.
#' Defaults to the calling environment.
#' @return An [bru()] object
#' @examples
#' \donttest{
#' if (bru_safe_inla() &&
#'   require(ggplot2, quietly = TRUE) &&
#'   require(fmesher, quietly = TRUE)) {
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
                 samplers = NULL,
                 domain = NULL,
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
  lik <- like(
    family = "cp",
    formula = formula, data = data, samplers = samplers,
    ips = ips, domain = domain,
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
#' @aliases predict.bru
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
#' @param include Character vector of component labels that are needed by the
#'   predictor expression; Default: the result of `[all.vars()]` on the
#'   predictor expression, unless the expression is not ".", in which case
#'   `include=NULL`, to include all components that are not
#'   explicitly excluded. The [bru_used()] methods are used
#'   to extract the variable names, followed by removal of non-component names
#'   when the components are available.
#' @param exclude Character vector of component labels that are not used by the
#'   predictor expression. The exclusion list is applied to the list
#'   as determined by the `include` parameter; Default: NULL (do not remove
#'   any components from the inclusion list)
#' @param used Either `NULL` or a [bru_used()] object, overriding `include` and
#'   `exclude`. Default `NULL`
#' @param drop logical; If `keep=FALSE`, `newdata` is a `Spatial*DataFrame`, and
#'   the prediciton summary has the same number of rows as `newdata`, then the
#'   output is a `Spatial*DataFrame` object. Default `FALSE`.
#' @param \dots Additional arguments passed on to `inla.posterior.sample()`
#' @param data `r lifecycle::badge("deprecated")` Use `newdata` instead.
#' @details
#' In addition to the component names (that give the effect of each component
#' evaluated for the input data), the suffix `_latent` variable name can be used
#' to directly access the latent state for a component, and the suffix function
#' `_eval` can be used to evaluate a component at other input values than the
#' expressions defined in the component definition itself, e.g.
#' `field_eval(cbind(x, y))` for a component that was defined with
#' `field(coordinates, ...)` (see also [component_eval()]).
#'
#' For "iid" models with `mapper = bru_mapper_index(n)`, `rnorm()` is used to
#' generate new realisations for indices greater than `n`.
#'
#' @return a `data.frame`, `sf`, or `Spatial*` object with predicted mean values and
#'   other summary statistics attached. Non-S4 object outputs have the class
#'   "bru_prediction" added at the front of the class list.
#' @example inst/examples/predict.bru.R

predict.bru <- function(object,
                        newdata = NULL,
                        formula = NULL,
                        n.samples = 100,
                        seed = 0L,
                        probs = c(0.025, 0.5, 0.975),
                        num.threads = NULL,
                        include = NULL,
                        exclude = NULL,
                        used = NULL,
                        drop = FALSE,
                        ...,
                        data = deprecated()) {
  object <- bru_check_object_bru(object)
  if (lifecycle::is_present(data)) {
    if (is.null(newdata)) {
      lifecycle::deprecate_warn("2.8.0", "predict(data)", "predict(newdata)",
        details = c("`data` provided but not `newdata`. Setting `newdata <- data`.")
      )
      newdata <- data
    } else {
      lifecycle::deprecate_warn("2.8.0", "predict(data)", "predict(newdata)",
        details = c("Both `newdata` and `data` provided. `data` will be ignored.")
      )
    }
    data <- NULL
  }

  # Convert data into list, data.frame or a Spatial object if not provided as such
  if (is.character(newdata)) {
    newdata <- as.list(setNames(newdata, newdata))
  } else if (inherits(newdata, c("fm_mesh_2d", "inla.mesh"))) {
    lifecycle::deprecate_warn(
      "2.8.0",
      "predict(newdata = 'should not be an fm_mesh_2d/inla.mesh object')",
      details = "Use 'newdata = fm_vertices(mesh, format = ...)' instead of 'newdata = mesh'"
    )
    newdata <- fm_vertices(newdata, format = "sp")
  } else if (inherits(newdata, "formula")) {
    stop("Formula supplied as data to predict.bru(). Please check your argument order/names.")
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
      warning("Some generated list elements are unnamed")
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
#' @aliases generate.bru
#' @export
#' @family sample generators
#' @param object A `bru` object obtained by calling [bru()].
#' @param newdata A data.frame or SpatialPointsDataFrame of covariates needed
#'   for sampling.
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
#' @param include Character vector of component labels that are needed by the
#'   predictor expression; Default: NULL (include all components that are not
#'   explicitly excluded) if `newdata` is provided, otherwise `character(0)`.
#' @param exclude Character vector of component labels that are not used by the
#'   predictor expression. The exclusion list is applied to the list
#'   as determined by the `include` parameter; Default: NULL (do not remove
#'   any components from the inclusion list)
#' @param used Either `NULL` or a [bru_used()] object, overriding `include` and
#'   `exclude`.
#' @param ... additional, unused arguments.
#' @param data Deprecated. Use `newdata` instead.
#' sampling.
#' @details
#' In addition to the component names (that give the effect of each component
#' evaluated for the input data), the suffix `_latent` variable name can be used
#' to directly access the latent state for a component, and the suffix function
#' `_eval` can be used to evaluate a component at other input values than the
#' expressions defined in the component definition itself, e.g.
#' `field_eval(cbind(x, y))` for a component that was defined with
#' `field(coordinates, ...)` (see also [component_eval()]).
#'
#' For "iid" models with `mapper = bru_mapper_index(n)`, `rnorm()` is used to
#' generate new realisations for indices greater than `n`.
#'
#' @return List of generated samples
#' @seealso [predict.bru]
#' @example inst/examples/generate.bru.R
#' @rdname generate

generate.bru <- function(object,
                         newdata = NULL,
                         formula = NULL,
                         n.samples = 100,
                         seed = 0L,
                         num.threads = NULL,
                         include = NULL,
                         exclude = NULL,
                         used = NULL,
                         ...,
                         data = deprecated()) {
  object <- bru_check_object_bru(object)
  if (lifecycle::is_present(data)) {
    if (is.null(newdata)) {
      lifecycle::deprecate_warn("2.8.0", "generate(data)", "generate(newdata)",
        details = c("Both `data` provided but not `newdata`. Setting `newdata <- data`.")
      )
      newdata <- data
    } else {
      lifecycle::deprecate_warn("2.8.0", "generate(data)", "generate(newdata)",
        details = c("Both `newdata` and `data` provided. `data` will be ignored.")
      )
    }
    data <- NULL
  }

  # Convert data into list, data.frame or a Spatial object if not provided as such
  if (is.character(newdata)) {
    newdata <- as.list(setNames(newdata, newdata))
  } else if (inherits(newdata, c("fm_mesh_2d", "inla.mesh"))) {
    lifecycle::deprecate_warn(
      "2.8.0",
      "predict(newdata = 'should not be an fm_mesh_2d/inla.mesh object')",
      details = "Use 'newdata = fm_vertices(mesh, format = ...)' instead of 'newdata = mesh'"
    )
    newdata <- fm_vertices(newdata, format = "sp")
  } else if (inherits(newdata, "formula")) {
    stop("Formula supplied as data to generate.bru(). Please check your argument order/names.")
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
      if (is.null(include)) {
        include <- bru_used_vars(formula)
      }
      used <-
        bru_used(
          effect = if (is.null(newdata) &&
            is.null(include)) {
            character(0)
          } else {
            include
          },
          effect_exclude = exclude,
          latent = NULL,
          labels = names(object$bru_info$model$effects)
        )
    }

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
      cat(paste0("hpd:", min(x), " ", max(x), ", err = ", err, ", n = ", n, "\n"))
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
#' @param x A `data.frame` of data columns that should be added to the summary data frame
#' @param cbind.only If TRUE, only `cbind` the samples and return a matrix where each column is a sample
#' @param max_moment integer, at least 2. Determines the largest moment
#'   order information to include in the output. If `max_moment > 2`,
#'   includes "skew" (skewness, `E[(x-m)^3/s^3]`), and
#'   if `max_moment > 3`, includes
#'   "ekurtosis" (excess kurtosis, `E[(x-m)^4/s^4] - 3`). Default 2.
#'   Note that the Monte Carlo variability of the `ekurtois` estimate may be large.
#' @return A `data.frame` or `Spatial[Points/Pixels]DataFrame` with summary statistics,
#' "mean", "sd", `paste0("q", probs)`, "mean.mc_std_err", "sd.mc_std_err"
#'
#' @examples
#' bru_summarise(matrix(rexp(10000), 10, 1000), max_moment = 4, probs = NULL)
#'
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
    smy[["mean.mc_std_err"]] <- smy[["sd"]] / sqrt(N)
    # Var(s) \approx (eK + 2) \sigma^2 / (4 n):
    # +2 replaced by 3-(n-3)/(n-1) = 2n/(n-1), from Rao 1973, p438
    smy[["sd.mc_std_err"]] <-
      sqrt(pmax(0, ekurtosis + 2 * N / (N - 1))) *
        smy[["sd"]] / sqrt(4 * N)
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
            input = param[["input"]][[lh_idx]],
            state = list(state),
            comp_simple = param[["comp_simple"]][[lh_idx]],
            predictor = bru_like_expr(
              param[["lhoods"]][[lh_idx]],
              param[["model"]][["effects"]]
            ),
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

line_search_optimisation_target_exact <- function(x, param, nonlin_param) {
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
    warning("NULL weights detected for line search. Using weights = 1 instead.",
      immediate. = TRUE
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
    warning(
      paste0(
        "Please notify the inlabru package developer:",
        "\nThe line search linear and nonlinear predictors have different lengths.",
        "\nThis should not happen!"
      ),
      immediate. = TRUE
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
        verbose = options$bru_verbose,
        verbose_store = options$bru_verbose_store,
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
        verbose = options$bru_verbose,
        verbose_store = options$bru_verbose_store,
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
        verbose = options$bru_verbose,
        verbose_store = options$bru_verbose_store,
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
        line_search_optimisation_target,
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
      verbose = options$bru_verbose,
      verbose_store = options$bru_verbose_store,
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
        verbose = options$bru_verbose,
        verbose_store = options$bru_verbose_store,
        verbosity = 3
      )
    }

    if ((norm1_opt > norm01) ||
      identical(options$bru_method$line_opt_method, "full")) {
      alpha <-
        optimise(
          line_search_optimisation_target_exact,
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
        verbose = options$bru_verbose,
        verbose_store = options$bru_verbose_store,
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
          verbose = options$bru_verbose,
          verbose_store = options$bru_verbose_store,
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
        verbose = options$bru_verbose,
        verbose_store = options$bru_verbose_store,
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
    verbose = options$bru_verbose,
    verbose_store = options$bru_verbose_store,
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
      verbose = options$bru_verbose,
      verbose_store = options$bru_verbose_store,
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
      verbose = options$bru_verbose,
      verbose_store = options$bru_verbose_store,
      verbosity = 2
    )
  }

  if (isTRUE(options[["bru_debug"]])) {
    df_debug$nonlinopt <- nonlin_pred

    requireNamespace("ggplot2")
    requireNamespace("patchwork")
    pl1 <- ggplot2::ggplot(df_debug) +
      ggplot2::geom_point(ggplot2::aes(.data$idx, (.data$lin0 - .data$lin1), col = "start")) +
      ggplot2::geom_point(ggplot2::aes(.data$idx, (.data$lin1 - .data$lin1), col = "inla")) +
      ggplot2::geom_point(ggplot2::aes(.data$idx, (.data$nonlin1 - .data$lin1), col = "full")) +
      ggplot2::geom_point(ggplot2::aes(.data$idx, (.data$nonlinopt - .data$lin1), col = "opt")) +
      ggplot2::geom_abline(slope = 0, intercept = 0) +
      ggplot2::scale_color_discrete(breaks = c("start", "inla", "full", "opt")) +
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
      ggplot2::scale_color_discrete(breaks = c("start", "inla", "full", "opt")) +
      ggplot2::ylab("(eta - lin1) * sqrt(w)")
    pl3 <- ggplot2::ggplot(df_debug) +
      ggplot2::geom_point(ggplot2::aes(.data$idx, .data$lin0, col = "start")) +
      ggplot2::geom_point(ggplot2::aes(.data$idx, .data$lin1, col = "inla")) +
      ggplot2::geom_point(ggplot2::aes(.data$idx, .data$nonlin1, col = "full")) +
      ggplot2::geom_point(ggplot2::aes(.data$idx, .data$nonlinopt, col = "opt")) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          as.numeric(.data$idx),
          ymin = .data$lin1 - 2 * .data$weights^-0.5,
          ymax = .data$lin1 + 2 * .data$weights^-0.5
        ),
        alpha = 0.1
      ) +
      ggplot2::geom_abline(slope = 0, intercept = 0) +
      ggplot2::scale_color_discrete(breaks = c("start", "inla", "full", "opt")) +
      ggplot2::ylab("eta")
    pl4 <- ggplot2::ggplot(data.frame(
      idx = seq_along(unlist(state0)),
      state0 = unlist(state0),
      state1 = unlist(state1)[names(unlist(state0))],
      state = unlist(state)[names(unlist(state0))]
    )) +
      ggplot2::geom_point(ggplot2::aes(.data$idx, .data$state0, col = "start")) +
      ggplot2::geom_point(ggplot2::aes(.data$idx, .data$state1, col = "inla")) +
      ggplot2::geom_point(ggplot2::aes(.data$idx, .data$state1, col = "full")) +
      ggplot2::geom_point(ggplot2::aes(.data$idx, .data$state, col = "opt")) +
      ggplot2::scale_color_discrete(breaks = c("start", "inla", "full", "opt")) +
      ggplot2::ylab("state")
    pl1234 <-
      ((pl1 | pl2) / (pl3 | pl4)) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "right")

    # Compute deviations in the delta = lin1-lin0 direction and the orthogonal direction
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
        line_search_optimisation_target_exact(
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
#' @aliases iinla
#' @export
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
  timing_start <- Sys.time()
  timing_setup <- Sys.time()
  timing_iterations <- Sys.time()

  inla.options <- bru_options_inla(options)

  bru_log_bookmark("iinla")
  original_timings <- NULL
  original_log <- character(0) # Updated further below
  # Local utility method for collecting information object:
  collect_misc_info <- function(...) {
    list(
      log = c(original_log, bru_log()["iinla"]),
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
      timings = {
        iteration_offset <- if (is.null(original_timings)) {
          0
        } else {
          max(c(0, original_timings$Iteration), na.rm = TRUE)
        }
        rbind(
          original_timings,
          data.frame(
            Task = c(
              "Setup",
              rep("Iteration", length(timing_iterations) - 1)
            ),
            Iteration = c(
              NA_integer_,
              iteration_offset +
                seq_len(length(timing_iterations) - 1)
            ),
            Time = c(
              timing_setup - timing_start,
              diff(timing_iterations)
            )
          )
        )
      },
      ...
    )
  }

  if (!is.null(inla.options[["offset"]])) {
    stop(paste0(
      "An offset option was specified which may interfere with the inlabru model construction.\n",
      "Please use an explicit offset component instead; e.g. ~ myoffset(value, model = 'offset')"
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
  family <- bru_like_inla_family(lhoods)

  # Extract the control.family information for each likelihood
  inla.options[["control.family"]] <-
    bru_like_control_family(lhoods, inla.options[["control.family"]])

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
    } else {
      old.result
    }
  )

  # Track variables
  track <- list()
  if (is.null(old.result[["bru_iinla"]][["track"]])) {
    original_track <- NULL
    track_size <- 0
  } else {
    original_track <- old.result[["bru_iinla"]][["track"]]
    track_size <- max(original_track[["iteration"]])
  }

  # Preserve old timings
  if (!is.null(old.result[["bru_iinla"]][["timings"]])) {
    original_timings <- old.result[["bru_iinla"]][["timings"]]
  }

  bru_log_message(
    "iinla: Evaluate component inputs",
    verbose = options$bru_verbose,
    verbose_store = options$bru_verbose_store,
    verbosity = 3
  )
  inputs <- evaluate_inputs(model, lhoods = lhoods, inla_f = TRUE)
  bru_log_message(
    "iinla: Evaluate component linearisations",
    verbose = options$bru_verbose,
    verbose_store = options$bru_verbose_store,
    verbosity = 3
  )
  comp_lin <- evaluate_comp_lin(model,
    input = inputs,
    state = states[[length(states)]],
    inla_f = TRUE
  )
  bru_log_message(
    "iinla: Evaluate component simplifications",
    verbose = options$bru_verbose,
    verbose_store = options$bru_verbose_store,
    verbosity = 3
  )
  comp_simple <- evaluate_comp_simple(model,
    input = inputs,
    inla_f = TRUE
  )
  bru_log_message(
    "iinla: Evaluate predictor linearisation",
    verbose = options$bru_verbose,
    verbose_store = options$bru_verbose_store,
    verbosity = 3
  )
  lin <- bru_compute_linearisation(
    model,
    lhoods = lhoods,
    input = inputs,
    state = states[[length(states)]],
    comp_simple = comp_simple
  )

  do_line_search <- (length(options[["bru_method"]][["search"]]) > 0)
  if (do_line_search) {
    # Always compute linearisation (above)
  }

  bru_log_message(
    "iinla: Construct inla stack",
    verbose = options$bru_verbose,
    verbose_store = options$bru_verbose_store,
    verbosity = 3
  )
  # Initial stack
  idx <- evaluate_index(model, lhoods)
  stk <- bru_make_stack(lhoods, lin, idx)

  stk.data <- INLA::inla.stack.data(stk)
  inla.options$control.predictor$A <- INLA::inla.stack.A(stk)
  bru_log_message(
    "iinla: Model initialisation completed",
    verbose = options$bru_verbose,
    verbose_store = options$bru_verbose_store,
    verbosity = 3
  )

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

  timing_setup <- Sys.time()
  timing_iterations <- Sys.time()

  track_df <- NULL
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

    result <- fm_try_callstack(
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
        index = 1,
        iteration = track_size + k,
        mode = result[["misc"]][["configs"]][["max.log.posterior"]],
        sd = Inf
      )
    track[[k]] <- do.call(rbind, track_df)

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
        }
        bru_log_message(
          "iinla: Evaluate component linearisations",
          verbose = options$bru_verbose,
          verbose_store = options$bru_verbose_store,
          verbosity = 3
        )
        comp_lin <- evaluate_comp_lin(model,
          input = inputs,
          state = state,
          inla_f = TRUE
        )
        bru_log_message(
          "iinla: Evaluate predictor linearisation",
          verbose = options$bru_verbose,
          verbose_store = options$bru_verbose_store,
          verbosity = 3
        )
        lin <- bru_compute_linearisation(
          model,
          lhoods = lhoods,
          input = inputs,
          state = state,
          comp_simple = comp_simple
        )
        stk <- bru_make_stack(lhoods, lin, idx)

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
      max.dev <- options$bru_method$rel_tol
      dev <- abs(track[[k - 1]]$mode - track[[k]]$mode) / track[[k]]$sd
      ## do.call(c, lapply(by(do.call(rbind, track),
      ##                           as.factor(do.call(rbind, track)$effect),
      ##                           identity),
      ##                        function(X) { abs(X$mean[k-1] - X$mean[k])/X$sd[k] }))
      bru_log_message(
        paste0(
          "iinla: Max deviation from previous: ",
          signif(100 * max(dev), 3),
          "% of SD, and line search is ",
          if (line_search[["active"]]) "active" else "inactive",
          " [stop if: <", 100 * max.dev, "% and line search inactive]"
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
    track[[k]][["new_linearisation"]] <- NA_real_
    for (label in names(states[[1]])) {
      track[[k]][["new_linearisation"]][track[[k]][["effect"]] == label] <-
        unlist(states[[length(states)]][[label]])
    }
    k <- k + 1

    timing_iterations <- c(timing_iterations, Sys.time())
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
  #  # elist <- elist[unlist(lapply(elist, function(x) !is.function(x)))]

  #  # Remove formulae. This can cause problems as well.
  #  # elist <- elist[unlist(lapply(elist, function(x) !inherits(x, "formula")))]

  # The formula expression is too general for this to be reliable:
  #  # Keep only purse formula variables
  #  # elist <- elist[names(elist) %in% all.vars(formula)]

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
#' [summary.component()].
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

  class(result) <- c("summary_bru", "list")
  return(result)
}


#' @export
#' @param x A `summary_bru` object
#' @rdname summary.bru

print.summary_bru <- function(x, ...) {
  print(x$bru_info, ...)
  print(x$inla, ...)
  invisible(x)
}
