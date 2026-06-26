#' Update used_component information objects
#'
#' Merge available component labels information with used components
#' information.
#'
#' @param x Object to be updated
#' @param labels character vector of component labels
#' @param \dots Unused
#' @returns An updated version of `x`
#' @keywords internal
#' @export
#' @family bru_used
bru_used_update <- function(x, labels, ...) {
  UseMethod("bru_used_update")
}

#' @rdname bru_used_update
#' @export
bru_used_update.bru_obs_list <- function(x, labels, ...) {
  for (k in seq_along(x)) {
    x[[k]] <- bru_used_update(x[[k]], labels = labels, ...)
  }
  x
}

#' @rdname bru_used_update
#' @export
bru_used_update.bru_obs <- function(x, labels, ...) {
  x[["pred_expr"]] <- bru_used_update(
    x[["pred_expr"]],
    labels = labels,
    ...
  )
  x <- bru_compat_pre_2_14_bru_obs(x)
  x
}

#' @rdname bru_used_update
#' @export
bru_used_update.bru_pred_expr <- function(x, labels, ...) {
  pre_used <- bru_used(x)
  used <- bru_used_update(pre_used, labels = labels, ...)
  if (isTRUE(x[["is_additive"]])) {
    if ((length(used$latent) > 0) ||
      (length(setdiff(pre_used$effect, used$effect)) > 0)) {
      x[["is_additive"]] <- FALSE
      x[["is_linear"]] <- FALSE

      if (is.null(x[["pred_expr"]])) {
        pred_text <- paste0(pre_used$effect, collapse = " + ")
        if (length(pre_used$latent) > 0) {
          x[["is_rowwise"]] <- FALSE
          pred_text <- paste(
            pred_text,
            paste0(pre_used$latent, "_latent", collapse = " + "),
            sep = " + "
          )
        }
        x[["pred_text"]] <- pred_text
        x[["pred_expr"]] <- rlang::parse_expr(x[["pred_text"]])
      }
    } else {
      if (is.null(x[["pred_expr"]])) {
        x[["pred_text"]] <- paste0(used$effect, collapse = " + ")
        x[["pred_expr"]] <- rlang::parse_expr(x[["pred_text"]])
      }
    }
  }
  x[["used"]] <- used
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
  if (is.null(used[["latent"]])) {
    used$latent <- character(0)
  }
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
#' @param labels character; component labels passed on to
#' [bru_used_update()]

#' @param join Whether to join list output into a single object; Default
#' may depend on the input object class
#' @param \dots Parameters passed on to the other methods
#' @returns A `bru_used` object (a list with elements `effect`
#' and `latent`), or a list of such objects
#' (for methods with `join = FALSE`)
#'
#' @examples
#' (used <- new_bru_used(~.))
#' bru_used(used, labels = c("a", "c"))
#' (used <- new_bru_used(~ a + b + c_latent + d_latent))
#' bru_used(used, labels = c("a", "c"))
#' (used <- new_bru_used(~ a_eval(.latent$c)))
#' bru_used(used, labels = c("a", "c"))
#'
#' @export
#' @keywords internal
#' @family bru_used
bru_used <- function(x = NULL, ...) {
  # Need to specify the dispatch object explicitly to handle the NULL case:
  UseMethod("bru_used", x)
}

#' @title Store information about used components
#'
#' @description Create a `bru_used` object from effect name character
#'   vectors.
#' Create information about which components are used by a model, or
#' its individual observation models. If a non-NULL `labels` argument is
#' supplied, also calls [bru_used_update()] on the `bru_used` object.
#'
#' @param x `NULL`, or an object representing an expression
#' @param effect character; components used as effects. When `NULL`, auto-detect
#' components to include all components in a predictor expression.
#' @param effect_exclude character; components to specifically exclude from
#' effect evaluation. When `NULL`, do not specifically exclude any components.
#' @param latent character; components used as `_latent` or `_eval()`. When
#' `NULL`, auto-detect components.
#' @param labels character; component labels passed on to
#' [bru_used_update()]

#' @param \dots Parameters passed on to the other methods
#' @returns A `bru_used` object (a list with elements `effect`
#' and `latent`)
#'
#' @details
#' The arguments `effect`, `effect_exclude`, and `latent` control what
#' components and effects are available for use in predictor expressions.
#' \describe{
#'   \item{`effect`}{
#'   Character vector of component labels that are used as effects
#'   by the predictor expression; If `NULL` (default), the names
#'   are extracted from the formula.
#'   }
#'   \item{`exclude`}{
#'   Character vector of component labels to be excluded from the effect list
#'   even if they have been auto-detected as being necessary.
#'   Default is `NULL`; do not remove any components from the inclusion list.
#'   }
#'   \item{`include_latent`}{Character vector.
#'   Specifies which latent state variables need to be directly available to the
#'   predictor expression, with a `_latent` suffix. This also makes evaluator
#'   functions with suffix `_eval` available, taking parameters `main`, `group`,
#'   and `replicate`, taking values for where to evaluate the component effect
#'   that are different than those defined in the component definition itself
#'   (see [bru_comp_eval()]). If `NULL`, the use of `_latent` and `_eval`
#'   in the predictor expression is detected automatically.
#'   }
#' }
#'
#' @examples
#' (used <- new_bru_used(~.))
#' bru_used(used, labels = c("a", "c"))
#' (used <- new_bru_used(~ a + b + c_latent + d_latent))
#' bru_used(used, labels = c("a", "c"))
#'
#' @export
#' @keywords internal
#' @family bru_used
#' @export
new_bru_used <- function(x = NULL,
                         ...,
                         effect = NULL,
                         effect_exclude = NULL,
                         latent = NULL,
                         labels = NULL) {
  # Need to specify the dispatch object explicitly to handle the NULL case:
  UseMethod("new_bru_used", x)
}

#' @describeIn new_bru_used Create a `bru_used` object from effect name
#'   character vectors.
#' @export
new_bru_used.NULL <- function(x = NULL,
                              ...,
                              effect = NULL,
                              effect_exclude = NULL,
                              latent = NULL,
                              labels = NULL) {
  used <- structure(
    list(
      effect = effect,
      latent = latent
    ),
    class = "bru_used"
  )
  used[["effect_exclude"]] <- effect_exclude

  if (!is.null(labels)) {
    used <- bru_used_update(used, labels = labels)
  }

  used
}

#' @describeIn bru_used Create a `bru_used` object by calling [new_bru_used()].
#' @export
bru_used.NULL <- function(x = NULL, ...) {
  new_bru_used(x = x, ...)
}


# Function from
# https://stackoverflow.com/questions/63580260/
#   is-there-a-way-to-stop-all-vars-returning-names-from-the-right-hand-side-of
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
      if (!is.null(expr[[i]])) {
        expr[[i]] <- replace_dollar(expr[[i]])
      }
    }
  }
  expr
}

#' @title Extract basic variable names from expression
#'
#' @description
#' Extracts the variable names and function names from an R expression by
#' traversing the expression structure. Internal helper function for
#' [new_bru_used()] and [bru_used()].
#'
#' @param x A `formula`, `expression`, or other supported class. For the
#'   `format` and `print` methods, a `bru_used_vars` object.
#' @param result A `bru_used_vars` object; a list with elements `vars`, `funs`,
#'   and `objects`, by default provided by [new_bru_used_vars()]
#'
#' @returns A `bru_used_vars` object with elements
#'   \describe{
#'   \item{vars}{character; names of directly accessed variables.}
#'   \item{funs}{character; names of functions called.}
#'   \item{objects}{named list; one character vector per container objects with
#'   variables names accessed via `$`, `[[`, or
#'   `[`. If
#'   the access is ambiguous, the container object name is stored in `vars`.
#'   }
#'   }
#'
#' @examples
#' bru_used_vars(~.)
#' bru_used_vars(~ a + b + c_latent + d_eval())
#'
#' # Ignores the LHS:
#' bru_used_vars(a ~ b)
#'
#' # Detects variables accessed via pronouns and objects,
#' # as well as function calls:
#' bru_used_vars(~ cos(x$z) + y_eval() + .latent$"q")
#'
#' @keywords internal
#' @export
#' @family bru_used
bru_used_vars <- function(x, result = new_bru_used_vars()) {
  UseMethod("bru_used_vars")
}
#' @export
#' @describeIn bru_used_vars Create a `bru_used_vars` object.
new_bru_used_vars <- function(x = list(
                                vars = character(0),
                                funs = character(0),
                                objects = list()
                              )) {
  stopifnot(is.list(x))
  stopifnot(all(c("vars", "funs", "objects") %in% names(x)))
  stopifnot(is.character(x$vars))
  stopifnot(is.character(x$funs))
  stopifnot(is.list(x$objects))
  stopifnot(length(x$objects) == length(names(x$objects)))
  structure(
    x,
    class = "bru_used_vars"
  )
}
#' @rdname bru_used_vars
#' @export
bru_used_vars.default <- function(x, result = new_bru_used_vars()) {
  if (inherits(x, "{") || inherits(x, "(")) {
    for (i in seq_along(x)[-1]) {
      if (!is.null(x[[i]])) {
        result <- bru_used_vars(x[[i]], result = result)
      }
    }
    return(result)
  }
  vars <- all.vars(x, functions = FALSE)
  result$vars <- union(result$vars, vars)
  if (!is.null(result$funs)) {
    funs <- all.vars(x, functions = TRUE)
    funs <- setdiff(funs, vars)
    result$funs <- union(result$funs, funs)
  }
  result
}
#' @rdname bru_used_vars
#' @export
`bru_used_vars.<-` <- function(x, result = new_bru_used_vars()) {
  if (!is.null(result$funs)) {
    result$funs <- union(result$funs, "<-")
  }
  # Ignore variables on the LHS on the assignment; they are outputs, not inputs.
  # Handle the RHS:
  for (i in seq_along(x)[-c(1, 2)]) {
    if (!is.null(x[[i]])) {
      result <- bru_used_vars(x[[i]], result = result)
    }
  }
  result
}
#' @rdname bru_used_vars
#' @export
bru_used_vars.call <- function(x, result = new_bru_used_vars()) {
  fun <- deparse1(x[[1]], collapse = "")
  if (fun %in% c("$", "[", "[[")) {
    # Unambiguous access to a variable via a pronoun or object,
    # e.g. .data$var or .data.[["var"]]; add `var` to result$objects$.data, etc.
    # Otherwise add the container object to result$vars, as well as recursively for
    # the accessor(s), as the access is ambiguous and might involve the entire object.
    obj <- as.character(x[[2]])
    if ((length(x) == 3) && is.symbol(x[[2]]) && (
      (fun == "$") ||
        ((fun %in% c("[", "[[")) && is.character(x[[3]]))
    )) {
      result$objects[[obj]] <- union(
        result$objects[[obj]],
        as.character(x[[3]])
      )
    } else if (is.symbol(x[[2]])) {
      result$vars <- union(result$vars, obj)
      for (i in seq_along(x)[-c(1, 2)]) {
        if (!is.null(x[[i]])) {
          result <- bru_used_vars(x[[i]], result = result)
        }
      }
    }
    return(result)
  }
  if (!is.null(result$funs) && is.symbol(x[[1]])) {
    result$funs <- union(result$funs, fun)
  }
  for (i in seq_along(x)[-1]) {
    if (!is.null(x[[i]])) {
      result <- bru_used_vars(x[[i]], result = result)
    }
  }
  result
}
#' @rdname bru_used_vars
#' @export
bru_used_vars.expression <- function(x, result = new_bru_used_vars()) {
  for (i in seq_along(x)) {
    if (!is.null(x[[i]])) {
      result <- bru_used_vars(x[[i]], result = result)
    }
  }
  result
}
#' @rdname bru_used_vars
#' @export
bru_used_vars.quosure <- function(x, result = new_bru_used_vars()) {
  bru_used_vars(rlang::quo_get_expr(x), result = result)
}
#' @rdname bru_used_vars
#' @export
bru_used_vars.formula <- function(x, result = new_bru_used_vars()) {
  bru_used_vars(rlang::as_quosure(x), result = result)
}
#' @rdname bru_used_vars
#' @export
format.bru_used_vars <- function(x, ...) {
  paste0(
    "vars: {", paste0(x$vars, collapse = ", "), "}, ",
    "funs: {", paste0(x$funs, collapse = ", "), "}, ",
    "objects[",
    if (length(x$objects) == 0) {
      ""
    } else {
      paste0(
        names(x$objects),
        ": {",
        vapply(x$objects, function(y) {
          paste0(y, collapse = ", ")
        }, character(1)),
        "}",
        collapse = ", "
      )
    },
    "]"
  )
}
#' @rdname bru_used_vars
#' @export
print.bru_used_vars <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}


#' @describeIn bru_used Create a `bru_used` object by calling [new_bru_used()]
#' @export
bru_used.default <- function(x, ...) {
  new_bru_used(x = x, ...)
}

#' @describeIn new_bru_used Create a `bru_used` object from an expression
#' object supported by [bru_used_vars()].
#' @export
new_bru_used.default <- function(x, ...,
                                 effect = NULL,
                                 effect_exclude = NULL,
                                 latent = NULL,
                                 labels = NULL) {
  if (is.null(effect) || is.null(latent)) {
    result <- bru_used_vars(x)
    if (is.null(effect)) {
      if ("." %in% result$vars) {
        effect <- NULL
      } else {
        effect <- union(
          result$vars[!grepl("^.*_latent$", result$vars)],
          c(
            result$objects[[".effect"]],
            result$objects[[".effect."]]
          )
        )
      }
    }
    if (is.null(latent)) {
      include_latent <- result$vars[grepl("^.*_latent$", result$vars)]
      include_latent <- gsub("_latent$", "", include_latent)
      include_latent <- union(
        include_latent,
        c(
          result$objects[[".latent"]],
          result$objects[[".latent."]]
        )
      )
      include_eval <- result$funs[grepl("^.*_eval$", result$funs)]
      include_eval <- gsub("_eval$", "", include_eval)
      latent <- union(include_latent, include_eval)
    }
  }

  new_bru_used(
    x = NULL,
    ...,
    effect = effect,
    effect_exclude = effect_exclude,
    latent = latent,
    labels = labels
  )
}


#' @describeIn new_bru_used Create a `bru_used` object from a string.
#' @export
new_bru_used.character <- function(x, ...,
                                   effect = NULL,
                                   effect_exclude = NULL,
                                   latent = NULL,
                                   labels = NULL) {
  new_bru_used(
    x = rlang::parse_expr(x),
    ...,
    effect = effect,
    effect_exclude = effect_exclude,
    latent = latent,
    labels = labels
  )
}

#' @describeIn bru_used Extract the `bru_used` information for the collection
#' of observation models used in a `bru` object.
#' @export
bru_used.bru <- function(x, ..., join = TRUE) {
  bru_used(x[["bru_info"]], ..., join = join)
}

#' @describeIn bru_used Extract the `bru_used` information for the collection
#' of observation models used in a `bru_info` object.
#' @export
bru_used.bru_info <- function(x, ..., join = TRUE) {
  bru_used(as_bru_obs_list(x), ..., join = join)
}

#' @describeIn bru_used Extract the `bru_used` information for each element
#'   of a list, and optionally join into a single `bru_used` object.
#' @export
bru_used.list <- function(x, ..., join = TRUE) {
  used <- lapply(x, function(y) bru_used(y, ...))
  if (join) {
    effect_exclude <- unique(unlist(lapply(used, function(y) {
      y[["effect_exclude"]]
    })))
    if (length(effect_exclude) > 0) {
      stop(paste0(
        "Cannot join 'bru_used' objects with non-null ",
        "'effect_exclude' information."
      ))
    }
    used <-
      bru_used(
        effect = unique(unlist(lapply(used, function(y) y[["effect"]]))),
        latent = unique(unlist(lapply(used, function(y) y[["latent"]])))
      )
  }
  used
}

#' @describeIn bru_used Extract the `bru_used` information for the observation
#'   model predictor used in a `bru` observation model `bru_obs` object.
#' @export
bru_used.bru_obs <- function(x, ...) {
  bru_used(x[["pred_expr"]], ...)
}

#' @describeIn bru_used Extract the `bru_used` information for an observation
#'   model predictor `bru_pred_expr` object.
#' @export
bru_used.bru_pred_expr <- function(x, ...) {
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


#' @describeIn bru_used Text formatting method for `bru_used` objects.
#' @export
format.bru_used <- function(x, ...) {
  if (is.null(x[["effect"]])) {
    s <- paste0("effect[<not yet initialised>]")
  } else {
    s <- paste0("effect[", paste0(x$effect, collapse = ", "), "]")
  }
  if (is.null(x[["latent"]])) {
    s <- paste0(s, ", latent[<not yet initialised>]")
  } else {
    s <- paste0(s, ", latent[", paste0(x$latent, collapse = ", "), "]")
  }
  if (!is.null(x[["effect_exclude"]])) {
    s <- paste0(s, ", exclude[", paste0(x$effect_exclude, collapse = ", "), "]")
  }
  s
}

#' @describeIn bru_used Summary method for `bru_used` objects.
#' @export
#' @method summary bru_used
summary.bru_used <- function(object, ...) {
  object
}

#' @describeIn bru_used Print method for `bru_used` objects.
#' @export
print.bru_used <- function(x, ...) {
  cat("Used ", format(x), "\n", sep = "")
  invisible(x)
}
