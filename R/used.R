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
#' (used <- new_bru_used(expression(a + b + c_latent + d_latent)))
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
#' (used <- new_bru_used(expression(a + b + c_latent + d_latent)))
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
  # Map '.' to NULL, to defer decision until later
  # Do not map character(0) to NULL, as it is a valid result
  if (identical(vars, ".")) {
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

#' @rdname bru_used_vars
#' @export
bru_used_vars.quosure <- function(x, functions = FALSE) {
  ex <- deparse1(rlang::quo_get_expr(x), collapse = "\n")
  bru_used_vars(ex, functions = functions)
}

#' @describeIn bru_used_vars Only the right-hand side is used.
#' @export
bru_used_vars.formula <- function(x, functions = FALSE) {
  form <- x[[length(x)]]
  ex <- deparse1(form, collapse = "\n")
  bru_used_vars(ex, functions = functions)
}


#' @describeIn bru_used Create a `bru_used` object by calling [new_bru_used()]
#' @export
bru_used.default <- function(x, ...) {
  new_bru_used(x = x, ...)
}

#' @describeIn new_bru_used Create a `bru_used` object from a `character`
#' representation of an expression.
#' @export
new_bru_used.character <- function(x, ...,
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

  new_bru_used(
    x = NULL,
    ...,
    effect = effect,
    effect_exclude = effect_exclude,
    latent = latent,
    labels = labels
  )
}


#' @describeIn new_bru_used Create a `bru_used` object from an expression
#'   object.
#' @export
new_bru_used.expression <- function(x, ...,
                                    effect = NULL,
                                    effect_exclude = NULL,
                                    latent = NULL,
                                    labels = NULL) {
  attributes(x) <- NULL
  y <- deparse1(x, collapse = "\n")
  new_bru_used(
    x = y,
    ...,
    effect = effect,
    effect_exclude = effect_exclude,
    latent = latent,
    labels = labels
  )
}

#' @describeIn new_bru_used Create a `bru_used` object from an `rlang` `expr`
#'   object.
#' @export
new_bru_used.quosure <- function(x, ...,
                                 effect = NULL,
                                 effect_exclude = NULL,
                                 latent = NULL,
                                 labels = NULL) {
  y <- deparse1(rlang::quo_get_expr(x), collapse = "\n")
  attributes(y) <- NULL
  new_bru_used(
    x = y,
    ...,
    effect = effect,
    effect_exclude = effect_exclude,
    latent = latent,
    labels = labels
  )
}


#' @describeIn new_bru_used Create a `bru_used` object from a formula (only the
#' right-hand side is used).
#' @export
new_bru_used.formula <- function(x, ...,
                                 effect = NULL,
                                 effect_exclude = NULL,
                                 latent = NULL,
                                 labels = NULL) {
  form <- x[[length(x)]]
  ex <- deparse1(form, collapse = "\n")
  new_bru_used(
    x = ex,
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
  bru_used(x[["lhoods"]], ..., join = join)
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
