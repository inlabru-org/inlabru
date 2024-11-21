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
#' @returns An updated version of `x`
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


#' @describeIn bru_used Text formatting method for `bru_used` objects.
#' @export
format.bru_used <- function(x, ...) {
  if (is.null(x[["effect"]])) {
    s <- paste0("effects[<not yet initialised>]")
  } else {
    s <- paste0("effects[", paste0(x$effect, collapse = ", "), "]")
  }
  if (is.null(x[["latent"]])) {
    s <- paste0(s, ", latent[<not yet initialised>]")
  } else {
    s <- paste0(s, ", latent[", paste0(x$latent, collapse = ", "), "]")
  }
  if (!is.null(x$effect_exclude)) {
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
