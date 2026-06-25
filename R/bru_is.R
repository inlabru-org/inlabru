#' @title Check for predictor linearity
#' @description Checks if a predictor expression is linear (or affine)
#' @param x A object containing a predictor definition
#' @param \dots Arguments passed on recursively.
#' @return `TRUE` if the expression is detected to be linear, `FALSE`
#'   otherwise.
#' @export
#' @examples
#' bru_is_linear(new_bru_pred_expr(~ x + y))
#' bru_is_linear(new_bru_pred_expr(~ x * y))
#' bru_is_linear(bm_scale())
#' bru_is_linear(bm_logsumexp())
#'
bru_is_linear <- function(x, ...) {
  UseMethod("bru_is_linear")
}

#' @rdname bru_is_linear
#' @export
bru_is_linear.bru <- function(x, ...) {
  all(bru_is_linear(x[["bru_info"]], ...))
}
#' @rdname bru_is_linear
#' @export
bru_is_linear.bru_info <- function(x, ...) {
  all(bru_is_linear(as_bru_obs_list(x), ...)) &&
    all(bru_is_linear(as_bru_comp_list(x), ...))
}
#' @rdname bru_is_linear
#' @export
bru_is_linear.bru_model <- function(x, ...) {
  bru_is_linear(as_bru_comp_list(x), ...)
}
#' @rdname bru_is_linear
#' @export
bru_is_linear.bru_obs <- function(x, ...) {
  bru_is_linear(x[["pred_expr"]])
}
#' @rdname bru_is_linear
#' @export
bru_is_linear.bru_obs_list <- function(x, ...) {
  vapply(x, function(lh) bru_is_linear(lh, ...), logical(1))
}
#' @rdname bru_is_linear
#' @export
bru_is_linear.bru_pred_expr <- function(x, ...) {
  isTRUE(x[["is_linear"]])
}
#' @rdname bru_is_linear
#' @export
bru_is_linear.bru_comp_list <- function(x, ...) {
  vapply(x, bru_is_linear, logical(1))
}
#' @rdname bru_is_linear
#' @export
bru_is_linear.bru_comp <- function(x, ...) {
  bru_is_linear(x[["mapper"]])
}
#' @rdname bru_is_linear
#' @export
bru_is_linear.bru_mapper <- function(x, ...) {
  ibm_is_linear(x)
}


#' @title Check for predictor expression additivity
#' @description Checks if a predictor expression is additive or not
#' @param x A predictor `expression`, `formula`, or parse information
#'   `data.frame`.
#' @param \dots Arguments passed on recursively.
#' @param verbose logical; if `TRUE`, print diagnostic parsing information.
#' @return `TRUE` if the expression is detected to be additive, `FALSE`
#'   otherwise.
#' @export
#' @examples
#' bru_is_additive(~ x + y)
#' bru_is_additive(~ x * y)
#'
bru_is_additive <- function(x, ...) {
  UseMethod("bru_is_additive")
}

#' @rdname bru_is_additive
#' @export
bru_is_additive.default <- function(x, ..., verbose = FALSE) {
  if (inherits(x, "{") || inherits(x, "(")) {
    if (length(x) < 2L) {
      if (verbose) {
        message("Empty '{}' or '()' found, assuming non-additive.")
      }
      return(FALSE)
    }
    return(bru_is_additive(x[[2]], ..., verbose = verbose))
  }
  if (verbose) {
    message(paste0("General class ",
                   paste0("'", class(x), "'", collapse = ", "),
                   " found, assuming non-additive."))
  }
  FALSE
}
#' @rdname bru_is_additive
#' @export
bru_is_additive.numeric <- function(x, ..., verbose = FALSE) {
  if (verbose) {
    message("Numeric found, assuming additive.")
  }
  TRUE
}
#' @rdname bru_is_additive
#' @export
bru_is_additive.name <- function(x, ..., verbose = FALSE) {
  if (verbose) {
    message("Symbol found, assuming additive.")
  }
  TRUE
}
#' @rdname bru_is_additive
#' @export
bru_is_additive.call <- function(x, ..., verbose = FALSE) {
  if (x[[1]] != as.name("+")) {
    if (verbose) {
      message("Non-'+' call found, assuming non-additive")
    }
    return(FALSE)
  }
  if (verbose) {
    message("'+' found, may be additive")
  }

  ok <- bru_is_additive(x[[2]], ..., verbose = verbose)
  if (length(x) == 2L) {
    return(ok)
  }

  ok && bru_is_additive(x[[3]], ..., verbose = verbose)
}
#' @rdname bru_is_additive
#' @export
bru_is_additive.expression <- function(x, ..., verbose = FALSE) {
  all(vapply(x, bru_is_additive, ..., verbose = verbose,
             FUN.VALUE = logical(1)))
}
#' @rdname bru_is_additive
#' @export
bru_is_additive.quosure <- function(x, ..., verbose = FALSE) {
  bru_is_additive(rlang::quo_get_expr(x), ..., verbose = verbose)
}
#' @rdname bru_is_additive
#' @export
bru_is_additive.formula <- function(x, ..., verbose = FALSE) {
  bru_is_additive(rlang::as_quosure(x), ..., verbose = verbose)
}
#' @rdname bru_is_additive
#' @export
bru_is_additive.character <- function(x, ..., verbose = FALSE) {
  bru_is_additive(rlang::parse_expr(x), ..., verbose = verbose)
}


#' @rdname bru_is_additive
#' @export
bru_is_additive.bru_pred_expr <- function(x, ...) {
  isTRUE(x[["is_additive"]])
}
#' @rdname bru_is_additive
#' @export
bru_is_additive.bru_obs <- function(x, ...) {
  bru_is_additive(x[["pred_expr"]], ...)
}
#' @rdname bru_is_additive
#' @export
bru_is_additive.bru_obs_list <- function(x, ...) {
  vapply(x, function(lh) bru_is_additive(lh, ...), logical(1))
}


#' @title Check for predictor rowwise evaluability
#' @description Checks if a predictor expression may be evaluated rowwise
#' @param x An object containing a predictor definition
#' @param \dots Arguments passed on to submethods.
#' @return `TRUE` if the expression is believed to be rowwise, `FALSE`
#'   otherwise.
#' @export
#' @examples
#' bru_is_rowwise(new_bru_pred_expr(~ x + y))
#'
bru_is_rowwise <- function(x, ...) {
  UseMethod("bru_is_rowwise")
}

#' @rdname bru_is_rowwise
#' @export
bru_is_rowwise.bru_pred_expr <- function(x, ...) {
  isTRUE(x[["is_rowwise"]])
}
#' @rdname bru_is_rowwise
#' @export
bru_is_rowwise.bru_obs <- function(x, ...) {
  bru_is_rowwise(x[["pred_expr"]], ...)
}
#' @rdname bru_is_rowwise
#' @export
bru_is_rowwise.bru_obs_list <- function(x, ...) {
  vapply(x, function(lh) bru_is_rowwise(lh, ...), logical(1))
}
#' @rdname bru_is_rowwise
#' @export
bru_is_rowwise.bru_comp_list <- function(x, ...) {
  vapply(x, bru_is_rowwise, logical(1))
}
#' @rdname bru_is_rowwise
#' @export
bru_is_rowwise.bru_comp <- function(x, ...) {
  bru_is_rowwise(x[["mapper"]])
}
#' @rdname bru_is_rowwise
#' @export
bru_is_rowwise.bru_mapper <- function(x, ...) {
  ibm_is_rowwise(x)
}
