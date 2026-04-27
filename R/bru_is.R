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
  all(bru_is_linear(x[["lhoods"]], ...)) &&
    all(bru_is_linear(x[["model"]], ...))
}
#' @rdname bru_is_linear
#' @export
bru_is_linear.bru_model <- function(x, ...) {
  bru_is_linear(x[["effects"]], ...)
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

bru_is_additive_data_frame <- function(x, root_id = 0, ..., verbose = FALSE) {
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
      add <- bru_is_additive_data_frame(
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
        bru_is_additive_data_frame(
          x = x,
          root_id = x_terms$id[1],
          verbose = verbose
        ),
        bru_is_additive_data_frame(
          x = x,
          root_id = x_terms$id[3],
          verbose = verbose
        )
      )
      return(all(add))
    }
    if ((x_terms$token[1] == "'('") && (x_terms$token[3] == "')'")) {
      if (x_terms$token[2] == "expr") {
        if (verbose) {
          message("(expr) found, may be additive")
        }
        add <- bru_is_additive_data_frame(
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

  FALSE
}

#' @rdname bru_is_additive
#' @export
bru_is_additive.character <- function(x, ..., verbose = FALSE) {
  bru_is_additive_data_frame(
    bru_get_parse_data(x),
    root_id = 0,
    ...
  )
}

#' @rdname bru_is_additive
#' @export
bru_is_additive.expression <- function(x, ..., verbose = FALSE) {
  bru_is_additive_data_frame(
    bru_get_parse_data(as.character(x)),
    root_id = 0,
    ...
  )
}

#' @rdname bru_is_additive
#' @export
bru_is_additive.formula <- function(x, ..., verbose = FALSE) {
  bru_is_additive_data_frame(
    bru_get_parse_data(as.character(x)[length(x)]),
    root_id = 0,
    ...
  )
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
  vapply(x, function(lh) bru_is_linear(lh, ...), logical(1))
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
