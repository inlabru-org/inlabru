#' @title Predictor expression handler object
#'
#' @description
#' Create and access predictor expression information, keeping track of what
#' latent variables are used, whether the predictor is additive/linear/rowwise.
#'
#' @param x For creation, a formula or character string giving the predictor
#' expression. If the character string is `"."`, it is taken to mean an
#' additive expression including all latent variables. For access, an object
#' containing a `bru_pred_expr` object.
#' @param \dots Passed on to submethods
#'
#' @keywords internal
#' @export
#' @rdname bru_pred_expr
#' @name bru_pred_expr
NULL

#' @describeIn bru_pred_expr Create a `bru_pred_expr` object from a formula or
#' character string.
#' @param used An optional `bru_used` object, overriding the automated
#'  detection. Default: `NULL`
#' @param is_rowwise logical; whether the predictor can be assumed to use the
#' input data rowwise, so that blockwise evaluation is possible. If `NULL`, it
#' needs to be set before use. Default: `NULL`
#' @param .envir The environment in which to evaluate the expression.
#' @export
#' @examples
#' (new_bru_pred_expr(~ x + z + Intercept))
#'
new_bru_pred_expr <- function(
  x,
  ...,
  used = NULL,
  is_rowwise = NULL,
  .envir = parent.frame()
) {
  if (is.character(x)) {
    x <- as.formula(glue::glue("~ {x}"), env = .envir)
  }
  pred_quo <- NULL
  resp_quo <- NULL
  if (is.null(x)) {
    is_additive <- NA
    is_additive_dot <- NA
  } else {
    is_additive <- bru_is_additive(x)
    if (rlang::is_quosure(x)) {
      pred_quo <- x
    } else if (inherits(x, "formula")) {
      .envir <- environment(x) %||% .envir
      pred_quo <- rlang::new_quosure(rlang::f_rhs(x), .envir)
      resp_expr <- rlang::f_lhs(x)
      if (!is.null(resp_expr)) {
        resp_quo <- rlang::new_quosure(resp_expr, .envir)
      }
    } else {
      pred_quo <- rlang::new_quosure(x, .envir)
    }
    is_additive_dot <- identical(rlang::get_expr(pred_quo), as.symbol("."))
    if (is_additive_dot) {
      pred_quo <- NULL
    }
  }
  structure(
    list(
      # TODO: Encode type of expression construction; full expression (expr),
      # additive (additive), list of effects and/or latent state vectors (list),
      # the latter also needing "hyper" information, so that generate.bru
      # can rely on this class for all use cases.
      type = "expr", # One of "expr", "additive", "list"
      pred_quo = pred_quo,
      is_additive = is_additive,
      # Also depends on component defs, so may be changed later:
      is_linear = is_additive,
      is_rowwise = is_rowwise,
      used = used %||% bru_used(pred_quo),
      .envir = .envir,
      resp_quo = resp_quo
    ),
    class = "bru_pred_expr"
  )
}

bru_compat_pre_2_14_bru_obs <- function(lh) {
  if (!isTRUE(bru_options_get("bru_compat_pre_2_14_enable"))) {
    return(lh)
  }
  lh$is_additive <- lh$pred_expr$is_additive
  # Also depends on component defs, so may be changed later:
  lh$linear <- lh$pred_expr$is_linear
  lh$expr_text <- lh$pred_expr$pred_text
  lh$expr <- lh$pred_expr$pred_expr
  lh$used <- lh$pred_expr$used
  lh$allow_combine <- !isTRUE(lh$pred_expr$is_rowwise)
  lh
}

expr_symbol_replace <- function(expr, sym, repl) {
  if (is.null(expr)) {
    return(NULL)
  }
  q <- rlang::is_quosure(expr)
  ex <- rlang::get_expr(expr)
  if (q) {
    env <- rlang::get_env(expr)
  }
  sym <- rlang::enexpr(sym)
  repl <- rlang::enexpr(repl)
  if (is.recursive(ex)) {
    for (k in seq_along(ex)) {
      ex[[k]] <- expr_symbol_replace(ex[[k]], !!sym, !!repl)
    }
  } else if (is.symbol(ex) && (ex == sym)) {
    ex <- repl
  }
  if (q) {
    return(rlang::new_quosure(ex, env))
  }
  ex
}

#' @describeIn bru_pred_expr Accessor generic for `bru_pred_expr` objects,
#'   including ones stored inside other objects.
#'
#' @keywords internal
#' @export
bru_pred_expr <- function(x, ...) {
  UseMethod("bru_pred_expr")
}


#' @describeIn bru_pred_expr Accessor for the `bru_pred_expr` object stored
#' inside a `bru_obs` object.
#' @export
bru_pred_expr.bru_obs <- function(x, ...) {
  bru_pred_expr(x[["pred_expr"]], ...)
}

#' @describeIn bru_pred_expr Accessor for the `bru_pred_expr` objects stored
#' inside a `bru_obs_list` object.
#' @export
bru_pred_expr.bru_obs_list <- function(x, ...) {
  lapply(x, function(y) bru_pred_expr(y[["pred_expr"]], ...))
}

#' @describeIn bru_pred_expr Accessor for the `bru_pred_expr` object stored
#' inside a `bru_info` object.
#' @export
bru_pred_expr.bru_info <- function(x, ...) {
  bru_pred_expr(as_bru_obs_list(x), ...)
}

#' @describeIn bru_pred_expr Accessor for the `bru_pred_expr` object stored
#' inside a `bru_model` object.
#' @export
bru_pred_expr.bru_model <- function(x, ...) {
  bru_pred_expr(as_bru_obs_list(x), ...)
}

#' @describeIn bru_pred_expr Accessor for the `bru_pred_expr` object stored
#' inside a `bru` object.
#' @export
bru_pred_expr.bru <- function(x, ...) {
  bru_pred_expr(as_bru_obs_list(x), ...)
}

#' @describeIn bru_pred_expr Access a `bru_pred_expr` object or convert it to
#' text or expression object.
#' @param format Character; one of
#' * `"object"` (the default),
#' * `"text"` (plain text),
#' * `"quo"` (an `rlang` quosure),
#' * `"expr"` (an `rlang` expression),
#' * `"formula"` (the original formula input if available, otherwise
#'   constructed from the expression),
#' * `"formula_text"` (a text version of the "formula"),
#' * `"resp_text"` (plain text of the response side of the formula).
#' @param raw Logical; whether to return the raw expression without
#'   substituting `BRU_EXPRESSION` with the used latent variables.
#'   Default: `FALSE`
#' @export
bru_pred_expr.bru_pred_expr <- function(
  x,
  ...,
  format = "object",
  raw = FALSE
) {
  format <- match.arg(
    format,
    c(
      "object",
      "text",
      "quo",
      "expr",
      "formula",
      "formula_text",
      "resp_quo",
      "resp_expr",
      "resp_text"
    )
  )
  if (identical(format, "object")) {
    return(x)
  }

  pred_quo <- x[["pred_quo"]]
  if (is.null(pred_quo)) {
    if (raw) {
      pred_quo <- rlang::new_quosure(as.symbol("BRU_EXPRESSION"), x[[".envir"]])
    } else {
      included <- bru_used(x)[["effect"]]
      if (length(included) == 0L) {
        pred_quo <- rlang::parse_quo(
          ".",
          env = x[[".envir"]]
        )
      } else {
        pred_quo <- rlang::parse_quo(
          paste0(included, collapse = " + "),
          env = x[[".envir"]]
        )
      }
    }
  } else if ((!raw) && ("BRU_EXPRESSION" %in% bru_used_vars(pred_quo)$vars)) {
    included <- bru_used(x)[["effect"]]
    pred_quo <- expr_symbol_replace(
      pred_quo,
      sym = !!as.symbol("BRU_EXPRESSION"),
      repl = !!rlang::parse_expr(paste0(included, collapse = " + "))
    )
  }

  switch(
    format,
    quo = pred_quo,
    expr = rlang::get_expr(pred_quo),
    text = deparse1(rlang::get_expr(pred_quo), collapse = "\n"),
    formula = rlang::new_formula(
      lhs = rlang::get_expr(x[["resp_quo"]]),
      rhs = rlang::get_expr(pred_quo),
      env = x[[".envir"]]
    ),
    formula_text = deparse1(
      bru_pred_expr(x, format = "formula", raw = raw),
      collapse = "\n"
    ),
    resp_quo = x[["resp_quo"]],
    resp_expr = rlang::get_expr(x[["resp_quo"]]),
    resp_text = deparse1(rlang::get_expr(x[["resp_quo"]]), collapse = "\n"),
    stop(glue::glue("Unknown format '{format}'"))
  )
}


#' @rdname bru_pred_expr
#' @export
#' @importFrom glue glue
format.bru_pred_expr <- function(x, ...) {
  glue::glue(
    "    Predictor: {predictor}\n",
    "    Additive/Linear/Rowwise: {is_additive}/{is_linear}/{is_rowwise}\n",
    "    Used components: {format(used)}",
    predictor = bru_pred_expr(x, format = "formula_text"),
    is_additive = bru_is_additive(x),
    is_linear = bru_is_linear(x),
    is_rowwise = bru_is_rowwise(x),
    used = bru_used(x)
  )
}

#' @rdname bru_pred_expr
#' @export
#' @importFrom glue glue
print.bru_pred_expr <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
