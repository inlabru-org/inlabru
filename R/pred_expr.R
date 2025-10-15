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
bru_pred_expr <- function(x, ...) {
  UseMethod("bru_pred_expr")
}

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
#' (bru_pred_expr(~ x + z + Intercept))
#'
bru_pred_expr.default <- function(x, ..., used = NULL, is_rowwise = NULL,
                                  .envir = parent.frame()) {
  resp_text <- NULL
  if (inherits(x, "formula")) {
    formula_text <- deparse1(x, collapse = "\n")
    x <- as.character(x)
    pred_text <- x[length(x)]
    if (length(x) > 2) {
      resp_text <- x[2]
    }
  } else {
    pred_text <- x
    formula_text <- glue::glue("~ {pred_text}")
  }
  is_additive <- bru_is_additive(pred_text)
  is_additive_dot <- identical(pred_text, ".")
  if (!is_additive_dot) {
    pred_expr <- parse(text = pred_text)
  } else {
    pred_text <- NULL
    pred_expr <- NULL
  }
  if (is.null(resp_text)) {
    resp_expr <- NULL
  } else {
    resp_expr <- parse(text = resp_text)
  }
  structure(
    list(
      formula_text = formula_text,
      pred_text = pred_text,
      pred_expr = pred_expr,
      is_additive = is_additive,
      # Also depends on component defs, so may be changed later:
      is_linear = is_additive,
      is_rowwise = is_rowwise,
      used = if (is.null(used)) bru_used(pred_text) else used,
      .envir = .envir,
      resp_text = resp_text,
      resp_expr = resp_expr
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



#' @describeIn bru_pred_expr Accessor for the `bru_pred_expr` object stored
#' inside a `bru_obs` object.
#' @export
bru_pred_expr.bru_obs <- function(x, ...) {
  bru_pred_expr(x[["pred_expr"]], ...)
}

#' @describeIn bru_pred_expr Access a `bru_pred_expr` object or convert it to
#' text or expression object.
#' @param format Character; one of
#' * `"object"` (the default),
#' * `"text"` (plain text),
#' * `"expression"` (a `base` `expression`),
#' * `"quo"` (an `rlang` quosure),
#' * `"expr"` (an `rlang` expression),
#' * `"formula"` (the original formula input if available, otherwise
#'   constructed from the expression),
#' * `"formula_text"` (a text version of the "formula"), or
#' * `"text_raw"` (plain text, without substituting `BRU_EXPRESSION`).
#' @export
bru_pred_expr.bru_pred_expr <- function(x, ..., format = "object") {
  format <- match.arg(
    format,
    c(
      "object",
      "text",
      "expression",
      "quo",
      "expr",
      "formula",
      "formula_text",
      "text_raw"
    )
  )
  if (identical(format, "object")) {
    return(x)
  }

  pred_text <- x[["pred_text"]]
  if (is.null(pred_text)) {
    pred_text <- "BRU_EXPRESSION"
  }
  if (!identical(format, "text_raw") &&
    grepl(
      pattern = "BRU_EXPRESSION",
      x = pred_text,
      fixed = TRUE
    )) {
    included <- bru_used(x)[["effect"]]
    pred_text <-
      gsub(
        pattern = "BRU_EXPRESSION",
        replacement = paste0(included, collapse = " + "),
        x = pred_text,
        fixed = TRUE
      )
  }

  switch(format,
    text = pred_text,
    expression = parse(text = pred_text),
    quo = rlang::parse_quo(pred_text, env = x[[".envir"]]),
    expr = rlang::parse_expr(pred_text),
    formula = if (is.null(x[["formula_text"]])) {
      as.formula(glue::glue("~ {pred_text}"),
        env = x[[".envir"]]
      )
    } else {
      as.formula(x[["formula_text"]], env = x[[".envir"]])
    },
    formula_text = if (is.null(x[["formula_text"]])) {
      glue::glue("~ {pred_text}")
    } else {
      x[["formula_text"]]
    },
    text_raw = pred_text,
    stop(glue::glue("Unknown format '{format}'"))
  )
}
