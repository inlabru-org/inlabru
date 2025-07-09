#' @title Extract predictor or component index information
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' Extract the index vector for a [bru_obs()] predictor,
#' or the whole or a subset of a full [bru()] predictor.
#'
#' @param object A [bru()] or [bru_obs()] output object
#' @param \dots Arguments passed on to sub-methods.
#' @export
bru_index <- function(object, ...) {
  UseMethod("bru_index")
}

#' @describeIn bru_index Extract the index vector for the predictor vector for a
#'   [bru_obs()] sub-model. The indices are relative to the sub-model, and need
#'   to be appropriately offset to be used in the full model predictor.
#' @param what `character` or `NULL`; One of `NULL`, "all", "observed", and
#'   "missing". If `NULL` (default) or "all", gives the index vector for the
#'   full sub-model predictor. If "observed", gives the index vector for the
#'   observed part (response is not `NA`). If "missing", gives the index vector
#'   for the missing part (response is `NA`) of the model.
#' @examplesIf bru_safe_inla()
#' fit <- bru(
#'   ~ 0 + x,
#'   bru_obs(
#'     y ~ .,
#'     data = data.frame(x = 1:3, y = 1:3 + rnorm(3)),
#'     tag = "A"
#'   ),
#'   bru_obs(
#'     y ~ .,
#'     data = data.frame(x = 1:4, y = c(NA, NA, 3:4) + rnorm(4)),
#'     tag = "B"
#'   )
#' )
#' bru_index(fit)
#' bru_index(fit, "A")
#' bru_index(fit, "B")
#' bru_index(fit, c("B", "A"))
#' bru_index(fit, what = "missing")
#' @export
#' @returns * `bru_index(bru_obs)`: An `integer` vector.
bru_index.bru_obs <- function(object, what = NULL, ...) {
  what <- match.arg(what, c("all", "observed", "missing"))
  size <- bru_response_size(object)
  idx <- seq_len(size)
  if (identical(what, "all")) {
    return(idx)
  }
  resp <- object[["response_data"]][[object[["response"]]]]
  if (is.vector(resp)) {
    miss <- is.na(resp)
  } else {
    # Should work for both inla.surv and inla.mdata objects:
    miss <- is.na(resp[[1]])
  }
  if (identical(what, "observed")) {
    return(idx[!miss])
  }
  return(idx[miss])
}

#' @describeIn bru_index Extract the index vector for "APredictor" for one or
#'   more specified observation [bru_obs()] sub-models. Accepts any combination
#'   of `tag` and `what`.
#' @param tag `character` or `integer`; Either a character vector identifying
#'   the tags of one or more of the [bru_obs()] observation models, or an
#'   integer vector identifying models by their [bru()] specification order. If
#'   `NULL` (default) computes indices for all sub-models.
#' @export
#' @returns * `bru_index(bru)`: An `integer` vector.
bru_index.bru <- function(object, tag = NULL, what = NULL, ...) {
  if (is.null(tag)) {
    tag <- seq_len(length(object[["bru_info"]][["lhoods"]]))
  }
  if (length(tag) == 0L) {
    return(integer(0))
  }
  size <- bru_response_size(object)
  off <- c(0L, cumsum(size))
  if (is.character(tag)) {
    if (!all(tag %in% names(object[["bru_info"]][["lhoods"]]))) {
      stop(paste0(
        "Invalid tag(s) '", paste(tag, collapse = ", "), "' for ",
        "bru object with tags '", paste(names(object[["bru_info"]][["lhoods"]]),
          collapse = ", "
        ), "'"
      ))
    }
    unlist(lapply(tag, function(x) {
      x_idx <- which(names(object$bru_info$lhoods) %in% x)
      off[x_idx] + bru_index(object$bru_info$lhoods[[x_idx]], what = what)
    }))
  } else {
    ok_tags <- (tag >= 1L) &
      (tag <= length(object[["bru_info"]][["lhoods"]]))
    if (any(!ok_tags)) {
      stop(paste0(
        "Invalid tag indices '", paste(tag[!ok_tags], collapse = ", "),
        "' for bru object with tag indices '1, ..., ",
        length(object[["bru_info"]][["lhoods"]]), "'"
      ))
    }
    unlist(lapply(tag, function(x) {
      off[x] + bru_index(object$bru_info$lhoods[[x]], what = what)
    }))
  }
}

#' @export
#' @keywords internal
#' @param object A [component].
#' @param inla_f logical; when `TRUE`, must result in
#' values compatible with `INLA::f(...)`
#' an specification and corresponding `INLA::inla.stack(...)` constructions.
#' @returns * `bru_index(bru_comp)`:
#'   A list of indices into the latent variables compatible with the
#'   component mapper.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com},
#'   Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname bru_index

bru_index.bru_comp <- function(object, inla_f, ...) {
  idx <- ibm_values(object[["mapper"]], inla_f = inla_f, multi = TRUE)
  names(idx) <- paste0(object[["label"]], c("", ".group", ".repl"))
  idx
}


#' @export
#' @rdname bru_index
#' @returns * `bru_index(bru_comp_list)`:
#'   A list of list of indices into the latent variables compatible with each
#'   component mapper.
bru_index.bru_comp_list <- function(object, inla_f, ...) {
  lapply(object, function(x) bru_index(x, inla_f = inla_f, ...))
}


#' @describeIn bru_index Compute all index values for a [bru_model()] object.
#' Computes the index value matrices for included components according to the
#' `used` argument.
#'
#' @param used A [bru_used()] object
#' @return * `bru_index(bru_model)`: A named list of `idx_full` and `idx_inla`,
#' named list of indices, and `inla_subset`, and `inla_subset`,
#' a named list of logical subset specifications for extracting the `INLA::f()`
#' compatible index subsets.
#' @export
bru_index.bru_model <- function(object, used, ...) {
  stopifnot(inherits(object, "bru_model"))
  included <- union(used[["effect"]], used[["latent"]])

  list(
    idx_full = bru_index(object[["effects"]][included], inla_f = FALSE),
    idx_inla = bru_index(object[["effects"]][included], inla_f = TRUE),
    inla_subset = inla_subset_eval(object[["effects"]][included])
  )
}


#' @export
#' @describeIn bru_index `r lifecycle::badge("deprecated")` since "2.12.0.9023".
#' Use `bru_index()` instead.
index_eval <- function(...) {
  lifecycle::deprecate_warn(
    when = "2.12.0.9023",
    "index_eval()",
    "bru_index()"
  )
  bru_index(...)
}

#' Obtain inla index subset information
#'
#' Subsets for `INLA::f()` compatible indexing
#'
#' @export
#' @param ... Passed on to submethods.
#' @rdname inla_subset_eval
inla_subset_eval <- function(...) {
  UseMethod("inla_subset_eval")
}

#' @export
#' @keywords internal
#' @param components A component list.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @rdname inla_subset_eval

inla_subset_eval.bru_comp_list <- function(components, ...) {
  lapply(components, function(x) ibm_inla_subset(x[["mapper"]], multi = TRUE))
}
