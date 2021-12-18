#' Build an inla data stack from linearisation information
#'
#' Combine linearisation for multiple likelihoods
#'
#' @param \dots Arguments passed on to other methods
#' @export
#' @rdname bru_make_stack
bru_make_stack <- function(...) {
  UseMethod("bru_make_stack")
}

#' @param lhood A `bru_like` object
#' @param lin Linearisation information
#' * For `.bru_like`, a linearisation information list with elements
#' `A` and `offset`
#' * For `.bru_like_list`, a list of linearisation information lists
#' @param idx Output from `evaluate_index(...)`
#' @export
#' @rdname bru_make_stack
bru_make_stack.bru_like <- function(lhood, lin, idx, ...) {
  if (is.null(lhood[["response_data"]])) {
    BRU.response <- lhood$data[[lhood$response]]
  } else {
    BRU.response <- lhood$response_data[[lhood$response]]
  }
  INLA::inla.stack(
    list(
      BRU.response = BRU.response,
      BRU.E = lhood[["E"]],
      BRU.Ntrials = lhood[["Ntrials"]],
      BRU.offset = as.vector(lin$offset)
    ),
    A = lapply(names(lin$A), function(nm) {
      lin$A[[nm]][, idx[["inla_subset"]][[nm]], drop = FALSE]
    }),
    effects = idx[["idx_inla"]][names(lin$A)],
    remove.unused = FALSE
  )
}

#' @param lhoods A `bru_like_list` object
#' @export
#' @rdname bru_make_stack
bru_make_stack.bru_like_list <- function(lhoods, lin, idx, ...) {
  stks <-
    lapply(
      seq_along(lhoods),
      function(lh_idx) {
        bru_make_stack(
          lhood = lhoods[[lh_idx]],
          lin = lin[[lh_idx]],
          idx = idx
        )
      }
    )

  # TODO: Either fix the check in inla that removes n&values information from
  # model="linear" input in INLA::f(), making it unable to handle all-NA linear
  # effect specifications (arguably, that check is wrong, since it doesn't need
  # that information for linear components), or remove all-but one unused value
  # above or in bru_make_stack.bru_like, and keep using remove.unused=FALSE here.
  stk <-
    do.call(
      inlabru::inla.stack.mjoin,
      c(stks, list(compress = TRUE, remove.unused = FALSE))
    )

  stk
}
