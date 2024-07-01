#' Build an inla data stack from linearisation information
#'
#' Combine linearisation for multiple likelihoods
#'
#' @param \dots Arguments passed on to other methods
#' @export
#' @keywords internal
#' @rdname bru_make_stack
bru_make_stack <- function(...) {
  UseMethod("bru_make_stack")
}

#' @param lhood A `bru_like` object
#' @param lin Linearisation information
#' * For `.bru_like`, a `bru_mapper_taylor` object
#' * For `.bru_like_list`, a list of `bru_mapper_taylor` objects
#' @param idx Output from `evaluate_index(...)`
#' @param family_index integer specifying the family sequence index of the observation model
#' @export
#' @rdname bru_make_stack
bru_make_stack.bru_like <- function(lhood, lin, idx, ..., family_index = 1L) {
  stopifnot(inherits(lin, "bru_mapper_taylor"))
  stopifnot(!is.null(lin[["offset"]]))
  stopifnot(is.null(lin[["jacobian"]]) || is.list(lin[["jacobian"]]))
  stopifnot(is.null(lin[["state0"]]))

  if (is.null(lhood[["response_data"]])) {
    BRU.response <- lhood$data[[lhood$response]]
  } else {
    BRU.response <- lhood$response_data[[lhood$response]]
  }
  nms <- names(lin$jacobian)
  if (utils::packageVersion("INLA") <= "24.06.02") {
    INLA::inla.stack(
      list(
        BRU.response = BRU.response,
        BRU.E = lhood[["E"]],
        BRU.Ntrials = lhood[["Ntrials"]],
        BRU.weights = lhood[["weights"]],
        BRU.scale = lhood[["scale"]],
        BRU.offset = as.vector(lin$offset),
        BRU.link = family_index
      ),
      A = lapply(nms, function(nm) {
        lin$jacobian[[nm]][, idx[["inla_subset"]][[nm]], drop = FALSE]
      }),
      effects = idx[["idx_inla"]][nms],
      remove.unused = FALSE
    )
  } else {
    INLA::inla.stack(
      list(
        BRU.E = lhood[["E"]],
        BRU.Ntrials = lhood[["Ntrials"]],
        BRU.weights = lhood[["weights"]],
        BRU.scale = lhood[["scale"]],
        BRU.offset = as.vector(lin$offset),
        BRU.link = family_index
      ),
      A = lapply(nms, function(nm) {
        lin$jacobian[[nm]][, idx[["inla_subset"]][[nm]], drop = FALSE]
      }),
      effects = idx[["idx_inla"]][nms],
      responses = list(BRU.response),
      remove.unused = FALSE
    )
  }
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
          idx = idx,
          family_index = lh_idx
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
      bru_inla.stack.mjoin,
      c(stks, list(compress = TRUE, remove.unused = FALSE))
    )

  stk
}
