# Create a stack from a model, data and prediction values
#
# # TODO: make sure data is allowed to be a list with unequal length vectors

make.stack <- function(data,
                       model,
                       y,
                       E = 1,
                       Ntrials = 1,
                       offset = 0,
                       expr = NULL,
                       state = NULL,
                       tag = "BRU.stack",
                       include = NULL,
                       exclude = NULL) {
  included <- parse_inclusion(names(model$effects),
    include = include,
    exclude = exclude
  )
  offsets <- names(model$effects)[vapply(
    model$effects,
    function(x) identical(x$type, "offset"),
    TRUE
  )]
  included <- setdiff(included, offsets)

  # Observations y
  if (length(y) == 1) {
    y <- rep(y, nrow(as.data.frame(data)))
  }

  # Expectation parameter E
  if (length(E) == 1) {
    E <- rep(E, nrow(as.data.frame(data)))
  }

  # Ntrials
  if (length(Ntrials) == 1) {
    Ntrials <- rep(Ntrials, nrow(as.data.frame(data)))
  }


  # Projection matrices (A) and inla-compatible indices
  A <- amatrix_eval(model$effects[included], data)
  effects <- index_eval(model$effects[included])


  # Taylor approximation
  if (!is.null(expr)) {
    rw <- nlinla.reweight(A, model, data, expr, state = state)
    A <- rw$A
    taylor.offset <- rw$const
  } else {
    taylor.offset <- 0
  }

  # Create and return stack
  stk <- INLA::inla.stack(
    data = list(
      BRU.response = y,
      BRU.E = E,
      BRU.Ntrials = Ntrials,
      BRU.offset = taylor.offset + offset
    ),
    A = A,
    tag = tag,
    effects = effects,
    # Make sure latent components with zero-derivatives aren't removed:
    remove.unused = FALSE
  )
  stk
}




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
#' @param idx Output from [evaluate_index()]
#' @export
#' @rdname bru_make_stack
bru_make_stack.bru_like <- function(lhood, lin, idx, ...) {
  INLA::inla.stack(
    list(
      BRU.response = lhood$data[[lhood$response]],
      BRU.E = lhood[["E"]],
      BRU.Ntrials = lhood[["Ntrials"]],
      BRU.offset = as.vector(lin$offset)
    ),
    A = lin$A,
    effects = idx[names(lin$A)],
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
          lhoods[[lh_idx]],
          lin = lin[[lh_idx]],
          idx
        )
      }
    )

  stk <-
    do.call(
      inlabru::inla.stack.mjoin,
      c(stks, list(compress = TRUE, remove.unused = FALSE))
    )

  stk
}
