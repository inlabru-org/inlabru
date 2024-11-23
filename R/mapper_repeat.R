#' @include deprecated.R
#' @include mappers.R

## _repeat ####

#' @title Mapper for repeating a mapper
#' @param mapper The mapper to be repeated.
#' @param n_rep The number of times to repeat the mapper. If a vector, the
#'   non-interleaved repeats are combined into a single repeat mapping,
#'   and combined with interleaved repeats via a [bru_mapper_sum()] of mappers.
#' @param interleaved logical; if `TRUE`, the repeated mapping columns are
#' interleaved; `(x1[1], x2[1], ..., x1[2], x2[2], ...)`.
#' If `FALSE` (default), the repeated
#' mapping columns are contiguous, `(x1[1], x1[2], ..., x2[1], x2[2], ...)`,
#' and the Jacobian is a `cbind()` of the Jacobians of the repeated mappers.
#'
#' If `n_rep` is a vector, `interleaved` should either be a single logical,
#' or a vector of the same length. Each element applies to the corresponding
#' `n_rep` repetition specification.
#' @export
#' @description Defines a repeated-space mapper that sums the contributions for
#'   each copy. The `ibm_n()` method returns `ibm_n(mapper) * n_rep`, and
#'   `ibm_values()` returns `seq_len(ibm_n(mapper))`.
#' @returns A `bru_mapper_repeat` or `bru_mapper_sum` object, or the original
#' input `mapper`.
#' @rdname bru_mapper_repeat
#' @inheritParams bru_mapper_generics
#' @inheritParams bru_mapper_scale
#' @inheritParams bru_mapper_multi
#' @seealso [bru_mapper], [bru_mapper_generics]
#' @family mappers
#' @examples
#' (m0 <- bru_mapper_index(3))
#' (m <- bru_mapper_repeat(m0, 5))
#' ibm_n(m)
#' ibm_values(m)
#' ibm_jacobian(m, 1:3)
#' ibm_eval(m, 1:3, seq_len(ibm_n(m)))
#'
#' # Interleaving and grouping
#' (m <- bru_mapper_repeat(m0 , c(2, 1, 2), c(TRUE, FALSE, FALSE)))
#' ibm_n(m)
#' ibm_values(m)
#' ibm_jacobian(m, 1:3)
#' ibm_eval(m, 1:3, seq_len(ibm_n(m)))
#'
bru_mapper_repeat <- function(mapper, n_rep, interleaved = FALSE) {
  stopifnot((length(n_rep) > 0L) && sum(n_rep) > 0L)
  if ((length(n_rep) == 1L) || !any(interleaved)) {
    n_rep <- sum(n_rep)
    if (n_rep == 1L) {
      return(mapper)
    }
    mapper <- list(
      mapper = mapper,
      n_rep = n_rep,
      is_linear = ibm_is_linear(mapper),
      interleaved = interleaved
    )
    return(bru_mapper_define(mapper, new_class = "bru_mapper_repeat"))
  }

  # Combine non-interleaved repeats
  if (length(interleaved) == 1L) {
    interleaved <- rep(interleaved, length(n_rep))
  }
  mappers <- list()
  m_idx <- 0L
  n_rep_cumulative <- 0L
  for (k in seq_along(n_rep)) {
    if (interleaved[k]) {
      if (n_rep_cumulative > 0L) {
        m_idx <- m_idx + 1L
        mappers[[m_idx]] <- bru_mapper_repeat(
          mapper,
          n_rep = n_rep_cumulative,
          interleaved = FALSE
        )
        n_rep_cumulative <- 0L
      }
      m_idx <- m_idx + 1L
      mappers[[m_idx]] <- bru_mapper_repeat(
        mapper,
        n_rep = n_rep[k],
        interleaved = TRUE
      )
    } else {
      n_rep_cumulative <- n_rep_cumulative + n_rep[k]
    }
  }
  if (n_rep_cumulative > 0L) {
    m_idx <- m_idx + 1L
    mappers[[m_idx]] <- bru_mapper_repeat(
      mapper,
      n_rep = n_rep_cumulative,
      interleaved = FALSE
    )
    n_rep_cumulative <- 0L
  }

  if (m_idx == 1L) {
    return(mappers[[1]])
  }
  mapper_ <- bru_mapper_sum(mappers, single_input = TRUE)
  return(mapper_)
}

#' @export
#' @rdname bru_mapper_repeat
ibm_n.bru_mapper_repeat <- function(mapper, ...) {
  ibm_n(mapper[["mapper"]], ...) * mapper[["n_rep"]]
}
#' @export
#' @rdname bru_mapper_repeat
ibm_n_output.bru_mapper_repeat <- function(mapper, ...) {
  ibm_n_output(mapper[["mapper"]], ...)
}

#' @export
#' @rdname bru_mapper_repeat
ibm_values.bru_mapper_repeat <- function(mapper, ...) {
  seq_len(ibm_n(mapper, ...))
}



bm_repeat_sub_lin <- function(mapper, input, state,
                              ...) {
  # We need all the sub_lin objects even for linear mappers
  n <- ibm_n(mapper[["mapper"]])
  if (mapper[["interleaved"]]) {
    n_offset <- seq_len(n) - 1L
  } else {
    n_offset <- n * (seq_len(mapper[["n_rep"]]) - 1L)
  }
  sub_lin <-
    lapply(
      seq_len(mapper[["n_rep"]]),
      function(x) {
        if (mapper[["interleaved"]]) {
          idx <- (seq_len(n) - 1L) * mapper[["n_rep"]] + 1L
        } else {
          idx <- seq_len(n)
        }
        state_subset <- state[n_offset[x] + idx]
        ibm_linear(
          mapper[["mapper"]],
          input = input,
          state = state_subset
        )
      }
    )
  sub_lin
}



#' @describeIn bru_mapper_repeat The input should take the format of the
#'   repeated submapper.
#' @export
ibm_jacobian.bru_mapper_repeat <- function(mapper, input, state = NULL,
                                           inla_f = FALSE,
                                           multi = FALSE,
                                           ...,
                                           sub_lin = NULL) {
  if (is.null(sub_lin)) {
    sub_lin <- bm_repeat_sub_lin(mapper, input, state)
  }
  A <- lapply(sub_lin, function(x) x[["jacobian"]])

  if (multi) {
    return(A)
  }

  # Combine the matrices (A1, A2, A3, ...) -> permuted cbind(A1, A2, A3, ...)
  A <- do.call(cbind, A)
  if (isTRUE(mapper[["interleaved"]])) {
    n_map <- ibm_n(mapper[["mapper"]])
    n_rep <- mapper[["n_rep"]]
    A <- A %*%
      Matrix::sparseMatrix(
        i = rep((seq_len(n_rep) - 1L) * n_map, times = n_map) +
          rep(seq_len(n_map), each = n_rep),
        j = seq_len(n_map * n_rep),
        x = 1,
        dims = c(1L, 1L) * n_map * n_rep
      )
  }
  return(A)
}


#' @export
#' @rdname bru_mapper_repeat
ibm_eval.bru_mapper_repeat <- function(mapper, input, state,
                                       multi = FALSE,
                                       ...,
                                       sub_lin = NULL) {
  if (is.null(sub_lin)) {
    sub_lin <- bm_repeat_sub_lin(mapper, input, state)
  }
  val <- lapply(sub_lin, function(x) x[["offset"]])
  if (multi) {
    return(val)
  }

  # Combine the vectors (b1, b2, b3, ...) -> b1 + b2 + b3 + ...
  val <- as.vector(Matrix::rowSums(do.call(cbind, val)))
  val
}


#' @export
#' @rdname bru_mapper_repeat
ibm_linear.bru_mapper_repeat <- function(mapper, input, state,
                                         ...) {
  sub_lin <-
    bm_repeat_sub_lin(
      mapper, input, state,
      ...
    )
  eval2 <- ibm_eval2(
    mapper,
    input = input,
    state = state,
    multi = FALSE,
    ...,
    sub_lin = sub_lin
  )
  bru_mapper_taylor(
    offset = eval2$offset,
    jacobian = eval2$jacobian,
    state0 = state,
    values_mapper = mapper
  )
}




#' @describeIn bru_mapper_repeat
#' Passes on the input to the corresponding method.
#' @export
ibm_invalid_output.bru_mapper_repeat <- function(mapper, input, state,
                                                 ...) {
  return(
    ibm_invalid_output(
      mapper[["mapper"]],
      input = input,
      state = state[seq_len(ibm_n(mapper[["mapper"]]))]
    )
  )
}
