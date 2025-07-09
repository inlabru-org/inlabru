#' @include deprecated.R
#' @include mappers.R

## _repeat ####

#' @title Mapper for repeating a mapper
#' @param mapper The mapper to be repeated.
#' @param n_rep The number of times to repeat the mapper. If a vector, the
#'   non-interleaved repeats are combined into a single repeat mapping,
#'   and combined with interleaved repeats via a [bm_sum()] of mappers.
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
#' @returns A `bm_repeat` or `bm_sum` object, or the original
#' input `mapper`.
#' @rdname bm_repeat
#' @inheritParams bru_mapper_generics
#' @inheritParams bm_scale
#' @inheritParams bm_multi
#' @seealso [bru_mapper], [bru_mapper_generics]
#' @family mappers
#' @examples
#' (m0 <- bm_index(3))
#' (m <- bm_repeat(m0, 5))
#' ibm_n(m)
#' ibm_values(m)
#' ibm_jacobian(m, 1:3)
#' ibm_eval(m, 1:3, seq_len(ibm_n(m)))
#'
#' # Interleaving and grouping
#' (m <- bm_repeat(m0, c(2, 1, 2), c(TRUE, FALSE, FALSE)))
#' ibm_n(m)
#' ibm_values(m)
#' ibm_jacobian(m, 1:3)
#' ibm_eval(m, 1:3, seq_len(ibm_n(m)))
#'
bm_repeat <- function(mapper, n_rep, interleaved = FALSE) {
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
    return(bru_mapper_define(mapper, new_class = "bm_repeat"))
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
        mappers[[m_idx]] <- bm_repeat(
          mapper,
          n_rep = n_rep_cumulative,
          interleaved = FALSE
        )
        n_rep_cumulative <- 0L
      }
      m_idx <- m_idx + 1L
      mappers[[m_idx]] <- bm_repeat(
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
    mappers[[m_idx]] <- bm_repeat(
      mapper,
      n_rep = n_rep_cumulative,
      interleaved = FALSE
    )
    n_rep_cumulative <- 0L
  }

  if (m_idx == 1L) {
    return(mappers[[1]])
  }
  mapper_ <- bm_sum(mappers, single_input = TRUE)
  return(mapper_)
}

#' @export
#' @rdname bm_repeat
bru_mapper_repeat <- function(...) {
  bm_repeat(...)
}

#' @export
#' @rdname bm_repeat
ibm_n.bm_repeat <- function(mapper, ...) {
  ibm_n(mapper[["mapper"]], ...) * mapper[["n_rep"]]
}
#' @export
#' @rdname bm_repeat
ibm_n_output.bm_repeat <- function(mapper, ...) {
  ibm_n_output(mapper[["mapper"]], ...)
}

#' @export
#' @rdname bm_repeat
ibm_values.bm_repeat <- function(mapper, ...) {
  seq_len(ibm_n(mapper, ...))
}



bm_repeat_sub_lin <- function(mapper, input, state,
                              ...) {
  # We need all the sub_lin objects even for linear mappers
  idx <- bm_repeat_indexing(
    n_map = ibm_n(mapper[["mapper"]]),
    n_rep = mapper[["n_rep"]],
    interleaved = mapper[["interleaved"]]
  )
  sub_lin <-
    lapply(
      seq_len(mapper[["n_rep"]]),
      function(x) {
        state_subset <- state[idx$offsets[x] + idx$index]
        ibm_linear(
          mapper[["mapper"]],
          input = input,
          state = state_subset
        )
      }
    )
  sub_lin
}



#' @describeIn bm_repeat The input should take the format of the
#'   repeated submapper.
#' @export
ibm_jacobian.bm_repeat <- function(mapper, input, state = NULL,
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
    A <- A %*% bm_repeat_indexing_matrix(
      n_map = ibm_n(mapper[["mapper"]]),
      n_rep = mapper[["n_rep"]]
    )
  }
  return(A)
}


#' @export
#' @rdname bm_repeat
ibm_eval.bm_repeat <- function(mapper, input, state,
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
#' @rdname bm_repeat
ibm_linear.bm_repeat <- function(mapper, input, state,
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
  bm_taylor(
    offset = eval2$offset,
    jacobian = eval2$jacobian,
    state0 = state,
    values_mapper = mapper
  )
}




#' @describeIn bm_repeat
#' Passes on the input to the corresponding method.
#' @export
ibm_invalid_output.bm_repeat <- function(mapper, input, state,
                                         ...) {
  idx <- bm_repeat_indexing(
    n_map = ibm_n(mapper[["mapper"]]),
    n_rep = mapper[["n_rep"]],
    interleaved = mapper[["interleaved"]]
  )
  ibm_invalid_output(
    mapper[["mapper"]],
    input = input,
    state = state[idx$offsets[1] + idx$index]
  )
}

# Interleaving ####

#' @title Re-indexing matrix for repeated mappers
#' @name bm_repeat_indexing
#' @rdname bm_repeat_indexing
#' @description Helper methods for regular and interleaved state vector
#'   indexing, meant for use in [bm_repeat] mappers.
#' @keywords internal
NULL

#' @describeIn bm_repeat_indexing Construct the offsets and within-block
#' index vectors for a repeated mapper, with or without interleaving.
#' @param n_map The block mapper size
#' @param n_rep The number of times the block is repeated
#' @param interleaved logical; if `TRUE`, the state vector indexing is
#'   interleaved. Default is `FALSE`, for blockwise indexing.
#' @returns `bm_repeat_indexing`: Returns a list with a vector of offsets,
#'   `offsets`, and a vector of relative index values, `index`; block `k`
#'   indexing is given by `offsets[k] + index`.
#' @export
#' @examples
#' (idx <- bm_repeat_indexing(3, 2, FALSE))
#' (idx <- bm_repeat_indexing(3, 2, TRUE))
bm_repeat_indexing <- function(n_map,
                               n_rep,
                               interleaved = FALSE) {
  if (isTRUE(interleaved)) {
    list(
      offsets = seq_len(n_rep) - 1L,
      index = n_rep * (seq_len(n_map) - 1L) + 1L
    )
  } else {
    list(
      offsets = n_map * (seq_len(n_rep) - 1L),
      index = seq_len(n_map)
    )
  }
}
#' @describeIn bm_repeat_indexing Creates a sparse matrix `A` such that
#'   `z <- A %*% x` constructs a blockwise version `z = c(x1, x2, ..., xn)` of
#'   an interleaved state vector `x = c(x1[1], x2[1], ..., x1[2], x2[2], ...)`.
#'   Each block is of size `n_map`, and there are `n_rep` blocks. The reverse
#'   operation, taking a blockwise `z = c(x1, x2, ..., xn)` to an interleaved
#'   vector is `x <- Matrix::t(A) %*% z`.
#' @returns `bm_repeat_indexing_matrix`: A `sparseMatrix` object.
#' @export
#' @examples
#' (A <- bm_repeat_indexing_matrix(3, 2))
#' (x_interleaved <- 1:6)
#' (x_blockwise <- as.vector(A %*% x_interleaved))
#' (x_recovered <- as.vector(Matrix::t(A) %*% x_blockwise))
#'
bm_repeat_indexing_matrix <- function(n_map, n_rep) {
  Matrix::sparseMatrix(
    i = rep((seq_len(n_rep) - 1L) * n_map, times = n_map) +
      rep(seq_len(n_map), each = n_rep),
    j = seq_len(n_map * n_rep),
    x = 1,
    dims = c(1L, 1L) * n_map * n_rep
  )
}
