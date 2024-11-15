#' @include deprecated.R
#' @include mappers.R

## _repeat ####

#' @title Mapper for repeating a mapper
#' @param mapper The mapper to be repeated.
#' @param n_rep The number of times to repeat the mapper.
#' @export
#' @description Defines a repeated-space mapper that sums the contributions for
#'   each copy. The `ibm_n()` method returns `ibm_n(mapper) * n_rep`, and
#'   `ibm_values()` returns `seq_len(ibm_n(mapper))`.
#' @returns A `bru_mapper_repeat` object.
#' @rdname bru_mapper_repeat
#' @inheritParams bru_mapper_generics
#' @inheritParams bru_mapper_scale
#' @inheritParams bru_mapper_multi
#' @seealso [bru_mapper], [bru_mapper_generics]
#' @family mappers
#' @examples
#' m <- bru_mapper_repeat(bru_mapper_index(3), 4)
#' ibm_n(m)
#' ibm_values(m)
#' ibm_jacobian(m, 1:3)
#' ibm_eval(m, 1:3, seq_len(ibm_n(m)))
#'
bru_mapper_repeat <- function(mapper, n_rep) {
  mapper <- list(
    mapper = mapper,
    n_rep = n_rep,
    is_linear = ibm_is_linear(mapper)
  )
  bru_mapper_define(mapper, new_class = "bru_mapper_repeat")
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
  n_offset <- n * (seq_len(mapper[["n_rep"]]) - 1L)
  sub_lin <-
    lapply(
      seq_len(mapper[["n_rep"]]),
      function(x) {
        state_subset <- state[n_offset[x] + seq_len(n)]
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
                                           ...,
                                           sub_lin = NULL) {
  if (is.null(sub_lin)) {
    sub_lin <- bm_repeat_sub_lin(mapper, input, state)
  }
  A <- lapply(sub_lin, function(x) x[["jacobian"]])
  # Combine the matrices (A1, A2, A3, ...) -> cbind(A1, A2, A3, ...)
  A <- do.call(cbind, A)
  return(A)
}


#' @export
#' @rdname bru_mapper_repeat
ibm_eval.bru_mapper_repeat <- function(mapper, input, state,
                                       ...,
                                       sub_lin = NULL) {
  if (is.null(sub_lin)) {
    sub_lin <- bm_repeat_sub_lin(mapper, input, state)
  }
  val <- lapply(sub_lin, function(x) x[["offset"]])

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
