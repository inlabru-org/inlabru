#' @include deprecated.R
#' @include mappers.R

## _sum ####

#' @title Mapper for adding multiple mappers
#' @param mappers A list of `bru_mapper` objects.
#' @param single_input logical. If `TRUE`, the input is passed to all
#' sub-mappers. Otherwise, the input should be a list, data.frame, or matrix.
#' If the `mappers` list has named entries, the `input` can reference their
#' corresponding sub-mapper using its name.
#' @export
#' @description Defines a mapper that adds the effects of each submapper.
#'   The `ibm_n()` method returns the sum of `ibm_n(mappers[[k]])`, and
#'   `ibm_values()` returns `seq_len(ibm_n(mapper))`.
#' @returns A `bm_sum` object.
#' @rdname bm_sum
#' @inheritParams bru_mapper_generics
#' @inheritParams bm_scale
#' @inheritParams bm_multi
#' @seealso [bru_mapper], [bru_mapper_generics]
#' @family mappers
#' @examples
#' (m <- bm_sum(list(a = bm_index(3), b = bm_index(2))))
#' ibm_n(m)
#' ibm_values(m)
#' ibm_jacobian(m, list(a = 1:3, b = c(1, 1, 2)))
#' ibm_eval(
#'   m,
#'   list(a = 1:3, b = c(1, 1, 2)),
#'   seq_len(ibm_n(m))
#' )
bm_sum <- function(mappers, single_input = FALSE) {
  mappers <- as_bm_list(mappers)
  if (is.null(names(mappers))) {
    names(mappers) <- as.character(seq_along(mappers))
  }
  mapper <- list(
    mappers = mappers,
    n_multi = lapply(mappers, ibm_n),
    values_multi = lapply(mappers, ibm_values),
    is_linear_multi = lapply(mappers, ibm_is_linear),
    single_input = single_input
  )
  mapper[["n"]] <- sum(unlist(mapper[["n_multi"]]))
  mapper[["values"]] <- seq_len(mapper[["n"]])
  mapper[["is_linear"]] <- all(unlist(mapper[["is_linear_multi"]]))
  bru_mapper_define(mapper, new_class = "bm_sum")
}

#' @export
#' @rdname bm_sum
bru_mapper_sum <- function(...) {
  bm_sum(...)
}

#' @export
#' @rdname bm_sum
ibm_n.bm_sum <- function(mapper,
                         inla_f = FALSE,
                         multi = FALSE,
                         ...) {
  if (multi) {
    mapper[["n_multi"]]
  } else {
    mapper[["n"]]
  }
}


bm_sum_prepare_input <- function(mapper, input) {
  if (mapper[["single_input"]]) {
    return(input)
  }
  if (is.matrix(input)) {
    input_names <- colnames(input)
    input <- as.data.frame(input)
  } else {
    input_names <- names(input)
  }
  if (!is.null(input_names)) {
    return(input)
  }
  if (is.data.frame(input)) {
    names(input) <- names(mapper[["mappers"]])[seq_len(ncol(input))]
  } else {
    if (is.null(input)) {
      input <- rep(list(NULL), length(mapper[["mappers"]]))
    }
    # input should now be a list
    names(input) <- names(mapper[["mappers"]])[seq_along(input)]
  }
  input
}



#' @export
#' @rdname bm_sum
ibm_n_output.bm_sum <- function(mapper, input, state = NULL, ...) {
  input <- bm_sum_prepare_input(mapper, input)
  # Assume that the first mapper fully handles the output size
  nm <- names(mapper[["mappers"]])[1]
  ibm_n_output(
    mapper[["mappers"]][[nm]],
    if (mapper[["single_input"]]) {
      input
    } else {
      input[[nm]]
    },
    state = state[seq_len(ibm_n(mapper[["mappers"]][[nm]]))],
    ...
  )
}

#' @export
#' @rdname bm_sum
ibm_values.bm_sum <- function(mapper,
                              inla_f = FALSE,
                              multi = FALSE,
                              ...) {
  if (multi) {
    mapper[["values_multi"]]
  } else {
    mapper[["values"]]
  }
}

#' @export
#' @rdname bm_sum
ibm_is_linear.bm_sum <- function(mapper,
                                 multi = FALSE,
                                 ...) {
  if (multi) {
    mapper[["is_linear_multi"]]
  } else {
    mapper[["is_linear"]]
  }
}


bm_sum_sub_lin <- function(mapper, input, state, ...) {
  input <- bm_sum_prepare_input(mapper, input)

  indexing <- names(mapper[["mappers"]])
  n_multi <- unlist(mapper[["n_multi"]])
  n_offset <- c(0, cumsum(n_multi))
  names(n_offset) <- c(names(n_multi), "THEEND")
  sub_lin <-
    lapply(
      indexing,
      function(x) {
        state_subset <- state[n_offset[x] + seq_len(mapper[["n_multi"]][[x]])]
        ibm_linear(
          mapper[["mappers"]][[x]],
          input =
            if (mapper[["single_input"]]) {
              input
            } else {
              input[[x]]
            },
          state = state_subset
        )
      }
    )
  sub_lin
}



#' @describeIn bm_sum
#' Accepts a list with
#' named entries, or a list with unnamed but ordered elements.
#' The names must match the sub-mappers, see [ibm_names.bm_sum()].
#' Each list element should take a format accepted by the corresponding
#' sub-mapper. In case each element is a vector, the input can be given as a
#' data.frame with named columns, a matrix with named columns, or a matrix
#' with unnamed but ordered columns.
#' @export
ibm_jacobian.bm_sum <- function(mapper, input, state = NULL,
                                inla_f = FALSE,
                                multi = FALSE,
                                ...,
                                sub_lin = NULL) {
  if (is.null(sub_lin)) {
    sub_lin <- bm_sum_sub_lin(mapper, input, state)
  }
  A <- lapply(sub_lin, function(x) x[["jacobian"]])

  if (multi) {
    return(A)
  }

  # Combine the matrices (A1, A2, A3, ...) -> cbind(A1, A2, A3, ...)
  A <- do.call(cbind, A)
  return(A)
}


#' @export
#' @rdname bm_sum
ibm_eval.bm_sum <- function(mapper, input, state,
                            multi = FALSE,
                            ...,
                            sub_lin = NULL) {
  if (is.null(sub_lin)) {
    sub_lin <- bm_sum_sub_lin(mapper, input, state)
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
#' @rdname bm_sum
ibm_linear.bm_sum <- function(mapper, input, state,
                              ...) {
  sub_lin <-
    bm_sum_sub_lin(
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




#' @describeIn bm_sum
#' Accepts a list with
#' named entries, or a list with unnamed but ordered elements.
#' The names must match the sub-mappers, see [ibm_names.bm_sum()].
#' Each list element should take a format accepted by the corresponding
#' sub-mapper. In case each element is a vector, the input can be given as a
#' data.frame with named columns, a matrix with named columns, or a matrix
#' with unnamed but ordered columns.
#' @export
ibm_invalid_output.bm_sum <- function(mapper, input, state,
                                      multi = FALSE,
                                      ...) {
  input <- bm_sum_prepare_input(mapper, input)
  invalid <-
    lapply(
      names(mapper[["mappers"]]),
      function(x) {
        if (mapper[["single_input"]]) {
          inp <- input
        } else {
          inp <- input[[x]]
        }
        ibm_invalid_output(
          mapper[["mappers"]][[x]],
          input = inp,
          multi = FALSE
        )
      }
    )

  if (multi) {
    return(invalid)
  }

  # Combine the vectors (v1, v2, v3) -> (v1 | v2 | v3)
  invalid_ <- invalid[[1]]
  for (k in seq_len(length(invalid) - 1L) + 1L) {
    invalid_ <- invalid_ | invalid[[k]]
  }
  return(invalid_)
}

#' @return
#' * `[`-indexing a `bm_sum` extracts a subset
#'   `bm_sum` object (for drop `FALSE`) or an individual sub-mapper
#'   (for drop `TRUE`, and `i` identifies a single element)
#' @export
#' @param x object from which to extract element(s)
#' @param i indices specifying element(s) to extract
#' @param drop logical;
#' For `[.bm_sum`, whether to extract an individual mapper when
#' `i` identifies a single element. If `FALSE`, a list of sub-mappers is
#' returned (suitable e.g. for creating a new `bm_sum` object).
#' Default: `TRUE`
#' @rdname bm_sum
`[.bm_sum` <- function(x, i, drop = TRUE) {
  if (is.logical(i)) {
    i <- which(i)
  }
  mapper <- x[["mappers"]][i]
  if (drop) {
    if (length(mapper) == 1) {
      mapper <- mapper[[1]]
    } else if (length(mapper) == 0) {
      mapper <- NULL
    }
  }
  mapper
}

#' @export
#' @rdname bm_sum
`[.bru_mapper_sum` <- function(x, i, drop = TRUE) {
  if (is.logical(i)) {
    i <- which(i)
  }
  mapper <- x[["mappers"]][i]
  if (drop) {
    if (length(mapper) == 1) {
      mapper <- mapper[[1]]
    } else if (length(mapper) == 0) {
      mapper <- NULL
    }
  }
  mapper
}

#' @return
#' * The `names()` method for `bm_sum` returns the names from the
#' sub-mappers list
#' @export
#' @rdname bm_sum
`ibm_names.bm_sum` <- function(mapper) {
  names(mapper[["mappers"]])
}

#' @export
#' @rdname bm_sum
`ibm_names<-.bm_sum` <- function(mapper, value) {
  names(mapper[["mappers"]]) <- value
  names(mapper[["n_multi"]]) <- value
  names(mapper[["values_multi"]]) <- value
  mapper
}

#' @export
#' @rdname bm_sum
`ibm_names<-.bru_mapper_sum` <- function(mapper, value) {
  names(mapper[["mappers"]]) <- value
  names(mapper[["n_multi"]]) <- value
  names(mapper[["values_multi"]]) <- value
  mapper
}
