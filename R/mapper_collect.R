#' @include deprecated.R
#' @include mappers.R

## _collect ####

#' @title Mapper for concatenated variables
#' @param mappers A list of `bru_mapper` objects
#' @param hidden `logical`, set to `TRUE` to flag that the mapper is to be used
#' as a first level input mapper for `INLA::f()` in a model that requires making
#' only the first mapper visible to `INLA::f()` and `INLA::inla.stack()`, such
#' as for "bym2" models, as activated by the `inla_f` argument to `ibm_n`,
#' `ibm_values`, and `ibm_jacobian`. Set to `FALSE` to always access the full
#' mapper, e.g. for `rgeneric` models
#' @description
#' Constructs a concatenated collection mapping
#' @export
#' @inheritParams bru_mapper_generics
#' @inheritParams bm_scale
#' @inheritParams bm_multi
#' @seealso [bru_mapper], [bru_mapper_generics]
#' @family mappers
#' @examples
#' (m <- bm_collect(list(
#'   a = bm_index(2),
#'   b = bm_index(3)
#' ), hidden = FALSE))
#' ibm_eval2(m, list(a = c(1, 2), b = c(1, 3, 2)), 1:5)
#'
#' @rdname bm_collect
bm_collect <- function(mappers, hidden = FALSE) {
  mapper <- list(
    mappers = as_bm_list(mappers),
    n_multi = lapply(mappers, ibm_n),
    values_multi = lapply(mappers, ibm_values),
    hidden = hidden,
    is_linear_multi = lapply(mappers, ibm_is_linear)
  )
  mapper[["n"]] <- sum(unlist(mapper[["n_multi"]]))
  mapper[["values"]] <- seq_len(mapper[["n"]])
  mapper[["is_linear"]] <- all(unlist(mapper[["is_linear_multi"]]))
  bru_mapper_define(mapper, new_class = "bm_collect")
}

#' @export
#' @rdname bm_collect
bru_mapper_collect <- function(...) {
  bm_collect(...)
}

#' @export
#' @rdname bm_collect
ibm_n.bm_collect <- function(mapper,
                             inla_f = FALSE,
                             multi = FALSE,
                             ...) {
  if (multi) {
    mapper[["n_multi"]]
  } else if (mapper[["hidden"]] && inla_f) {
    mapper[["n_multi"]][[1]]
  } else {
    mapper[["n"]]
  }
}


# Only for inla_f = FALSE or not hidden
bm_collect_indexing <- function(mapper, input) {
  if (is.matrix(input)) {
    nms <- colnames(input)
    len <- ncol(input)
  } else {
    stopifnot(is.list(input))
    nms <- names(input)
    len <- length(input)
  }

  nms_mapper <- names(mapper[["mappers"]])
  if (is.null(nms)) {
    indexing <- seq_len(min(length(nms_mapper), length(input)))
  } else {
    indexing <- intersect(nms_mapper, nms)
  }
  names(indexing) <- nms_mapper[indexing]
  indexing
}


#' @export
#' @rdname bm_collect
ibm_n_output.bm_collect <- function(mapper, input,
                                    state = NULL,
                                    inla_f = FALSE,
                                    multi = FALSE, ...) {
  if (mapper[["hidden"]] && inla_f) {
    return(ibm_n_output(mapper[["mappers"]][[1]], input = input))
  }

  indexing <- bm_collect_indexing(mapper, input)
  if (is.matrix(input)) {
    n <- vapply(
      indexing,
      function(x) {
        as.integer(
          ibm_n_output(mapper[["mapper"]][[x]], input[, x], ...)
        )
      },
      0L
    )
  } else {
    n <- vapply(
      indexing,
      function(x) {
        as.integer(
          ibm_n_output(mapper[["mapper"]][[x]], input[[x]], ...)
        )
      },
      0L
    )
  }

  if (!multi) {
    n <- sum(n)
  }

  n
}


#' @export
#' @rdname bm_collect
ibm_values.bm_collect <- function(mapper,
                                  inla_f = FALSE,
                                  multi = FALSE,
                                  ...) {
  if (multi) {
    mapper[["values_multi"]]
  } else if (mapper[["hidden"]] && inla_f) {
    mapper[["values_multi"]][[1]]
  } else {
    mapper[["values"]]
  }
}

#' @export
#' @rdname bm_collect
ibm_is_linear.bm_collect <- function(mapper,
                                     inla_f = FALSE,
                                     multi = FALSE,
                                     ...) {
  if (mapper[["hidden"]] && inla_f && !multi) {
    ibm_is_linear(mapper[["mappers"]][[1]])
  } else if (multi) {
    mapper[["is_linear_multi"]]
  } else {
    mapper[["is_linear"]]
  }
}



bm_collect_sub_lin <- function(mapper, input, state,
                               inla_f = FALSE,
                               ...) {
  if (mapper[["hidden"]] && inla_f) {
    input <- list(input)
  }
  indexing <- bm_collect_indexing(mapper, input)
  if (is.matrix(input)) {
    nms <- colnames(input)
    input <- as.data.frame(input)
  } else {
    nms <- names(input)
  }
  if (is.null(nms)) {
    names(input) <- names(indexing)
  }

  # We need all the sub_lin objects even if some input is NULL
  indexing <- names(mapper[["mappers"]])
  n_multi <- unlist(mapper[["n_multi"]])
  n_offset <- c(0, cumsum(n_multi))
  names(n_offset) <- c(names(n_multi), "THEEND")
  sub_lin <-
    lapply(
      indexing,
      function(x) {
        state_subset <- state[n_offset[x] + seq_len(n_multi[x])]
        ibm_linear(
          mapper[["mappers"]][[x]],
          input =
            if (is.numeric(x) && (x > length(input))) {
              NULL
            } else {
              input[[x]]
            },
          state = state_subset
        )
      }
    )
  sub_lin
}



#' @describeIn bm_collect
#' Accepts a list with
#' named entries, or a list with unnamed but ordered elements.
#' The names must match the sub-mappers, see [ibm_names.bm_collect()].
#' Each list element should take a format accepted by the corresponding
#' sub-mapper. In case each element is a vector, the input can be given as a
#' data.frame with named columns, a matrix with named columns, or a matrix
#' with unnamed but ordered columns. When `inla_f=TRUE` and `hidden=TRUE` in
#' the mapper definition, the input format should instead match that of
#' the first, non-hidden, sub-mapper.
#' @export
ibm_jacobian.bm_collect <- function(mapper, input, state = NULL,
                                    inla_f = FALSE, multi = FALSE,
                                    ...,
                                    sub_lin = NULL) {
  if (is.null(sub_lin)) {
    sub_lin <- bm_collect_sub_lin(mapper, input, state, inla_f = inla_f)
  }
  A <- lapply(sub_lin, function(x) x[["jacobian"]])

  if (multi) {
    return(A)
  }

  # Combine the matrices (A1, A2, A3, ...) -> bdiag(A1, A2, A3, ...)
  A <- Matrix::.bdiag(A)
  return(A)
}


#' @export
#' @rdname bm_collect
ibm_eval.bm_collect <- function(mapper, input, state,
                                inla_f = FALSE, multi = FALSE,
                                ...,
                                sub_lin = NULL) {
  if (is.null(sub_lin)) {
    sub_lin <- bm_collect_sub_lin(mapper, input, state, inla_f = inla_f)
  }
  val <- lapply(sub_lin, function(x) x[["offset"]])
  if (multi) {
    return(val)
  }

  # Combine the vectors (b1, b2, b3, ...) -> c(b1, b2, b3, ...)
  val <- do.call(c, val)
  val
}


#' @export
#' @rdname bm_collect
ibm_linear.bm_collect <- function(mapper, input, state,
                                  inla_f = FALSE,
                                  ...) {
  if (mapper[["hidden"]] && inla_f) {
    input <- list(input)
  }
  sub_lin <-
    bm_collect_sub_lin(mapper, input, state,
      inla_f = FALSE,
      ...
    )
  eval2 <- ibm_eval2(
    mapper,
    input = input,
    state = state,
    inla_f = FALSE,
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




#' @describeIn bm_collect
#' Accepts a list with
#' named entries, or a list with unnamed but ordered elements.
#' The names must match the sub-mappers, see [ibm_names.bm_collect()].
#' Each list element should take a format accepted by the corresponding
#' sub-mapper. In case each element is a vector, the input can be given as a
#' data.frame with named columns, a matrix with named columns, or a matrix
#' with unnamed but ordered columns.
#' @export
ibm_invalid_output.bm_collect <- function(mapper, input, state,
                                          inla_f = FALSE,
                                          multi = FALSE, ...) {
  if (mapper[["hidden"]] && inla_f) {
    return(
      ibm_invalid_output(
        mapper[["mappers"]][[1]],
        input = input,
        multi = FALSE
      )
    )
  }

  indexing <- bm_collect_indexing(mapper, input)
  if (is.matrix(input)) {
    invalid <-
      lapply(
        indexing,
        function(x) {
          ibm_invalid_output(
            mapper[["mappers"]][[x]],
            input = input[, x],
            multi = FALSE
          )
        }
      )
  } else if (is.list(input)) {
    invalid <-
      lapply(
        indexing,
        function(x) {
          ibm_invalid_output(
            mapper[["mappers"]][[x]],
            input = input[[x]],
            multi = FALSE
          )
        }
      )
  } else {
    # TODO: check and fix this!
    invalid <- as.list(rep(TRUE, length(mapper[["mappers"]])))
  }

  if (multi) {
    return(invalid)
  }

  # Combine the vectors (v1, v2, v3) -> c(v1, v2, v3)
  invalid_ <- do.call(c, invalid)
  return(invalid_)
}

#' @return
#' * `[`-indexing a `bm_collect` extracts a subset
#'   `bm_collect` object (for drop `FALSE`) or an individual sub-mapper
#'   (for drop `TRUE`, and `i` identifies a single element)
#' @export
#' @param x object from which to extract element(s)
#' @param i indices specifying element(s) to extract
#' @param drop logical;
#' For `[.bm_collect`, whether to extract an individual mapper when
#' `i` identifies a single element. If `FALSE`, a list of sub-mappers is
#' returned (suitable e.g. for creating a new `bm_collect` object).
#' Default: `TRUE`
#' @rdname bm_collect
`[.bm_collect` <- function(x, i, drop = TRUE) {
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
#' @rdname bm_collect
`[.bru_mapper_collect` <- function(x, i, drop = TRUE) {
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
#' * The `names()` method for `bm_collect` returns the names from the
#' sub-mappers list
#' @export
#' @rdname bm_collect
`ibm_names.bm_collect` <- function(mapper) {
  names(mapper[["mappers"]])
}

#' @export
#' @rdname bm_collect
`ibm_names<-.bm_collect` <- function(mapper, value) {
  names(mapper[["mappers"]]) <- value
  names(mapper[["n_multi"]]) <- value
  names(mapper[["values_multi"]]) <- value
  mapper
}

#' @export
#' @rdname bm_collect
`ibm_names<-.bru_mapper_collect` <- function(mapper, value) {
  names(mapper[["mappers"]]) <- value
  names(mapper[["n_multi"]]) <- value
  names(mapper[["values_multi"]]) <- value
  mapper
}
