#' Methods for `bru_mapper` objects
#'
#' @details
#' * `bru_mapper` Generic mapper S3 constructor. See below for details of the
#' default constructor that can be used to define new mappers in user code.
#' @param \dots Arguments passed on to other methods
#' @param mapper A mapper S3 object, normally inheriting from `bru_mapper`.
#' For the default `bru_mapper` method, a list that will be converted to a
#' `bru_mapper` object by adding class information and (optional) methods.
#' @export
#' @seealso [bru_mapper_methods] for specific method implementations.
#' @rdname bru_mapper
#' @examples
#' mapper <- bru_mapper_index(5)
#' ibm_amatrix(mapper, c(1, 3, 4, 5, 2))
bru_mapper <- function(...) {
  UseMethod("bru_mapper")
}
#' @details
#' * `ibm_n` Generic. Implementations must return the size of the latent vector
#' being mapped to.
#' @param inla_f logical; when `TRUE` in `ibm_n`, `ibm_values`, and
#' `ibm_amatrix` methods, these must result in values compatible with `INLA::f(...)`
#' an specification and corresponding `INLA::inla.stack(...)` constructions.
#' Implementations do not normally need to do anything different, except
#' for mappers of the type needed for hidden multicomponent models such
#' as "bym2", which can be handled by `bru_mapper_collect`.
#'
#' @export
#' @rdname bru_mapper
ibm_n <- function(mapper, inla_f = FALSE, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_n"]])) {
    mapper[[".envir"]][["ibm_n"]](mapper, inla_f = inla_f, ...)
  } else {
    UseMethod("ibm_n")
  }
}
#' @details
#' * `ibm_values` Generic. Implementations must return a vector that
#' would be interpretable by an `INLA::f(..., values = ...)` specification.
#' The exception is the method for `bru_mapper_multi`, that returns a
#' multi-column data frame
#' @export
#' @rdname bru_mapper
ibm_values <- function(mapper, inla_f = FALSE, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_values"]])) {
    mapper[[".envir"]][["ibm_values"]](mapper, inla_f = inla_f, ...)
  } else {
    UseMethod("ibm_values")
  }
}
#' @details
#' * `ibm_amatrix` Generic.
#' Implementations must return a (sparse) matrix of size `NROW(input)`
#' (except for the `bru_mapper_multi` and `bru_mapper_collect` methods,
#' that require `list()` inputs, and the input size is determined by the
#' combined inputs)
#' by `ibm_n(mapper, inla_f = FALSE)`. The `inla_f=TRUE` argument should only affect
#' the allowed type of input format.
#' @param input The values for which to produce a mapping matrix
#' @export
#' @rdname bru_mapper
ibm_amatrix <- function(mapper, input, inla_f = FALSE, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_amatrix"]])) {
    mapper[[".envir"]][["ibm_amatrix"]](
      mapper, input, inla_f = inla_f, ...)
  } else {
    UseMethod("ibm_amatrix")
  }
}
#' @details
#' * `ibm_inla_subset` Generic.
#' Implementations must return a logical vector of `TRUE/FALSE` for
#' the subset such that, given the full A matrix and values output,
#' `A[, subset, drop = FALSE]` and `values[subset]`
#' (or `values[subset, , drop = FALSE]` for data.frame values) are equal
#' to the `inla_f = TRUE` version of A and values. The default method uses
#' the `ibm_values` output to construct the subset indexing.
#' @export
#' @rdname bru_mapper
ibm_inla_subset <- function(mapper, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_inla_subset"]])) {
    mapper[[".envir"]][["ibm_inla_subset"]](
      mapper, ...)
  } else {
    UseMethod("ibm_inla_subset")
  }
}
#' @details
#' * `ibm_valid_input` Generic.
#' Implementations must return a logical vector of length `NROW(input)` (
#' or for `bru_mapper_multi` and `bru_mapper_collect` a list of such
#' vectors)
#' @param input The values for which to produce validity information
#' @export
#' @rdname bru_mapper
ibm_valid_input <- function(mapper, input, inla_f = FALSE, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_valid_input"]])) {
    mapper[[".envir"]][["ibm_valid_input"]](
      mapper, input, inla_f = inla_f, ...)
  } else {
    UseMethod("ibm_valid_input")
  }
}


# MAPPERS ----

#' @details * `bru_mapper.default` adds the "bru_mapper" class and `new_class`
#' to an object. If provided, mapper method functions are added to an environment
#' `.envir` in the object.  The generic methods `ibm_n`, `ibm_n_inla`,
#' `ibm_values`, `ibm_values_inla`,
#' `ibm_amatrix`, `ibm_amatrix_inla`,
#' `ibm_valid_input`, and `ibm_valid_input_inla` look for these
#' functions first,
#' and otherwise call `UseMethod()`.  This is an alternative to using `.S3method()`
#' to register the methods, e.g.
#' `.S3method("ibm_amatrix", "my_mapper_class", ibm_amatrix.my_mapper_class)`.
#' @param new_class If non-`NULL`, this is added at the front of the class definition
#' @param methods, optional `list` of named method definitions; See Details.
#'
#' @export
#' @rdname bru_mapper
bru_mapper.default <- function(mapper,
                               new_class = NULL,
                               methods = NULL,
                               ...) {
  if (any(c(
    "ibm_n", "ibm_n_inla",
    "ibm_values", "ibm_values_inla",
    "ibm_amatrix", "ibm_amatrix_inla",
    "ibm_valid_input", "ibm_valid_input_inla"
  ) %in%
    names(list(...)))) {
    warning(
      paste0(
        "Deprecated use of named method arguments for 'bru_mapper'.\n",
        "Use methods = list(methodname = ..., ...) instead."
      )
    )
    method_names <- intersect(
      c(
        "ibm_n",
        "ibm_values",
        "ibm_amatrix",
        "ibm_inla_subset",
        "ibm_valid_input"
      ),
      names(list(...))
    )
    methods <- list(...)[method_names]
  }
  if (!inherits(mapper, "bru_mapper")) {
    class(mapper) <- c("bru_mapper", class(mapper))
  }
  if (!is.null(new_class) && !inherits(mapper, new_class)) {
    class(mapper) <- c(new_class, class(mapper))
  }
  if (is.null(mapper[[".envir"]])) {
    mapper$.envir <- new.env()
  }
  if (!is.null(methods)) {
    for (method in names(methods)) {
      if (!is.null(methods[[method]])) {
        assign(method, methods[[method]], envir = mapper[[".envir"]])
      }
    }
  }
  mapper
}


#' Implementation methods for mapper objects
#'
#' A `bru_mapper` sub-class implementation must provide an
#' `ibm_matrix()` method. If the model size 'n' and definition
#' values 'values' are stored in the object itself, default methods are
#' available (see Details). Otherwise the
#' `ibm_n()` and `ibm_values()` methods also need to be provided.
#'
#' @param \dots Arguments passed on to other methods
#' @param mapper A mapper S3 object, normally inheriting from `bru_mapper`
#' @param inla_f logical; when `TRUE` in `ibm_n` and `ibm_values`, these must result in values compatible with `INLA::f(...)`
#' an specification and corresponding `INLA::inla.stack(...)` constructions.
#' For `ibm_amatrix` methods, it may influence how the input data is interpreted.
#' Implementations do not normally need to do anything different, except
#' for mappers of the type needed for hidden multicomponent models such
#' as "bym2", which can be handled by `bru_mapper_collect`.
#' @seealso [bru_mapper()]
#' @name bru_mapper_methods
#'
#' @details
#' * The default `ibm_n()` method returns a non-null element 'n' from the
#' mapper object, and gives an error if it doesn't exist. If `inla_f=TRUE`,
#' first checks for a 'n_inla' element.
#' @export
#' @rdname bru_mapper_methods
ibm_n.default <- function(mapper, inla_f = FALSE, ...) {
  if (inla_f && !is.null(mapper[["n_inla"]])) {
    mapper[["n_inla"]]
  } else if (!is.null(mapper[["n"]])) {
    mapper[["n"]]
  } else {
    stop("Default 'ibm_n()' method called but mapper doesn't have an 'n' element.")
  }
}


#' @details
#' * The default `ibm_values()` method returns a non-null element
#' 'values' from the mapper object, and `seq_len(ibm_n(mapper))` if
#' it doesn't exist.
#' @export
#' @rdname bru_mapper_methods
ibm_values.default <- function(mapper, inla_f = FALSE, ...) {
  if (inla_f && !is.null(mapper[["values_inla"]])) {
    mapper[["values_inla"]]
  } else if (!inla_f && !is.null(mapper[["values"]])) {
    mapper[["values"]]
  } else {
    seq_len(ibm_n(mapper, inla_f = inla_f))
  }
}

#' @details
#' * The default `ibm_amatrix()` gives an error message.
#' Mapper classes must implement their own `ibm_amatrix` method.
#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.default <- function(mapper, input, inla_f = FALSE, ...) {
  stop(paste0(
    "Missing implementation of 'ibm_amatrix()' for mapper of class '",
    class(mapper)[1], "'."
  ))
}

#' @details
#' * The default `ibm_inla_subset` method uses
#' the `ibm_values` output to construct the inla subset indexing, passing
#' extra arguments such as `multi` on to the methods (this means it supports
#' both regular vector values and `multi=1` data.frame values).
#' @export
#' @rdname bru_mapper
ibm_inla_subset.default <- function(mapper, ...) {
  values_full <- ibm_values(mapper, inla_f = FALSE, ...)
  values_inla <- ibm_values(mapper, inla_f = TRUE, ...)
  if (is.data.frame(values_full)) {
    subset <- logical(NROW(values_full))
    if (length(subset) > 0) {
      subset[
        plyr::match_df(
          cbind(
            .inla_subset = seq_len(NROW(values_full)),
            values_full
          ),
          values_inla,
          on = names(values_full)
        )[[".inla_subset"]]
      ] <- TRUE
    }
  } else {
    subset <- base::match(values_full, values_inla, nomatch = 0) > 0
  }
  subset
}


#' @details
#' * The default `ibm_valid_input()` method returns an all-TRUE logical vector.
#' @export
#' @rdname bru_mapper_methods
ibm_valid_input.default <- function(mapper, input, ...) {
  rep(TRUE, NROW(input))
}





## inla.mesh ####

#' @param mesh An `inla.mesh.1d` or `inla.mesh.2d` object to use as a mapper
#' @param indexed logical; If `TRUE`, the `ibm_values()` output will be the
#' integer indexing sequence for the latent variables. If `FALSE`, the knot
#' locations are returned (useful as an interpolator for `rw2` models
#' and similar).
#' Default: `TRUE`
#' @export
#' @rdname bru_mapper
bru_mapper.inla.mesh <- function(mesh, ...) {
  mapper <- list(mesh = mesh)
  class(mapper) <- c("bru_mapper_inla_mesh_2d", "list")
  bru_mapper.default(mapper)
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_inla_mesh_2d <- function(mapper, ...) {
  mapper[["mesh"]]$n
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_inla_mesh_2d <- function(mapper, ...) {
  seq_len(mapper[["mesh"]]$n)
}
#' @param input The values for which to produce a mapping matrix
#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_inla_mesh_2d <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  }
  if (!is.matrix(input) && !inherits(input, "Spatial")) {
    input <- as.matrix(input)
  }
  INLA::inla.spde.make.A(mapper[["mesh"]], loc = input)
}

## inla.mesh.1d ####

#' @param indexed logical; If `TRUE`, the `ibm_values()` output will be the
#' integer indexing sequence for the latent variables (needed for `spde` models).
#' If `FALSE`, the knot
#' locations are returned (useful as an interpolator for `rw2` models
#' and similar).
#' Default: `NULL`, to force user specification of this parameter
#' @export
#' @rdname bru_mapper
bru_mapper.inla.mesh.1d <- function(mesh, indexed = NULL, ...) {
  if (is.null(indexed)) {
    stop("indexed=TRUE/FALSE needs to be specified to convert inla.mesh.1d to a bru_mapper")
  }
  mapper <- list(
    mesh = mesh,
    indexed = indexed
  )
  class(mapper) <- c("bru_mapper_inla_mesh_1d", "list")
  bru_mapper.default(mapper)
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_inla_mesh_1d <- function(mapper, ...) {
  mapper[["mesh"]]$m
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_inla_mesh_1d <- function(mapper, ...) {
  if (mapper[["indexed"]]) {
    seq_len(mapper[["mesh"]]$m)
  } else {
    mapper[["mesh"]]$loc
  }
}
#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_inla_mesh_1d <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  }
  INLA::inla.spde.make.A(mapper[["mesh"]], loc = input)
}

## _index ####

#' @param n Size of a model for `bru_mapper_index`
#' @export
#' @rdname bru_mapper
bru_mapper_index <- function(n = 1L, ...) {
  mapper <- list(
    n = n
  )
  class(mapper) <- c("bru_mapper_index", "bru_mapper", "list")
  mapper
}

#' @export
#' @rdname bru_mapper_methods
ibm_valid_input.bru_mapper_index <- function(mapper, input, ...) {
  ok <- !is.na(input)
  ok[ok] <- (input[ok] >= 1) & (input[ok] <= ibm_n(mapper))
  ok
}

#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_index <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  }
  ok <- which(ibm_valid_input(mapper, input, ...))
  Matrix::sparseMatrix(
    i = ok,
    j = input[ok],
    x = rep(1, length(ok)),
    dims = c(NROW(input), ibm_n(mapper))
  )
}

## _linear ####

#' @export
#' @rdname bru_mapper
bru_mapper_linear <- function(...) {
  mapper <- list()
  class(mapper) <- c("bru_mapper_linear", "bru_mapper", "list")
  mapper
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_linear <- function(mapper, ...) {
  1L
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_linear <- function(mapper, ...) {
  1.0
}

#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_linear <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  }
  ok <- !is.na(input)
  A <- Matrix::sparseMatrix(
    i = which(ok),
    j = rep(1, sum(ok)),
    x = input[ok],
    dims = c(NROW(input), ibm_n(mapper))
  )
  A
}

## _matrix ####

#' @param labels Column labels for matrix mappings
#' @export
#' @rdname bru_mapper
bru_mapper_matrix <- function(labels, ...) {
  if (is.factor(labels)) {
    mapper <- list(
      labels = levels(labels)
    )
  } else {
    mapper <- list(
      labels = as.character(labels)
    )
  }
  class(mapper) <- c("bru_mapper_matrix", "bru_mapper", "list")
  mapper
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_matrix <- function(mapper, ...) {
  length(mapper$labels)
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_matrix <- function(mapper, ...) {
  factor(x = mapper$labels, levels = mapper$labels)
}

#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_matrix <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  } else if (is.matrix(input)) {
    A <- as(input, "Matrix")
  } else if (inherits(input, "Matrix")) {
    A <- input
  } else if (inherits(input, "Spatial")) {
    A <- sp::coordinates(input)
    A <- as(A, "Matrix")
  } else {
    A <- as(input, "Matrix")
  }
  if (ncol(A) != ibm_n(mapper)) {
    stop(paste0(
      "Input to matrix mapper has ", ncol(A),
      " columns but should have ", ibm_n(mapper),
      " columns."
    ))
  }
  A
}



## _factor ####

#' @param values Input values calculated by [input_eval.bru_input()]
#' @param factor_mapping character; selects the type of factor mapping.
#' * `'contrast'` for leaving out the first factor level.
#' * `'full'` for keeping all levels.
#' @export
#' @rdname bru_mapper
bru_mapper_factor <- function(values, factor_mapping, ...) {
  factor_mapping <- match.arg(factor_mapping, c("full", "contrast"))
  if (is.factor(values)) {
    mapper <- list(
      levels = levels(values),
      factor_mapping = factor_mapping
    )
  } else if (is.character(values)) {
    mapper <- list(
      levels = unique(values),
      factor_mapping = factor_mapping
    )
  } else {
    mapper <- list(
      levels = as.character(sort(unique(values))),
      factor_mapping = factor_mapping
    )
  }
  class(mapper) <- c("bru_mapper_factor", "bru_mapper", "list")
  mapper
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_factor <- function(mapper, ...) {
  length(ibm_values(mapper))
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_factor <- function(mapper, ...) {
  if (identical(mapper$factor_mapping, "contrast")) {
    mapper$levels[-1L]
  } else {
    mapper$levels
  }
}

#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_factor <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  }
  if (is.factor(input)) {
    if (!identical(levels(input), mapper$levels)) {
      input <- factor(as.character(input), levels = mapper$levels)
    }
  } else {
    input <- factor(input, levels = mapper$levels)
  }
  if (mapper$factor_mapping == "full") {
    input <- as.numeric(input)
  } else if (mapper$factor_mapping == "contrast") {
    input <- as.numeric(input) - 1L
  } else {
    stop("Unknown factor mapping '", mapper$factor_mapping, "'.")
  }
  ok <- !is.na(input)
  ok[which(ok)] <- (input[ok] > 0L)
  A <- Matrix::sparseMatrix(
    i = which(ok),
    j = input[ok],
    x = rep(1.0, sum(ok)),
    dims = c(NROW(input), ibm_n(mapper))
  )
  A
}



## _offset ####

#' @param values Input values calculated by [input_eval.bru_input()]
#' @export
#' @rdname bru_mapper
bru_mapper_offset <- function(...) {
  mapper <- list()
  class(mapper) <- c("bru_mapper_offset", "bru_mapper", "list")
  mapper
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_offset <- function(mapper, ...) {
  0L
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_offset <- function(mapper, ...) {
  NULL
}

#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_offset <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, 1L))
  }
  ok <- !is.na(input)
  ok[which(ok)] <- (input[ok] > 0L)
  A <- Matrix::sparseMatrix(
    i = which(ok),
    j = rep(1L, sum(ok)),
    x = input[ok],
    dims = c(NROW(input), 1L)
  )
  A
}


## _multi ####

#' @param values Input values calculated by [input_eval.bru_input()]
#' @param mappers A list of `bru_mapper` objects
#' @details * `bru_mapper_multi` constructs a kronecker product mapping
#' @export
#' @rdname bru_mapper
bru_mapper_multi <- function(mappers, ...) {
  mapper <- list(
    mappers = mappers,
    n_multi = lapply(mappers, function(x) {
      ibm_n(x)
    }),
    n_inla_multi = lapply(mappers, function(x) {
      ibm_n(x, inla_f = TRUE)
    }),
    values_multi = lapply(mappers, function(x) {
      ibm_values(x)
    }),
    values_inla_multi = lapply(mappers, function(x) {
      ibm_values(x, inla_f = TRUE)
    })
  )
  mapper[["n"]] <- prod(unlist(mapper[["n_multi"]]))
  mapper[["n_inla"]] <- prod(unlist(mapper[["n_inla_multi"]]))
  class(mapper) <- c("bru_mapper_multi", "bru_mapper", "list")
  mapper
}

#' @param multi integer or logical;
#' If positive, the number of levels to recurse in a `bru_multi_mapper`.
#' If `TRUE`, equivalent to `1L`. If `FALSE`, equivalent to `0L`.
#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_multi <- function(mapper, inla_f = FALSE, multi = 0L, ...) {
  if (multi > 1) {
    lapply(mapper[["mappers"]], function(x) {
      ibm_n(x, multi = multi - 1)
    })
  } else if (multi == 1) {
    if (inla_f) {
      mapper[["n_inla_multi"]]
    } else {
      mapper[["n_multi"]]
    }
  } else {
    if (inla_f) {
      mapper[["n_inla"]]
    } else {
      mapper[["n"]]
    }
  }
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_multi <- function(mapper, inla_f = FALSE, multi = 0L, ...) {
  if (multi > 1) {
    lapply(mapper[["mappers"]], function(x) {
      ibm_values(x, multi = multi - 1)
    })
  } else if (multi == 1) {
    # Expand indices/values. First sub-mapper varies fastest
    as.data.frame(
      do.call(
        expand.grid,
        c(
          if (inla_f) {
            mapper[["values_inla_multi"]]
          } else {
            mapper[["values_multi"]]
          },
          list(
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
          )
        )
      )
    )
  } else {
    if (inla_f) {
      seq_len(mapper[["n_inla"]])
    } else {
      seq_len(mapper[["n"]])
    }
  }
}


#' @details
#' * `ibm_amatrix` for `bru_mapper_multi` accepts a list with
#' named entries, or a list with unnamed but ordered elements.
#' The names must match the sub-mappers, see [names.bru_mapper_multi()].
#' Each list element should take a format accepted by the corresponding
#' sub-mapper. In case each element is a vector, the input can be given as a
#' data.frame with named columns, a matrix with named columns, or a matrix
#' with unnamed but ordered columns.
#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_multi <- function(mapper, input,
                                         inla_f = FALSE, multi = 0L, ...) {
  by_number <- FALSE
  if (is.matrix(input)) {
    if (is.null(colnames(input))) {
      by_number <- TRUE
    }
    input <- as.data.frame(input)
  }
  if (by_number || is.null(names(input))) {
    indexing <- seq_along(mapper[["mappers"]])
  } else {
    indexing <- names(mapper[["mappers"]])
  }
  A <-
    lapply(
      indexing,
      function(x) {
        ibm_amatrix(mapper[["mappers"]][[x]],
          input = input[[x]],
          inla_f = inla_f,
          multi = multi - 1
        )
      }
    )
  if (multi < 1) {
    # Combine the matrices (A1, A2, A3) -> rowkron(A3, rowkron(A2, A1))
    A_ <- A[[1]]
    for (k in seq_len(length(mapper[["mappers"]]) - 1)) {
      A_ <- row_kron(A[[k + 1]], A_)
    }
    return(A_)
  }
  A
}

#' @details
#' * `ibm_valid_input` for `bru_mapper_multi` accepts a list with
#' named entries, or a list with unnamed but ordered elements.
#' The names must match the sub-mappers, see [names.bru_mapper_multi()].
#' Each list element should take a format accepted by the corresponding
#' sub-mapper. In case each element is a vector, the input can be given as a
#' data.frame with named columns, a matrix with named columns, or a matrix
#' with unnamed but ordered columns.
#' @export
#' @rdname bru_mapper_methods
ibm_valid_input.bru_mapper_multi <- function(mapper, input,
                                             inla_f = FALSE, multi = 0L, ...) {
  by_number <- FALSE
  if (is.matrix(input)) {
    if (is.null(colnames(input))) {
      by_number <- TRUE
    }
    input <- as.data.frame(input)
  }
  if (!is.list(input)) {
    validity <- as.list(rep(FALSE, length(mapper[["mappers"]])))
  } else {
    if (by_number || is.null(names(input))) {
      indexing <- seq_along(mapper[["mappers"]])
    } else {
      indexing <- names(mapper[["mappers"]])
    }
    validity <-
      lapply(
        indexing,
        function(x) {
          ibm_valid_input(mapper[["mappers"]][[x]],
            input = input[[x]],
            inla_f = inla_f,
            multi = 0
          )
        }
      )
  }
  if (multi < 1) {
    # Combine the vectors (v1, v2, v3) -> v1 & v2 & v3
    validity_ <- validity[[1]]
    for (k in seq_len(length(mapper[["mappers"]]) - 1)) {
      validity_ <- validity_ & validity[[k + 1]]
    }
    return(validity_)
  }
  validity
}

#' @return
#' * `[`-indexing a `bru_mapper_multi` extracts a subset
#'   `bru_mapper_multi` object (for drop `FALSE`) or an individual sub-mapper
#'   (for drop `TRUE`, and `i` identifies a single element)
#' @export
#' @param x object from which to extract element(s)
#' @param i indices specifying element(s) to extract
#' @param drop logical;
#' For `[.bru_mapper_multi`, whether to extract an individual mapper when
#' `i` identifies a single element. If `FALSE`, a list of sub-mappers is
#' returned (suitable e.g. for creating a new `bru_mapper_multi` object).
#' Default: `TRUE`
#' @rdname bru_mapper_methods
`[.bru_mapper_multi` <- function(x, i, drop = TRUE) {
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
#' * The `names()` method for `bru_mapper_multi` returns the names from the
#' sub-mappers list
#' @export
#' @rdname bru_mapper_methods
`names.bru_mapper_multi` <- function(x) {
  names(x[["mappers"]])
}

#' @param value a character vector of up to the same length as the number
#' of mappers in the multi-mapper x
#' @export
#' @rdname bru_mapper_methods
`names<-.bru_mapper_multi` <- function(x, value) {
  names(x[["mappers"]]) <- value
  names(x[["n_multi"]]) <- value
  names(x[["n_inla_multi"]]) <- value
  names(x[["values_multi"]]) <- value
  names(x[["values_inla_multi"]]) <- value
  x
}



## _collect ####

#' @param values Input values calculated by [input_eval.bru_input()]
#' @param mappers A list of `bru_mapper` objects
#' @param hidden `logical`, set to `TRUE` to flag that the mapper is to be used
#' as a first level input mapper for `INLA::f()` in a model that requires making
#' only the first mapper visible to `INLA::f()` and `INLA::inla.stack()`, such
#' as for "bym2" models, as activated by the `inla_f` argument to `ibm_n`,
#' `ibm_values`, and `ibm_amatrix`. Set to `FALSE` to always access the full
#' mapper, e.g. for `rgeneric` models
#' @details * `bru_mapper_collect` constructs concatenated collection mapping
#' @export
#' @rdname bru_mapper
bru_mapper_collect <- function(mappers, hidden = FALSE, ...) {
  mapper <- list(
    mappers = mappers,
    n_multi = lapply(mappers, function(x) {
      ibm_n(x)
    }),
    values_multi = lapply(mappers, function(x) {
      ibm_values(x)
    }),
    hidden = hidden
  )
  mapper[["n"]] <- sum(unlist(mapper[["n_multi"]]))
  mapper[["values"]] <- seq_len(mapper[["n"]])
  class(mapper) <- c("bru_mapper_collect", "bru_mapper", "list")
  mapper
}

#' @param multi integer or logical;
#' If positive, the number of levels to recurse in a `bru_collect_mapper`.
#' If `TRUE`, equivalent to `1L`. If `FALSE`, equivalent to `0L`.
#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_collect <- function(mapper, inla_f = FALSE, multi = 0L, ...) {
  if (mapper[["hidden"]] && inla_f) {
    if (multi > 0) {
      ibm_n(mapper, multi = multi, ...)
    } else {
      ibm_n(mapper[["mappers"]][[1]])
    }
  } else {
    if (multi > 1) {
      lapply(mapper[["mappers"]], function(x) {
        ibm_n(x, multi = multi - 1)
      })
    } else if (multi == 1) {
      mapper[["n_multi"]]
    } else {
      mapper[["n"]]
    }
  }
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_collect <- function(mapper, inla_f = FALSE, multi = 0L, ...) {
  if (mapper[["hidden"]] && inla_f) {
    if (multi > 0) {
      ibm_values(mapper, multi = multi, ...)
    } else {
      ibm_values(mapper[["mappers"]][[1]])
    }
  } else {
    if (multi > 1) {
      lapply(mapper[["mappers"]], function(x) {
        ibm_values(x, multi = multi - 1)
      })
    } else if (multi == 1) {
      mapper[["values_multi"]]
    } else {
      mapper[["values"]]
    }
  }
}
#' @details
#' * `ibm_amatrix` for `bru_mapper_collect` accepts a list with
#' named entries, or a list with unnamed but ordered elements.
#' The names must match the sub-mappers, see [names.bru_mapper_collect()].
#' Each list element should take a format accepted by the corresponding
#' sub-mapper. In case each element is a vector, the input can be given as a
#' data.frame with named columns, a matrix with named columns, or a matrix
#' with unnamed but ordered columns. When `inla_f=TRUE` and `hidden=TRUE` in
#' the mapper definition, the input format should instead match that of
#' the first, non-hidden, sub-mapper.
#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_collect <- function(mapper, input, inla_f = FALSE, multi = 0L, ...) {
  if (mapper[["hidden"]] && inla_f) {
    input <- list(input)
  }
  if (is.null(names(input))) {
    indexing <- seq_len(length(mapper[["mappers"]]))
    A <-
      lapply(
        indexing,
        function(x) {
          ibm_amatrix(
            mapper[["mappers"]][[x]],
            input = if (x <= length(input)) {
              input[[x]]
            } else {
              NULL
            },
            multi = multi - 1
          )
        }
      )
  } else {
    indexing <- names(mapper[["mappers"]])
    A <-
      lapply(
        indexing,
        function(x) {
          ibm_amatrix(
            mapper[["mappers"]][[x]],
            input = input[[x]],
            multi = multi - 1
          )
        }
      )
  }
  if (multi < 1) {
    # Combine the matrices (A1, A2, A3, ...) -> bdiag(A1, A2, A3, ...)
    A_ <- Matrix::.bdiag(A)
    return(A_)
  }
  A
}

#' @details
#' * `ibm_valid_input` for `bru_mapper_collect` accepts a list with
#' named entries, or a list with unnamed but ordered elements.
#' The names must match the sub-mappers, see [names.bru_mapper_collect()].
#' Each list element should take a format accepted by the corresponding
#' sub-mapper. In case each element is a vector, the input can be given as a
#' data.frame with named columns, a matrix with named columns, or a matrix
#' with unnamed but ordered columns.
#' @export
#' @rdname bru_mapper_methods
ibm_valid_input.bru_mapper_collect <- function(mapper, input, inla_f = FALSE, multi = 0L, ...) {
  if (mapper[["hidden"]] && inla_f) {
    validity <- ibm_valid_input(mapper[["mappers"]][[1]],
      input = input,
      multi = multi - 1
    )
  } else {
    if (!is.list(input)) {
      validity <- as.list(rep(FALSE, mapper[["mappers"]]))
    } else {
      if (is.null(names(input))) {
        indexing <- seq_along(mapper[["mappers"]])
      } else {
        indexing <- names(mapper[["mappers"]])
      }
      validity <-
        lapply(
          indexing,
          function(x) {
            ibm_valid_input(mapper[["mappers"]][[x]],
              input = input[[x]],
              multi = 0L
            )
          }
        )
    }
    if (multi < 1) {
      # Combine the vectors (v1, v2, v3) -> c(v1, v2, v3)
      validity_ <- do.call(c, validity)
      return(validity_)
    }
  }
  validity
}

#' @return
#' * `[`-indexing a `bru_mapper_collect` extracts a subset
#'   `bru_mapper_collect` object (for drop `FALSE`) or an individual sub-mapper
#'   (for drop `TRUE`, and `i` identifies a single element)
#' @export
#' @param x object from which to extract element(s)
#' @param i indices specifying element(s) to extract
#' @param drop logical;
#' For `[.bru_mapper_collect`, whether to extract an individual mapper when
#' `i` identifies a single element. If `FALSE`, a list of sub-mappers is
#' returned (suitable e.g. for creating a new `bru_mapper_collect` object).
#' Default: `TRUE`
#' @rdname bru_mapper_methods
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
#' * The `names()` method for `bru_mapper_collect` returns the names from the
#' sub-mappers list
#' @export
#' @rdname bru_mapper_methods
`names.bru_mapper_collect` <- function(x) {
  names(x[["mappers"]])
}

#' @param value a character vector of up to the same length as x
#' @export
#' @rdname bru_mapper_methods
`names<-.bru_mapper_collect` <- function(x, value) {
  names(x[["mappers"]]) <- value
  names(x[["n_multi"]]) <- value
  names(x[["values_multi"]]) <- value
  x
}








## _harmonics ####


#' @param order For `bru_mapper_harmonics`, specifies the maximum `cos`/`sin`
#' order. (Default 1)
#' @param scaling For `bru_mapper_harmonics`, specifies an optional vector of
#'  scaling factors of length `intercept + order`, or a common single scalar.
#' @param intercept logical; For `bru_mapper_harmonics`, if `TRUE`, the first
#'   basis function is a constant. (Default `TRUE`)
#' @param interval numeric length-2 vector specifying a domain interval.
#'   Default `c(0, 1)`.
#' @details * `bru_mapper_harmonics` constructs a mapper for `cos`/`sin` functions
#'   of orders 1 (if `intercept` is `TRUE`, otherwise 0) through `order`. The total
#'   number of basis functions is `intercept + 2 * order`.
#'
#'   Optionally, each order can be given a non-unit scaling, via the `scaling`
#'   vector, of length `intercept + order`. This can be used to
#'   give an effective spectral prior. For example, let
#'   ```
#'   scaling = 1 / (1 + (0:4)^2)
#'   A1 = bru_mapper_harmonics(order = 4)
#'   u1 <- A1 %*% rnorm(9, sd = scaling)
#'   ```
#'   Then, with
#'   ```
#'   A2 = bru_mapper_harmonics(order = 4, scaling = scaling)
#'   u2 = A2 %*% rnorm(9)
#'   ```
#'   the stochastic properties of `u1` and `u2` will be the same, with `scaling^2`
#'   determining the variance for each frequency contribution.
#'
#'   The period for the first order harmonics is shifted and scaled to match
#'   `interval`.
#' @export
#' @rdname bru_mapper
bru_mapper_harmonics <- function(order = 1,
                                 scaling = 1,
                                 intercept = TRUE,
                                 interval = c(0, 1),
                                 ...) {
  if (length(scaling) == 1) {
    scaling <- rep(scaling, intercept + order)
  } else {
    stopifnot(length(scaling) == intercept + order)
  }
  mapper <- list(
    order = order,
    scaling = scaling,
    intercept = intercept,
    interval = interval
  )
  bru_mapper(mapper, new_class = "bru_mapper_harmonics")
}


#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_harmonics <- function(mapper, inla_f = FALSE, ...) {
  mapper[["order"]] * 2 + mapper[["intercept"]]
}

#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.bru_mapper_harmonics <- function(mapper, input, inla_f = FALSE, ...) {
  A <- Matrix::Matrix(0.0, NROW(input), ibm_n(mapper))
  off <- 0
  if (mapper[["intercept"]]) {
    A[, 1] <- 1.0 * mapper[["scaling"]][1]
    off <- off + 1
  }
  if (mapper[["order"]] > 0) {
    input <- (input - mapper[["interval"]][1]) / diff(mapper[["interval"]])
    for (ord in seq_len(mapper[["order"]])) {
      scale <- mapper[["scaling"]][mapper[["intercept"]] + ord]
      A[, off + 1] <- cos(2 * pi * input * ord) * scale
      A[, off + 2] <- sin(2 * pi * input * ord) * scale
      off <- off + 2
    }
  }
  A
}
