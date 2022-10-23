#' Constructors for `bru_mapper` objects
#'
#' @details
#' * `bru_mapper` Generic mapper S3 constructor, used for constructing
#' mappers for special objects. See below for details of the
#' default constructor [bru_mapper_define()] that can be used to define
#' new mappers in user code.
#' @param \dots Arguments passed on to other methods
#' @export
#' @seealso [bru_mapper_methods] for specific method implementations, and
#' [bru_get_mapper] for hooks to extract mappers from latent model object
#' class objects.
#' @rdname bru_mapper
#' @examples
#' mapper <- bru_mapper_index(5)
#' ibm_jacobian(mapper, input = c(1, 3, 4, 5, 2))
bru_mapper <- function(...) {
  UseMethod("bru_mapper")
}

#' Methods for bru_mapper objects
#'
#' @details
#' * `ibm_n` Generic. Implementations must return the size of the latent vector
#' being mapped to.
#' @param mapper A mapper S3 object, inheriting from `bru_mapper`.
#' For the `bru_mapper_define` method, instead a
#' list that will be converted to a `bru_mapper` object by adding
#' class information and (optional) methods.
#' @param inla_f logical; when `TRUE` for `ibm_n()` and `ibm_values()`, the
#' result must be compatible with the `INLA::f(...)` and corresponding
#' `INLA::inla.stack(...)` constructions.  For `ibm_{eval,offset,jacobian,linear,amatrix}`,
#' the `input` interpretation may be different.
#' Implementations do not normally need to do anything different, except
#' for mappers of the type needed for hidden multicomponent models such
#' as "bym2", which can be handled by `bru_mapper_collect`.
#'
#' @export
#' @rdname bru_mapper_methods
#' @name bru_mapper_methods
ibm_n <- function(mapper, inla_f = FALSE, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_n"]])) {
    mapper[[".envir"]][["ibm_n"]](mapper, inla_f = inla_f, ...)
  } else {
    UseMethod("ibm_n")
  }
}

#' @details
#' * `ibm_n_output` Generic.
#' Implementations must return an integer denoting the
#' mapper output length, for valid inputs as determined by
#' `ibm_valid_input(mapper, input, inla_f = inla_f)`.
#' The default implementation returns `NROW(input)`.
#' Mappers such as `bru_mapper_multi` and `bru_mapper_collect`,
#' that can accept `list()` inputs require their own methods implementations.
#' @export
#' @rdname bru_mapper_methods
ibm_n_output <- function(mapper, input, inla_f = FALSE, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_n"]])) {
    mapper[[".envir"]][["ibm_n"]](mapper, input, inla_f = inla_f, ...)
  } else {
    UseMethod("ibm_n_output")
  }
}



#' @details
#' * `ibm_values` Generic. When `inla_f=TRUE`, implementations must return a vector that
#' would be interpretable by an `INLA::f(..., values = ...)` specification.
#' The exception is the method for `bru_mapper_multi`, that returns a
#' multi-column data frame.
#' @export
#' @rdname bru_mapper_methods
ibm_values <- function(mapper, inla_f = FALSE, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_values"]])) {
    mapper[[".envir"]][["ibm_values"]](mapper, inla_f = inla_f, ...)
  } else {
    UseMethod("ibm_values")
  }
}
#' @details
#' * `ibm_amatrix` Generic, will become deprecated in 2.7.0. Use `ibm_jacobian`
#' instead.
#' Implementations must return a (sparse) matrix of size `ibm_n_output(...)`
#' by `ibm_n(...)`. The `inla_f=TRUE` argument should only affect
#' the allowed type of input format.
#' @export
#' @rdname bru_mapper_methods
ibm_amatrix <- function(mapper, input, state = NULL, inla_f = FALSE, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_amatrix"]])) {
    mapper[[".envir"]][["ibm_amatrix"]](
      mapper, input, state = state, inla_f = inla_f, ...)
  } else {
    UseMethod("ibm_amatrix")
  }
}

#' @details
#' * `ibm_is_linear` Generic.
#' Implementations must return `TRUE` or `FALSE`.
#' If `TRUE` (returned by the default method unless the mapper
#' contains an `is_linear` variable), users of the mapper
#' may assume the mapper is linear.
#' @export
#' @rdname bru_mapper_methods
ibm_is_linear <- function(mapper, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_is_linear"]])) {
    mapper[[".envir"]][["ibm_is_linear"]](
      mapper = mapper, ...)
  } else {
    UseMethod("ibm_is_linear")
  }
}

#' @details
#' * `ibm_jacobian` Generic.
#' Implementations must return a (sparse) matrix of size
#' `ibm_n_output(mapper, input, inla_f)`
#' by `ibm_n(mapper, inla_f = FALSE)`. The `inla_f=TRUE` argument should
#' only affect the allowed type of input format.
#' @param input Data input for the mapper.
#' @export
#' @rdname bru_mapper_methods
ibm_jacobian <- function(mapper, input, state = NULL, inla_f = FALSE, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_jacobian"]])) {
    mapper[[".envir"]][["ibm_jacobian"]](
      mapper = mapper, input = input, state = state, inla_f = inla_f, ...)
  } else {
    UseMethod("ibm_jacobian")
  }
}


#' @details
#' * `ibm_linear()` Generic.
#' Implementations must return a [bru_mapper_taylor] object
#' The linearisation information includes `offset`, `jacobian`, and `state0`.
#' The state information indicates for which state the `offset` was evaluated,
#' with `NULL` meaning all-zero.
#' The linearised mapper output is defined as
#' `effect(input, state) = offset(input, state0) + jacobian(input, state0) %*% (state - state0)`.
#' The default method calls `ibm_eval()` and `ibm_jacobian()` to generate
#' the needed information.
#' @export
#' @rdname bru_mapper_methods
ibm_linear <- function(mapper, input, state = NULL, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_linear"]])) {
    mapper[[".envir"]][["ibm_linear"]](
      mapper, input, state = state, ...)
  } else {
    UseMethod("ibm_linear")
  }
}

#' @details
#' * `ibm_eval` Generic.
#' Implementations must return a vector of length `ibm_n_output(...)`.
#' The `input` contents must
#' be in a format accepted by `ibm_jacobian(...)`
#' for the mapper.
#' @param state A vector of latent state values for the mapping,
#' of length `ibm_n(mapper, inla_f = FALSE)`
#' @export
#' @rdname bru_mapper_methods
ibm_eval <- function(mapper, input, state = NULL, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_eval"]])) {
    mapper[[".envir"]][["ibm_eval"]](
      mapper, input, state = state, ...)
  } else {
    UseMethod("ibm_eval")
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
#' @rdname bru_mapper_methods
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
#' Implementations should return a logical vector of length `NROW(input)`,
#' or for `bru_mapper_multi` and `bru_mapper_collect` a list of such
#' vectors.
#' @export
#' @rdname bru_mapper_methods
ibm_valid_input <- function(mapper, input, inla_f = FALSE, ...) {
  if (!is.null(mapper[[".envir"]][["ibm_valid_input"]])) {
    mapper[[".envir"]][["ibm_valid_input"]](
      mapper, input, inla_f = inla_f, ...)
  } else {
    UseMethod("ibm_valid_input")
  }
}


# MAPPERS ----

#' @details * `bru_mapper_define` adds the `new_class` and "bru_mapper" class
#' names to the inheritance list for the input `mapper` object, unless the object
#' already inherits from these.
#' If provided, mapper method functions are added to an environment
#' `.envir` in the object.  The generic methods look for these
#' functions first,
#' and otherwise call `UseMethod()`.  This is an alternative to using `.S3method()`
#' to register the methods, e.g.
#' `.S3method("ibm_jacobian", "my_mapper_class", ibm_jacobian.my_mapper_class)`.
#' @param mapper For `bru_mapper_define`, a prototype mapper object, see Details.
#' For `bru_mapper_scale`, a mapper to be scaled.
#' @param new_class If non-`NULL`, this is added at the front of the class definition
#' @param methods Optional `list` of named method definitions; See Details.
#' @param \dots Deprecated, alternative way to supply optional method definitions.
#'
#' @export
#' @rdname bru_mapper
bru_mapper_define <- function(mapper,
                              new_class = NULL,
                              methods = NULL,
                              ...) {
  valid_method_names <-
    c(
      "ibm_n",
      "ibm_n_output",
      "ibm_values",
      "ibm_amatrix",
      "ibm_is_linear",
      "ibm_jacobian",
      "ibm_linear",
      "ibm_eval",
      "ibm_inla_subset",
      "ibm_valid_input"
    )
  if (any(c(
    valid_method_names,
    "ibm_n_inla",
    "ibm_values_inla",
    "ibm_amatrix_inla",
    "ibm_valid_input_inla"
  ) %in%
    names(list(...)))) {
    warning(
      paste0(
        "Deprecated use of named method arguments for 'bru_mapper'.\n",
        "Use methods = list(methodname = ..., ...) instead."
      )
    )
    method_names <-
      setdiff(
        intersect(
          valid_method_names,
          names(list(...))
        ),
        names(methods)
      )
    methods <- c(list(...)[method_names], methods)
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
    deprecated_methods <- setdiff(names(methods),
                                  valid_method_names)
    if (length(deprecated_methods) > 0) {
      warning(
        paste0(
          "Unknown bru_mapper method names detected, ignoring methods\n",
          paste0(deprecated_methods, collapse = ", "),
          "."
        )
      )
    }
    for (method in setdiff(names(methods), deprecated_methods)) {
      if (!is.null(methods[[method]])) {
        assign(method,
               methods[[method]],
               envir = mapper[[".envir"]])
      }
    }
  }
  mapper
}

#' @details * `bru_mapper.default` calls `bru_mapper_define`, passing all
#' arguments along. Mapper implementations should call [bru_mapper_define()]
#' instead, and supply at least a `new_class` class name.
#' Use of the `bru_mapper.default` will be deprecated from version 2.7.0.
#' @export
#' @rdname bru_mapper
bru_mapper.default <- function(...) {
  # TODO: Mark deprecated from version 2.7.0
  # .Deprecated("bru_mapper_define")
  bru_mapper_define(...)
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
#' @param inla_f logical; when `TRUE` in `ibm_n` and `ibm_values`,
#' these must result in values compatible with `INLA::f(...)`
#' an specification and corresponding `INLA::inla.stack(...)` constructions.
#' For the `ibm_eval` and `ibm_jacobian` methods, it may influence how the
#' input data is interpreted.
#' Implementations do not normally need to do anything different, except
#' for mappers of the type needed for hidden multicomponent models such
#' as "bym2", which can be handled by `bru_mapper_collect`.
#' @seealso [bru_mapper] for constructor methods, and
#' [bru_get_mapper] for hooks to extract mappers from latent model object
#' class objects.
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

#' @export
#' @rdname bru_mapper_methods
ibm_n_output.default <- function(mapper, input, inla_f = FALSE, ...) {
  NROW(input)
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
#' * The default `ibm_amatrix()` method gives an error message.
#' Mapper classes must implement their own `ibm_jacobian` or
#' `ibm_amatrix` methods. New implementations should use
#' a `ibm_jacobian` method. `ibm_amatrix` may become deprecated
#' in a future version.
#' @export
#' @rdname bru_mapper_methods
ibm_amatrix.default <- function(mapper, ...) {
  stop(paste0(
    "Missing implementation of 'ibm_jacobian()/ibm_amatrix()' for mapper of class '",
    paste0(class(mapper), collapse = ", "), "'.\n",
    "New implementations should implement a 'ibm_jacobian()' method."
  ))
}

#' @details
#' * The default `ibm_is_linear()` method returns logical
#' 'is_linear' from the mapper object if it exists, and otherwise `TRUE`.
#' @export
#' @rdname bru_mapper_methods
ibm_is_linear.default <- function(mapper, ...) {
  if (!is.null(mapper[["is_linear"]])) {
    mapper[["is_linear"]]
  } else {
    TRUE
  }
}

#' @details
#' * The default `ibm_jacobian()` calls `ibm_amatrix`, which
#' by default gives an error.
#' Mapper classes should implement their own `ibm_jacobian` method.
#' @export
#' @rdname bru_mapper_methods
ibm_jacobian.default <- function(mapper, input, state, ...) {
  if (!ibm_is_linear(mapper)) {
    stop(paste0(
      "Non-linear mappers must implement their own 'ibm_jacobian()' method.",
      " Missing method for class '",
      class(mapper)[1], "'."
    ))
  }
  ibm_amatrix(mapper, input = input, state = state, ...)
}

#' @export
#' @rdname bru_mapper_methods
ibm_offset.default <- function(mapper, input, state = NULL, ...) {
  stop("ibm_offset.default method has been removed")
}

#' @details
#' * The default `ibm_linear()` method calls `ibm_eval()` and `ibm_jacobian()`
#' and returns a `bru_mapper_taylor` object.
#' The `state0` information in the affine mapper indicates for which state
#' the `offset` was evaluated; The affine mapper output is defined as
#' `effect(input, state) = offset(input, state0) + jacobian(input, state0) %*% (state - state0)`
#' @export
#' @rdname bru_mapper_methods
ibm_linear.default <- function(mapper, input, state, ...) {
  bru_mapper_taylor(
    offset = ibm_eval(mapper, input, state, ...),
    jacobian = ibm_jacobian(mapper, input, state, ...),
    state0 = state,
    values_mapper = mapper
  )
}

#' @details
#' * The default `ibm_eval()` method verifies that the mapper is linear
#' with `ibm_is_linear()`, and then computes a linear mapping
#' as `ibm_jacobian(...) %*% state`.  When `state` is `NULL`,
#' a zero vector of length `ibm_n_output(...)` is returned.
#' @export
#' @rdname bru_mapper_methods
ibm_eval.default <- function(mapper, input, state = NULL, ...) {
  if (!ibm_is_linear(mapper)) {
    stop("Non-linear mappers must implement their own ibm_eval() method.")
  }

  val <- numeric(ibm_n_output(mapper, input, ...))

  if ((ibm_n(mapper) > 0) && !is.null(state)) {
    A <- ibm_jacobian(mapper, input = input, state = state, ...)
    val + A %*% state
  }

  val
}


#' @details
#' * The default `ibm_inla_subset` method uses
#' the `ibm_values` output to construct the inla subset indexing, passing
#' extra arguments such as `multi` on to the methods (this means it supports
#' both regular vector values and `multi=1` data.frame values).
#' @export
#' @rdname bru_mapper_methods
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
  bru_mapper_define(mapper, new_class = "bru_mapper_inla_mesh_2d")
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
#' @export
#' @rdname bru_mapper_methods
ibm_jacobian.bru_mapper_inla_mesh_2d <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  }
  if (inherits(input, "sfc_POINT")) {
    # TODO: Add direct sf support to inla.spde.make.A,
    # to support crs passthrough.
    A <- sf::st_coordinates(input)
    nm <- intersect(colnames(A), c("X", "Y", "Z"))
    input <- as.matrix(A[, nm, drop = FALSE])
  } else if (!is.matrix(input) && !inherits(input, "Spatial")) {
    input <- as.matrix(input)
  }
  INLA::inla.spde.make.A(mapper[["mesh"]], loc = input)
}

#' @export
#' @rdname bru_mapper_methods
ibm_jacobian.bru_mapper_inla_mesh_2d <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  }
  if (!is.matrix(input) && !inherits(input, "Spatial")) {
    input <- as.matrix(input)
  }
  INLA::inla.spde.make.A(mapper[["mesh"]], loc = input)
}

#' @export
#' @describeIn inlabru-deprecated Replaced by [ibm_jacobian()]
ibm_amatrix.bru_mapper_inla_mesh_2d <- function(...) {
  .Deprecated("ibm_jacobian")
  ibm_jacobian(...)
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
  bru_mapper_define(mapper, new_class = "bru_mapper_inla_mesh_1d")
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
ibm_jacobian.bru_mapper_inla_mesh_1d <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  }
  INLA::inla.spde.make.A(mapper[["mesh"]], loc = input)
}

#' @export
#' @describeIn inlabru-deprecated Replaced by [ibm_jacobian()]
ibm_amatrix.bru_mapper_inla_mesh_1d <- function(...) {
  .Deprecated("ibm_jacobian")
  ibm_jacobian(...)
}

## _index ####

#' @param n Size of a model for `bru_mapper_index`
#' @export
#' @rdname bru_mapper
bru_mapper_index <- function(n = 1L, ...) {
  bru_mapper_define(list(n = n), new_class = "bru_mapper_index")
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
ibm_jacobian.bru_mapper_index <- function(mapper, input, ...) {
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

#' @export
#' @describeIn inlabru-deprecated Replaced by [ibm_jacobian()]
ibm_amatrix.bru_mapper_index <- function(...) {
  .Deprecated("ibm_jacobian")
  ibm_jacobian(...)
}

## _taylor ####

#' @details `bru_mapper_taylor` provides a pre-computed affine mapping,
#' internally used to represent and evaluate linearisation information.
#' The `state0` information indicates for which state the `offset` was evaluated;
#' The affine mapper output is defined as
#' `effect(state) = offset + jacobian %*% (state - state0)`
#' @param offset For `bru_mapper_taylor`, an offset vector evaluated
#' at `state0`
#' @param jacobian For `bru_mapper_taylor`, the Jacobian matrix,
#' evaluated at `state0`, or, a named list of such matrices.
#' May be `NULL` or an empty list, for a constant mapping.
#' @param state0 For `bru_mapper_taylor`, the state the linearisation
#' was evaluated at, or a list of length matching the `jacobian` list.
#' `NULL` is interpreted as 0.
#' @param values_mapper mapper object to be used for `ibm_n` and
#' `ibm_values` for `inla_f=TRUE` (experimental, currently unused)
#' @export
#' @rdname bru_mapper
bru_mapper_taylor <- function(offset, jacobian, state0, ...,
                              values_mapper = NULL) {
  if (is.null(state0)) {
    n_state <- 0
  } else if (is.list(state0)) {
      n_state <- vapply(state0,
                        function(x) {
                          ifelse(is.null(x),
                                 0L,
                                 length(x))
                        },
                        0L)
  } else {
    n_state <- length(state0)
  }

  if (is.null(jacobian)) {
    n_jacobian <- 0
  } else if (is.list(jacobian)) {
    n_jacobian <- vapply(jacobian,
                         function(x) {
                           ifelse(is.null(x),
                                  0L,
                                  ncol(x))
                         },
                         0L)
  } else {
    n_jacobian <- ncol(jacobian)
  }

  if (!is.null(state0)) {
    if ((sum(n_jacobian) == 0) && (sum(n_state) > 0)) {
      stop("Jacobian information indicates an empty state, but the state information has length > 0")
    }
    stopifnot(sum(n_jacobian) == sum(n_state))
  }

  n_multi <- n_jacobian
  if (is.null(jacobian)) {
    state0 <- NULL
    n_state <- 0
  } else if (!is.list(jacobian) && is.list(state0)) {
    state0 <- unlist(state0)
    n_state <- length(state0)
  } else if (is.list(jacobian) && !is.null(state0) && !is.list(state0)) {
    # Split vector
    state0 <- lapply(
      seq_along(n_multi),
      function(k) {
        idx <- sum(n_multi[seq_len(k - 1)])
        idx <- idx + seq_len(n_multi[k]) - 1
        state0[idx]
      })
    n_state <- n_multi
  }
  if (!is.null(state0)) {
    stopifnot(all(n_jacobian == n_state))
  }

  bru_mapper_define(list(offset = offset,
                         jacobian = jacobian,
                         state0 = state0,
                         n_multi = n_multi,
                         n = sum(n_multi),
                         noutput = length(offset),
                         values_mapper = NULL),
                    # TODO: maybe allow values_mapper
                    new_class = "bru_mapper_taylor")
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_taylor <- function(mapper, inla_f = FALSE, multi = FALSE, ...) {
  if (inla_f && !is.null(mapper[["values_mapper"]])) {
    stop("ibm_n.bru_mapper_taylor should not be used with inla_f = TRUE")

    ibm_n(mapper[["values_inla"]], inla_f = inla_f, multi = multi)
  } else if (multi) {
    mapper[["n_multi"]]
  } else {
    mapper[["n"]]
  }
}

#' @export
#' @rdname bru_mapper_methods
ibm_n_output.bru_mapper_taylor <- function(mapper, input, inla_f = FALSE, ...) {
  mapper[["n_output"]]
}


#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_taylor <- function(mapper,
                                         inla_f = FALSE,
                                         multi = FALSE,
                                         ...) {
  if (inla_f && !is.null(mapper[["values_mapper"]])) {
    stop("ibm_values.bru_mapper_taylor should not be used with inla_f = TRUE")

    ibm_values(mapper[["values_inla"]], inla_f = inla_f, multi = multi)
  } else if (multi) {
    lapply(mapper[["n_multi"]], function(k) seq_len(k))
  } else {
    seq_len(mapper[["n"]])
  }
}

#' @rdname bru_mapper_methods
ibm_jacobian.bru_mapper_taylor <- function(mapper, ..., multi = FALSE) {
  if (is.null(mapper[["jacobian"]])) {
    return(Matrix::Matrix(0, mapper[["n_output"]], 0))
  }
  if (multi) {
    mapper[["jacobian"]]
  } else if (is.list(mapper[["jacobian"]])) {
    do.call(cbind, mapper[["jacobian"]])
  } else {
    mapper[["jacobian"]]
  }
}


#' @details
#' * The `ibm_eval.bru_mapper_taylor()` evaluates linearised
#' mapper information at the given `state`. The `input` argument is ignored,
#' so that the usual argument order
#' `ibm_eval(mapper, input, state)` syntax can be used, but also
#' `ibm_eval(mapper, state = state)`.  For a mapper with a named jacobian list,
#' the `state` argument must also be a named list.  If `state` is `NULL`,
#' all-zero is assumed.
#' @export
#' @rdname bru_mapper_methods
ibm_eval.bru_mapper_taylor <- function(mapper, input = NULL, state = NULL, ...) {
  if (is.null(mapper[["jacobian"]]) ||
      (mapper[["n"]] == 0) ||
      (is.null(state) && is.null(mapper[["state0"]]))) {
    mapper[["offset"]]
  } else if (is.list(mapper[["jacobian"]])) {
    stopifnot(is.null(state) || is.list(state))
    val <- mapper[["offset"]]
    if (is.null(mapper[["state0"]])) {
      for (nm in names(mapper[["jacobian"]])) {
        val <- val + mapper[["jacobian"]][[nm]] %*% state[[nm]]
      }
    } else if (is.null(state)) {
      for (nm in names(mapper[["jacobian"]])) {
        val <- val - mapper[["jacobian"]][[nm]] %*% mapper[["state0"]][[nm]]
      }
    } else {
      for (nm in names(mapper[["jacobian"]])) {
        val <- val + mapper[["jacobian"]][[nm]] %*%
          (state[[nm]] - mapper[["state0"]][[nm]])
      }
    }
    val
  } else {
    stopifnot(is.null(state) || !is.list(state))
    if (is.null(mapper[["state0"]])) {
      mapper[["offset"]] + mapper[["jacobian"]] %*% state
    } else if (is.null(state)) {
      mapper[["offset"]] - mapper[["jacobian"]] %*% mapper[["state0"]]
    } else {
      mapper[["offset"]] + mapper[["jacobian"]] %*% (state - mapper[["state0"]])
    }
  }
}


## _linear ####

#' @export
#' @rdname bru_mapper
bru_mapper_linear <- function(...) {
  bru_mapper_define(list(), new_class = "bru_mapper_linear")
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
ibm_jacobian.bru_mapper_linear <- function(mapper, input, ...) {
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

#' @export
#' @describeIn inlabru-deprecated Replaced by [ibm_jacobian()]
ibm_amatrix.bru_mapper_linear <- function(...) {
  .Deprecated("bru_jacobian")
  ibm_jacobian(...)
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
  bru_mapper_define(mapper, new_class = "bru_mapper_matrix")
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
ibm_jacobian.bru_mapper_matrix <- function(mapper, input, state = NULL,
                                           inla_f = FALSE, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  } else if (is.matrix(input)) {
    A <- as(input, "Matrix")
  } else if (inherits(input, "Matrix")) {
    A <- input
  } else if (inherits(input, "Spatial")) {
    A <- sp::coordinates(input)
    A <- as(A, "Matrix")
  } else if (inherits(input, "sfc_POINT")) {
    A <- sf::st_coordinates(input)
    nm <- intersect(colnames(A), c("X", "Y", "Z", "M"))
    A <- as(A[, nm, drop = FALSE], "Matrix")
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
  colnames(A) <- mapper$labels
  A
}

#' @export
#' @describeIn inlabru-deprecated Replaced by [ibm_jacobian()]
ibm_amatrix.bru_mapper_matrix <- function(...) {
  .Deprecated("bru_jacobian")
  ibm_jacobian(...)
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
  bru_mapper_define(mapper, new_class = "bru_mapper_factor")
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
ibm_jacobian.bru_mapper_factor <- function(mapper, input, ...) {
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


#' @export
#' @describeIn inlabru-deprecated Replaced by [ibm_jacobian()]
ibm_amatrix.bru_mapper_factor <- function(...) {
  .Deprecated("bru_jacobian")
  ibm_jacobian(...)
}



## _const ####

#' @export
#' @rdname bru_mapper
bru_mapper_const <- function(...) {
  bru_mapper_define(list(), new_class = "bru_mapper_const")
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_const <- function(mapper, ...) {
  0L
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_const <- function(mapper, ...) {
  NULL
}

#' @export
#' @rdname bru_mapper_methods
ibm_jacobian.bru_mapper_const <- function(mapper, input, ...) {
  Matrix::Matrix(0, ibm_n_output(mapper, input, ...), 0L)
}

#' @export
#' @rdname bru_mapper_methods
ibm_eval.bru_mapper_const <- function(mapper, input, state = NULL, ...) {
  if (is.null(input)) {
    return(numeric(0))
  }
  input <- as.vector(input)
  ok <- !is.na(input)
  input[!ok] <- 0
  input
}


## _offset ####

#' @export
#' @rdname bru_mapper
bru_mapper_offset <- function(...) {
#  .Deprecated("bru_mapper_const")
  bru_mapper_define(bru_mapper_const(), new_class = "bru_mapper_offset")
}

#' @export
#' @describeIn inlabru-deprecated Replaced by [bru_mapper_const] methods
ibm_n.bru_mapper_offset <- function(...) {
  .Deprecated("bru_mapper_const")
  NextMethod()
}
#' @export
#' @describeIn inlabru-deprecated Replaced by [bru_mapper_const] methods
ibm_values.bru_mapper_offset <- function(...) {
  .Deprecated("bru_mapper_const")
  NextMethod()
}

#' @export
#' @describeIn inlabru-deprecated Replaced by [bru_mapper_const] methods
ibm_amatrix.bru_mapper_offset <- function(...) {
  .Deprecated("bru_mapper_const")
  NextMethod()
}



## _scale ####

#' @export
#' @details For `bru_mapper_scale()`, `mapper` is a mapper to be scaled.
#' The `input` format for the `ibm_eval` and `ibm_jacobian` methods is
#' `list(mapper = input_to_the_inner_mapper, scale = scaling_weights)`
#' @rdname bru_mapper
bru_mapper_scale <- function(mapper, ...) {
  bru_mapper_define(list(mapper = mapper,
                         is_linear = ibm_is_linear(mapper)),
                    new_class = "bru_mapper_scale")
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_scale <- function(mapper, ...) {
  ibm_n(mapper[["mapper"]], ...)
}
#' @export
#' @rdname bru_mapper_methods
ibm_n_output.bru_mapper_scale <- function(mapper, input, ...) {
  ibm_n_output(mapper[["mapper"]], input, ...)
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_scale <- function(mapper, ...) {
  ibm_values(mapper[["mapper"]], ...)
}

#' @export
#' @param sub_lin Internal, optional pre-computed sub-mapper information
#' @details For `bru_mapper_scale`, `input` values without a `scale` element
#' are interpreted as no scaling.
#' @rdname bru_mapper_methods
ibm_jacobian.bru_mapper_scale <- function(mapper, input, state = NULL, ...,
                                         sub_lin = NULL) {
  if (is.null(sub_lin)) {
    A <- ibm_jacobian(mapper[["mapper"]],
                      input = input[["mapper"]],
                      state = state, ...)
  } else {
    A <- sub_lin$jacobian
  }
  if ((NROW(A) > 0) && !is.null(input[["scale"]])) {
    scale <- as.vector(input[["scale"]])
    ok <- !is.na(scale)
    scale[!ok] <- 0
    # Scale each row of A
    scale * A
  } else {
    A
  }
}

#' @export
#' @param sub_lin Internal, optional pre-computed sub-mapper information
#' @rdname bru_mapper_methods
ibm_offset.bru_mapper_scale <- function(mapper, input, ..., sub_lin = NULL) {
  off <- ibm_eval(mapper[["mapper"]], input = input[["mapper"]], ...,
                  lin = sub_lin)
  if ((NROW(off) > 0) && !is.null(input[["scale"]])) {
    scale <- as.vector(input[["scale"]])
    ok <- !is.na(scale)
    scale[!ok] <- 0
    # Scale each row of A
    scale * off
  } else {
    off
  }
}


#' @export
#' @rdname bru_mapper_methods
ibm_linear.bru_mapper_scale <- function(mapper, input, state, ...) {
  sub_lin <- ibm_linear(mapper[["mapper"]],
                        input[["mapper"]],
                        state = state, ...)
  bru_mapper_taylor(
    offset = ibm_eval(mapper, input, state, ..., sub_lin = sub_lin),
    jacobian = ibm_jacobian(mapper, input, state, ..., sub_lin = sub_lin),
    state0 = state,
    values_mapper = mapper
  )
}


#' @export
#' @rdname bru_mapper_methods
ibm_eval.bru_mapper_scale <- function(mapper, input, state = NULL, ...,
                                      sub_lin = NULL) {
  if (!is.null(sub_lin)) {
    values <- ibm_eval(sub_lin, input = NULL, state = state)
  } else {
    values <- ibm_eval(mapper[["mapper"]],
                       input = input[["mapper"]],
                       state = state, ...)
  }
  if ((NROW(values) > 0) && !is.null(input[["scale"]])) {
    scale <- as.vector(input[["scale"]])
    ok <- !is.na(scale)
    scale[!ok] <- 0
    scale * values
  } else {
    values
  }
}


#' @details
#' * `ibm_valid_input` for `bru_mapper_scale` accepts a list with
#' named entries `mapper` and `scale`. The contents of the `mapper`
#' element is checked for validity for the submapper
#' with `ibm_valid_input(mapper$mapper, input$mapper, ...)`
#' @export
#' @rdname bru_mapper_methods
ibm_valid_input.bru_mapper_scale <- function(mapper, input, ...) {
  ibm_valid_input(mapper[["mapper"]], input[["mapper"]], ...)
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
    n_multi = lapply(mappers, ibm_n),
    n_inla_multi = lapply(mappers, function(x) {
      ibm_n(x, inla_f = TRUE)
    }),
    values_multi = lapply(mappers, ibm_values),
    values_inla_multi = lapply(mappers, function(x) {
      ibm_values(x, inla_f = TRUE)
    }),
    is_linear_multi = lapply(mappers, ibm_is_linear)
  )
  mapper[["n"]] <- prod(unlist(mapper[["n_multi"]]))
  mapper[["n_inla"]] <- prod(unlist(mapper[["n_inla_multi"]]))
  mapper[["is_linear"]] <- all(unlist(mapper[["is_linear_multi"]]))

  if (!mapper[["is_linear"]]) {
    stop("bru_mapper_multi sub-mappers must be linear mappers")
  }
  bru_mapper_define(mapper, new_class = "bru_mapper_multi")
}

#' @param multi logical;
#' If `TRUE` (or positive), recurse one level into sub-mappers
#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_multi <- function(mapper, inla_f = FALSE, multi = FALSE, ...) {
  if (multi) {
    if (inla_f) {
      mapper[["n_inla_multi"]]
    } else {
      mapper[["n_multi"]]
    }
  } else if (inla_f) {
    mapper[["n_inla"]]
  } else {
    mapper[["n"]]
  }
}

#' @export
#' @rdname bru_mapper_methods
ibm_n_output.bru_mapper_multi <- function(mapper, input, ...) {
  # Assume that the first mapper fully handles the output size
  if (is.matrix(input)) {
    ibm_n_output(mapper[["mapper"]][[1]], input[, 1, drop = TRUE], ...)
  } else {
    stopifnot(is.list(input)) # data.frame or list
    ibm_n_output(mapper[["mapper"]][[1]], input[[1]], ...)
  }
}


#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_multi <- function(mapper, inla_f = FALSE, multi = FALSE, ...) {
  if (multi) {
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
  } else if (inla_f) {
    seq_len(mapper[["n_inla"]])
  } else {
    seq_len(mapper[["n"]])
  }
}

#' @export
#' @rdname bru_mapper_methods
ibm_is_linear.bru_mapper_multi <- function(mapper, multi = FALSE, ...) {
  if (multi) {
    mapper[["is_linear_multi"]]
  } else {
    mapper[["is_linear"]]
  }
}

bru_mapper_multi_prepare_input <- function(mapper, input) {
  if (is.matrix(input)) {
    input_names <- colnames(input)
    input <- as.data.frame(input)
  } else {
    input_names <- names(input)
  }
  if (is.null(input_names)) {
    if (is.data.frame(input)) {
      names(input) <- names(mapper[["mappers"]])[seq_len(ncol(input))]
    } else {
      # input should now be a list
      names(input) <- names(mapper[["mappers"]])[seq_along(input)]
    }
  }
  input
}

#' @details
#' * `ibm_jacobian` for `bru_mapper_multi` accepts a list with
#' named entries, or a list with unnamed but ordered elements.
#' The names must match the sub-mappers, see [names.bru_mapper_multi()].
#' Each list element should take a format accepted by the corresponding
#' sub-mapper. In case each element is a vector, the input can be given as a
#' data.frame with named columns, a matrix with named columns, or a matrix
#' with unnamed but ordered columns.
#' @param sub_A Internal; precomputed Jacobian matrices.
#' @export
#' @rdname bru_mapper_methods
ibm_jacobian.bru_mapper_multi <- function(mapper,
                                          input,
                                          state = NULL,
                                          inla_f = FALSE,
                                          multi = FALSE,
                                          ...,
                                          sub_A = NULL) {
  input <- bru_mapper_multi_prepare_input(mapper, input)

  # Note: the multi-mapper only handles linear sub-mappers.
  if (is.null(sub_A)) {
    sub_A <-
      lapply(
        names(mapper[["mappers"]]),
        function(x) {
          ibm_jacobian(mapper[["mappers"]][[x]],
                       input = input[[x]],
                       state = NULL,
                       inla_f = inla_f,
                       multi = FALSE
          )
        }
      )
  }
  if (multi) {
    return(sub_A)
  }
  # Combine the matrices
  # (A1, A2, A3) -> rowkron(A3, rowkron(A2, A1))
  A_ <- sub_A[[1]]
  for (k in seq_len(length(mapper[["mappers"]]) - 1)) {
    A_ <- row_kron(sub_A[[k + 1]], A_)
  }
  return(A_)
}

#' @export
#' @describeIn inlabru-deprecated Replaced by [ibm_jacobian()]
ibm_amatrix.bru_mapper_multi <- function(...) {
  .Deprecated("bru_jacobian")
  ibm_jacobian(...)
}



#' @param pre_A Internal; precomputed Jacobian matrix
#' @export
#' @rdname bru_mapper_methods
ibm_linear.bru_mapper_multi <- function(mapper, input, state,
                                        inla_f = FALSE,
                                        ...) {
  input <- bru_mapper_multi_prepare_input(mapper, input)
  A <- ibm_jacobian(mapper, input, state = state,
                    inla_f = inla_f, multi = FALSE,
                    ...)
  bru_mapper_taylor(
    offset = ibm_eval(mapper, input, state,
                      inla_f = inla_f, multi = FALSE, ...,
                      pre_A = A),
    jacobian = A,
    state0 = state,
    values_mapper = mapper
  )
}


#' @export
#' @rdname bru_mapper_methods
ibm_eval.bru_mapper_multi <- function(mapper, input, state = NULL,
                                      inla_f = FALSE, ..., pre_A = NULL) {
  input <- bru_mapper_multi_prepare_input(mapper, input)
  if ((ibm_n(mapper) == 0) || is.null(state)) {
    # Handle the case when the mapper is a _const mapper
    val <- ibm_eval(mapper[["mappers"]][[1]],
                    input = input[[1]],
                    state = NULL,
                    inla_f = inla_f,
                    multi = FALSE,
                    ...)
  } else {
    if (is.null(pre_A)) {
      pre_A <- ibm_jacobian(mapper, input = input, state = state,
                            inla_f = inla_f, multi = FALSE, ...)
    }
    val <- pre_A %*% state
  }
  val
}


bm_multi_indexing <- function(mapper, input) {
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
                                             inla_f = FALSE, multi = FALSE, ...) {
  indexing <- bm_multi_indexing(mapper, input)
  if (is.matrix(input)) {
    validity <-
      lapply(
        indexing,
        function(x) {
          ibm_valid_input(mapper[["mappers"]][[x]],
                          input = input[, x],
                          inla_f = inla_f,
                          multi = FALSE
          )
        }
      )
  } else if (is.list(input)) {
    validity <-
      lapply(
        indexing,
        function(x) {
          ibm_valid_input(mapper[["mappers"]][[x]],
                          input = input[[x]],
                          inla_f = inla_f,
                          multi = FALSE
          )
        }
      )
  } else {
    validity <- as.list(rep(FALSE, length(mapper[["mappers"]])))
    names(validity) <- names(mapper[["mappers"]])
  }
  if (!multi) {
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
#' `ibm_values`, and `ibm_jacobian`. Set to `FALSE` to always access the full
#' mapper, e.g. for `rgeneric` models
#' @details * `bru_mapper_collect` constructs concatenated collection mapping
#' @export
#' @rdname bru_mapper
bru_mapper_collect <- function(mappers, hidden = FALSE, ...) {
  mapper <- list(
    mappers = mappers,
    n_multi = lapply(mappers, ibm_n),
    values_multi = lapply(mappers, ibm_values),
    hidden = hidden,
    is_linear_multi = lapply(mappers, ibm_is_linear)
  )
  mapper[["n"]] <- sum(unlist(mapper[["n_multi"]]))
  mapper[["values"]] <- seq_len(mapper[["n"]])
  mapper[["is_linear"]] <- all(unlist(mapper[["is_linear_multi"]]))
  bru_mapper_define(mapper, new_class = "bru_mapper_collect")
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_collect <- function(mapper, inla_f = FALSE, multi = FALSE, ...) {
  if (multi) {
    mapper[["n_multi"]]
  } else if (mapper[["hidden"]] && inla_f) {
    mapper[["n_multi"]][[1]]
  } else {
    mapper[["n"]]
  }
}


# Only for inla_f = FALSE
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
#' @rdname bru_mapper_methods
ibm_n_output.bru_mapper_collect <- function(mapper, input,
                                            inla_f = FALSE,
                                            multi = FALSE, ...) {
  if (mapper[["hidden"]] && inla_f) {
    return(ibm_n_output(mapper[["mappers"]][[1]], input = input))
  }

  indexing <- bm_collect_indexing(mapper, input)
  if (is.matrix(input)) {
    n <- vapply(indexing,
                function(x) {
                  ibm_n_output(mapper[["mapper"]][[x]], input[, x], ...)
                },
                0L)
  } else {
    n <- vapply(indexing,
                function(x) {
                  ibm_n_output(mapper[["mapper"]][[x]], input[[x]], ...)
                },
                0L)
  }

  if (!multi) {
    n <- sum(n)
  }

  n
}


#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_collect <- function(mapper, inla_f = FALSE, multi = FALSE, ...) {
  if (multi) {
    mapper[["values_multi"]]
  } else if (mapper[["hidden"]] && inla_f) {
    mapper[["values_multi"]][[1]]
  } else {
    mapper[["values"]]
  }
}

#' @export
#' @rdname bru_mapper_methods
ibm_is_linear.bru_mapper_collect <- function(mapper,
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



#' @details
#' * `ibm_jacobian` for `bru_mapper_collect` accepts a list with
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
ibm_jacobian.bru_mapper_collect <- function(mapper, input, state = NULL,
                                           inla_f = FALSE, multi = FALSE,
                                           ...,
                                           sub_lin = NULL) {
  if (is.null(sub_lin)) {
    sub_lin <- bm_collect_sub_lin(mapper, input, state, inla_f = inla_f)
  }
  A <- lapply(sub_lin, function(x) x[["jacobian"]])

  if (!multi) {
    # Combine the matrices (A1, A2, A3, ...) -> bdiag(A1, A2, A3, ...)
    A_ <- Matrix::.bdiag(A)
    return(A_)
  }
  A
}

#' @export
#' @describeIn inlabru-deprecated Replaced by [ibm_jacobian()]
ibm_amatrix.bru_mapper_collect <- function(...) {
  .Deprecated("bru_jacobian")
  ibm_jacobian(...)
}


#' @export
#' @rdname bru_mapper_methods
ibm_eval.bru_mapper_collect <- function(mapper, input, state,
                                        inla_f = FALSE, multi = FALSE,
                                        ...,
                                        sub_lin = NULL) {
  if (is.null(sub_lin)) {
    sub_lin <- bm_collect_sub_lin(mapper, input, state, inla_f = inla_f)
  }
  val <- lapply(sub_lin, function(x) x[["offset"]])

  if (!multi) {
    # Combine the vectors (b1, b2, b3, ...) -> c(b1, b2, b3, ...)
    off <- do.call(c, off)
  }
  off
}


#' @export
#' @rdname bru_mapper_methods
ibm_linear.bru_mapper_collect <- function(mapper, input, state,
                                          inla_f = FALSE,
                                          ...) {
  if (mapper[["hidden"]] && inla_f) {
    input <- list(input)
  }
  sub_lin <-
    bm_collect_sub_lin(mapper, input, state,
                       inla_f = FALSE,
                       ...)
  A <- ibm_jacobian(mapper, input, state,
                    inla_f = FALSE, multi = FALSE, ...,
                    sub_lin = sub_lin)
  bru_mapper_taylor(
    offset = ibm_eval(mapper, input, state,
                      inla_f = FALSE, multi = FALSE, ...,
                      sub_lin = sub_lin),
    jacobian = A,
    state0 = state,
    values_mapper = mapper
  )
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
ibm_valid_input.bru_mapper_collect <- function(mapper, input, inla_f = FALSE, multi = FALSE, ...) {
  if (mapper[["hidden"]] && inla_f) {
    return(ibm_valid_input(mapper[["mappers"]][[1]],
      input = input,
      multi = FALSE
    ))
  }

  indexing <- bm_collect_indexing(mapper, input)
  if (is.matrix(input)) {
    validity <-
      lapply(
        indexing,
        function(x) {
          ibm_valid_input(mapper[["mappers"]][[x]],
                          input = input[, x],
                          multi = FALSE
          )
        }
      )
  } else if (is.list(input)) {
    validity <-
      lapply(
        indexing,
        function(x) {
          ibm_valid_input(mapper[["mappers"]][[x]],
                          input = input[[x]],
                          multi = FALSE
          )
        }
      )
  } else {
    validity <- as.list(rep(FALSE, length(mapper[["mappers"]])))
  }

  if (!multi) {
    # Combine the vectors (v1, v2, v3) -> c(v1, v2, v3)
    validity_ <- do.call(c, validity)
    return(validity_)
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
  bru_mapper_define(mapper, new_class = "bru_mapper_harmonics")
}


#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_harmonics <- function(mapper, inla_f = FALSE, ...) {
  mapper[["order"]] * 2 + mapper[["intercept"]]
}

#' @export
#' @rdname bru_mapper_methods
ibm_jacobian.bru_mapper_harmonics <- function(mapper, input, state = NULL, inla_f = FALSE, ...) {
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
