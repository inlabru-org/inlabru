# Generics ----

#' Constructors for `bru_mapper` objects
#' @name bru_mapper
NULL

#' @param \dots Arguments passed on to other methods
#' @export
#' @seealso [bru_mapper_generics] for generic methods,
#' [bru_mapper_methods] for specific method implementations, and
#' [bru_get_mapper] for hooks to extract mappers from latent model object
#' class objects.
#' @describeIn bru_mapper
#' Generic mapper S3 constructor, used for constructing
#' mappers for special objects. See below for details of the
#' default constructor [bru_mapper_define()] that can be used to define
#' new mappers in user code.
#'
#' @examples
#' mapper <- bru_mapper_index(5)
#' ibm_jacobian(mapper, input = c(1, 3, 4, 5, 2))
bru_mapper <- function(...) {
  UseMethod("bru_mapper")
}

#' Generic methods for bru_mapper objects
#'
#' @description
#' A `bru_mapper` sub-class implementation must provide an
#' `ibm_jacobian()` method. If the model size 'n' and definition
#' values 'values' are stored in the object itself, default methods are
#' available (see Details). Otherwise the
#' `ibm_n()` and `ibm_values()` methods also need to be provided.
#'
#' @seealso [bru_mapper] for constructor methods, and
#' [bru_get_mapper] for hooks to extract mappers from latent model object
#' class objects.
#'
#' @name bru_mapper_generics
#' @seealso [bru_mapper], [bru_mapper_methods]
NULL

#' Methods for bru_mapper objects
#'
#' @description
#' A `bru_mapper` sub-class implementation must provide an
#' `ibm_jacobian()` method. If the model size 'n' and definition
#' values 'values' are stored in the object itself, default methods are
#' available (see Details). Otherwise the
#' `ibm_n()` and `ibm_values()` methods also need to be provided.
#'
#' @seealso [bru_mapper] for constructor methods, and
#' [bru_get_mapper] for hooks to extract mappers from latent model object
#' class objects.
#'
#' @name bru_mapper_methods
#' @inheritParams bru_mapper_generics
#' @seealso [bru_mapper], [bru_mapper_generics]
NULL

#' @param mapper A mapper S3 object, inheriting from `bru_mapper`.
#' @param inla_f logical; when `TRUE` for `ibm_n()` and `ibm_values()`, the
#' result must be compatible with the `INLA::f(...)` and corresponding
#' `INLA::inla.stack(...)` constructions.  For `ibm_{eval,jacobian,linear}`,
#' the `input` interpretation may be different.
#' Implementations do not normally need to do anything different, except
#' for mappers of the type needed for hidden multicomponent models such
#' as "bym2", which can be handled by `bru_mapper_collect`.
#' @param \dots Arguments passed on to other methods
#' @describeIn bru_mapper_generics
#' Implementations must return the size of the latent vector
#' being mapped to.
#' @export
ibm_n <- function(mapper, inla_f = FALSE, ...) {
  UseMethod("ibm_n")
}

#' @describeIn bru_mapper_generics
#' Implementations must return an integer denoting the
#' mapper output length.
#' The default implementation returns `NROW(input)`.
#' Mappers such as `bru_mapper_multi` and `bru_mapper_collect`,
#' that can accept `list()` inputs require their own methods implementations.
#' @export
ibm_n_output <- function(mapper, input, state = NULL, inla_f = FALSE, ...) {
  UseMethod("ibm_n_output")
}



#' @describeIn bru_mapper_generics
#' When `inla_f=TRUE`, implementations must return a vector that
#' would be interpretable by an `INLA::f(..., values = ...)` specification.
#' The exception is the method for `bru_mapper_multi`, that returns a
#' multi-column data frame.
#' @export
ibm_values <- function(mapper, inla_f = FALSE, ...) {
  UseMethod("ibm_values")
}

#' @describeIn bru_mapper_generics
#' will become deprecated in 2.7.0. Use `ibm_jacobian`
#' instead.
#' Implementations must return a (sparse) matrix of size `ibm_n_output(...)`
#' by `ibm_n(...)`. The `inla_f=TRUE` argument should only affect
#' the allowed type of input format.
#' @export
ibm_amatrix <- function(mapper, input, state = NULL, inla_f = FALSE, ...) {
  lifecycle::deprecate_soft(
    "2.7.0",
    "ibm_amatrix()",
    "ibm_jacobian()"
  )
  UseMethod("ibm_amatrix")
}

#' @describeIn bru_mapper_generics
#' Implementations must return `TRUE` or `FALSE`.
#' If `TRUE` (returned by the default method unless the mapper
#' contains an `is_linear` variable), users of the mapper
#' may assume the mapper is linear.
#' @export
ibm_is_linear <- function(mapper, ...) {
  UseMethod("ibm_is_linear")
}

#' @param input Data input for the mapper.
#' @describeIn bru_mapper_generics
#' Implementations must return a (sparse) matrix of size
#' `ibm_n_output(mapper, input, inla_f)`
#' by `ibm_n(mapper, inla_f = FALSE)`. The `inla_f=TRUE` argument should
#' only affect the allowed type of input format.
#' @export
ibm_jacobian <- function(mapper, input, state = NULL, inla_f = FALSE, ...) {
  UseMethod("ibm_jacobian")
}


#' @describeIn bru_mapper_generics
#' Implementations must return a [bru_mapper_taylor] object
#' The linearisation information includes `offset`, `jacobian`, and `state0`.
#' The state information indicates for which state the `offset` was evaluated,
#' with `NULL` meaning all-zero.
#' The linearised mapper output is defined as
#' `effect(input, state) = offset(input, state0) + jacobian(input, state0) %*% (state - state0)`.
#' The default method calls `ibm_eval()` and `ibm_jacobian()` to generate
#' the needed information.
#' @export
ibm_linear <- function(mapper, input, state = NULL, ...) {
  UseMethod("ibm_linear")
}

#' @param state A vector of latent state values for the mapping,
#' of length `ibm_n(mapper, inla_f = FALSE)`
#' @describeIn bru_mapper_generics
#' Implementations must return a vector of length `ibm_n_output(...)`.
#' The `input` contents must
#' be in a format accepted by `ibm_jacobian(...)`
#' for the mapper.
#' @export
ibm_eval <- function(mapper, input, state = NULL, ...) {
  UseMethod("ibm_eval")
}

#' @describeIn bru_mapper_generics
#' Implementations must return a character vector of sub-mapper names, or `NULL`.
#' Intended for providing information about multi-mappers and mapper collections.
#' @export
#' @examples
#' # ibm_names
#' mapper <- bru_mapper_multi(list(
#'   A = bru_mapper_index(2),
#'   B = bru_mapper_index(2)
#' ))
#' ibm_names(mapper)
#' ibm_names(mapper) <- c("new", "names")
#' ibm_names(mapper)
ibm_names <- function(mapper) {
  UseMethod("ibm_names")
}

#' @param value a character vector of the same length as the number
#' of sub-mappers in the mapper
#' @describeIn bru_mapper_generics Set mapper names.
#' @export
`ibm_names<-` <- function(mapper, value) {
  UseMethod("ibm_names<-")
}


#' @describeIn bru_mapper_generics
#' Implementations must return a logical vector of `TRUE/FALSE` for
#' the subset such that, given the full A matrix and values output,
#' `A[, subset, drop = FALSE]` and `values[subset]`
#' (or `values[subset, , drop = FALSE]` for data.frame values) are equal
#' to the `inla_f = TRUE` version of A and values. The default method uses
#' the `ibm_values` output to construct the subset indexing.
#' @export
ibm_inla_subset <- function(mapper, ...) {
  UseMethod("ibm_inla_subset")
}

#' @describeIn bru_mapper_generics
#' Implementations should return a logical vector of length
#' `ibm_n_output(mapper, input, state, ...)` indicating which, if any,
#' output elements of `ibm_eval(mapper, input, state, ...)` are known to be
#' invalid.
#' For for multi/collect mappers, a list, when given a `multi=TRUE` argument.
#' @export
ibm_invalid_output <- function(mapper, input, state, ...) {
  UseMethod("ibm_invalid_output")
}


#' @describeIn inlabru-deprecated Use case changed to [ibm_invalid_output()]
#' @export
ibm_valid_input <- function(...) {
  lifecycle::deprecate_stop(
    "2.7.0",
    "ibm_valid_input()"
  )
  UseMethod("ibm_valid_input")
}


#' Methods for mapper lists
#'
#' `bru_mapper` lists can be combined into `bm_list` lists.
#' @name bm_list
#' @param \dots Objects to be combined.
NULL

#' @export
#' @describeIn bm_list The `...` arguments should be `bru_mapper`
#' objects.
#' @examples
#' m <- c(A = bru_mapper_const(), B = bru_mapper_scale())
#' str(m)
#' str(m[2])
`c.bru_mapper` <- function(...) {
  stopifnot(all(vapply(
    list(...),
    function(x) inherits(x, "bru_mapper"),
    TRUE
  )))
  mappers <- list(...)
  class(mappers) <- c("bm_list", "list")
  mappers
}

#' @export
#' @describeIn bm_list The `...` arguments should be `bm_list`
#' objects.
`c.bm_list` <- function(...) {
  stopifnot(all(vapply(
    list(...),
    function(x) inherits(x, "bm_list"),
    TRUE
  )))
  object <- NextMethod()
  class(object) <- c("bm_list", "list")
  object
}

#' @export
#' @param x `bm_list` object from which to extract element(s)
#' @param i indices specifying elements to extract
#' @describeIn bm_list Extract sub-list
`[.bm_list` <- function(x, i) {
  object <- NextMethod()
  class(object) <- c("bm_list", "list")
  object
}




# Summaries ----

ibm_shortname <- function(mapper, ...) {
  shortname <- sub("^bru_mapper_", "", class(mapper)[1])
  if (nchar(shortname) == 0) {
    shortname <- "default"
  }
  shortname
}


#' @title mapper object summaries
#'
#' @param object `bru_mapper` object to summarise
#' @param \dots Unused arguments
#' @param prefix character prefix for each line. Default `""`.
#' @param initial character prefix for the first line. Default `initial=prefix`.
#' @param depth The recursion depth for multi/collection/pipe mappers.
#' Default 1, to only show the collection, and not the contents of the sub-mappers.
#' @export
#' @method summary bru_mapper
#' @rdname bru_mapper_summary
summary.bru_mapper <- function(object, ...,
                               prefix = "",
                               initial = prefix,
                               depth = 1) {
  summary_object <- paste0(initial, ibm_shortname(object))
  class(summary_object) <- c(
    "summary_bru_mapper",
    class(summary_object)
  )
  summary_object
}

#' @export
#' @method summary bru_mapper_multi
#' @rdname bru_mapper_summary
#' @examples
#' mapper <-
#'   bru_mapper_pipe(
#'     list(
#'       bru_mapper_multi(list(
#'         A = bru_mapper_index(2),
#'         B = bru_mapper_index(3)
#'       )),
#'       bru_mapper_index(2)
#'     )
#'   )
#' summary(mapper, depth = 2)
summary.bru_mapper_multi <- function(object, ...,
                                     prefix = "",
                                     initial = prefix,
                                     depth = 1) {
  summary_object <- NextMethod()
  if (depth <= 0) {
    class(summary_object) <- c("summary_bru_mapper", class(summary_object))
    return(summary_object)
  }
  sub_prefix <- paste0(prefix, "      ")
  summary_object <-
    paste0(
      summary_object,
      "(",
      paste0(
        vapply(
          names(object[["mappers"]]),
          function(nm) {
            paste0(
              nm,
              " = ",
              summary(
                object[["mappers"]][[nm]],
                prefix = sub_prefix,
                initial = "",
                depth = depth - 1
              )
            )
          },
          ""
        ),
        collapse = ", "
      ),
      ")"
    )
  class(summary_object) <- c("summary_bru_mapper", class(summary_object))
  summary_object
}

#' @export
#' @method summary bru_mapper_pipe
#' @rdname bru_mapper_summary
summary.bru_mapper_pipe <- function(object, ...,
                                    prefix = "",
                                    initial = prefix,
                                    depth = 1) {
  summary_object <- NextMethod()
  if (depth <= 0) {
    class(summary_object) <- c("summary_bru_mapper", class(summary_object))
    return(summary_object)
  }
  sub_prefix <- paste0(prefix, "      ")
  summary_object <-
    paste0(
      summary_object,
      " = ",
      paste0(
        vapply(
          object[["mappers"]],
          function(x) {
            summary(x,
              prefix = sub_prefix,
              initial = "",
              depth = depth - 1
            )
          },
          ""
        ),
        collapse = paste0(" -> ")
      )
    )
  class(summary_object) <- c("summary_bru_mapper", class(summary_object))
  summary_object
}

#' @export
#' @method summary bru_mapper_collect
#' @rdname bru_mapper_summary
summary.bru_mapper_collect <- function(object, ...,
                                       prefix = "",
                                       initial = prefix,
                                       depth = 1) {
  summary_object <- NextMethod()
  if (depth <= 0) {
    class(summary_object) <- c("summary_bru_mapper", class(summary_object))
    return(summary_object)
  }
  sub_prefix <- paste0(prefix, "      ")
  summary_object <-
    paste0(
      summary_object,
      "(",
      paste0(
        vapply(
          seq_along(object[["mappers"]]),
          function(k) {
            nm <- names(object[["mappers"]])[k]
            if (is.null(nm)) {
              nm <- as.character(k)
            }
            paste0(
              nm,
              " = ",
              summary(
                object[["mappers"]][[k]],
                prefix = sub_prefix,
                initial = "",
                depth = depth - 1
              ),
              if (object[["hidden"]] && (k > 1)) {
                "(hidden)"
              } else {
                NULL
              }
            )
          },
          ""
        ),
        collapse = ", "
      ),
      ")"
    )
  class(summary_object) <- c("summary_bru_mapper", class(summary_object))
  summary_object
}


#' @param x Object to be printed
#' @export
#' @method print summary_bru_mapper
#' @rdname bru_mapper_summary
print.summary_bru_mapper <- function(x, ...) {
  cat(x)
  invisible(x)
}


# MAPPERS ----
## Constructor ----

#' @param mapper For `bru_mapper_define`, a prototype mapper object, see Details.
#' For `bru_mapper_scale`, a mapper to be scaled.
#' @param new_class If non-`NULL`, this is added at the front of the class definition
#' @param \dots Deprecated, alternative way to supply optional method definitions.
#' @param methods Deprecated.
#'
#' @describeIn bru_mapper Adds the `new_class` and "bru_mapper" class
#' names to the inheritance list for the input `mapper` object, unless the object
#' already inherits from these.
#'
#' To register mapper classes and methods in scripts, use `.S3method()`
#' to register the methods, e.g.
#' `.S3method("ibm_jacobian", "my_mapper_class", ibm_jacobian.my_mapper_class)`.
#'
#' In packages with `Suggests: inlabru`, add method information for delayed
#' registration, e.g.:
#' ```
#' #' @rawNamespace S3method(inlabru::bru_get_mapper, inla_rspde)
#' #' @rawNamespace S3method(inlabru::ibm_n, bru_mapper_inla_rspde)
#' #' @rawNamespace S3method(inlabru::ibm_values, bru_mapper_inla_rspde)
#' #' @rawNamespace S3method(inlabru::ibm_jacobian, bru_mapper_inla_rspde)
#' ```
#' or before each method, use `@exportS3Method`:
#' ```
#' #' @exportS3Method inlabru::bru_get_mapper
#' ```
#' etc., which semi-automates it.
#'
#' @export
bru_mapper_define <- function(mapper,
                              new_class = NULL,
                              ...,
                              methods = NULL) {
  if (!is.null(methods)) {
    txt <-
      c(
        "In packages with Suggests: inlabru, add method information, e.g.:",
        "  #' @rawNamespace S3method(inlabru::bru_get_mapper, inla_rspde)",
        "  #' @rawNamespace S3method(inlabru::ibm_n, bru_mapper_inla_rspde)",
        "  #' @rawNamespace S3method(inlabru::ibm_values, bru_mapper_inla_rspde)",
        "  #' @rawNamespace S3method(inlabru::ibm_jacobian, bru_mapper_inla_rspde)",
        "or before each method, use @exportS3Method:",
        "  #' @exportS3Method inlabru::bru_get_mapper",
        "etc., which semi-automates it."
      )
    lifecycle::deprecate_warn(
      when = "2.7.0",
      what = "bru_mapper_define(methods)",
      details =
        c(
          "In scripts, use '.S3method()' to register each method.",
          paste(txt, sep = "\n"),
          "In packages with Imports: inlabru, `#' @export` should be sufficient."
        )
    )
  }
  if (!inherits(mapper, "bru_mapper")) {
    class(mapper) <- c("bru_mapper", class(mapper))
  }
  if (!is.null(new_class) && !inherits(mapper, new_class)) {
    class(mapper) <- c(new_class, class(mapper))
  }
  mapper
}

## Default methods ----

#' @describeIn bru_mapper Calls `bru_mapper_define`, passing all
#' arguments along. Mapper implementations should call [bru_mapper_define()]
#' instead, and supply at least a `new_class` class name.
#' Use of the `bru_mapper.default` method will be deprecated from version 2.7.0.
#' @export
bru_mapper.default <- function(...) {
  lifecycle::deprecate_warn(
    "2.7.0",
    "bru_mapper.default()",
    "bru_mapper_define()"
  )
  bru_mapper_define(...)
}


#' @describeIn bru_mapper_generics
#' Returns a non-null element 'n' from the
#' mapper object, and gives an error if it doesn't exist. If `inla_f=TRUE`,
#' first checks for a 'n_inla' element.
#' @export
ibm_n.default <- function(mapper, inla_f = FALSE, ...) {
  if (inla_f && !is.null(mapper[["n_inla"]])) {
    mapper[["n_inla"]]
  } else if (!is.null(mapper[["n"]])) {
    mapper[["n"]]
  } else {
    stop(paste0(
      "Default 'ibm_n()' method called but mapper doesn't have an 'n' element;",
      "\nClass: ",
      paste0(class(mapper), collapse = ", ")
    ))
  }
}

#' @export
#' @describeIn bru_mapper_generics Returns `NROW(input)`
ibm_n_output.default <- function(mapper, input, state = NULL, inla_f = FALSE, ...) {
  NROW(input)
}


#' @describeIn bru_mapper_generics
#' Returns a non-null element
#' 'values' from the mapper object, and `seq_len(ibm_n(mapper))` if
#' it doesn't exist.
#' @export
ibm_values.default <- function(mapper, inla_f = FALSE, ...) {
  if (inla_f && !is.null(mapper[["values_inla"]])) {
    mapper[["values_inla"]]
  } else if (!inla_f && !is.null(mapper[["values"]])) {
    mapper[["values"]]
  } else {
    seq_len(ibm_n(mapper, inla_f = inla_f))
  }
}

#' @describeIn bru_mapper_generics
#' Gives an error message.
#' Mapper classes must implement their own `ibm_jacobian` or
#' `ibm_amatrix` methods. New implementations should use
#' a `ibm_jacobian` method. `ibm_amatrix` may become deprecated
#' in a future version.
#' @export
ibm_amatrix.default <- function(mapper, ...) {
  stop(paste0(
    "Missing implementation of 'ibm_jacobian()/ibm_amatrix()' for mapper of class '",
    paste0(class(mapper), collapse = ", "), "'.\n",
    "New implementations should implement a 'ibm_jacobian()' method."
  ))
}

#' @describeIn bru_mapper_generics
#' Returns logical
#' `is_linear` from the mapper object if it exists, and otherwise `TRUE`.
#' @export
ibm_is_linear.default <- function(mapper, ...) {
  if (!is.null(mapper[["is_linear"]])) {
    mapper[["is_linear"]]
  } else {
    TRUE
  }
}

#' @describeIn bru_mapper_generics
#' Calls `ibm_amatrix`, which
#' by default gives an error.
#' Mapper classes should implement their own `ibm_jacobian` method.
#' @export
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

#' @describeIn bru_mapper_generics
#' Calls `ibm_eval()` and `ibm_jacobian()`
#' and returns a `bru_mapper_taylor` object.
#' The `state0` information in the affine mapper indicates for which state
#' the `offset` was evaluated; The affine mapper output is defined as
#' `effect(input, state) = offset(input, state0) + jacobian(input, state0) %*% (state - state0)`
#' @export
ibm_linear.default <- function(mapper, input, state, ...) {
  bru_mapper_taylor(
    offset = ibm_eval(mapper, input, state, ...),
    jacobian = ibm_jacobian(mapper, input, state, ...),
    state0 = state,
    values_mapper = mapper
  )
}

#' @describeIn bru_mapper_generics
#' Verifies that the mapper is linear
#' with `ibm_is_linear()`, and then computes a linear mapping
#' as `ibm_jacobian(...) %*% state`.  When `state` is `NULL`,
#' a zero vector of length `ibm_n_output(...)` is returned.
#' @export
ibm_eval.default <- function(mapper, input, state = NULL, ...) {
  if (!ibm_is_linear(mapper)) {
    stop("Non-linear mappers must implement their own ibm_eval() method.")
  }

  val <- numeric(ibm_n_output(mapper, input, state = state, ...))

  if ((ibm_n(mapper) > 0) && !is.null(state)) {
    A <- ibm_jacobian(mapper, input = input, state = state, ...)
    val <- val + as.vector(A %*% state)
  }

  val
}


#' @export
#' @describeIn bru_mapper_generics Returns `NULL`
ibm_names.default <- function(mapper, ...) {
  NULL
}


#' @describeIn bru_mapper_generics
#' Uses
#' the `ibm_values` output to construct the inla subset indexing, passing
#' extra arguments such as `multi` on to the methods (this means it supports
#' both regular vector values and `multi=1` data.frame values).
#' @export
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


#' @describeIn bru_mapper_generics
#' Returns an all-`FALSE` logical vector.
#' @export
ibm_invalid_output.default <- function(mapper, input, state, ...) {
  rep(FALSE, ibm_n_output(mapper, input = input, state = state, ...))
}





## inla.mesh ####

#' @param mesh An `inla.mesh.1d` or `inla.mesh.2d` object to use as a mapper
#' @param indexed logical; If `TRUE`, the `ibm_values()` output will be the
#' integer indexing sequence for the latent variables. If `FALSE`, the knot
#' locations are returned (useful as an interpolator for `rw2` models
#' and similar).
#' Default: `TRUE`
#' @export
#' @describeIn bru_mapper Creates a mapper for 2D `inla.mesh` objects
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
    input <- fm_transform(input,
      crs = fm_crs(mapper[["mesh"]]),
      passthrough = TRUE
    )
    A <- sf::st_coordinates(input)
    nm <- intersect(colnames(A), c("X", "Y", "Z"))
    input <- as.matrix(A[, nm, drop = FALSE])
  } else if (inherits(input, "Spatial")) {
    input <- fm_transform(input,
      crs = fm_crs(mapper[["mesh"]]),
      passthrough = TRUE
    )
    input <- sp::coordinates(input)
  } else if (!is.matrix(input)) {
    input <- as.matrix(input)
  }
  INLA::inla.spde.make.A(mapper[["mesh"]], loc = input)
}


#' @export
#' @describeIn inlabru-deprecated Replaced by [ibm_jacobian()]
ibm_amatrix.bru_mapper_inla_mesh_2d <- function(...) {
  lifecycle::deprecate_warn("2.6.0", "ibm_amatrix()", "ibm_jacobian()")
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
#' @describeIn bru_mapper Create mapper for an `inla.mesh.1d` object
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
  ok <- !is.na(input)
  if (all(ok)) {
    A <- INLA::inla.spde.make.A(mapper[["mesh"]], loc = input)
  } else {
    A <- Matrix::Matrix(0, length(input), ibm_n(mapper))
    A[ok, ] <- INLA::inla.spde.make.A(mapper[["mesh"]], loc = input[ok])
  }
  A
}

#' @export
#' @describeIn inlabru-deprecated Replaced by [ibm_jacobian()]
ibm_amatrix.bru_mapper_inla_mesh_1d <- function(...) {
  lifecycle::deprecate_warn("2.6.0", "ibm_amatrix()", "ibm_jacobian()")
  ibm_jacobian(...)
}

## _index ####

#' @param n Size of a model for `bru_mapper_index`
#' @export
#' @describeIn bru_mapper Create a an indexing mapper
bru_mapper_index <- function(n = 1L, ...) {
  bru_mapper_define(list(n = n), new_class = "bru_mapper_index")
}

#' @export
#' @rdname bru_mapper_methods
ibm_invalid_output.bru_mapper_index <- function(mapper, input, state, ...) {
  # Allow NA to give effect zero
  nok <- !is.na(input)
  nok[nok] <-
    !((input[nok] >= 1) &
      (input[nok] <= ibm_n_output(mapper, input, state, ...)))
  nok
}


#' @export
#' @rdname bru_mapper_methods
ibm_jacobian.bru_mapper_index <- function(mapper, input, state, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, ibm_n(mapper)))
  }
  ok <- which((1L <= input) & (input <= ibm_n(mapper)))
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
  lifecycle::deprecate_warn("2.6.0", "ibm_amatrix()", "ibm_jacobian()")
  ibm_jacobian(...)
}

## _taylor ####

#' @param offset For `bru_mapper_taylor`, an offset vector evaluated
#' at `state0`.
#' May be `NULL`, interpreted as an all-zero vector of length determined by
#' a non-null Jacobian.
#' @param jacobian For `bru_mapper_taylor`, the Jacobian matrix,
#' evaluated at `state0`, or, a named list of such matrices.
#' May be `NULL` or an empty list, for a constant mapping.
#' @param state0 For `bru_mapper_taylor`, the state the linearisation
#' was evaluated at, or a list of length matching the `jacobian` list.
#' `NULL` is interpreted as 0.
#' @param values_mapper mapper object to be used for `ibm_n` and
#' `ibm_values` for `inla_f=TRUE` (experimental, currently unused)
#' @describeIn bru_mapper
#' Provides a pre-computed affine mapping,
#' internally used to represent and evaluate linearisation information.
#' The `state0` information indicates for which state the `offset` was evaluated;
#' The affine mapper output is defined as
#' `effect(state) = offset + jacobian %*% (state - state0)`
#' @export
bru_mapper_taylor <- function(offset = NULL, jacobian = NULL, state0 = NULL, ...,
                              values_mapper = NULL) {
  stopifnot(!is.null(offset) || !is.null(jacobian))
  if (is.null(state0)) {
    n_state <- 0
  } else if (is.list(state0)) {
    n_state <- vapply(
      state0,
      function(x) {
        ifelse(is.null(x),
          0L,
          length(x)
        )
      },
      0L
    )
  } else {
    n_state <- length(state0)
  }

  if (is.null(jacobian)) {
    n_jacobian <- 0
  } else if (is.list(jacobian)) {
    n_jacobian <- vapply(
      jacobian,
      function(x) {
        ifelse(is.null(x),
          0L,
          ncol(x)
        )
      },
      0L
    )
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
      }
    )
    n_state <- n_multi
  }
  if (!is.null(state0)) {
    stopifnot(all(n_jacobian == n_state))
  }

  if (is.null(offset)) {
    offset <- numeric(nrow(jacobian))
  }
  bru_mapper_define(
    list(
      offset = as.vector(offset),
      jacobian = jacobian,
      state0 = state0,
      n_multi = n_multi,
      n = sum(n_multi),
      n_output = length(offset),
      values_mapper = NULL
    ),
    # TODO: maybe allow values_mapper
    new_class = "bru_mapper_taylor"
  )
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
ibm_n_output.bru_mapper_taylor <- function(mapper, input, ...) {
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

#' @export
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
    val <- mapper[["offset"]]
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
  } else {
    stopifnot(is.null(state) || !is.list(state))
    if (is.null(mapper[["state0"]])) {
      val <- mapper[["offset"]] + mapper[["jacobian"]] %*% state
    } else if (is.null(state)) {
      val <- mapper[["offset"]] - mapper[["jacobian"]] %*% mapper[["state0"]]
    } else {
      val <- mapper[["offset"]] + mapper[["jacobian"]] %*% (state - mapper[["state0"]])
    }
  }
  as.vector(val)
}


## _linear ####

#' @export
#' @describeIn bru_mapper Create a mapper for linear effects
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
  lifecycle::deprecate_warn("2.6.0", "ibm_amatrix()", "ibm_jacobian()")
  ibm_jacobian(...)
}


## _matrix ####

#' @param labels Column labels for matrix mappings
#' @export
#' @describeIn bru_mapper Create a matrix mapper, for a given number of columns
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
  lifecycle::deprecate_warn("2.6.0", "ibm_amatrix()", "ibm_jacobian()")
  ibm_jacobian(...)
}



## _factor ####

#' @param values Input values calculated by [input_eval.bru_input()]
#' @param factor_mapping character; selects the type of factor mapping.
#' * `'contrast'` for leaving out the first factor level.
#' * `'full'` for keeping all levels.
#' @param indexed logical; if `TRUE`, the `ibm_values()` method
#' will return an integer vector instead of the factor levels.
#' This is needed e.g. for `group` and `replicate` mappers, since
#' `INLA::f()` doesn't accept factor values. Default: `FALSE`, which
#' works for the main input mappers. The default mapper constructions
#' will set it the required setting.
#' @export
#' @describeIn bru_mapper Create a factor mapper
bru_mapper_factor <- function(values, factor_mapping, indexed = FALSE, ...) {
  factor_mapping <- match.arg(factor_mapping, c("full", "contrast"))
  if (is.factor(values)) {
    mapper <- list(
      levels = levels(values),
      factor_mapping = factor_mapping,
      indexed = indexed
    )
  } else if (is.character(values)) {
    mapper <- list(
      levels = unique(values),
      factor_mapping = factor_mapping,
      indexed = indexed
    )
  } else {
    mapper <- list(
      levels = as.character(sort(unique(values))),
      factor_mapping = factor_mapping,
      indexed = indexed
    )
  }
  mapper$n <- length(mapper$levels) - identical(factor_mapping, "contrast")
  if (indexed) {
    bru_mapper_define(mapper,
      new_class = c(
        "bru_mapper_factor_index",
        "bru_mapper_factor"
      )
    )
  } else {
    bru_mapper_define(mapper, new_class = "bru_mapper_factor")
  }
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_factor <- function(mapper, ...) {
  length(mapper[["levels"]]) - identical(mapper[["factor_mapping"]], "contrast")
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_factor <- function(mapper, ...) {
  if (is.null(mapper[["indexed"]]) || !mapper[["indexed"]]) {
    if (identical(mapper[["factor_mapping"]], "contrast")) {
      mapper$levels[-1L]
    } else {
      mapper$levels
    }
  } else {
    seq_len(ibm_n(mapper))
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
  lifecycle::deprecate_warn("2.6.0", "ibm_amatrix()", "ibm_jacobian()")
  ibm_jacobian(...)
}



## _const ####

#' @export
#' @describeIn bru_mapper Create a constant mapper
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
#' @describeIn inlabru-deprecated Creates a [bru_mapper_const()] mapper.
bru_mapper_offset <- function(...) {
  lifecycle::deprecate_soft("2.6.0", "bru_mapper_offset()", "bru_mapper_const()")
  bru_mapper_define(bru_mapper_const(), new_class = "bru_mapper_offset")
}

#' @export
#' @describeIn inlabru-deprecated Replaced by [bru_mapper_const] methods
ibm_n.bru_mapper_offset <- function(...) {
  lifecycle::deprecate_soft("2.6.0", "bru_mapper_offset()", "bru_mapper_const()")
  NextMethod()
}
#' @export
#' @describeIn inlabru-deprecated Replaced by [bru_mapper_const] methods
ibm_values.bru_mapper_offset <- function(...) {
  lifecycle::deprecate_soft("2.6.0", "bru_mapper_offset()", "bru_mapper_const()")
  NextMethod()
}

#' @export
#' @describeIn inlabru-deprecated Replaced by [bru_mapper_const] methods
ibm_amatrix.bru_mapper_offset <- function(...) {
  lifecycle::deprecate_soft("2.6.0", "bru_mapper_offset()", "bru_mapper_const()")
  lifecycle::deprecate_warn("2.6.0", "ibm_amatrix()", "ibm_jacobian()")
  NextMethod()
}



## _scale ####

#' @export
#' @describeIn bru_mapper Create a standalone
#' scaling mapper that can be used as part of a `bru_mapper_pipe`.
#' If `mapper` is non-null, the `bru_mapper_scale()` constructor
#' returns
#' `bru_mapper_pipe(list(mapper = mapper, scale = bru_mapper_scale()))`
bru_mapper_scale <- function(mapper = NULL, ...) {
  # TODO:
  # 1. If mapper is NULL, implement a plain scaling mapper. Done!
  # 2. First, if mapper is non-null, keep current behaviour. Still allowed
  #    in the methods, but the constructor doesn't do it anymore, but
  #    instead returns a pipe mapper.
  # 3. Later, if mapper is non-null, make a pipe mapper, that pipes
  #   the mapper into a plain null-mapper scaling mapper. Done!
  # 4. Remove the old non-null mapper support from the methods. Future.
  if (is.null(mapper)) {
    bru_mapper_define(
      list(
        is_linear = ibm_is_linear(mapper)
      ),
      new_class = "bru_mapper_scale"
    )
  } else {
    bru_mapper_pipe(list(mapper = mapper, scale = bru_mapper_scale()))
  }
}

#' @param n_state integer giving the length of the state vector for mappers
#' that have state dependent output size.
#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_scale <- function(mapper, ..., state = NULL, n_state = NULL) {
  if (is.null(mapper[["mapper"]])) {
    # Output size depends on the state size
    if (!is.null(state)) {
      length(state)
    } else if (is.null(n_state)) {
      NA_integer_
    } else {
      n_state
    }
  } else {
    ibm_n(mapper[["mapper"]], ...)
  }
}
#' @export
#' @rdname bru_mapper_methods
ibm_n_output.bru_mapper_scale <- function(mapper, input, state = NULL, ...,
                                          n_state = NULL) {
  if (is.null(mapper[["mapper"]])) {
    if (!is.null(state)) {
      length(state)
    } else if (is.null(n_state)) {
      # To allow scalar input weights, do not assume NROW(input) size
      NA_integer_
    } else {
      n_state
    }
  } else {
    ibm_n_output(mapper[["mapper"]], input, state = state, ..., n_state = n_state)
  }
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_scale <- function(mapper, ...,
                                        state = NULL, n_state = NULL) {
  if (is.null(mapper[["mapper"]])) {
    n_state <- ibm_n(mapper, state = state, n_state = n_state)
    if (is.na(n_state)) {
      # Don't know how big the mapper will be
      NULL
    } else {
      seq_len(n_state)
    }
  } else {
    ibm_values(mapper[["mapper"]], ...)
  }
}

#' @export
#' @param sub_lin Internal, optional pre-computed sub-mapper information
#' @details For `bru_mapper_scale`, `input` values without a `scale` element
#' are interpreted as no scaling.
#' @rdname bru_mapper_methods
ibm_jacobian.bru_mapper_scale <- function(mapper, input, state = NULL, ...,
                                          sub_lin = NULL) {
  if (is.null(mapper[["mapper"]])) {
    stopifnot(!is.null(state))
    if (is.null(input)) {
      # No scaling
      return(Matrix::Diagonal(n = length(state), 1.0))
    } else {
      scale <- as.vector(input)
      ok <- !is.na(scale)
      scale[!ok] <- 0
      return(Matrix::Diagonal(n = length(state), scale))
    }
  }

  # Behaviour for v2.6.0, non-null mapper
  if (is.null(sub_lin)) {
    A <- ibm_jacobian(mapper[["mapper"]],
      input = input[["mapper"]],
      state = state, ...
    )
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
#' @rdname bru_mapper_methods
ibm_linear.bru_mapper_scale <- function(mapper, input, state, ...) {
  if (is.null(mapper[["mapper"]])) {
    stopifnot(!is.null(state))
    return(
      bru_mapper_taylor(
        offset = ibm_eval(mapper, input, state, ...),
        jacobian = ibm_jacobian(mapper, input, state, ...),
        state0 = state,
        values_mapper = mapper
      )
    )
  }

  # Behaviour for v2.6.0, non-null mapper
  sub_lin <- ibm_linear(mapper[["mapper"]],
    input[["mapper"]],
    state = state, ...
  )
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
  if (is.null(mapper[["mapper"]])) {
    stopifnot(!is.null(state))
    if (is.null(input)) {
      return(state)
    }
    scale <- as.vector(input)
    ok <- !is.na(scale)
    scale[!ok] <- 0
    return(scale * state)
  }

  # Behaviour for v2.6.0, non-null mapper
  if (!is.null(sub_lin)) {
    values <- ibm_eval(sub_lin, input = NULL, state = state)
  } else {
    values <- ibm_eval(mapper[["mapper"]],
      input = input[["mapper"]],
      state = state, ...
    )
  }
  if ((NROW(values) > 0) && !is.null(input[["scale"]])) {
    scale <- as.vector(input[["scale"]])
    ok <- !is.na(scale)
    scale[!ok] <- 0
    values <- scale * values
  }
  as.vector(values)
}




## _aggregate ####

#' @export
#' @param rescale
#' logical; For `bru_mapper_aggregate` and `bru_mapper_logsumexp`,
#' specifies if the blockwise sums should be normalised by the blockwise weight
#' sums or not:
#' * `FALSE`: (default) Straight weighted sum, no rescaling.
#' * `TRUE`: Divide by the sum of the weight values within each block.
#'   This is useful for integration averages, when the given weights are plain
#'   integration weights. If the weights are `NULL` or all ones, this is
#'   the same as dividing by the number of entries in each block.
#' @param n_block Predetermined number of output blocks. If `NULL`, overrides
#' the maximum block index in the inputs.
#' @describeIn bru_mapper
#' Constructs a mapper
#' that aggregates elements of the input state, so it can be used e.g.
#' for weighted summation or integration over blocks of values.
bru_mapper_aggregate <- function(rescale = FALSE,
                                 n_block = NULL,
                                 ...) {
  bru_mapper_define(
    list(
      rescale = rescale,
      n_block = n_block
    ),
    new_class = "bru_mapper_aggregate"
  )
}

#' @param n_state integer giving the length of the state vector for mappers
#' that have state dependent output size.
#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_aggregate <- function(mapper, ...,
                                       input = NULL,
                                       state = NULL,
                                       n_state = NULL) {
  # Output size depends on the state size
  if (!is.null(state)) {
    length(state)
  } else if (!is.null(n_state)) {
    n_state
  } else if (!is.null(input) &&
    (!is.null(input[["block"]]) ||
      !is.null(input[["weights"]]) ||
      !is.null(input[["log_weights"]]))) {
    max(
      length(input[["block"]]),
      length(input[["weights"]]),
      length(input[["log_weights"]])
    )
    length(input[["block"]])
  } else {
    NA_integer_
  }
}
#' @export
#' @rdname bru_mapper_methods
ibm_n_output.bru_mapper_aggregate <- function(mapper,
                                              input = NULL, ...) {
  if (!is.null(mapper[["n_block"]])) {
    mapper[["n_block"]]
  } else if (!is.null(input)) {
    max(input[["block"]])
  } else {
    NA_integer_
  }
}
#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_aggregate <- function(mapper, ...,
                                            state = NULL,
                                            n_state = NULL) {
  n_state <- ibm_n(mapper, state = state, n_state = n_state)
  if (is.na(n_state)) {
    # Don't know how big the mapper will be
    NULL
  } else {
    seq_len(n_state)
  }
}

bm_aggregate_input <- function(input,
                               state = NULL,
                               n_state = NULL,
                               allow_log = TRUE,
                               force_log = FALSE) {
  if (is.null(input)) {
    block <- NULL
    weights <- NULL
    log_weights <- NULL
  } else {
    block <- input[["block"]]
    weights <- input[["weights"]]
    if (allow_log) {
      log_weights <- input[["log_weights"]]
    } else {
      log_weights <- NULL
    }
  }
  if (is.null(state) && is.null(n_state)) {
    n_state <- max(length(block), length(weights), length(log_weights))
  } else if (!is.null(state)) {
    n_state <- length(state)
  }
  if (is.null(block)) {
    block <- rep(1L, n_state)
  } else if (length(block) == 1) {
    block <- rep(block, n_state)
  }
  if (allow_log && force_log) {
    if (is.null(log_weights)) {
      if (is.null(weights)) {
        log_weights <- rep(0.0, n_state)
      } else {
        log_weights <- log(weights)
        weights <- NULL
      }
    } else if (!is.null(weights)) {
      warning("Both weights and log_weights supplied. Using log_weights.",
        immediate. = TRUE
      )
      weights <- NULL
    }
    if (length(log_weights) == 1) {
      log_weights <- rep(log_weights, n_state)
    }
  } else if (allow_log && !force_log) {
    if (is.null(log_weights) && is.null(weights)) {
      weights <- rep(1.0, n_state)
    } else if (!is.null(weights)) {
      # log_weights is non-NULL
      if (!is.null(log_weights)) {
        warning("Both weights and log_weights supplied. Using weights.",
          immediate. = TRUE
        )
        log_weights <- NULL
      }
      if (length(weights) == 1) {
        weights <- rep(weights, n_state)
      }
    } else {
      # log_weights is non-NULL, weights is null
      if (length(log_weights) == 1) {
        log_weights <- rep(log_weights, n_state)
      }
    }
  } else {
    if (is.null(weights)) {
      weights <- rep(1.0, n_state)
    } else if (length(weights) == 1) {
      weights <- rep(weights, n_state)
    }
  }
  list(block = block, weights = weights, log_weights = log_weights)
}



bm_calc_weights <- function(block, log_weights, weights,
                            n_out, n_state, rescale) {
  if (rescale) {
    # Compute blockwise normalised weights
    if (!is.null(log_weights)) {
      log_weights <- bm_calc_log_weights(
        log_weights = log_weights,
        block = block,
        n_out = n_out,
        n_state = n_state
      )
      weights <- exp(log_weights)
    } else {
      scale <- as.vector(
        Matrix::sparseMatrix(
          i = block,
          j = rep(1L, n_state),
          x = weights,
          dims = c(n_out, 1)
        )
      )
      weights <- weights / scale[block]
    }
  } else {
    if (!is.null(log_weights)) {
      weights <- exp(log_weights)
    }
  }
  weights
}


#' @export
#' @details
#' * For `bru_mapper_aggregate`, `input` should be a list with elements `block`
#' and `weights`. `block`
#' should be a vector of the same length as the `state`, or `NULL`, with `NULL`
#' equivalent to all-1.
#' If `weights` is `NULL`, it's interpreted as all-1.
#' @rdname bru_mapper_methods
ibm_jacobian.bru_mapper_aggregate <- function(mapper, input, state = NULL, ...) {
  input <- bm_aggregate_input(input,
    state = state,
    allow_log = TRUE,
    force_log = FALSE
  )
  n_state <- ibm_n(mapper, input = input, state = state)
  n_out <- ibm_n_output(mapper, input = input)

  weights <-
    bm_calc_weights(
      block = input[["block"]],
      log_weights = input[["log_weights"]],
      weights = input[["weights"]],
      n_out = n_out,
      n_state = n_state,
      rescale = mapper[["rescale"]]
    )

  Matrix::sparseMatrix(
    i = input[["block"]],
    j = seq_len(n_state),
    x = weights,
    dims = c(n_out, n_state)
  )
}


#' @export
#' @rdname bru_mapper_methods
ibm_eval.bru_mapper_aggregate <- function(mapper, input, state = NULL, ...,
                                          sub_lin = NULL) {
  input <- bm_aggregate_input(input,
    state = state,
    allow_log = TRUE,
    force_log = FALSE
  )
  n_state <- ibm_n(mapper, input = input, state = state)
  n_out <- ibm_n_output(mapper, input = input)

  weights <-
    bm_calc_weights(
      block = input[["block"]],
      log_weights = input[["log_weights"]],
      weights = input[["weights"]],
      n_out = n_out,
      n_state = n_state,
      rescale = mapper[["rescale"]]
    )

  values <-
    Matrix::sparseMatrix(
      i = input[["block"]],
      j = rep(1L, n_state),
      x = state * weights,
      dims = c(n_out, 1)
    )
  as.vector(values)
}


#' @export
#' @rdname bru_mapper_methods
ibm_linear.bru_mapper_aggregate <- function(mapper, input, state, ...) {
  bru_mapper_taylor(
    offset = ibm_eval(mapper, input, state, ...),
    jacobian = ibm_jacobian(mapper, input, state, ...),
    state0 = state,
    values_mapper = mapper
  )
}





## _logsumexp ####

#' @export
#' @describeIn bru_mapper
#' Constructs a mapper
#' that aggregates elements of `exp(state)`, with optional non-negative weighting,
#' and then takes the `log()`, so it can be used e.g.
#' for  \eqn{v_k=\log[\sum_{i\in I_k} w_i \exp(u_i)]}{log(blocksum(w*exp(state)))}
#' and \eqn{v_k=\log[\sum_{i\in I_k} w_i \exp(u_i) / \sum_{i\in I_k} w_i]}{log(blocksum(w*exp(state)/blocksum(w)))}
#' calculations.  Relies on the input handling methods for `bru_mapper_aggregate`,
#' but also allows the weights to be supplied on a logarithmic scale as `log_weights`.
#' To avoid numerical overflow, it uses the common method of internally
#' shifting the state blockwise with
#' `(state-log_weights)[block] - max((state-log_weights)[block])`,
#' and shifting the result back afterwards.
bru_mapper_logsumexp <- function(rescale = FALSE,
                                 n_block = NULL,
                                 ...) {
  # Arguments documented for bru_mapper_aggregate
  # Inherit class bru_mapper_aggregate to reuse common methods
  bru_mapper_define(
    list(
      rescale = rescale,
      n_block = n_block
    ),
    new_class = c("bru_mapper_logsumexp", "bru_mapper_aggregate")
  )
}



bm_calc_log_shift <- function(block, log_weights, n_out, n_state) {
  block_k <- sort(unique(block))
  shift <- numeric(n_out)
  shift[block_k] <-
    vapply(
      block_k,
      function(k) {
        max(log_weights[block == k])
      },
      0.0
    )
  shift
}

bm_calc_log_weights <- function(block, log_weights, weights = NULL,
                                n_out, n_state,
                                rescale) {
  if (is.null(log_weights)) {
    log_weights <- log(weights)
  }
  if (rescale) {
    shift <- bm_calc_log_shift(
      block = block, log_weights = log_weights,
      n_out = n_out, n_state = n_state
    )
    log_rescale <- as.vector(
      Matrix::sparseMatrix(
        i = block,
        j = rep(1L, n_state),
        x = exp(log_weights - shift[block]),
        dims = c(n_out, 1)
      )
    )
    log_rescale <- (log(log_rescale) + shift)[block]
    log_weights <- log_weights - log_rescale
  }
  log_weights
}



#' @export
#' @details
#' * For `bru_mapper_logsumexp`, `input` should be a list with elements `block`
#' and `weights`. `block`
#' should be a vector of the same length as the `state`, or `NULL`, with `NULL`
#' equivalent to all-1.
#' If `weights` is `NULL`, it's interpreted as all-1.
#' @rdname bru_mapper_methods
ibm_jacobian.bru_mapper_logsumexp <- function(mapper, input, state = NULL, ...) {
  input <- bm_aggregate_input(input,
    state = state,
    allow_log = TRUE, force_log = TRUE
  )
  n_state <- ibm_n(mapper, input = input, state = state)
  n_out <- ibm_n_output(mapper, input = input)

  log_weights <- bm_calc_log_weights(
    log_weights = input[["log_weights"]],
    block = input[["block"]],
    n_out = n_out,
    n_state = n_state,
    rescale = mapper[["rescale"]]
  )

  # Compute shift for stable log-sum-exp
  w_state <- state + log_weights
  shift <- bm_calc_log_shift(
    log_weights = w_state,
    block = input[["block"]],
    n_out = n_out,
    n_state = n_state
  )

  sum_values <-
    as.vector(
      Matrix::sparseMatrix(
        i = input[["block"]],
        j = rep(1L, n_state),
        x = exp(w_state - shift[input[["block"]]]),
        dims = c(n_out, 1)
      )
    )
  scale <- sum_values[input[["block"]]]

  Matrix::sparseMatrix(
    i = input[["block"]],
    j = seq_len(n_state),
    x = exp(w_state - shift[input[["block"]]]) / scale,
    dims = c(n_out, n_state)
  )
}



#' @export
#' @rdname bru_mapper_methods
ibm_eval.bru_mapper_logsumexp <- function(mapper, input, state = NULL, ...,
                                          sub_lin = NULL) {
  input <- bm_aggregate_input(input,
    state = state,
    allow_log = TRUE, force_log = TRUE
  )
  n_state <- ibm_n(mapper, input = input, state = state)
  n_out <- ibm_n_output(mapper, input = input)

  log_weights <- bm_calc_log_weights(
    log_weights = input[["log_weights"]],
    block = input[["block"]],
    n_out = n_out,
    n_state = n_state,
    rescale = mapper[["rescale"]]
  )

  # Compute shift for stable log-sum-exp
  w_state <- state + log_weights
  shift <- bm_calc_log_shift(
    log_weights = w_state,
    block = input[["block"]],
    n_out = n_out,
    n_state = n_state
  )

  values <-
    Matrix::sparseMatrix(
      i = input[["block"]],
      j = rep(1L, n_state),
      x = exp(w_state - shift[input[["block"]]]),
      dims = c(n_out, 1)
    )
  log(as.vector(values)) + shift
}


#' @export
#' @rdname bru_mapper_methods
ibm_linear.bru_mapper_logsumexp <- function(mapper, input, state, ...) {
  bru_mapper_taylor(
    offset = ibm_eval(mapper, input, state, ...),
    jacobian = ibm_jacobian(mapper, input, state, ...),
    state0 = state,
    values_mapper = mapper
  )
}





## _pipe ####

#' @export
#' @describeIn bru_mapper
#' Create a pipe mapper, where, `mappers` is a list of mappers,
#' where the evaluated output of each mapper is handed as the state to the next
#' mapper..
#' The `input` format for the `ibm_eval` and `ibm_jacobian` methods is
#' a list of inputs, one for each mapper.
bru_mapper_pipe <- function(mappers, ...) {
  is_linear_multi <- vapply(mappers, function(x) ibm_is_linear(x), TRUE)
  n_multi <- vapply(mappers, function(x) as.integer(ibm_n(x)), 0L)
  n <- ibm_n(mappers[[1]])
  is_linear <- all(is_linear_multi)
  if (is.null(names(mappers))) {
    names(mappers) <- as.character(seq_along(mappers))
  } else if ("" %in% names(mappers)) {
    warning("Either all or none of the pipe sub-mappers should be named.", immediate. = TRUE)
    names(mappers) <- as.character(seq_along(mappers))
  }
  bru_mapper_define(
    list(
      mappers = mappers,
      is_linear_multi,
      is_linear = is_linear,
      n_multi = n_multi,
      n = n
    ),
    new_class = "bru_mapper_pipe"
  )
}

#' @export
#' @rdname bru_mapper_methods
ibm_n.bru_mapper_pipe <- function(mapper, ..., state = NULL) {
  if (is.null(mapper[["mappers"]][[1]]) && is.null(state)) {
    return(NA_integer_)
  }
  ibm_n(mapper[["mappers"]][[1]], ..., state = state)
}

#' @export
#' @rdname bru_mapper_methods
ibm_n_output.bru_mapper_pipe <- function(mapper, input, state = NULL, ...) {
  final <- length(mapper[["mappers"]])
  ibm_n_output(mapper[["mappers"]][[final]], input[[final]], state = state, ..., )
}

#' @export
#' @rdname bru_mapper_methods
ibm_values.bru_mapper_pipe <- function(mapper, ...) {
  ibm_values(mapper[["mappers"]][[1]], ...)
}

#' @export
#' @rdname bru_mapper_methods
ibm_jacobian.bru_mapper_pipe <- function(mapper, input, state = NULL, ...) {
  state_k <- state
  if (!is.null(input)) {
    if (is.null(names(input))) {
      names(input) <- names(mapper[["mappers"]])[seq_along(input)]
    }
  }
  first <- names(mapper[["mappers"]])[1]
  for (k in names(mapper[["mappers"]])) {
    # TODO: Introduce an "eval and Jacobian" method to avoid
    # double-computing the evaluations when needing both.
    A_k <- ibm_jacobian(mapper[["mappers"]][[k]],
      input = input[[k]],
      state = state_k, ...
    )
    state_k <- ibm_eval(mapper[["mappers"]][[k]],
      input = input[[k]],
      state = state_k,
      ...
    )
    if (k == first) {
      A <- A_k
    } else {
      A <- A_k %*% A
    }
  }
  A
}



#' @export
#' @rdname bru_mapper_methods
ibm_linear.bru_mapper_pipe <- function(mapper, input, state, ...) {
  bru_mapper_taylor(
    offset = ibm_eval(mapper, input, state, ...),
    jacobian = ibm_jacobian(mapper, input, state, ...),
    state0 = state,
    values_mapper = mapper
  )
}


#' @export
#' @rdname bru_mapper_methods
ibm_eval.bru_mapper_pipe <- function(mapper, input, state = NULL, ...) {
  state_k <- state
  for (k in seq_along(mapper[["mappers"]])) {
    # TODO: Introduce an "eval and Jacobian" method...
    state_k <- ibm_eval(mapper[["mappers"]][[k]],
      input = input[[k]],
      state = state_k,
      ...
    )
  }
  state_k
}



## _multi ####

#' @param values Input values calculated by [input_eval.bru_input()]
#' @param mappers A list of `bru_mapper` objects
#' @describeIn bru_mapper Constructs a rowwise Kronecker product mapping
#' @export
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
#' The names must match the sub-mappers, see [ibm_names.bru_mapper_multi()].
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
  lifecycle::deprecate_warn("2.6.0", "ibm_amatrix()", "ibm_jacobian()")
  ibm_jacobian(...)
}



#' @param pre_A Internal; precomputed Jacobian matrix
#' @export
#' @rdname bru_mapper_methods
ibm_linear.bru_mapper_multi <- function(mapper, input, state,
                                        inla_f = FALSE,
                                        ...) {
  input <- bru_mapper_multi_prepare_input(mapper, input)
  A <- ibm_jacobian(mapper, input,
    state = state,
    inla_f = inla_f, multi = FALSE,
    ...
  )
  bru_mapper_taylor(
    offset = ibm_eval(mapper, input, state,
      inla_f = inla_f, multi = FALSE, ...,
      pre_A = A
    ),
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
      ...
    )
  } else {
    if (is.null(pre_A)) {
      pre_A <- ibm_jacobian(mapper,
        input = input, state = state,
        inla_f = inla_f, multi = FALSE, ...
      )
    }
    val <- pre_A %*% state
  }
  as.vector(val)
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
#' * `ibm_invalid_output` for `bru_mapper_multi` accepts a list with
#' named entries, or a list with unnamed but ordered elements.
#' The names must match the sub-mappers, see [ibm_names.bru_mapper_multi()].
#' Each list element should take a format accepted by the corresponding
#' sub-mapper. In case each element is a vector, the input can be given as a
#' data.frame with named columns, a matrix with named columns, or a matrix
#' with unnamed but ordered columns.
#' @export
#' @rdname bru_mapper_methods
ibm_invalid_output.bru_mapper_multi <- function(mapper,
                                                input,
                                                state,
                                                inla_f = FALSE,
                                                multi = FALSE,
                                                ...) {
  indexing <- bm_multi_indexing(mapper, input)
  if (is.matrix(input)) {
    invalid <-
      lapply(
        indexing,
        function(x) {
          ibm_invalid_output(
            mapper[["mappers"]][[x]],
            input = input[, x],
            inla_f = inla_f,
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
            inla_f = inla_f,
            multi = FALSE
          )
        }
      )
  } else {
    invalid <- as.list(rep(TRUE, length(mapper[["mappers"]])))
    names(invalid) <- names(mapper[["mappers"]])
  }
  if (!multi) {
    # Combine the vectors (v1, v2, v3) -> v1 | v2 | v3
    invalid_ <- invalid[[1]]
    for (k in seq_len(length(mapper[["mappers"]]) - 1)) {
      invalid_ <- invalid_ | invalid[[k + 1]]
    }
    return(invalid_)
  }
  invalid
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

#' @describeIn bru_mapper_methods
#' Returns the names from the sub-mappers list
#' @export
`ibm_names.bru_mapper_multi` <- function(mapper) {
  names(mapper[["mappers"]])
}

#' @param value a character vector of up to the same length as the number
#' of mappers in the multi-mapper x
#' @export
#' @rdname bru_mapper_methods
`ibm_names<-.bru_mapper_multi` <- function(mapper, value) {
  names(mapper[["mappers"]]) <- value
  names(mapper[["n_multi"]]) <- value
  names(mapper[["n_inla_multi"]]) <- value
  names(mapper[["values_multi"]]) <- value
  names(mapper[["values_inla_multi"]]) <- value
  mapper
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
#' @describeIn bru_mapper
#' Constructs a concatenated collection mapping
#' @export
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
#' The names must match the sub-mappers, see [ibm_names.bru_mapper_collect()].
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
  lifecycle::deprecate_warn("2.6.0", "ibm_amatrix()", "ibm_jacobian()")
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
    val <- do.call(c, val)
  }
  val
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
      ...
    )
  A <- ibm_jacobian(mapper, input, state,
    inla_f = FALSE, multi = FALSE, ...,
    sub_lin = sub_lin
  )
  bru_mapper_taylor(
    offset = ibm_eval(mapper, input, state,
      inla_f = FALSE, multi = FALSE, ...,
      sub_lin = sub_lin
    ),
    jacobian = A,
    state0 = state,
    values_mapper = mapper
  )
}




#' @details
#' * `ibm_invalid_output` for `bru_mapper_collect` accepts a list with
#' named entries, or a list with unnamed but ordered elements.
#' The names must match the sub-mappers, see [ibm_names.bru_mapper_collect()].
#' Each list element should take a format accepted by the corresponding
#' sub-mapper. In case each element is a vector, the input can be given as a
#' data.frame with named columns, a matrix with named columns, or a matrix
#' with unnamed but ordered columns.
#' @export
#' @rdname bru_mapper_methods
ibm_invalid_output.bru_mapper_collect <- function(mapper, input, state,
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
    invalid <- as.list(rep(TRUE, length(mapper[["mappers"]])))
  }

  if (!multi) {
    # Combine the vectors (v1, v2, v3) -> c(v1, v2, v3)
    invalid_ <- do.call(c, invalid)
    return(invalid_)
  }

  invalid
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
`ibm_names.bru_mapper_collect` <- function(mapper) {
  names(mapper[["mappers"]])
}

#' @export
#' @rdname bru_mapper_methods
`ibm_names<-.bru_mapper_collect` <- function(mapper, value) {
  names(mapper[["mappers"]]) <- value
  names(mapper[["n_multi"]]) <- value
  names(mapper[["values_multi"]]) <- value
  mapper
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
#' @describeIn bru_mapper Constructs a mapper for `cos`/`sin` functions
#'   of orders 1 (if `intercept` is `TRUE`, otherwise 0) through `order`. The total
#'   number of basis functions is `intercept + 2 * order`.
#'
#'   Optionally, each order can be given a non-unit scaling, via the `scaling`
#'   vector, of length `intercept + order`. This can be used to
#'   give an effective spectral prior. For example, let
#'   ```
#'   scaling = 1 / (1 + (0:4)^2)
#'   x <- seq(0, 1, length.out = 11)
#'   bmh1 = bru_mapper_harmonics(order = 4, interval = c(0, 1))
#'   u1 <- ibm_eval(
#'     bmh1,
#'     input = x,
#'     state = rnorm(9, sd = rep(scaling, c(1, 2, 2, 2, 2)))
#'   )
#'   ```
#'   Then, with
#'   ```
#'   bmh2 = bru_mapper_harmonics(order = 4, scaling = scaling)
#'   u2 = ibm_eval(bmh2, input = x, state = rnorm(9))
#'   ```
#'   the stochastic properties of `u1` and `u2` will be the same, with `scaling^2`
#'   determining the variance for each frequency contribution.
#'
#'   The period for the first order harmonics is shifted and scaled to match
#'   `interval`.
#' @export
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
