#' @include deprecated.R
#' @include mappers.R

## _expr ####

#' @title Mapper for general expressions
#' @description
#' Constructs an expression mapping
#' @export
#' @param expr An expression object, typically created using `rlang::quo()`,
#'   that can be evaluated in a data mask containing the root variables, derived
#'   variables, and data variables.
#' @param labels list of `root`, `derived`, and `suffix`.
#'   The `root` and `derived` elements specify the pronoun
#'   labels for the data mask for root and derived variables. Default "root" and
#'   "derived", respectively. For example, if `root = "latent"`, then the
#'   root variables will be accessible in the expression as `.latent$varname`.
#'   If `suffix` is non-NULL, a character string specifying a suffix to add
#'   to the root variable names when making them directly available to the
#'   expression. Default is NULL, equivalent to `""`. For example, if
#'   `suffix = "_latent"`,
#'   then a root variable named `"x"` will be available as `"x_latent"` in the
#'   expression, in addition to the data mask pronoun version.
#' @param .envir The environment for the expression evaluation. By default, this
#'   is set to the caller environment.
#' @returns A `bm_expr` mapper object
#' @details
#' The `input` is currently ignored.
#'
#' `data` should be a list with data objects, with the main object called
#' `data`.
#' If `is_rowwise == TRUE`, the number of rows in the `data` data.frame
#' determines the number of rows in the output, and the columns can be used as
#' constants in the expression, accessed via `.data$colname`,
#' `.data.$colname`, or `.data.[["colname"]]`.
#'
#' @seealso [bru_mapper], [bru_mapper_generics]
#' @family mappers
#' @examples
#' # Basic expression with only root variables ("x").
#' (m <- bm_expr(rlang::quo(cos(x))))
#' ibm_eval(m, list(), list(x = 1:5))
#' ibm_eval2(m, list(), list(x = 1:5))
#'
#' # Expression with data
#' (m <- bm_expr(rlang::quo(cos(x) * .data$z)))
#' the_data <- list(data = data.frame(z = 11:15))
#' ibm_eval(m, list(), list(x = 1:5), data = the_data)
#' ibm_eval2(m, list(), list(x = 1:5), data = the_data)
#'
#' # Expression with data, root variables, and derived variables.
#' (m <- bm_expr(
#'   rlang::quo(sin(x_latent) + cos(.effect$y) * .data$z),
#'   labels = list(root = "latent", derived = "effect", suffix = "_latent")
#' )
#' )
#' ibm_eval(
#'   m,
#'   list(),
#'   derived = list(y = 2:6), # y = x + 1
#'   jacobians = list(y = list(x = Matrix::Diagonal(1.0, 5))),
#'   data = the_data,
#'   state = list(x = 1:5)
#' )
#'
#' @rdname bm_expr
#'
bm_expr <- function(
  expr,
  #' @param assume character vector listing valid assumptions for the
  #'   combination of the expression and data, derived variables, and root
  #'   variables. This can be used to specify assumptions that allow for more
  #'   efficient Jacobian calculations. For example, if the expression
  #'   is row-wise in the data and derived variables, then
  #'   `assume = "rowwise"` can be used to indicate that the Jacobian with
  #'   respect to derived variables can be calculated simultaneously for all
  #'   rows. Supported assumptions are "rowwise", "linear", "additive", and
  #'   "no_root".
  assume = character(0),
  labels = NULL,
  .envir = rlang::caller_env()
) {
  if (is.null(labels)) {
    labels <- list()
  } else if (is.null(names(labels))) {
    names(labels) <- c("root", "derived", "suffix")[seq_along(labels)]
  } else if (!all(names(labels) %in% c("root", "derived", "suffix"))) {
    stop(
      glue::glue(
        "`bm_expr` `labels` must be a list with elements ",
        "'root', 'derived', and 'suffix'."
      )
    )
  }
  labels <- modifyList(
    list(root = "root", derived = "derived", suffix = NULL),
    labels,
    keep.null = TRUE
  )

  mapper <- list(
    expr = rlang::as_quosure(expr, env = .envir),
    assume = assume,
    labels = labels
  )
  bru_mapper_define(mapper, new_class = "bm_expr")
}

#' @export
#' @rdname ibm_n
#' @param data should be a list with data objects, with the main object called
#' `data`; see [bm_expr()] for details.
#'
ibm_n.bm_expr <- function(
  mapper,
  ...,
  input = NULL,
  state = NULL,
  multi = FALSE,
  data = NULL
) {
  if (is.null(state)) {
    return(NA_integer_)
  }

  n <- lengths(state)
  if (!multi) {
    n <- sum(n)
  }

  n
}


#' @export
#' @rdname ibm_n_output
#' @inheritParams ibm_jacobian
ibm_n_output.bm_expr <- function(
  mapper,
  input,
  state = NULL,
  inla_f = FALSE,
  ...,
  derived = NULL,
  n_state = NULL,
  data = NULL
) {
  if (!("rowwise" %in% mapper[["assume"]])) {
    return(NA_integer_)
  }
  if (!is.null(data[["data"]])) {
    return(NROW(data[["data"]][[1]]))
  }
  if (!is.null(derived[[1]])) {
    return(NROW(derived[[1]]))
  }
  if (!is.null(n_state)) {
    return(n_state)
  }
  NA_integer_
}

#' @export
#' @rdname ibm_values
#'
ibm_values.bm_expr <- function(mapper, inla_f = FALSE, ...) {
  NULL
}

#' @export
#' @rdname ibm_is_linear
#'
ibm_is_linear.bm_expr <- function(mapper, ...) {
  "linear" %in% mapper[["assume"]]
}


# ibm_jacobian_bm_expr <- function(
#   mapper,
#   input,
#   state,
#   derived = NULL,
#   jacobians = NULL,
#   data = NULL,
#   ...,
#   var,
#   offset,
#   n_output = NROW(offset),
#   eps = 1e-6,
#   env = rlang::caller_env()
# ) {
#   if (length(state[[var]]) == 0L) {
#     return(Matrix::sparseMatrix(
#       i = c(),
#       j = c(),
#       x = c(1),
#       dims = c(NROW(offset), 0)
#     ))
#   }
#
#   # TODO: Store and access adjusted bru_used info for the expression!
#   # Or is it enough to precompute "assume_rowwise"?
#   used <- NULL # TODO!
#   allow_root <- var %in% used[[mapper[["labels"]][["root"]]]]
#
#   if (is.null(comp_simple)) { # TODO: Should this be here or elsewhere?
#     A <- NULL
#     assume_rowwise <- FALSE
#   } else {
#     # Jacobian of each derived variable with respect to the root variable.
#     # Each missing entry implies an all-zero matrix.
#     A <- lapply(jacobians, function(x) x[[var]])
#
#     assume_rowwise <- !allow_root &&
#       is_rowwise &&
#       is.data.frame(data[["data"]])
#     if (assume_rowwise) {
#       if (!is.null(n_output) && (NROW(offset) != n_output)) {
#         stop(
#           "Number of rows (",
#           NROW(offset),
#           ") in the predictor for component '",
#           var,
#           "' does not match the length implied by the response data (",
#           n_output,
#           ")."
#         )
#       }
#       if (NROW(A) == 1L) {
#         A <- Matrix::kronecker(rep(1, NROW(offset)), A)
#       }
#     }
#   }
#
#   triplets <- list(
#     i = integer(0),
#     j = integer(0),
#     x = numeric(0)
#   )
#
#   if (!all(is.finite(offset))) {
#     warning(
#       "Non-finite (-Inf/Inf/NaN) entries detected in predictor.\n",
#       immediate. = TRUE
#     )
#   }
#
#   symmetric_diffs <- FALSE
#   for (k in seq_len(NROW(state[[var]]))) {
#     if (is.null(A)) {
#       row_subset <- seq_len(NROW(offset))
#     } else {
#       Ak <- lapply(A, function(x) x[, k, drop = TRUE])
#       row_subset <- which(Ak != 0.0)
#     }
#     if (length(row_subset) > 0) {
#       if (symmetric_diffs) {
#         state_eps <- list(state, state)
#         state_eps[[1]][[label]][k] <- state[[label]][k] - eps
#         state_eps[[2]][[label]][k] <- state[[label]][k] + eps
#       } else {
#         state_eps <- state
#         state_eps[[label]][k] <- state[[label]][k] + eps
#       }
#       # TODO:
#       # Option: filter out the data and effect rows for which
#       # the rows of A have some non-zeros, or all if !is_rowwise
#       # Option: compute predictor for multiple different states. This requires
#       # constructing multiple states and corresponding effects before calling
#       # evaluate_predictor
#
#       if (symmetric_diffs) {
#         effects_eps <- list(effects, effects)
#       } else {
#         effects_eps <- effects
#       }
#       if (!is.null(A)) {
#         if (assume_rowwise) {
#           if (symmetric_diffs) {
#             for (label_loop in names(effects)) {
#               if (NROW(effects[[label_loop]]) == 1) {
#                 effects_eps[[1]][[label_loop]] <-
#                   rep(effects[[label_loop]], length(row_subset))
#                 effects_eps[[2]][[label_loop]] <-
#                   rep(effects[[label_loop]], length(row_subset))
#               } else {
#                 effects_eps[[1]][[label_loop]] <-
#                   effects[[label_loop]][row_subset]
#                 effects_eps[[2]][[label_loop]] <-
#                   effects[[label_loop]][row_subset]
#               }
#             }
#             effects_eps[[1]][[label]] <-
#               effects_eps[[1]][[label]] - Ak[row_subset] * eps
#             effects_eps[[2]][[label]] <-
#               effects_eps[[2]][[label]] + Ak[row_subset] * eps
#           } else {
#             for (label_loop in names(effects)) {
#               if (NROW(effects[[label_loop]]) == 1) {
#                 effects_eps[[label_loop]] <-
#                   rep(effects[[label_loop]], length(row_subset))
#               } else {
#                 effects_eps[[label_loop]] <- effects[[label_loop]][row_subset]
#               }
#             }
#             effects_eps[[label]] <- effects_eps[[label]] +
#               Ak[row_subset] * eps
#           }
#         } else {
#           if (symmetric_diffs) {
#             effects_eps <- list(effects, effects)
#             effects_eps[[1]][[label]] <- effects_eps[[1]][[label]] - Ak * eps
#             effects_eps[[2]][[label]] <- effects_eps[[2]][[label]] + Ak * eps
#           } else {
#             effects_eps <- effects
#             effects_eps[[label]] <- effects_eps[[label]] + Ak * eps
#           }
#         }
#       }
#       pred_eps <- evaluate_predictor(
#         model,
#         state = if (symmetric_diffs) {
#           state_eps
#         } else {
#           list(state_eps)
#         },
#         data =
#           if (assume_rowwise) {
#             data[row_subset, , drop = FALSE]
#           } else {
#             data
#           },
#         data_extra = data_extra,
#         effects =
#           if (symmetric_diffs) {
#             effects_eps
#           } else {
#             list(effects_eps)
#           },
#         predictor = lhood_expr,
#         used = used,
#         format = "matrix",
#         n_pred =
#           if (assume_rowwise) {
#             length(row_subset)
#           } else {
#             n_pred
#           }
#       )
#       # Store sparse triplet information
#       if (symmetric_diffs) {
#         if (assume_rowwise) {
#           values <- (pred_eps[, 2] - pred_eps[, 1]) / 2
#         } else {
#           values <- (pred_eps[, 2] - pred_eps[, 1]) / 2
#         }
#       } else {
#         if (!all(is.finite(pred_eps))) {
#           warning(
#             "Non-finite (-Inf/Inf/NaN) entries detected in predictor '",
#             label,
#             "' plus eps.\n",
#             immediate. = TRUE
#           )
#         }
#         if (assume_rowwise) {
#           values <- (pred_eps - pred0[row_subset])
#         } else {
#           values <- (pred_eps - pred0)
#         }
#       }
#       nonzero <- is.finite(values)
#       if (!all(nonzero)) {
#         warning(
#           "Non-finite (-Inf/Inf/NaN) entries detected in predictor ",
#           "derivatives for '",
#           label,
#           "'; treated as 0.0.\n",
#           immediate. = TRUE
#         )
#       }
#       nonzero[nonzero] <- (values[nonzero] != 0.0) # Detect exact (non)zeros
#       if (assume_rowwise) {
#         triplets$i <- c(triplets$i, row_subset[nonzero])
#       } else {
#         triplets$i <- c(triplets$i, which(nonzero))
#       }
#       triplets$j <- c(triplets$j, rep(k, sum(nonzero)))
#       triplets$x <- c(triplets$x, values[nonzero] / eps)
#     }
#   }
#   B <- Matrix::sparseMatrix(
#     i = triplets$i,
#     j = triplets$j,
#     x = triplets$x,
#     dims = c(NROW(pred0), NROW(state[[label]]))
#   )
#   if (NROW(B) != NROW(pred0)) {
#     stop(
#       "Jacobian matrix for component '",
#       label,
#       "' has ",
#       NROW(B),
#       " rows, but expected ",
#       NROW(pred0),
#       " rows based on the predictor length."
#     )
#   }
#   B
#
#   NULL
# }

ibm_jacobian_bm_expr_var <- function(
  mapper,
  input,
  state,
  derived = NULL,
  jacobians = NULL,
  data = NULL,
  var,
  offset,
  n_output = NROW(offset),
  eps = 1e-6,
  .envir = rlang::caller_env(),
  ...
) {
  n_offset <- NROW(offset)
  if (is.null(state[[var]]) || (length(state[[var]]) == 0L)) {
    # Zero-column matrix, since there are no latent state variables.
    return(Matrix::sparseMatrix(
      i = c(),
      j = c(),
      x = c(1),
      dims = c(n_offset, 0)
    ))
  }

  affected_nms <- names(jacobians)[
    vapply(
      jacobians,
      function(x) !is.null(x[[var]]),
      logical(1)
    )
  ]
  affected_nms <- intersect(affected_nms, names(derived))
  A_affected <- lapply(jacobians[affected_nms], function(x) x[[var]])

  N <- length(state[[var]])

  # TODO: Take advantage of "additive" being in `mapper$assume`
  if (
    ("rowwise" %in% mapper$assume) &&
      ("no_root" %in% mapper$assume)
  ) {
    # Sum of dE/dv * dv/du for all affected derived variables v and root
    # variable u, computed element-wise.
    B <- Matrix::Matrix(0.0, n_offset, N)
    derived_eps <- derived
    for (nm in affected_nms) {
      derived_eps[[nm]] <- derived[[nm]] + eps
      offset_eps <- ibm_eval(
        mapper,
        input,
        state,
        derived = derived_eps,
        data = data,
        .envir = .envir,
        ...
      )
      the_diff <- (offset_eps - offset) / eps
      nonzero <- the_diff != 0.0
      if ((n_offset > 1L) && (nrow(A_affected[[nm]]) == 1L)) {
        B <- B + Matrix::Matrix(the_diff, n_offset, 1L) %*% A_affected[[nm]]
      } else if (nrow(A_affected[[nm]]) != n_offset) {
        stop(
          "The number of rows (",
          nrow(A_affected[[nm]]),
          ") in the Jacobian for derived variable '",
          nm,
          "' with respect to\nroot variable '",
          var,
          "' does not match the number of rows (",
          n_offset,
          ") in the expression result;\n",
          "This indicates the expression was incorrectly inferred",
          " to be rowwise."
        )
      } else {
        B <- B +
          Matrix::Diagonal(n = n_offset, x = the_diff) %*% A_affected[[nm]]
      }
      # Restore the original derived variable values for the next iteration.
      derived_eps[[nm]] <- derived[[nm]]
    }

    return(B)
  }

  ii <- list()
  jj <- list()
  xx <- list()
  for (k in seq_len(N)) {
    Ak <- lapply(A_affected, function(AA) AA[, k, drop = TRUE])

    state_eps <- state
    state_eps[[var]][k] <- state_eps[[var]][k] + eps
    derived_eps <- derived
    for (nm in affected_nms) {
      # row_subset <- (Ak[[nm]] != 0.0)
      # row_subset <- which(row_subset)
      # derived_eps[[nm]][row_subset] <-
      #   derived_eps[[nm]][row_subset] + Ak[[nm]][row_subset] * eps
      derived_eps[[nm]] <-
        derived_eps[[nm]] + Ak[[nm]] * eps
    }
    offset_eps <- ibm_eval(
      mapper,
      input,
      state_eps,
      derived = derived_eps,
      data = data,
      .envir = .envir,
      ...
    )

    the_diff <- (offset_eps - offset) / eps

    nonzero <- is.finite(offset_eps)
    if (!all(nonzero)) {
      warning(
        "Non-finite (-Inf/Inf/NaN) entries detected in predictor ",
        "derivatives for '",
        var,
        "'; treated as 0.0.\n",
        immediate. = TRUE
      )
    }
    nonzero[nonzero] <- (the_diff[nonzero] != 0.0) # Detect exact (non)zeros
    ii[[k]] <- which(nonzero)
    jj[[k]] <- rep(k, sum(nonzero))
    xx[[k]] <- the_diff[ii[[k]]]
  }
  ii <- unlist(ii)
  jj <- unlist(jj)
  xx <- unlist(xx)
  B <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n_offset, N))

  B
}


#' @describeIn ibm_jacobian
#' Accepts a `state` list with named entries, one for each variable.
#' The `input` format should match the description given for [bm_expr()].
#' @param offset The offset value, pre-calculated by [ibm_eval.bm_expr()].
#' @param derived The state vectors of variables derived from the root
#'   variables. If `NULL` or missing, defaults to an empty list.
#' @param jacobians The state vectors of variables derived from the root
#' If `jacobians` is `NULL` or missing, defaults to an empty list.
#' If not `NULL`, should be a list with named entries, one for variable derived
#' from the root variables. Each list element should be a named list of Jacobian
#' matrices, with names matching the root variables.
#' Missing entries are treated as all-zero matrices.
#' @param data A list with data objects, with the main object called `data`;
#' see [bm_expr()]
#' @param eps The finite difference step size to use for numerical
#' differentiation. Default is `1e-6`.
#' @inheritParams bm_expr .envir
#' @export
#'
ibm_jacobian.bm_expr <- function(
  mapper,
  input,
  state = NULL,
  inla_f = FALSE,
  ...,
  derived = NULL,
  jacobians = NULL,
  data = NULL,
  multi = FALSE,
  offset = NULL,
  eps = 1e-6,
  .envir = rlang::caller_env()
) {
  if (is.null(offset)) {
    offset <- ibm_eval(
      mapper,
      input,
      state,
      derived = derived,
      data = data,
      .envir = .envir,
      ...
    )
  }
  n_offset <- NROW(offset)
  assume_rowwise <-
    ("no_root" %in% mapper[["assume"]]) &&
    ("rowwise" %in% mapper[["assume"]]) &&
    is.data.frame(data[["data"]])
  if (assume_rowwise) {
    if (!is.null(n_offset) && (NROW(offset) != n_offset)) {
      stop(
        "The number of values (",
        NROW(offset),
        ") in the expression '",
        format(mapper[["expr"]]),
        "' does not match the expected length (",
        n_offset,
        ")."
      )
    }
    for (nm in names(derived)) {
      if (NROW(derived[[nm]]) == 1L) {
        derived[[nm]] <- rep(derived[[nm]], n_offset)
      }
    }
    for (nm1 in names(jacobians)) {
      for (nm2 in names(jacobians[[nm1]])) {
        if (NROW(jacobians[[nm1]][[nm2]]) == 1L) {
          jacobians[[nm1]][[nm2]] <- Matrix::kronecker(
            rep(1, n_offset),
            jacobians[[nm1]][[nm2]]
          )
        }
      }
    }
  }

  A <- lapply(
    setNames(nm = names(state)),
    function(var) {
      ibm_jacobian_bm_expr_var(
        mapper,
        input,
        state,
        derived = derived,
        jacobians = jacobians,
        data = data,
        var = var,
        offset = offset,
        eps = eps,
        .envir = .envir
      )
    }
  )

  if (multi) {
    return(A)
  }

  # Combine the matrices (A1, A2, A3, ...) -> cbind(A1, A2, A3, ...)
  A <- do.call(cbind, A)
  A
}

bm_expr_data_mask <- function(
  mapper,
  input,
  state = NULL,
  derived = NULL,
  data = data,
  .envir = rlang::caller_env()
) {
  suffix <- mapper[["labels"]][["suffix"]]
  if (is.null(suffix) || (is.character(suffix) && identical(suffix, ""))) {
    state_with_suffix <- NULL
  } else {
    state_names_with_suffix <-
      expand_labels(names(state), names(state), mapper[["labels"]][["suffix"]])
    state_with_suffix <- stats::setNames(state, state_names_with_suffix)
  }
  data_mask <- bru_data_mask(
    stats::setNames(
      c(
        list(
          derived = derived,
          state_with_suffix,
          root = state
        ),
        data
      ),
      c(
        mapper[["labels"]][["derived"]],
        "",
        mapper[["labels"]][["root"]],
        names(data)
      )
    )
  )
}

#' @export
#' @describeIn ibm_eval_methods
#' Accepts a `state` list with named entries, one for each variable.
#' The `input` and `data` formats should match the description given for
#' [bm_expr()].
#' @param data_mask A data mask object to use for evaluating the expression. If
#'   `NULL` or missing, a data mask will be constructed with [bru_data_mask()]
#'   from the `data`, `state`, `derived`, and `.envir` arguments. This can be
#'   used to avoid redundant construction of the data mask when evaluating
#'   different expressions multiple times with the same input and state.
#' @param .envir The environment in which to evaluate the expression. By
#'   default, this is set to the caller environment.
#' @param derived The state vectors of variables derived from the root
#'   variables. If `NULL` or missing, defaults to an empty list.
ibm_eval.bm_expr <- function(
  mapper,
  input,
  state = NULL,
  ...,
  derived = NULL,
  data = NULL,
  data_mask = NULL,
  .envir = rlang::caller_env()
) {
  if (is.null(data_mask)) {
    data_mask <- bm_expr_data_mask(
      mapper,
      input,
      state = state,
      derived = derived,
      data = data,
      .envir = .envir
    )
  }
  val <- rlang::eval_tidy(
    mapper[["expr"]],
    data = data_mask,
    env = .envir
  )
  val
}

#' @export
#' @rdname ibm_as_taylor
#' @inheritParams ibm_jacobian
#'
ibm_as_taylor.bm_expr <- function(
  mapper,
  input,
  state,
  ...,
  derived = NULL,
  jacobians = NULL,
  data = data,
  inla_f = FALSE
) {
  eval2 <- ibm_eval2(
    mapper,
    input = input,
    state = state,
    derived = derived,
    jacobians = jacobians,
    inla_f = FALSE,
    multi = TRUE,
    ...
  )
  bm_taylor(
    offset = eval2$offset,
    jacobian = eval2$jacobian,
    state0 = state,
    values_mapper = mapper
  )
}

#' @export
#' @rdname ibm_eval2
#'
ibm_eval2.bm_expr <- function(mapper, input, state = NULL, ..., data = NULL) {
  offset <- ibm_eval(
    mapper,
    input,
    state,
    ...,
    data = data
  )
  jacobian <- ibm_jacobian(
    mapper,
    input,
    state,
    ...,
    offset = offset,
    data = data
  )
  list(offset = offset, jacobian = jacobian)
}


#' @export
#' @method format bm_expr
#' @rdname bm_summary
#' @examples
#' mapper <- bm_expr(~ cos(x))
#' summary(mapper)
#' summary(mapper, depth = 1)
format.bm_expr <- function(x, ..., prefix = "", initial = prefix, depth = 1) {
  txt <- NextMethod()
  if (depth <= 0) {
    return(txt)
  }
  sub_prefix <- paste0(prefix, "      ")
  txt <-
    paste0(
      txt,
      "(",
      format(x[["expr"]]),
      ")"
    )
  txt
}


# bru_obs ####

na_or_value <- function(value, default) {
  if (is.null(value) || anyNA(value)) {
    default
  } else {
    value
  }
}


#' @title Obtain bru component evaluation functions
#' @keywords internal
#' @description Internal generic function and methods for generating evaluation
#'   functions for `bru_comp` objects, that can be used in [generate.bru()] and
#'   [predict.bru()] expressions.
#' @returns A component evaluation function or a list of functions
#' @export
bru_eval_fun <- function(x, ...) {
  UseMethod("bru_eval_fun")
}

#' @describeIn bru_eval_fun Returns function that takes the arguments `main`,
#'   `group`, `replicate`, `weights`, and `.state` and can be evaluated in the
#'   context of a data mask from [bru_data_mask()].
#' @export
bru_eval_fun.bru_comp <- function(x, ...) {
  .comp <- x
  .is_offset <- .comp$main$type %in% c("offset", "const")
  .is_iid <- .comp$main$type %in% c("iid")
  .mapper <- .comp$mapper
  .label <- paste0(.comp$label, "_latent")
  .iid_precision <- paste0("Precision_for_", .comp$label)
  .iid_cache <- list()
  .iid_cache_index <- NULL
  .get_mask <- function(frame = parent.frame(2L)) {
    rlang::eval_tidy(rlang::parse_expr(".mask."), env = frame)
  }
  .eval_tidy <- function(expr, mask, frame = parent.frame(2L)) {
    rlang::eval_tidy(expr, data = mask, env = frame)
  }
  eval_fun <- function(
    main,
    group = NULL,
    replicate = NULL,
    weights = NULL,
    .state = NULL
  ) {
    .mask <- .get_mask()
    n_input <- ibm_n_output(
      .mapper[["mappers"]][["core"]][["mappers"]][["main"]],
      input = main
    )
    if (is.null(group)) {
      group <- rep(1, n_input)
    }
    if (is.null(replicate)) {
      replicate <- rep(1, n_input)
    }
    if (!.is_offset && is.null(.state)) {
      .state <- .eval_tidy(
        rlang::parse_expr(.label),
        mask = .mask
      )
    }
    .input <- list(
      core = list(
        main = main,
        group = group,
        replicate = replicate
      ),
      scale = weights
    )

    .values <- ibm_eval(
      .mapper,
      input = .input,
      state = .state
    )
    if (.is_iid) {
      # .is_iid, invalid indices give new samples
      # Check for known invalid output elements, based on the
      # initial mapper (subsequent mappers in the component pipe
      # are assumed to keep the same length and validity)
      not_ok <- ibm_invalid_output(
        .mapper[["mappers"]][[1]],
        input = .input[[1]],
        state = .state
      )
      if (any(not_ok)) {
        .cache_state_index <- .eval_tidy(
          rlang::parse_expr(".cache_state_index"),
          mask = .mask
        )
        if (!identical(.cache_state_index, .iid_cache_index)) {
          .iid_cache_index <<- .cache_state_index
          .iid_cache <<- list()
        }
        key <- as.character(main[not_ok])
        not_cached <- !(key %in% names(.iid_cache))
        if (any(not_cached)) {
          .prec <- .eval_tidy(
            rlang::parse_expr(.iid_precision),
            mask = .mask
          )
          for (k in unique(key[not_cached])) {
            .iid_cache[k] <<- rnorm(1, mean = 0, sd = .prec^-0.5)
          }
        }
        .values[not_ok] <- vapply(
          key,
          function(k) .iid_cache[[k]],
          0.0
        )
      }
    } else {
      not_ok <- ibm_invalid_output(
        .mapper[["mappers"]][[1]],
        input = .input[[1]],
        state = .state
      )
      if (any(not_ok)) {
        warning(
          "Inputs for `ibm_eval()` for '",
          .comp[["label"]],
          '" give some invalid outputs.'
        )
      }
    }

    as.matrix(.values)
  }
  eval_fun
}

#' @describeIn bru_eval_fun Returns a list of component evaluation functions
#'   named `<label>_eval` for each component with label `<label>`.
#' @export
bru_eval_fun.bru_comp_list <- function(x, ...) {
  comp_lst <- x
  eval_list <- lapply(comp_lst, bru_eval_fun, ...)
  names(eval_list) <- paste0(names(eval_list), "_eval")
  eval_list
}

#' @describeIn bru_eval_fun Returns a list of component evaluation functions
#'   named `<label>_eval` for each model component with label `<label>`.
#' @export
bru_eval_fun.bru_model <- function(x, ...) {
  bru_eval_fun(as_bru_comp_list(x, ...))
}


bru_expr_as_mapper <- function(pred_expr) {
  expr_mapper <- bm_expr(
    bru_pred_expr(pred_expr, format = "quo"),
    labels = list(
      root = "latent",
      derived = "effects",
      suffix = "_latent"
    ),
    assume = c(
      if (isTRUE(bru_is_rowwise(pred_expr))) "rowwise" else NULL,
      if (isTRUE(bru_is_linear(pred_expr))) "linear" else NULL,
      if (isTRUE(bru_is_additive(pred_expr))) "additive" else NULL,
      if (isTRUE(length(bru_used(pred_expr)[["latent"]]) == 0L)) {
        "no_root"
      } else {
        NULL
      }
    )
  )
}

#' @rdname ibm_eval2
#' @export
#' @param comp_mappers A list of mappers, typically from
#'   `as_bm_list<bru_comp_list>`.
#' @param eval_fun A list of functions, typically from
#'   [bru_eval_fun()].
ibm_eval2.bru_obs <- function(
  mapper,
  input,
  state,
  ...,
  multi = FALSE,
  comp_mappers,
  eval_fun = NULL
) {
  used <- bru_used(mapper)
  nms <- unique(c(used$effect, used$latent))
  param_nms <- setdiff(names(state), nms)
  state_nms <- setdiff(names(state), param_nms)
  param <- state[param_nms]

  # TODO: check consistency; nms should be equal to names(state) or a subset of
  # it, and should include nothing from param_nms.

  der_nms <- intersect(nms, used$effect)
  res <- ibm_eval2(
    comp_mappers[der_nms],
    input = input[["comp"]][der_nms],
    state = state[der_nms],
    ...
  )
  derived <- lapply(res, function(x) x$offset)
  jacobians <- lapply(
    setNames(nm = der_nms),
    function(x) {
      setNames(list(res[[x]]$jacobian), nm = x)
    }
  )

  pred_expr <- mapper[["pred_expr"]]
  expr_mapper <- bru_expr_as_mapper(pred_expr)

  res2 <- ibm_eval2(
    expr_mapper,
    input = list(),
    state = state[state_nms],
    derived = derived,
    jacobians = jacobians,
    multi = TRUE,
    data = list(
      param = param,
      data = mapper[["data"]],
      response_data = mapper[["response_data"]],
      data_extra = mapper[["data_extra"]],
      fun = eval_fun
    )
  )

  # Feed forward into optional transformation mapper.
  # TODO: work in progress... Should perhaps be part of the pred_expr object,
  # or in a wrapper pipeline.
  post_mapper <- mapper[["aggregate"]]
  if (is.null(post_mapper)) {
    n_resp <- bru_response_size(mapper)
    if ((length(res2$offset) == 1L) && (n_resp > 1L)) {
      res2$offset <- rep(res2$offset, n_resp)
      extended_scalar <- TRUE
    } else if (length(res2$offset) != n_resp) {
      stop(
        glue(
          "The number of values ({length(res2$offset)})",
          " in the predictor for observation model ",
          "{na_or_value(mapper$tag, '<unknown>')} does not match",
          " the expected length ({n_resp})."
        )
      )
    } else {
      extended_scalar <- FALSE
    }

    B <- res2$jacobian
    if (extended_scalar) {
      B <- lapply(B, function(x) Matrix::Matrix(1.0, n_resp, 1L) %*% x)
    }
    if (!multi) {
      B <- do.call(cbind, B)
    }
    return(list(offset = res2$offset, jacobian = B))
  }

  # The input is constructed as part of the overall bru_obs input evaluation,
  # so that we have a this pre-computed (input for ibm_eval2.bru_obs_list):
  # input = list(`obs1` = list(comp = input for components,
  #                            post = post_mapper input), `obs2` = ...)

  n_post <- ibm_n(
    post_mapper,
    input = input[["post"]],
    state = res2$offset,
    multi = FALSE
  )
  if ((length(res2$offset) == 1L) && is.finite(n_post) && (n_post > 1L)) {
    res2$offset <- rep(res2$offset, n_post)
    extended_scalar <- TRUE
  } else {
    extended_scalar <- FALSE
  }

  res3 <- ibm_eval2(
    post_mapper,
    input = input[["post"]],
    state = res2$offset,
    multi = FALSE
  )

  # Combine jacobians from the expression mapper and the aggregate mapper.
  offset <- res3$offset
  B <- list()
  for (nm in names(state)) {
    B[[nm]] <- res2$jacobian[[nm]]
    if (is.null(res3$jacobian)) {
      if (extended_scalar) {
        B[[nm]] <- Matrix::Matrix(1.0, n_post, 1L) %*% B[[nm]]
      }
    } else {
      if (extended_scalar) {
        B[[nm]] <-
          (res3$jacobian %*% Matrix::Matrix(1.0, n_post, 1L)) %*% B[[nm]]
      } else {
        B[[nm]] <- res3$jacobian %*% B[[nm]]
      }
    }
  }

  n_resp <- bru_response_size(mapper)
  if ((length(res3$offset) == 1L) && (n_resp > 1L)) {
    res3$offset <- rep(res3$offset, n_resp)
    extended_scalar <- TRUE
  } else if (length(res3$offset) != n_resp) {
    stop(
      glue(
        "The number of values ({length(res3$offset)})",
        " in the predictor for {na_or_value(mapper$tag, '<unknown>')}",
        " does not match",
        " the expected length ({n_resp})."
      )
    )
  } else {
    extended_scalar <- FALSE
  }

  if (extended_scalar) {
    B <- lapply(B, function(x) Matrix::Matrix(1.0, n_resp, 1L) %*% x)
  }

  if (!multi) {
    B <- do.call(cbind, B)
  }

  list(offset = offset, jacobian = B)
}


#' @rdname ibm_eval
#' @inheritParams ibm_eval2
#' @export
ibm_eval.bru_obs <- function(
  mapper,
  input,
  state,
  ...,
  multi = FALSE,
  comp_mappers,
  eval_fun = NULL
) {
  used <- bru_used(mapper)
  nms <- unique(c(used$effect, used$latent))
  param_nms <- setdiff(names(state), nms)
  state_nms <- setdiff(names(state), param_nms)
  param <- state[param_nms]

  # TODO: check consistency; nms should be equal to names(state) or a subset of
  # it, and should include nothing from param_nms.

  derived <- list()
  der_nms <- intersect(nms, used$effect)
  derived <- ibm_eval(
    comp_mappers[der_nms],
    input = input[["comp"]][der_nms],
    state = state[der_nms],
    ...
  )

  pred_expr <- mapper[["pred_expr"]]
  expr_mapper <- bru_expr_as_mapper(pred_expr)

  val <- ibm_eval(
    expr_mapper,
    input = list(),
    state = state[state_nms],
    derived = derived,
    multi = TRUE,
    data = list(
      param = param,
      data = mapper[["data"]],
      response_data = mapper[["response_data"]],
      data_extra = mapper[["data_extra"]],
      fun = eval_fun
    )
  )

  # Feed forward into optional transformation mapper.
  # TODO: work in progress... Should be part of the pred_expr object
  post_mapper <- mapper[["aggregate"]]
  # If there is no post_mapper, we can return the expression result directly.
  if (!is.null(post_mapper)) {
    # The input is constructed as part of the overall bru_obs input evaluation,
    # so that we have a this pre-computed (input for ibm_eval2.bru_obs_list):
    # input = list(`obs1` = list(comp = input for components,
    #                            post = post_mapper input), `obs2` = ...)

    n_post <- ibm_n(
      post_mapper,
      input = input[["post"]],
      state = val,
      multi = FALSE
    )
    if ((length(val) == 1L) && is.finite(n_post) && (n_post > 1L)) {
      val <- rep(val, n_post)
    }

    val <- ibm_eval(
      post_mapper,
      input = input[["post"]],
      state = val,
      multi = FALSE
    )
  }

  n_resp <- bru_response_size(mapper)
  if ((length(val) == 1L) && (n_resp > 1L)) {
    val <- rep(val, n_resp)
  } else if (length(val) != n_resp) {
    stop(
      glue(
        "The number of values ({length(val)})",
        " in the predictor for {na_or_value(mapper$tag, '<unknown>')}",
        " does not match",
        " the expected length ({n_resp})."
      )
    )
  }

  val
}

#' @rdname ibm_jacobian
#' @inheritParams ibm_eval2
#' @export
ibm_jacobian.bru_obs <- function(
  mapper,
  input,
  state,
  ...,
  multi = FALSE,
  comp_mappers,
  eval_fun = NULL
) {
  result <- ibm_eval2(
    mapper,
    input = input,
    state = state,
    ...,
    multi = multi,
    comp_mappers = comp_mappers,
    eval_fun = eval_fun
  )

  B <- result$jacobian

  B
}

#' @rdname ibm_as_taylor
#' @export
ibm_as_taylor.bru_obs <- function(
  mapper,
  input,
  state,
  ...,
  multi = FALSE,
  comp_mappers
) {
  stopifnot(isTRUE(multi))
  eval2 <- ibm_eval2(
    mapper,
    input = input,
    state = state,
    multi = TRUE,
    comp_mappers = comp_mappers
  )

  for (nm in names(state)) {
    if (nm %in% names(eval2$jacobian)) {
      eval2$offset <- eval2$offset - eval2$jacobian[[nm]] %*% state[[nm]]
    }
  }

  mappers <- bm_taylor(
    offset = eval2$offset,
    jacobian = eval2$jacobian,
    state0 = NULL
  )

  mappers
}

# bru_obs_list ####

#' @rdname ibm_eval2
#' @export
ibm_eval2.bru_obs_list <- function(
  mapper,
  input,
  state,
  ...,
  multi = FALSE,
  comp_mappers,
  eval_fun = NULL
) {
  results <- lapply(
    setNames(seq_along(mapper), names(mapper)),
    function(m) {
      ibm_eval2.bru_obs(
        mapper[[m]],
        input[[m]],
        state,
        ...,
        multi = TRUE,
        comp_mappers = if (inherits(comp_mappers, "bm_list")) {
          comp_mappers
        } else {
          comp_mappers[[m]]
        },
        eval_fun = eval_fun
      )
    }
  )

  # Combine jacobians from the expression mapper and the aggregate mapper.
  offset <- lapply(results, function(r) r[["offset"]])
  jacobian <- lapply(results, function(r) r[["jacobian"]])

  if (!multi) {
    jac_names <- lapply(jacobian, function(r) names(r))
    jac_names <- unique(unlist(jac_names))
    names(jac_names) <- jac_names
    jac_ncol <- vapply(
      jac_names,
      function(nm) {
        max(vapply(
          jacobian,
          function(r) {
            if (nm %in% names(r)) {
              ncol(r[[nm]])
            } else {
              0L
            }
          },
          0L
        ))
      },
      0L
    )
    jacobian <- lapply(
      jac_names,
      function(nm) {
        do.call(
          rbind,
          lapply(seq_along(jacobian), function(k) {
            if (nm %in% names(jacobian[[k]])) {
              jacobian[[k]][[nm]]
            } else {
              Matrix::sparseMatrix(
                i = c(),
                j = c(),
                x = c(1),
                dims = c(length(offset[[k]]), jac_ncol[nm])
              )
            }
          })
        )
      }
    )

    offset <- unlist(offset)
  }

  list(offset = offset, jacobian = jacobian)
}

#' @rdname ibm_eval
#' @export
ibm_eval.bru_obs_list <- function(
  mapper,
  input,
  state,
  ...,
  multi = FALSE,
  comp_mappers,
  eval_fun = NULL
) {
  offset <- lapply(
    setNames(seq_along(mapper), names(mapper)),
    function(m) {
      ibm_eval.bru_obs(
        mapper[[m]],
        input[[m]],
        state,
        ...,
        multi = TRUE,
        comp_mappers = if (inherits(comp_mappers, "bm_list")) {
          comp_mappers
        } else {
          comp_mappers[[m]]
        },
        eval_fun = eval_fun
      )
    }
  )

  if (!multi) {
    offset <- unlist(offset)
  }

  offset
}

#' @rdname ibm_jacobian
#' @export
ibm_jacobian.bru_obs_list <- function(
  mapper,
  input,
  state,
  ...,
  multi = FALSE,
  comp_mappers,
  eval_fun = NULL
) {
  result <- ibm_eval2(
    mapper,
    input = input,
    state = state,
    ...,
    multi = multi,
    comp_mappers = comp_mappers,
    eval_fun = eval_fun
  )

  B <- result$jacobian

  B
}

#' @rdname ibm_as_taylor
#' @export
ibm_as_taylor.bru_obs_list <- function(
  mapper,
  input,
  state,
  ...,
  multi = FALSE,
  comp_mappers,
  eval_fun = NULL
) {
  stopifnot(isTRUE(multi))

  lapply(
    setNames(seq_along(mapper), names(mapper)),
    function(k) {
      ibm_as_taylor(
        mapper[[k]],
        input[[k]],
        state = state,
        multi = TRUE,
        comp_mappers = if (inherits(comp_mappers, "bm_list")) {
          comp_mappers
        } else {
          comp_mappers[[k]]
        },
        eval_fun = eval_fun
      )
    }
  )
}
