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
#' @param root_label,derived_label Character strings specifying the pronoun
#'   labels for the data mask for root and derived variables. Default "root" and
#'   "derived", respectively. For example, if `root_label = "latent"`, then the
#'   root variables will be accessible in the expression as `.latent$varname`.
#' @param root_suffix If non-NULL, aharacter string specifying a suffix to add
#'   to the root variable names when making them directly available to the
#'   expression. Default is `""`. For example, if `root_suffix = "_latent"`,
#'   then a root variable named `"x"` will be available as `"x_latent"` in the
#'   expression, in addition to the data mask pronoun version.
#' @inheritParams bru_mapper_generics
#' @details
#' The `input` should be a list with elements
#' \describe{
#' \item{data}{Defaults to `NULL`. If not `NULL`, should be a list of
#' data.frame or similar objects, with `data$data` being the main data
#' container.
#' If `is_rowwise == TRUE`, the number of rows in the `data$data` data.frame
#' determines the number of rows in the output, and the columns can be used as
#' constants in the expression, accessed via `.data$colname`,
#' `.data.$colname`, or `.data.[["colname"]]`.}
#' \item{roots}{If `NULL` or missing, defaults to `names(state)`, and must
#' otherwise be a subset of `names(state)`.
#' If not `NULL`, should be a character vector of variable names.
#' By definition, the pre-Jacobian for each root variable is an identity
#' matrix.}
#' \item{derived}{The state vectors of variables derived from the root
#'   variables. If `NULL` or missing, defaults to an empty list.}
#' \item{jacobians}{If `NULL` or missing, defaults to an empty list.
#' If not `NULL`, should be a list with named entries, one for variable derived
#' from the root variables. Each list element should be a named list of Jacobian
#' matrices, with names matching the root variables.
#' Missing entries are treated as all-zero matrices.
#' }
#' }
#'
#' @seealso [bru_mapper], [bru_mapper_generics]
#' @family mappers
#' @examples
#' # Basic expression with only root variables ("x"). Implies
#' # input$roots = names(state) = "x" and identity Jacobian for "x".
#' (m <- bm_expr(rlang::quo(cos(x))))
#' ibm_eval(m, list(), list(x = 1:5))
#' ibm_eval2(m, list(), list(x = 1:5))
#'
#' # Expression with data
#' (m <- bm_expr(rlang::quo(cos(x) * .data$z)))
#' ibm_eval(m, list(data = list(data = data.frame(z = 11:15))), list(x = 1:5))
#' ibm_eval2(m, list(data = data.frame(z = 11:15)), list(x = 1:5))
#'
#' # Expression with data, root variables, and derived variables.
#' (m <- bm_expr(
#'   rlang::quo(sin(x_latent) + cos(.effect$y) * .data$z),
#'   root_label = "latent",
#'   derived_label = "effect",
#'   root_suffix = "_latent"
#' )
#' )
#' ibm_eval(
#'   m,
#'   list(
#'     data = list(data = data.frame(z = 11:15)),
#'     derived = list(y = 2:6), # y = x + 1
#'     jacobians = list(y = list(x = Matrix::Diagonal(1.0, 5)))
#'   ),
#'   state = list(x = 1:5)
#' )
#'
#' @rdname bm_expr
#'
bm_expr <- function(
  expr,
  ...,
  root_label = "root",
  derived_label = "derived",
  root_suffix = NULL,
  .envir = parent.frame()
) {
  mapper <- list(
    expr = rlang::as_quosure(expr, env = .envir),
    is_linear = FALSE,
    is_additive = FALSE,
    is_rowwise = FALSE,
    root_label = root_label,
    derived_label = derived_label,
    root_suffix = root_suffix
  )
  bru_mapper_define(mapper, new_class = "bm_expr")
}

#' @export
#' @rdname ibm_n
#'
ibm_n.bm_expr <- function(
  mapper,
  ...,
  input = NULL,
  state = NULL,
  multi = FALSE
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
#'
ibm_n_output.bm_expr <- function(
  mapper,
  input,
  state = NULL,
  inla_f = FALSE,
  ...,
  n_state = NULL
) {
  if (!mapper[["is_rowwise"]]) {
    return(NA_integer_)
  }
  if (!is.null(input[["data"]][["data"]])) {
    return(NROW(input[["data"]][["data"]][[1]]))
  }
  if (!is.null(input[["derived"]][[1]])) {
    return(length(input[["derived"]][[1]]))
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
  mapper[["is_linear"]]
}


ibm_jacobian_bm_expr <- function(
    mapper,
    input,
    state,
    ...,
    var,
    offset,
    n_output = NROW(offset),
    eps = 1e-6,
    env = rlang::caller_env()
) {
  if (length(state[[var]]) == 0L) {
    return(Matrix::sparseMatrix(
      i = c(),
      j = c(),
      x = c(1),
      dims = c(NROW(offset), 0)
    ))
  }

  # TODO: Store and access adjusted bru_used info for the expression!
  # Or is it enough to precompute "assume_rowwise"?
  used <- NULL # TODO!
  allow_root <- var %in% used[[mapper[["root_label"]]]]

  if (is.null(comp_simple)) { # TODO: Should this be here or elsewhere?
    A <- NULL
    assume_rowwise <- FALSE
  } else {
    # Jacobian of each derived variable with respect to the root variable.
    # Each missing entry implies an all-zero matrix.
    A <- lapply(input[["jacobians"]], function(x) x[[var]])

    assume_rowwise <- !allow_root &&
      is_rowwise &&
      is.data.frame(input[["data"]][["data"]])
    if (assume_rowwise) {
      if (!is.null(n_output) && (NROW(offset) != n_output)) {
        stop(
          "Number of rows (",
          NROW(offset),
          ") in the predictor for component '",
          var,
          "' does not match the length implied by the response data (",
          n_output,
          ")."
        )
      }
      if (NROW(A) == 1L) {
        A <- Matrix::kronecker(rep(1, NROW(offset)), A)
      }
    }
  }

  triplets <- list(
    i = integer(0),
    j = integer(0),
    x = numeric(0)
  )

  if (!all(is.finite(offset))) {
    warning(
      "Non-finite (-Inf/Inf/NaN) entries detected in predictor.\n",
      immediate. = TRUE
    )
  }

  symmetric_diffs <- FALSE
  for (k in seq_len(NROW(state[[var]]))) {
    if (is.null(A)) {
      row_subset <- seq_len(NROW(offset))
    } else {
      Ak <- lapply(A, function(x) x[, k, drop = TRUE])
      row_subset <- which(Ak != 0.0)
    }
    if (length(row_subset) > 0) {
      if (symmetric_diffs) {
        state_eps <- list(state, state)
        state_eps[[1]][[label]][k] <- state[[label]][k] - eps
        state_eps[[2]][[label]][k] <- state[[label]][k] + eps
      } else {
        state_eps <- state
        state_eps[[label]][k] <- state[[label]][k] + eps
      }
      # TODO:
      # Option: filter out the data and effect rows for which
      # the rows of A have some non-zeros, or all if !is_rowwise
      # Option: compute predictor for multiple different states. This requires
      # constructing multiple states and corresponding effects before calling
      # evaluate_predictor

      if (symmetric_diffs) {
        effects_eps <- list(effects, effects)
      } else {
        effects_eps <- effects
      }
      if (!is.null(A)) {
        if (assume_rowwise) {
          if (symmetric_diffs) {
            for (label_loop in names(effects)) {
              if (NROW(effects[[label_loop]]) == 1) {
                effects_eps[[1]][[label_loop]] <-
                  rep(effects[[label_loop]], length(row_subset))
                effects_eps[[2]][[label_loop]] <-
                  rep(effects[[label_loop]], length(row_subset))
              } else {
                effects_eps[[1]][[label_loop]] <-
                  effects[[label_loop]][row_subset]
                effects_eps[[2]][[label_loop]] <-
                  effects[[label_loop]][row_subset]
              }
            }
            effects_eps[[1]][[label]] <-
              effects_eps[[1]][[label]] - Ak[row_subset] * eps
            effects_eps[[2]][[label]] <-
              effects_eps[[2]][[label]] + Ak[row_subset] * eps
          } else {
            for (label_loop in names(effects)) {
              if (NROW(effects[[label_loop]]) == 1) {
                effects_eps[[label_loop]] <-
                  rep(effects[[label_loop]], length(row_subset))
              } else {
                effects_eps[[label_loop]] <- effects[[label_loop]][row_subset]
              }
            }
            effects_eps[[label]] <- effects_eps[[label]] + Ak[row_subset] * eps
          }
        } else {
          if (symmetric_diffs) {
            effects_eps <- list(effects, effects)
            effects_eps[[1]][[label]] <- effects_eps[[1]][[label]] - Ak * eps
            effects_eps[[2]][[label]] <- effects_eps[[2]][[label]] + Ak * eps
          } else {
            effects_eps <- effects
            effects_eps[[label]] <- effects_eps[[label]] + Ak * eps
          }
        }
      }
      pred_eps <- evaluate_predictor(
        model,
        state = if (symmetric_diffs) {
          state_eps
        } else {
          list(state_eps)
        },
        data =
          if (assume_rowwise) {
            data[row_subset, , drop = FALSE]
          } else {
            data
          },
        data_extra = data_extra,
        effects =
          if (symmetric_diffs) {
            effects_eps
          } else {
            list(effects_eps)
          },
        predictor = lhood_expr,
        used = used,
        format = "matrix",
        n_pred =
          if (assume_rowwise) {
            length(row_subset)
          } else {
            n_pred
          }
      )
      # Store sparse triplet information
      if (symmetric_diffs) {
        if (assume_rowwise) {
          values <- (pred_eps[, 2] - pred_eps[, 1]) / 2
        } else {
          values <- (pred_eps[, 2] - pred_eps[, 1]) / 2
        }
      } else {
        if (!all(is.finite(pred_eps))) {
          warning(
            "Non-finite (-Inf/Inf/NaN) entries detected in predictor '",
            label,
            "' plus eps.\n",
            immediate. = TRUE
          )
        }
        if (assume_rowwise) {
          values <- (pred_eps - pred0[row_subset])
        } else {
          values <- (pred_eps - pred0)
        }
      }
      nonzero <- is.finite(values)
      if (!all(nonzero)) {
        warning(
          "Non-finite (-Inf/Inf/NaN) entries detected in predictor ",
          "derivatives for '",
          label,
          "'; treated as 0.0.\n",
          immediate. = TRUE
        )
      }
      nonzero[nonzero] <- (values[nonzero] != 0.0) # Detect exact (non)zeros
      if (assume_rowwise) {
        triplets$i <- c(triplets$i, row_subset[nonzero])
      } else {
        triplets$i <- c(triplets$i, which(nonzero))
      }
      triplets$j <- c(triplets$j, rep(k, sum(nonzero)))
      triplets$x <- c(triplets$x, values[nonzero] / eps)
    }
  }
  B <- Matrix::sparseMatrix(
    i = triplets$i,
    j = triplets$j,
    x = triplets$x,
    dims = c(NROW(pred0), NROW(state[[label]]))
  )
  if (NROW(B) != NROW(pred0)) {
    stop(
      "Jacobian matrix for component '",
      label,
      "' has ",
      NROW(B),
      " rows, but expected ",
      NROW(pred0),
      " rows based on the predictor length."
    )
  }
  B

  NULL
}

#' @describeIn ibm_jacobian
#' Accepts a `state` list with named entries, one for each variable.
#' The `input` format should match the description given for [bm_expr()].
#' @param offset The offset value, pre-calculated by [ibm_eval.bm_expr()].
#' @export
#'
ibm_jacobian.bm_expr <- function(
  mapper,
  input,
  state = NULL,
  inla_f = FALSE,
  multi = FALSE,
  ...,
  offset = NULL,
  env = rlang::caller_env()
) {
  if (is.null(offset)) {
    offset <- ibm_eval(
      mapper,
      input,
      state,
      env = env,
      ...
    )
  }
  A <- lapply(
    setNames(nm = names(state)),
    function(var) {
      ibm_jacobian_bm_expr_var(
        mapper,
        input,
        state,
        var = var,
        offset = offset,
        env = env
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
  env = rlang::caller_env()
) {
  if (
    is.null(mapper[["root_suffix"]]) ||
      (is.character(mapper[["root_suffix"]]) &&
        identical(mapper[["root_suffix"]], ""))
  ) {
    state_with_suffix <- NULL
  } else {
    state_names_with_suffix <-
      expand_labels(names(state), names(state), mapper[["root_suffix"]])
    state_with_suffix <- stats::setNames(state, state_names_with_suffix)
  }
  data_mask <- bru_data_mask(
    stats::setNames(
      c(
        list(
          derived = input[["derived"]],
          state_with_suffix,
          root = state
        ),
        input[["data"]]
      ),
      c(
        mapper[["derived_label"]],
        "",
        mapper[["root_label"]],
        names(input[["data"]])
      )
    )
  )
}

#' @export
#' @describeIn ibm_eval
#' Accepts a `state` list with named entries, one for each variable.
#' The `input` format should match the description given for [bm_expr()].
#'
ibm_eval.bm_expr <- function(
  mapper,
  input,
  state = NULL,
  ...,
  data_mask = NULL,
  env = rlang::caller_env()
) {
  if (is.null(data_mask)) {
    data_mask <- bm_expr_data_mask(mapper, input, state = state, env = env)
  }
  val <- rlang::eval_tidy(
    mapper[["expr"]],
    data = data_mask,
    env = env
  )
  val
}

#' @export
#' @rdname ibm_linear
#'
ibm_linear.bm_expr <- function(mapper, input, state, inla_f = FALSE, ...) {
  eval2 <- ibm_eval2(
    mapper,
    input = input,
    state = state,
    inla_f = FALSE,
    multi = TRUE,
    ...,
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
ibm_eval2.bm_expr <- function(mapper, input, state = NULL, ...) {
  offset <- ibm_eval(
    mapper,
    input,
    state,
    ...
  )
  jacobian <- ibm_jacobian(mapper, input, state, ..., offset = offset)
  list(offset = offset, jacobian = jacobian)
}
