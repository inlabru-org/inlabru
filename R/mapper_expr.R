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
#' @inheritParams bru_mapper_generics
#' @details
#' The `input` should be a list with data objects, with the main object called
#' `data`.
#' If `is_rowwise == TRUE`, the number of rows in the `data` data.frame
#' determines the number of rows in the output, and the columns can be used as
#' constants in the expression, accessed via `.data$colname`,
#' `.data.$colname`, or `.data.[["colname"]]`.
#'
#' @param derived,jacobians The state vectors of variables derived from the root
#'   variables. If `NULL` or missing, defaults to an empty list.
#' If `jacobians` is `NULL` or missing, defaults to an empty list.
#' If not `NULL`, should be a list with named entries, one for variable derived
#' from the root variables. Each list element should be a named list of Jacobian
#' matrices, with names matching the root variables.
#' Missing entries are treated as all-zero matrices.
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
#' ibm_eval(m, list(data = data.frame(z = 11:15)), list(x = 1:5))
#' ibm_eval2(m, list(data = data.frame(z = 11:15)), list(x = 1:5))
#'
#' # Expression with data, root variables, and derived variables.
#' (m <- bm_expr(
#'   rlang::quo(sin(x_latent) + cos(.effect$y) * .data$z),
#'   labels = list(root = "latent", derived = "effect", suffix = "_latent")
#' )
#' )
#' ibm_eval(
#'   m,
#'   list(data = data.frame(z = 11:15)),
#'   derived = list(y = 2:6), # y = x + 1
#'   jacobians = list(y = list(x = Matrix::Diagonal(1.0, 5))),
#'   state = list(x = 1:5)
#' )
#'
#' @rdname bm_expr
#'
bm_expr <- function(
  expr,
  ...,
  #' @param assume character vector listing valid assumptions for the
  #'   combination of the expression and input data, derived variables, and root
  #'   variables. This can be used to specify assumptions that allow for more
  #'   efficient Jacobian calculations. For example, if the expression
  #'   is row-wise in the data and derived variables, then
  #'   `assume = "rowwise"` can be used to indicate that the Jacobian with
  #'   respect to derived variables can be calculated simultaneously for all
  #'   rows. Supported assumptions are "rowwise", "linear", and "additive".
  assume = character(0),
  labels = NULL,
  .envir = parent.frame()
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
  derived = NULL,
  n_state = NULL
) {
  if (!("rowwise" %in% mapper[["assume"]])) {
    return(NA_integer_)
  }
  if (!is.null(input[["data"]])) {
    return(NROW(input[["data"]][[1]]))
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


ibm_jacobian_bm_expr <- function(
  mapper,
  input,
  state,
  derived = NULL,
  jacobians = NULL,
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
  allow_root <- var %in% used[[mapper[["labels"]][["root"]]]]

  if (is.null(comp_simple)) { # TODO: Should this be here or elsewhere?
    A <- NULL
    assume_rowwise <- FALSE
  } else {
    # Jacobian of each derived variable with respect to the root variable.
    # Each missing entry implies an all-zero matrix.
    A <- lapply(input[["jacobians"]], function(x) x[[var]])

    assume_rowwise <- !allow_root &&
      is_rowwise &&
      is.data.frame(input[["data"]])
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


ibm_jacobian_bm_expr_var <- function(
  mapper,
  input,
  state,
  derived = NULL,
  jacobians = NULL,
  var,
  offset,
  n_output = NROW(offset),
  eps = 1e-6,
  assume_rowwise,
  allow_root,
  env = rlang::caller_env(),
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

  eps <- 1e-6
  affected_nms <- names(jacobians)[
    vapply(
      jacobians,
      function(x) !is.null(x[[var]]),
      logical(1)
    )
  ]
  affected_nms <- intersect(affected_nms, names(derived))
  A_affected <- lapply(input[["jacobians"]][affected_nms], function(x) x[[var]])

  N <- length(state[[var]])

  if (assume_rowwise && !allow_root) {
    # Sum of dE/dv * dv/du for all affected derived variables v and root
    # variable u, computed element-wise.
    B <- 0.0
    derived_eps <- derived
    for (nm in affected_nms) {
      derived_eps[[nm]] <- derived[[nm]] + eps
      offset_eps <- ibm_eval(
        mapper,
        input,
        state,
        derived = derived_eps,
        env = env,
        ...
      )
      the_diff <- (offset_eps - offset) / eps
      nonzero <- the_diff != 0.0
      B <- B + Matrix::Diagonal(n = n_offset, x = the_diff) %*% A_affected[[nm]]
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
      row_subset <- (Ak[[nm]] != 0.0)
      row_subset <- which(row_subset)
      derived_eps[[nm]][row_subset] <-
        derived_eps[[nm]][row_subset] + Ak[[nm]][row_subset] * eps
    }
    offset_eps <- ibm_eval(
      mapper,
      input,
      state_eps,
      derived = derived_eps,
      env = env,
      ...
    )

    the_diff <- (offset_eps - offset) / eps
    nonzero <- the_diff != 0.0
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
  multi = FALSE,
  offset = NULL,
  env = rlang::caller_env()
) {
  if (is.null(offset)) {
    offset <- ibm_eval(
      mapper,
      input,
      state,
      derived = derived,
      env = env,
      ...
    )
  }
  n_offset <- NROW(offset)
  # !allow_latent && is_rowwise && is.data.frame(data)
  assume_rowwise <- ("rowwise" %in% mapper[["labels"]]) &&
    is.data.frame(input[["data"]])
  if (assume_rowwise) {
    if (!is.null(n_offset) && (NROW(offset) != n_offset)) {
      stop(
        "Number of values (",
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
        var = var,
        offset = offset,
        assume_rowwise = assume_rowwise,
        allow_root = FALSE,
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
  derived = NULL,
  env = rlang::caller_env()
) {
  suffix <- mapper[["labels"]][["suffix"]]
  if (
    is.null(suffix) || (is.character(suffix) && identical(suffix, ""))
  ) {
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
        input
      ),
      c(
        mapper[["labels"]][["derived"]],
        "",
        mapper[["labels"]][["root"]],
        names(input)
      )
    )
  )
}

#' @export
#' @describeIn ibm_eval_methods
#' Accepts a `state` list with named entries, one for each variable.
#' The `input` format should match the description given for [bm_expr()].
#' @param data_mask A data mask object to use for evaluating the expression. If
#'   `NULL` or missing, a data mask will be constructed from the `input`,
#'   `state`, `derived`, and `env` arguments. This can be used to avoid
#'   redundant construction of the data mask when evaluating different
#'   expressions multiple times with the same input and state.
#' @param env The environment in which to evaluate the expression. By default,
#'   this is set to the caller environment.
ibm_eval.bm_expr <- function(
  mapper,
  input,
  state = NULL,
  derived = NULL,
  ...,
  data_mask = NULL,
  env = rlang::caller_env()
) {
  if (is.null(data_mask)) {
    data_mask <- bm_expr_data_mask(
      mapper, input,
      state = state, derived = derived, env = env
    )
  }
  val <- rlang::eval_tidy(
    mapper[["expr"]],
    data = data_mask,
    env = env
  )
  val
}

#' @export
#' @rdname ibm_as_taylor
#'
ibm_as_taylor.bm_expr <- function(
  mapper,
  input,
  state,
  derived = NULL,
  jacobians = NULL,
  inla_f = FALSE,
  ...
) {
  eval2 <- ibm_eval2(
    mapper,
    input = input,
    state = state,
    derived = derived,
    jacobians = jacobians,
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


#' @export
#' @method format bm_expr
#' @rdname bm_summary
#' @examples
#' mapper <- bm_expr(~ cos(x))
#' summary(mapper)
#' summary(mapper, depth = 1)
format.bm_expr <- function(x, ...,
                           prefix = "",
                           initial = prefix,
                           depth = 1) {
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





#' @rdname ibm_eval2
#' @export
ibm_eval2.bru_obs <- function(
    mapper, input, state, ..., multi = FALSE, comp_list
) {
  derived <- list()
  jacobians <- list()
  for (nm in names(comp_list)) {
    res <- ibm_eval2(
      comp_list[[nm]],
      input = input[["comp"]][[nm]],
      state = state[[nm]],
      ...
    )
    derived[[nm]] <- res$offset
    jacobians[[nm]] <- setNames(list(res$jacobian), nm = nm)
  }

  expr_mapper <- bm_expr(
    bru_pred_expr(mapper, format = "quo"),
    labels = list(
      root = "latent",
      derived = "effects",
      suffix = "_latent"
    ),
    assume = c("rowwise", "additive")
  )

  res2 <- ibm_eval2(
    expr_mapper,
    input = list(
      data = mapper[["data"]],
      response_data = mapper[["response_data"]],
      data_extra = mapper[["data_extra"]]
    ),
    state = state,
    derived = derived,
    jacobians = jacobians,
    multi = TRUE,
    ...
  )

  # Feed forward into optional transformation mapper.
  # TODO: work in progress...
  post_mapper <- mapper[["aggregate"]]

  # This should be constructed as part of the overall bru_obs input evaluation,
  # so that we have a this pre-computed (input for ibm_eval2.bru_obs_list):
  # input = list(`obs1` = list(comp = input for components,
  #                            post = post_mapper input), `obs2` = ...)

  res3 <- ibm_eval2(
    post_mapper,
    input = input[["post"]],
    state = res2$offset,
    multi = FALSE,
    ...
  )

  # Combine jacobians from the expression mapper and the aggregate mapper.
  offset <- res3$offset
  B <- list()
  for (nm in names(state)) {
    B[[nm]] <- res2$jacobian[[nm]]
    if (!is.null(res3$jacobian)) {
      B[[nm]] <- res3$jacobian %*% B[[nm]]
    }
  }

  if (!multi) {
    B <- do.call(cbind, B)
  }

  list(offset = offset, jacobian = B)
}


#' @rdname ibm_eval2
#' @export
ibm_eval2.bru_obs_list <- function(
    mapper, input, state, ..., multi = FALSE, comp_list
) {
  results <- lapply(
    setNames(seq_along(mapper), names(mapper)),
    function(m) ibm_eval2.bru_obs(mapper[[m]], input[[m]], state, ...,
                                  multi = TRUE, comp_list = comp_list)
  )

  # Combine jacobians from the expression mapper and the aggregate mapper.
  offset <- lapply(results, function(r) r[["offset"]])
  jacobian <- lapply(results, function(r) r[["jacobian"]])

  if (!multi) {
    jac_names <- lapply(jacobian, function(r) names(r))
    jac_names <- unique(unlist(jac_names))
    names(jac_names) <- jac_names
    jac_ncol <- vapply(jac_names, function(nm) {
      max(vapply(jacobian, function(r) if (nm %in% names(r)) ncol(r[[nm]]) else 0L, 0L))
    }, 0L)
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
          }
          )
        )
      }
    )

    offset <- unlist(offset)
  }

  list(offset = offset, jacobian = jacobian)
}

#' @rdname ibm_jacobian
#' @export
ibm_jacobian.bru_obs <- function(
    mapper, input, state, ..., multi = FALSE, comp_list
) {
  result <- ibm_eval2(
    mapper,
    input = input,
    state = state,
    ...,
    multi = TRUE,
    comp_list = comp_list
  )

  B <- result$jacobian

  if (!multi) {
    B <- do.call(cbind, B)
  }

  B
}

#' @rdname ibm_eval
#' @export
ibm_eval.bru_obs <- function(
    mapper, input, state, ..., multi = FALSE, comp_list
) {
  derived <- list()
  for (nm in names(comp_list)) {
    res <- ibm_eval(
      comp_list[[nm]],
      input = input[["comp"]][[nm]],
      state = state[[nm]],
      ...
    )
    derived[[nm]] <- res
  }

  expr_mapper <- bm_expr(
    bru_pred_expr(mapper, format = "quo"),
    labels = list(
      root = "latent",
      derived = "effects",
      suffix = "_latent"
    ),
    assume = c("rowwise", "additive")
  )

  res2 <- ibm_eval(
    expr_mapper,
    input = list(data = mapper[["data"]], data_extra = mapper[["data_extra"]]),
    state = state,
    derived = derived,
    multi = FALSE,
    ...
  )

  # Feed forward into optional transformation mapper.
  # TODO: work in progress...
  post_mapper <- mapper[["aggregate_mapper"]]
  post_input <- bru_input(ibm_input_get(post_mapper),
                          data = mapper[["data"]],
                          data_extra = mapper[["data_extra"]],
                          comp = input[["comp"]]
  )

  res3 <- ibm_eval(
    post_mapper,
    input = post_input,
    state = res2$offset,
    multi = FALSE,
    ...
  )

  res3
}
