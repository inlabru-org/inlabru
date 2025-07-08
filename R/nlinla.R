# Linearisation ----

#' Compute inlabru model linearisation information
#'
#' @param \dots Parameters passed on to other methods
#' @export
#' @rdname bru_compute_linearisation
#' @keywords internal
bru_compute_linearisation <- function(...) {
  UseMethod("bru_compute_linearisation")
}

#' @param cmp A [bru_comp] object
#' @param model A [bru_model] object
#' @param lhood_expr A predictor expression
#' @param data Input data
#' @param data_extra Additional data for the predictor
#' @param input Precomputed component inputs from `bru_input()`
#' @param state The state information, as a list of named vectors
#' @param comp_simple Component evaluation information
#' * For `bru_comp`: A [bm_taylor] object
#' * For `bru_obs`: A [bm_list] object
#'   for the components in the likelihood
#' * For `bru_obs_list`: A list of [bm_list] objects
#' @param effects
#' * For `bru_comp`:
#' Precomputed effect list for all components involved in the likelihood
#' expression
#' @param pred0 Precomputed predictor for the given state
#' @param used A [bru_used()] object for the predictor expression
#' @param allow_combine logical; If `TRUE`, the predictor expression may
#' involve several rows of the input data to influence the same row.
#' @param eps The finite difference step size
#' @param n_pred The length of the predictor expression. If not `NULL`, scalar
#' predictor evaluations are expanded to vectors of length `n_pred`.
#'
#' @export
#' @rdname bru_compute_linearisation
bru_compute_linearisation.bru_comp <- function(cmp,
                                               model,
                                               lhood_expr,
                                               data,
                                               data_extra,
                                               input,
                                               state,
                                               comp_simple,
                                               effects,
                                               pred0,
                                               used,
                                               allow_combine,
                                               eps,
                                               n_pred = NULL,
                                               ...) {
  label <- cmp[["label"]]
  bru_log_message(
    paste0("Linearise with respect to component '", label, "'"),
    verbosity = 5
  )

  if (cmp[["main"]][["type"]] %in% c("offset", "const")) {
    # Zero-column matrix, since a const/offset has no latent state variables.
    return(Matrix::sparseMatrix(
      i = c(),
      j = c(),
      x = c(1),
      dims = c(NROW(pred0), 0)
    ))
  }

  allow_latent <- label %in% used[["latent"]]

  if (is.null(comp_simple)) {
    A <- NULL
    assume_rowwise <- FALSE
  } else {
    A <- ibm_jacobian(
      comp_simple,
      input = input[[label]],
      state = state[[label]]
    )

    assume_rowwise <- !allow_latent && !allow_combine && is.data.frame(data)

    if (assume_rowwise) {
      if (!is.null(n_pred) && (NROW(pred0) != n_pred)) {
        stop(
          "Number of rows (",
          NROW(pred0),
          ") in the predictor for component '",
          label,
          "' does not match the length implied by the response data (",
          n_pred,
          ")."
        )
      }
      if (NROW(A) == 1L) {
        A <- Matrix::kronecker(rep(1, NROW(pred0)), A)
      }
    }
  }

  triplets <- list(
    i = integer(0),
    j = integer(0),
    x = numeric(0)
  )

  if (any(!is.finite(pred0))) {
    warning(
      "Non-finite (-Inf/Inf/NaN) entries detected in predictor.\n",
      immediate. = TRUE
    )
  }

  symmetric_diffs <- FALSE
  for (k in seq_len(NROW(state[[label]]))) {
    if (is.null(A)) {
      row_subset <- seq_len(NROW(pred0))
    } else {
      Ak <- A[, k, drop = TRUE]
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
      # the rows of A have some non-zeros, or all if allow_combine
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
        if (any(!is.finite(pred_eps))) {
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
      if (any(!nonzero)) {
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
}

#' @param lhood A `bru_obs` object
#' @param model A `bru_model` object
#' @export
#' @rdname bru_compute_linearisation
bru_compute_linearisation.bru_obs <- function(lhood,
                                              model,
                                              data,
                                              input,
                                              state,
                                              comp_simple,
                                              eps,
                                              ...) {
  used <- bru_used(lhood)
  allow_combine <- lhood[["allow_combine"]]
  effects <- evaluate_effect_single_state(
    comp_simple[used[["effect"]]],
    input = input[used[["effect"]]],
    state = state[used[["effect"]]]
  )

  lhood_expr <- bru_obs_expr(lhood, model[["effects"]])
  n_pred <- bru_response_size(lhood)

  pred0 <- evaluate_predictor(
    model,
    state = list(state),
    data = data,
    data_extra = lhood[["data_extra"]],
    effects = list(effects),
    predictor = lhood_expr,
    used = used,
    format = "matrix",
    n_pred = n_pred
  )

  # Compute derivatives for each non-const/offset component
  B <- list()
  offset <- pred0
  # Either this loop or the internal bru_comp specific loop
  # can in principle be parallelised.
  for (label in union(used[["effect"]], used[["latent"]])) {
    if (ibm_n(model[["effects"]][[label]][["mapper"]]) > 0) {
      if (lhood[["is_additive"]] && !lhood[["allow_combine"]]) {
        # If additive and no combinations allowed, just need to copy the
        # non-offset A matrix, and possibly expand to full size
        A <- ibm_jacobian(
          comp_simple[[label]],
          input[[label]],
          state[[label]]
        )
        if (NROW(A) == 1) {
          if (NROW(offset) > 1) {
            B[[label]] <- Matrix::kronecker(rep(1, NROW(offset)), A)
          } else {
            B[[label]] <- A
          }
        } else {
          B[[label]] <- A
        }
      } else {
        B[[label]] <-
          bru_compute_linearisation(
            model[["effects"]][[label]],
            model = model,
            lhood_expr = lhood_expr,
            data = data,
            data_extra = lhood[["data_extra"]],
            input = input,
            state = state,
            comp_simple = comp_simple[[label]],
            effects = effects,
            pred0 = pred0,
            used = used,
            allow_combine = lhood[["allow_combine"]],
            eps = eps,
            n_pred = n_pred,
            ...
          )
      }
      if ((NROW(offset) == 1L) && (NROW(B[[label]]) > 1L)) {
        offset <- matrix(offset, NROW(B[[label]]), 1)
      }
      offset <- offset - B[[label]] %*% state[[label]]
    }
  }

  bm_taylor(offset = offset, jacobian = B, state0 = NULL)
}

#' @param lhoods A `bru_obs_list` object
#' @export
#' @rdname bru_compute_linearisation
bru_compute_linearisation.bru_obs_list <- function(lhoods,
                                                   model,
                                                   input,
                                                   state,
                                                   comp_simple,
                                                   eps = 1e-5,
                                                   ...) {
  # TODO: set the eps default more intelligently
  lapply(seq_along(lhoods), function(idx) {
    x <- lhoods[[idx]]
    bru_compute_linearisation(
      x,
      model = model,
      data = x[["data"]],
      input = input[[idx]],
      state = state,
      comp_simple = comp_simple[[idx]],
      eps = eps,
      ...
    )
  })
}

#' @export
#' @rdname bru_compute_linearisation
bru_compute_linearisation.bru_model <- function(model, lhoods,
                                                input, state,
                                                comp_simple, ...) {
  bru_compute_linearisation(lhoods,
    model = model,
    input = input, state = state,
    comp_simple = comp_simple, ...
  )
}
