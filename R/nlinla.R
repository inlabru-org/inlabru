
# Linearisation ----

#' Compute inlabru model linearisation information
#'
#' @param \dots Parameters passed on to other methods
#' @export
#' @rdname bru_compute_linearisation
bru_compute_linearisation <- function(...) {
  UseMethod("bru_compute_linearisation")
}

#' @param cmp A [bru_component] object
#' @param model A [bru_model] object
#' @param lhood_expr A predictor expression
#' @param data Input data
#' @param input Precomputed component inputs from `evaluate_inputs()`
#' @param state The state information, as a list of named vectors
#' @param comp_simple Component evaluation information
#' * For `bru_component`: `bru_mapper_taylor` object
#' * For `bru_like`: A `comp_simple_list` object
#'   for the components in the likelihood
#' * For `bru_like_list`: A `comp_simple_list_list` object
#' @param effects
#' * For `bru_component`:
#' Precomputed effect list for all components involved in the likelihood
#' expression
#' @param pred0 Precomputed predictor for the given state
#' @param allow_latent logical. If `TRUE`, the latent state of each component is
#' directly available to the predictor expression, with a `_latent` suffix.
#' @param allow_combine logical; If `TRUE`, the predictor expression may
#' involve several rows of the input data to influence the same row.
#' @export
#' @rdname bru_compute_linearisation
bru_compute_linearisation.component <- function(cmp,
                                                model,
                                                lhood_expr,
                                                data,
                                                input,
                                                state,
                                                comp_simple,
                                                effects,
                                                pred0,
                                                allow_latent,
                                                allow_combine,
                                                eps,
                                                ...) {
  label <- cmp[["label"]]

  if (!inherits(comp_simple, "bru_mapper_taylor")) {
    warning("Non-linear component mappers not fully supported!", immediate. = TRUE)
  }
  A <- ibm_jacobian(comp_simple,
    input = input[[label]],
    state = state[[label]]
  )

  assume_rowwise <- !allow_latent && !allow_combine && is.data.frame(data)

  if (assume_rowwise) {
    if (NROW(A) == 1) {
      if (NROW(pred0) > 1) {
        A <- Matrix::kronecker(rep(1, NROW(pred0)), A)
      } else if (is.data.frame(data)) {
        A <- Matrix::kronecker(rep(1, NROW(data)), A)
        pred0 <- rep(pred0, NROW(A))
      }
    }
  }

  if (cmp[["main"]][["type"]] %in% c("offset", "const")) {
    # Zero-column matrix, since a const/offset has no latent state variables.
    return(Matrix::sparseMatrix(
      i = c(),
      j = c(),
      x = c(1),
      dims = c(NROW(pred0), 0)
    ))
  }

  triplets <- list(
    i = integer(0),
    j = integer(0),
    x = numeric(0)
  )
  for (k in seq_len(NROW(state[[label]]))) {
    row_subset <- which(A[, k] != 0.0)
    if (length(row_subset) > 0) {
      state_eps <- state
      state_eps[[label]][k] <- state[[label]][k] + eps
      # TODO:
      # Option: filter out the data and effect rows for which
      # the rows of A have some non-zeros, or all if allow_combine
      # Option: compute predictor for multiple different states. This requires
      # constructing multiple states and corresponding effects before calling
      # evaluate_predictor

      if (assume_rowwise) {
        effects_eps <- list()
        for (label_loop in names(effects)) {
          if (NROW(effects[[label_loop]]) == 1) {
            effects_eps[[label_loop]] <-
              rep(effects[[label_loop]], length(row_subset))
          } else {
            effects_eps[[label_loop]] <- effects[[label_loop]][row_subset]
          }
        }
        effects_eps[[label]] <- effects_eps[[label]] + A[row_subset, k] * eps
      } else {
        effects_eps <- effects
        effects_eps[[label]] <- effects_eps[[label]] + A[, k] * eps
      }
      pred_eps <- evaluate_predictor(
        model,
        state = list(state_eps),
        data =
          if (assume_rowwise) {
            data[row_subset, , drop = FALSE]
          } else {
            data
          },
        effects = list(effects_eps),
        predictor = lhood_expr,
        format = "matrix"
      )
      # Store sparse triplet information
      if (assume_rowwise) {
        values <- (pred_eps - pred0[row_subset])
      } else {
        values <- (pred_eps - pred0)
      }
      nonzero <- is.finite(values)
      if (any(!nonzero)) {
        warning("Non-finite (-Inf/Inf/NaN) entries detected in predictor derivatives; treated as 0.0",
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
}

#' @param lhood A `bru_like` object
#' @param model A `bru_model` object
#' @export
#' @rdname bru_compute_linearisation
bru_compute_linearisation.bru_like <- function(lhood,
                                               model,
                                               data,
                                               input,
                                               state,
                                               comp_simple,
                                               eps,
                                               ...) {
  allow_latent <- lhood[["allow_latent"]]
  allow_combine <- lhood[["allow_combine"]]
  included <- parse_inclusion(
    names(model[["effects"]]),
    lhood[["include_components"]],
    lhood[["exclude_components"]]
  )
  effects <- evaluate_effect_single_state(
    comp_simple[included],
    input = input[included],
    state = state[included],
  )

  lhood_expr <- bru_like_expr(lhood, model[["effects"]])

  pred0 <- evaluate_predictor(
    model,
    state = list(state),
    data = data,
    effects = list(effects),
    predictor = lhood_expr,
    format = "matrix"
  )
  if (lhood[["linear"]]) {
    # If linear, can check if the predictor is a scalar or vector,
    # and possibly expand to full size
    if (length(pred0) == 1) {
      if (is.data.frame(data) ||
        inherits(data, c(
          "SpatialPointsDataFrame",
          "SpatialPolygonsDataFrame",
          "SpatialLinesDataFrame"
        ))) {
        pred0 <- rep(pred0, NROW(data))
      }
    }
  }

  # Compute derivatives for each non-const/offset component
  B <- list()
  offset <- pred0
  # Either this loop or the internal bru_component specific loop
  # can in principle be parallelised.
  for (label in included) {
    if (ibm_n(model[["effects"]][[label]][["mapper"]]) > 0) {
      if (lhood[["linear"]] && !lhood[["allow_combine"]]) {
        # If linear and no combinations allowed, just need to copy the
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
            input = input,
            state = state,
            comp_simple = comp_simple[[label]],
            effects = effects,
            pred0 = pred0,
            allow_latent = lhood[["allow_latent"]],
            allow_combine = lhood[["allow_combine"]],
            eps = eps,
            ...
          )
      }
      offset <- offset - B[[label]] %*% state[[label]]
    }
  }

  bru_mapper_taylor(offset = offset, jacobian = B, state0 = NULL)
}

#' @param lhoods A `bru_like_list` object
#' @param eps The finite difference step size
#' @export
#' @rdname bru_compute_linearisation
bru_compute_linearisation.bru_like_list <- function(lhoods,
                                                    model,
                                                    input,
                                                    state,
                                                    comp_simple,
                                                    eps = 1e-5, # TODO: set more intelligently
                                                    ...) {
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
