
nlinla.taylor <- function(expr, epunkt, data, env, offsets = c()) {
  if (nrow(epunkt) <= 1000) {
    effects <- colnames(epunkt)
    df <- as.data.frame(data)
    df <- df[, setdiff(names(df), names(epunkt)), drop = FALSE]
    wh <- cbind(df, epunkt)
    myenv <- new.env()
    invisible(lapply(colnames(wh), function(x) myenv[[x]] <- wh[, x]))
    invisible(lapply(setdiff(names(env), names(myenv)), function(x) myenv[[x]] <- env[[x]]))
    active <- setdiff(effects, offsets)
    tmp <- numericDeriv(expr[[1]], active, rho = myenv)
    gra <- attr(tmp, "gradient")
    nr <- NROW(gra)
    ngrd <- matrix(0.0, nrow = nr, ncol = length(active))
    colnames(ngrd) <- active
    for (k in seq_along(active)) {
      # as.matrix required since diag() will not work if gra has single entry
      ngrd[, k] <- diag(as.matrix(
        gra[, ((k - 1) * nr + 1):(k * nr), drop = FALSE]
      ))
    }
    nconst <- as.vector(tmp) - Matrix::rowSums(
      ngrd * as.matrix(epunkt[, active, drop = FALSE])
    )
    ngrd <- as.data.frame(ngrd)
    return(list(gradient = ngrd, const = nconst))
  } else {
    blk <- floor((seq_len(nrow(epunkt)) - 1L) / 1000)
    qq <- by(seq_len(nrow(epunkt)), blk, function(idx) {
      nlinla.taylor(
        expr,
        epunkt[idx, , drop = FALSE],
        data[idx, , drop = FALSE],
        env
      )
    })
    nconst <- do.call(c, lapply(qq, function(x) {
      x$const
    }))
    ngrd <- do.call(rbind, lapply(qq, function(x) {
      x$grad
    }))
    return(list(gradient = ngrd, const = nconst))
  }
}

nlinla.epunkt <- function(model, data, state = NULL) {
  # This function determines the current point around which
  # to perform the taylor approximation
  # (1) If result is NULL set all all effects to 0
  # (2) If result is a data.frame, use the entries as to where to approximate
  # (3) if result is an inla object, use these estimates as to where to approximate

  dfdata <- as.data.frame(data) # data as data.frame (may have been supplied as Spatial* object)
  if (is.null(state)) {
    df <- data.frame(matrix(0, nrow = nrow(dfdata), ncol = length(model$effects)))
    colnames(df) <- names(model$effects)
    df
  } else {
    evaluate_model(model, state = state, data = data)[[1]]
  }
}

nlinla.reweight <- function(A, model, data, expr, state) {
  offsets <- names(model$effects)[vapply(
    model$effects,
    function(x) identical(x$type, "offset"),
    TRUE
  )]
  epkt <- nlinla.epunkt(model, data, state = state)
  ae <- nlinla.taylor(expr, epkt, data, environment(model$formula),
    offsets = offsets
  )
  for (nm in setdiff(names(A), offsets)) {
    A[[nm]] <- A[[nm]] * ae$gradient[[nm]]
  }
  return(list(A = A, const = ae$const))
}


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
#' @param state The state information, as a list of named vectors
#' @param A A-matrix information:
#' * For `bru_component`: Precomputed A-matrix for the component
#' * For `bru_like`: A list of named A-matrices for the components in the
#' likelihood for the component
#' * For `bru_like_list`: A list, where each element is a list of named
#' A-matrices.
#' @param effects
#' * For `bru_component`:
#' Precomputed effect data.frame for all components involved in the likelihood
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
                                                state,
                                                A,
                                                effects,
                                                pred0,
                                                allow_latent,
                                                allow_combine,
                                                eps,
                                                ...) {
  label <- cmp[["label"]]
  if (identical(cmp[["main"]][["type"]], "offset")) {
    # Zero-column matrix, since an offset has no latent state variables.
    return(Matrix::sparseMatrix(
      i = c(),
      j = c(),
      x = c(1),
      dims = c(NROW(pred0), 0)
    ))
  }

  label <- cmp[["label"]]
  triplets <- list(
    i = integer(0),
    j = integer(0),
    x = numeric(0)
  )
  for (k in seq_len(NROW(state[[label]]))) {
    row_subset <- which(A[, k] != 0.0)
    if (allow_latent || (length(row_subset) > 0)) {
      state_eps <- state
      state_eps[[label]][k] <- state[[label]][k] + eps
      # TODO:
      # Option: filter out the data and effect rows for which
      # the rows of A have some non-zeros, or all if allow_combine
      # Option: compute predictor for multiple different states. This requires
      # constructing multiple states and corresponding effects before calling
      # evaluate_predictor
      if (allow_latent || allow_combine) {
        # TODO: Allow some grouping specification to allow subsetting even
        # when allow_combine is TRUE
        row_subset <- seq_len(NROW(A))
      }
      effects_eps <- effects[row_subset, , drop = FALSE]
      effects_eps[, label] <- effects_eps[, label] + A[row_subset, k] * eps

      pred_eps <- evaluate_predictor(
        model,
        state = list(state_eps),
        data = data[row_subset, , drop = FALSE],
        effects = list(effects_eps),
        predictor = lhood_expr,
        format = "matrix"
      )
      # Store sparse triplet information
      values <- (pred_eps - pred0[row_subset])
      nonzero <- (values != 0.0) # Detect exact (non)zeros
      triplets$i <- c(triplets$i, row_subset[nonzero])
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
                                               state,
                                               A,
                                               eps,
                                               ...) {
  allow_latent <- lhood[["allow_latent"]]
  allow_combine <- lhood[["allow_combine"]]
  included <- parse_inclusion(
    names(model[["effects"]]),
    lhood[["include_components"]],
    lhood[["exclude_components"]]
  )
  effects <- evaluate_effect_single(
    model[["effects"]][included],
    state = state,
    data = NULL,
    A = A
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
  # Compute derivatives for each noon-offset component
  B <- list()
  offset <- pred0
  # Either this loop or the internal bru_component specific loop
  # can in principle be parallelised.
  for (label in included) {
    if (!identical(model[["effects"]][[label]][["main"]][["type"]], "offset")) {
      # If linear, just need to copy the non-offset A matrix
      if (lhood[["linear"]]) {
        B[[label]] <- A[[label]]
      } else {
        B[[label]] <-
          bru_compute_linearisation(
            model[["effects"]][[label]],
            model = model,
            lhood_expr = lhood_expr,
            data = data,
            state = state,
            A = A[[label]],
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

  list(A = B, offset = offset)
}

#' @param lhoods A `bru_like_list` object
#' @param eps The finite difference step size
#' @export
#' @rdname bru_compute_linearisation
bru_compute_linearisation.bru_like_list <- function(
                                                    lhoods,
                                                    model,
                                                    state,
                                                    A,
                                                    eps = 1e-5, # TODO: set more intelligently
                                                    ...) {
  lapply(seq_along(lhoods), function(idx) {
    x <- lhoods[[idx]]
    bru_compute_linearisation(
      x,
      model = model,
      data = x[["data"]],
      state = state,
      A = A[[idx]],
      eps = eps,
      ...
    )
  })
}

#' @export
#' @rdname bru_compute_linearisation
bru_compute_linearisation.bru_model <- function(model, lhoods, state, A, ...) {
  bru_compute_linearisation(lhoods, model = model, state = state, A = A, ...)
}
