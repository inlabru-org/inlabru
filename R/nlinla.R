
nlinla.taylor <- function(expr, epunkt, data, env) {
  if ("offset" %in% names(epunkt)) {
    stop("One of your model components is an offset. 
          However, you are using a non-linear predictor (formula),
          which would set the respective term to zero. 
          Please remove the offset component and add its value
          to the predictor formula.")
  }

  if (nrow(epunkt) <= 1000) {
    effects <- colnames(epunkt)
    df <- as.data.frame(data)
    df <- df[, setdiff(names(df), names(epunkt)), drop = FALSE]
    wh <- cbind(df, epunkt)
    myenv <- new.env()
    invisible(lapply(colnames(wh), function(x) myenv[[x]] <- wh[, x]))
    invisible(lapply(setdiff(names(env), names(myenv)), function(x) myenv[[x]] <- env[[x]]))
    tmp <- numericDeriv(expr[[1]], effects, rho = myenv)
    gra <- attr(tmp, "gradient")
    nr <- nrow(gra)
    ngrd <- matrix(NA, nrow = nr, ncol = length(effects))
    for (k in 1:length(effects)) {
      # as.matrix required since diag() will not work if gra has single entry
      ngrd[, k] <- diag(as.matrix(gra[, ((k - 1) * nr + 1):(k * nr), drop = FALSE]))
    }
    nconst <- as.vector(tmp) - rowSums(ngrd * epunkt)
    ngrd <- data.frame(ngrd)
    colnames(ngrd) <- effects
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

nlinla.epunkt <- function(model, data, result = NULL) {
  # This function determines the current point around which
  # to perform the taylor approximation
  # (1) If result is NULL set all all effects to 0
  # (2) If result is a data.frame, use the entries as to where to approximate
  # (3) if result is an inla object, use these estimates as to where to approximate

  dfdata <- data.frame(data) # data as data.frame (may have been supplied as Spatial* object)
  if (is.null(result)) {
    df <- data.frame(matrix(0, nrow = nrow(dfdata), ncol = length(model$effects)))
    colnames(df) <- names(model$effects)
    df
  } else if (!inherits(result, "inla") & is.data.frame(result)) {
    # If result contains only a single row data frame repeat it to match the data
    if ((nrow(result) == 1) & (nrow(dfdata) > 1)) {
      result <- result[rep(1, nrow(dfdata)), , drop = FALSE]
    }
    # Check if all variables have been supplied. Those that aren't are set to 0
    for (eff in setdiff(names(model$effects), names(result))) {
      result[[eff]] <- 0
    }
    return(result)
  } else {
    evaluate_model(model, result = result, data = data, property = "mean")[[1]]
  }
}

nlinla.reweight <- function(A, model, data, expr, result) {
  epkt <- nlinla.epunkt(model, data, result = result)
  ae <- nlinla.taylor(expr, epkt, data, environment(model$formula))
  for (k in 1:length(A)) {
    nm <- names(A)[k]
    if (!(is.null(nm) || nm == "")) {
      A[[k]] <- A[[k]] * ae$gradient[[nm]]
    }
  }
  return(list(A = A, const = ae$const))
}


# Linearisation ----

#' Compute inlabru model linearisation information
#' 
#' @export
#' @rdname compute_linearisation
bru_compute_linearisation <- function(...) {
  UseMethod("compute_linearisation")
}
#' @export
#' @rdname compute_linearisation
bru_compute_linearisation.bru_component <- function(cmp,
                                               data,
                                               state,
                                               A,
                                               ...) {
}
#' @export
#' @rdname compute_linearisation
bru_compute_linearisation.bru_like <- function(lhood,
                                               model,
                                               data,
                                               state,
                                               A,
                                               ...) {
  allow_latent <- lhood[["allow_latent"]]
  allow_combine <- lhood[["allow_combine"]]
  # TODO: filter out the unused effects/latent for the likelihood
  # TODO: if linear, just need to copy the A matrices and return
  effects <- evaluate_effect_single(
    model[["effects"]],
    state = state,
    data = NULL,
    A = A
  )
  stopifnot(!is.null(lhood$expr)) # Check details for purely linear models
  pred0 <- evaluate_predictor(
    model,
    state = list(state),
    data = data,
    effects = list(effects),
    predictor = lhood$expr,
    format = "matrix"
  )
  # Compute derivatives for each component
  B <- list()
  eps <- 1e-5 # TODO: set more intelligently
  # Both of these loops could be parallelised, in principle.
  for (cmp in names(model[["effects"]])) {
    triplets <- list(
      i = integer(0),
      j = integer(0),
      x = numeric(0)
    )
    for (k in seq_len(NROW(state[[cmp]]))) {
      row_subset <- which(A[[cmp]][, k] != 0.0)
      if (allow_latent || (length(row_subset) > 0)) {
        state_eps <- state
        state_eps[[cmp]][k] <- state[[cmp]][k] + eps
        # TODO:
        # Option: filter out the data and effect rows for which
        # the rows of A have some non-zeros, or all if allow_combine
        # Option: compute predictor for multiple different states. This requires
        # constructing multiple states and corresponding effects before calling
        # evaluate_predictor
        if (allow_latent || allow_combine) {
          # TODO: Allow some grouping specification to allow subsetting even
          # when allow_combine is TRUE
          row_subset <- seq_len(NROW(A[[cmp]]))
        }
        effects_eps <- effects[row_subset, , drop = FALSE]
        effects_eps[row_subset, k] <- evaluate_effect_single(
          model[["effects"]][[cmp]],
          state = state_eps[[cmp]],
          data = NULL,
          A = A[[cmp]][row_subset, , drop = FALSE]
        )
        pred_eps <- evaluate_predictor(
          model,
          state = list(state_eps),
          data = data[row_subset, , drop = FALSE],
          effects = list(effects_eps),
          predictor = lhood$expr,
          format = "matrix"
        )
        # Store sparse triplet information
        values <- (pred_eps - pred0)
        nonzero <- (values != 0.0) # Detect exact (non)zeros
        triplets$i <- c(triplets$i, row_subset[nonzero])
        triplets$j <- c(triplets$j, rep(k, sum(nonzero)))
        triplets$x <- c(triplets$x, values[nonzero] / eps)
      }
    }
    B[[cmp]] <- Matrix::sparseMatrix(
      i = triplets$i,
      j = triplets$j,
      x = triplets$x,
      dims = c(NROW(pred0), NROW(state[[cmp]]))
    )
  }
  # TODO: compute offset
  list(A = B, offset = offset)
}
#' @export
#' @rdname compute_linearisation
bru_compute_linearisation.bru_like_list <- function(lhoods, A, ...) {
  lapply(lhoods, function(x) bru_compute_linearisation(x,
                                                       A = x[["A"]],
                                                       ...))
}
#' @export
#' @rdname compute_linearisation
bru_compute_linearisation.bru_model <- function(model, lhoods, ...) {
  bru_compute_linearisation(lhoods, model = model, ...)
}
