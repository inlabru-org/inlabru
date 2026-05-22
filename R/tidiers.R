#' @importFrom generics tidy glance augment
#' @importFrom tibble tibble
#' @importFrom dplyr bind_cols
#' @importFrom rlang %||%
#' @importFrom stats predict
NULL

.with_seed <- function(seed, expr) {
  if (is.null(seed)) expr else withr::with_seed(seed, expr)
}

#' Tidy a bru model fit
#'
#' @param x A fitted `bru` object.
#' @param effects `"fixed"` (default) or `"hyperpar"`.
#' @param ... Unused.
#' @return A tibble with one row per term.
#' @method tidy bru
#' @export
tidy.bru <- function(x, effects = "fixed", ...) {
  if (effects == "fixed") {
    fe <- x$summary.fixed
    if (is.null(fe)) {
      stop("No fixed effects found on bru fit.")
    }
    tibble(
      term = rownames(fe),
      estimate = fe$mean,
      std.error = fe$sd,
      conf.low = fe$`0.025quant`,
      conf.high = fe$`0.975quant`
    )
  } else if (effects == "hyperpar") {
    hp <- x$summary.hyperpar
    if (is.null(hp)) {
      stop("No hyperparameters found on bru fit.")
    }
    tibble(
      term = rownames(hp),
      estimate = hp$mean,
      std.error = hp$sd,
      conf.low = hp$`0.025quant`,
      conf.high = hp$`0.975quant`
    )
  } else {
    stop("`effects` must be \"fixed\" or \"hyperpar\".")
  }
}

#' Glance at a bru model fit
#'
#' @param x A fitted `bru` object.
#' @param ... Unused.
#' @return A one-row tibble of model-level summaries. `elapsed` is wall-clock
#'   seconds summed across the phases in `fit$bru_timings`, not CPU time.
#'   `nobs` is the total row count summed across all likelihoods, returned as
#'   `NA` for joint models that combine likelihoods of different families
#'   (where the sum is not statistically meaningful).
#' @method glance bru
#' @export
glance.bru <- function(x, ...) {
  dic_val <- x$dic$dic %||% NA_real_
  waic_val <- x$waic$waic %||% NA_real_
  mlik_val <- tryCatch(x$mlik[1, 1], error = function(e) NA_real_) %||% NA_real_

  nobs_val <- tryCatch(
    {
      if (any(bru_obs_family(x) %in% "cp")) {
        NA_integer_
      } else {
        sum(bru_response_size(x))
      }
    },
    error = function(e) NA_integer_
  )

  elapsed_val <- tryCatch(
    {
      bt <- bru_timings(x)
      if (!is.null(bt) && "Elapsed" %in% names(bt)) {
        sum(as.numeric(bt$Elapsed), na.rm = TRUE)
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_
  )

  tibble(
    dic = dic_val,
    waic = waic_val,
    marginal_loglik = mlik_val,
    nobs = nobs_val,
    elapsed = elapsed_val
  )
}

#' Augment a bru model fit with fitted values
#'
#' Adds posterior mean and credible interval columns to `data`. The user must
#' supply `pred_formula` because inlabru prediction expressions are arbitrary R
#' and cannot be recovered from the fit object alone.
#'
#' @param x A fitted `bru` object.
#' @param data A data frame of covariate values.
#' @param pred_formula A formula passed to `inlabru::predict()`.
#' @param n_samples Number of posterior samples (default 500).
#' @param seed Optional integer. If non-`NULL`, the session RNG is reset to
#'   this value before `inlabru::predict()` is called so that two identical
#'   calls return identical posterior summaries. `NULL` leaves the session
#'   RNG state unchanged (predictions then vary run-to-run because
#'   `inlabru::predict()` samples from the posterior).
#' @param ... Passed to `inlabru::predict()`.
#' @return `data` with additional columns `.fitted`, `.fitted_low`,
#'   `.fitted_high`, `.fitted_sd`.
#' @method augment bru
#' @export
augment.bru <- function(
  x,
  data,
  pred_formula,
  n_samples = 500L,
  seed = NULL,
  ...
) {
  preds <- .with_seed(
    seed,
    predict(
      x,
      newdata = data,
      formula = pred_formula,
      n.samples = n_samples,
      ...
    )
  )
  # Drop any pre-existing .fitted* columns so augment(augment(...)) does not
  # produce duplicate-name warnings or stale numbers.
  data <- tibble::as_tibble(data)
  data[, intersect(
    names(data),
    c(".fitted", ".fitted_low", ".fitted_high", ".fitted_sd")
  )] <- NULL
  bind_cols(
    data,
    tibble(
      .fitted = preds$mean,
      .fitted_low = preds$q0.025,
      .fitted_high = preds$q0.975,
      .fitted_sd = preds$sd
    )
  )
}
