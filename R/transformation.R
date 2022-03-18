#' @title Transformation tools
#' @description Tools for transforming between N(0,1) variables and other
#' distributions in predictor expressions
#' @param qfun A quantile function object, such as `qexp`
#' @param x Values to be transformed
#' @param ... Distribution parameters passed on to the `qfun` and `pfun` functions
#' @param tail.split. For x-values larger than `tail.split.`, upper quantile calculations
#' are used internally, and for smaller values lower quantile calculations are used. This
#' can avoid lack of accuracy in the distribution tails. If `NULL`, forward calculations split at 0,
#' and inverse calculations use lower tails only, potentially losing accuracy in the upper tails.
#' @return * For `bru_forward_transformation`, a numeric vector
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @export
#' @rdname bru_transformation

bru_forward_transformation <- function(qfun, x, ..., tail.split. = 0) {
  if (is.null(tail.split.)) {
    # By default, split at 0
    upper <- x >= 0
  } else {
    upper <- x >= tail.split.
  }
  res <- numeric(length(x))
  if (sum(upper) > 0) {
    res[upper] <-
      qfun(pnorm(x[upper],
        lower.tail = FALSE,
        log.p = TRUE
      ),
      ...,
      lower.tail = FALSE,
      log.p = TRUE
      )
  }
  if (sum(!upper) > 0) {
    res[!upper] <-
      qfun(pnorm(x[!upper],
        lower.tail = TRUE,
        log.p = TRUE
      ),
      ...,
      lower.tail = TRUE,
      log.p = TRUE
      )
  }
  res
}

#' @param pfun A CDF function object, such as `pexp`
#' @param x Values to be transformed
#' @return * For `bru_inverse_transformation`, a numeric vector
#' @export
#' @rdname bru_transformation

bru_inverse_transformation <- function(pfun, x, ..., tail.split. = NULL) {
  if (is.null(tail.split.)) {
    # By default, upper = FALSE
    upper <- logical(length(x))
  } else {
    upper <- x >= tail.split.
  }
  res <- numeric(length(x))
  if (sum(upper) > 0) {
    res[upper] <-
      qnorm(pfun(x[upper],
        ...,
        lower.tail = FALSE,
        log.p = TRUE
      ),
      lower.tail = FALSE,
      log.p = TRUE
      )
  }
  if (sum(!upper) > 0) {
    res[!upper] <-
      qnorm(pfun(x[!upper],
        ...,
        lower.tail = TRUE,
        log.p = TRUE
      ),
      lower.tail = TRUE,
      log.p = TRUE
      )
  }
  res
}




# p = 0.5 + 0.5 * sign(q) * (1 - exp(-abs(q)))
#   = 0.5 + 0.5 * sign(q) - 0.5 * sign(q) * exp(-abs(q))
# 2p = 1 + sign(q) - sign(q) * exp(-abs(q))
# 2p - sign(2p-1) -1 = -sign(2p-1) * exp(-abs(q))
# -(2p-1-sign(2p-1))/sign(2p-1) = exp(-abs(q))
# -(2p-1)/sign(2p-1)+1 = exp(-abs(q))
# 1 - |2p-1| = exp(-abs(q))
# -log(1 - |2p-1|) = abs(q)
# -sign(2p-1) * log(1 - |2p-1|) = q
# -sign(2p-1) * log1p(- |2p-1|) = q
#
# p <= 0.5: log(2p) = log(2) + log(p)
# p >= 0.5: -log(2-2p) = log(2) + log(1-p) = q

qlaplace <- function(p, lower.tail = TRUE, log.p = FALSE) {
  q <- numeric(length(p))
  if (lower.tail) {
    if (log.p) {
      upper <- p >= log(1 / 2)
      q[upper] <- log(2) + p[upper]
      q[!upper] <- -log(2) - log1p(-exp(p[!upper]))
    } else {
      upper <- p >= 0.5
      q[upper] <- log(2) + log(p[upper])
      q[!upper] <- -log(2) - log1p(-p[!upper])
    }
  } else {
    if (log.p) {
      upper <- p >= log(1 / 2)
      q[upper] <- -log(2) - log1p(-exp(p[upper]))
      q[!upper] <- log(2) + p[!upper]
    } else {
      upper <- p >= 0.5
      q[upper] <- -log(2) - log1p(-p[upper])
      q[!upper] <- log(2) + log(p[!upper])
    }
  }
  q
}
