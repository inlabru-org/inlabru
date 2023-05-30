qlaplace <- function(p, mean = 0, rate = 1, lower.tail = TRUE, log.p = FALSE) {
  q <- numeric(length(p))
  if (log.p) {
    pos <- p >= log(0.5)
    if (lower.tail) {
      q[!pos] <- mean - qexp(log(2) + p[!pos], lower.tail = FALSE, log.p = TRUE) / rate
      q[pos] <- mean + qexp(log(2) + log1p(-exp(p[pos])), lower.tail = FALSE, log.p = TRUE) / rate
    } else {
      q[pos] <- mean - qexp(log(2) + log1p(-exp(p[pos])), lower.tail = FALSE, log.p = TRUE) / rate
      q[!pos] <- mean + qexp(log(2) + p[!pos], lower.tail = FALSE, log.p = TRUE) / rate
    }
  } else {
    outer <- abs(p - 0.5) > 0.25
    if (lower.tail) {
      upper <- p >= 0.5
      q[!outer] <- mean - sign(p[!outer] - 0.5) * log1p(-2 * abs(p[!outer] - 0.5)) / rate
      q[outer & upper] <-
        mean - sign(p[outer & upper] - 0.5) *
        (log(2) + log1p(-p[outer & upper])) / rate
      q[outer & !upper] <-
        mean - sign(p[outer & !upper] - 0.5) *
        (log(2) + log(p[outer & !upper])) / rate
    } else {
      upper <- p >= 0.5
      q[!outer] <- mean + sign(p[!outer] - 0.5) * log1p(-2 * abs(p[!outer] - 0.5)) / rate
      q[outer & upper] <-
        mean + sign(p[outer & upper] - 0.5) *
        (log(2) + log1p(-p[outer & upper])) / rate
      q[outer & !upper] <-
        mean + sign(p[outer & !upper] - 0.5) *
        (log(2) + log(p[outer & !upper])) / rate
    }
  }
  q
}

dlaplace <- function(x, mean = 0, rate = 1, log = FALSE) {
  if (log) {
    (log(rate) - log(2)) - abs(x - mean) * rate
  } else {
    (rate / 2) * exp(-abs(x - mean) * rate)
  }
}

plaplace <- function(q, mean = 0, rate = 1, lower.tail = TRUE, log.p = FALSE) {
  p <- numeric(length(q))
  if (log.p) {
    upper <- q > mean
    if (lower.tail) {
      p[upper] <- log1p(-exp(-(q[upper] - mean) * rate) / 2)
      p[!upper] <- ((q[!upper] - mean) * rate) - log(2)
    } else {
      p[upper] <- (-(q[upper] - mean) * rate) - log(2)
      p[!upper] <- log1p(-exp((q[!upper] - mean) * rate) / 2)
    }
  } else {
    upper <- q > mean
    if (lower.tail) {
      p[upper] <- 1 - exp(-(q[upper] - mean) * rate) / 2
      p[!upper] <- exp((q[!upper] - mean) * rate) / 2
    } else {
      p[upper] <- exp(-(q[upper] - mean) * rate) / 2
      p[!upper] <- 1 - exp((q[!upper] - mean) * rate) / 2
    }
  }
  p
}

rlaplace <- function(n, mean = 0, rate = 1) {
  mean + rexp(n, rate = rate) * sign(runif(n) * 2 - 1)
}
