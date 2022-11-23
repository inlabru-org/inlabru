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

