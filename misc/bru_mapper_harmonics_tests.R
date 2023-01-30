ibm_jacobian.bru_mapper_harmonics.sparseMatrix <- function(mapper, input, state = NULL, inla_f = FALSE, ...) {
  # Indexing into sparseMatrix is slower than into a dense Matrix,
  # so make sure we create a dense Matrix
  A <- Matrix::Matrix(0.0, NROW(input), ibm_n(mapper))
  off <- 0
  if (mapper[["intercept"]]) {
    A[, 1] <- 1.0 * mapper[["scaling"]][1]
    off <- off + 1
  }
  if (mapper[["order"]] > 0) {
    input <- (input - mapper[["interval"]][1]) / diff(mapper[["interval"]])
    for (ord in seq_len(mapper[["order"]])) {
      scale <- mapper[["scaling"]][mapper[["intercept"]] + ord]
      angle <- (2 * pi * ord) * input
      A[, off + 1] <- cos(angle) * scale
      A[, off + 2] <- sin(angle) * scale
      off <- off + 2
    }
  }
  A
}
ibm_jacobian.bru_mapper_harmonics.Matrix <- function(mapper, input, state = NULL, inla_f = FALSE, ...) {
  # Indexing into sparseMatrix is slower than into a dense Matrix,
  # so make sure we create a dense Matrix
  A <- Matrix::Matrix(1.0, NROW(input), ibm_n(mapper))
  off <- 0
  if (mapper[["intercept"]]) {
    A[, 1] <- 1.0 * mapper[["scaling"]][1]
    off <- off + 1
  }
  if (mapper[["order"]] > 0) {
    input <- (input - mapper[["interval"]][1]) / diff(mapper[["interval"]])
    for (ord in seq_len(mapper[["order"]])) {
      scale <- mapper[["scaling"]][mapper[["intercept"]] + ord]
      angle <- (2 * pi * ord) * input
      A[, off + 1] <- cos(angle) * scale
      A[, off + 2] <- sin(angle) * scale
      off <- off + 2
    }
  }
  A
}
ibm_jacobian.bru_mapper_harmonics.matrix <- function(mapper, input, state = NULL, inla_f = FALSE, ...) {
  # Indexing into sparseMatrix is slower than into a dense Matrix,
  # so make sure we create a dense Matrix
  A <- matrix(0.0, NROW(input), ibm_n(mapper))
  off <- 0
  if (mapper[["intercept"]]) {
    A[, 1] <- 1.0 * mapper[["scaling"]][1]
    off <- off + 1
  }
  if (mapper[["order"]] > 0) {
    input <- (input - mapper[["interval"]][1]) / diff(mapper[["interval"]])
    for (ord in seq_len(mapper[["order"]])) {
      scale <- mapper[["scaling"]][mapper[["intercept"]] + ord]
      angle <- (2 * pi * ord) * input
      A[, off + 1] <- cos(angle) * scale
      A[, off + 2] <- sin(angle) * scale
      off <- off + 2
    }
  }
  as(A, "Matrix")
}

library(inlabru)
NN <- c(10, 100, 1000, 10000, 10000, 100000)
timings <- list()
for (ord in 1:5) {
  timings[[ord]] <- list()
  m <- bru_mapper_harmonics(ord)
  for (Ni in seq_along(NN)) {
    N <- NN[Ni]
    input <- seq(0, 1, length.out = N)
    timings[[ord]][[Ni]] <-
      cbind(
      as.data.frame(
        bench::mark(
        "Matrix(0)" = ibm_jacobian.bru_mapper_harmonics.sparseMatrix(m, input = input),
        "Matrix(1)" = ibm_jacobian.bru_mapper_harmonics.Matrix(m, input = input),
        "matrix" = ibm_jacobian.bru_mapper_harmonics.matrix(m, input = input),
        check = FALSE
      ))[, c("expression", "itr/sec")],
      ord = ord,
      N = N)[, c("ord", "N", "expression", "itr/sec")]
    timings[[ord]][[Ni]]["expression"] <-
      attr(timings[[ord]][[Ni]][["expression"]], "description")
  }
}

timings <- do.call(rbind, lapply(timings, function(x) do.call(rbind, x)))
timings$order <- as.character(timings$ord)

library(ggplot2)
ggplot(timings, aes(N, 1/`itr/sec`, col = expression, shape = order)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  ylab("sec/itr")
