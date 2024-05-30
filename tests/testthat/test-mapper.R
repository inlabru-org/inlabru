test_that("Linear mapper", {
  mapper <- bru_mapper_linear()

  input <- seq_len(4)
  state <- 10
  expect_equal(ibm_n(mapper), 1)
  expect_equal(ibm_n_output(mapper, input = input), length(input))
  expect_equal(ibm_values(mapper), 1)
  expect_equal(ibm_values(mapper, inla_f = TRUE), 1)
  expect_equal(
    ibm_eval(mapper, input = input, state = state),
    input * state
  )
  expect_equal(
    ibm_eval(
      ibm_linear(mapper, input = input, state = state),
      state = state
    ),
    input * state
  )
})


test_that("Index mapper", {
  values <- seq_len(4)
  state <- c(10, 20, 30, 40)
  mapper <- bru_mapper_index(max(values))
  input <- c(2, 2, 1, 3)
  val <- state[input]

  expect_equal(ibm_n(mapper), length(values))
  expect_equal(ibm_n_output(mapper, input = input), length(input))
  expect_equal(ibm_values(mapper), values)
  expect_equal(ibm_values(mapper, inla_f = TRUE), values)
  expect_equal(
    ibm_eval(mapper, input = input, state = state),
    val
  )
  expect_equal(
    ibm_eval(
      ibm_linear(mapper, input = input, state = state),
      state = state
    ),
    val
  )
})


test_that("Factor mapper", {
  all_values <-
    list(
      full = c("a", "b", "c"),
      contrast = c("b", "c")
    )
  all_state <-
    list(
      full = c(a = 10, b = 30, c = 60),
      contrast = c(b = 20, c = 50)
    )
  for (factor_mapping in rev(c("full", "contrast"))) {
    values <- all_values[[factor_mapping]]
    mapper <- bru_mapper_factor(all_values[["full"]],
      factor_mapping = factor_mapping
    )

    input <- c("b", "b", "a", "c")
    state <- all_state[[factor_mapping]]
    val <- {
      val <- rep(0, length(input))
      ok <- input %in% names(state)
      val[ok] <- state[input[ok]]
      val
    }

    expect_equal(ibm_n(mapper), length(values))
    expect_equal(ibm_n_output(mapper, input = input), length(input))
    expect_equal(ibm_values(mapper), values)
    expect_equal(ibm_values(mapper, inla_f = TRUE), values)
    expect_equal(
      ibm_eval(mapper, input = input, state = state),
      val
    )
    expect_equal(
      ibm_eval(
        ibm_linear(mapper, input = input, state = state),
        state = state
      ),
      val
    )
  }
})


test_that("Pipe mapper", {
  mapper <-
    bru_mapper_pipe(
      list(
        bru_mapper_multi(list(
          space = bru_mapper_index(4),
          time = bru_mapper_index(3)
        )),
        bru_mapper_scale()
      )
    )

  expect_equal(ibm_n(mapper), 12)
  expect_equal(
    ibm_n(mapper, multi = TRUE),
    list(space = 4, time = 3)
  )
  expect_equal(ibm_values(mapper), seq_len(12))
  expect_equal(
    as.data.frame(ibm_values(mapper, multi = TRUE)),
    expand.grid(space = seq_len(4), time = seq_len(3)),
    ignore_attr = TRUE
  )

  olist_data <-
    list(
      list(space = 2:4, time = 1:3),
      2
    )
  A <- Matrix::sparseMatrix(
    i = 1:3,
    j = c(2, 7, 12),
    x = 1 * 2,
    dims = c(3, 12)
  )
  expect_equal(ibm_jacobian(mapper, olist_data), A)
})


test_that("Aggregate mapper", {
  mapper <- bru_mapper_aggregate(rescale = FALSE)

  expect_equal(ibm_n(mapper), NA_integer_)
  expect_equal(
    ibm_n(mapper, state = c(1, 1)),
    2
  )
  expect_equal(
    ibm_n_output(mapper, input = list(block = c(1, 2, 2))),
    2
  )
  expect_equal(ibm_values(mapper, state = c(1, 1)), seq_len(2))

  input <-
    list(block = c(1, 2, 2, 1, 3), weights = c(1, 1, 1, 2, 3))
  state <- c(1, 2, 3, 4, 5)

  val <- c(9, 5, 15)
  A <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2, 3),
    j = c(1, 4, 2, 3, 5),
    x = c(1, 2, 1, 1, 3),
    dims = c(3, 5)
  )
  expect_equal(ibm_eval(mapper, input = input, state = state), val)
  expect_equal(ibm_jacobian(mapper, input = input, state = state), A)


  mapper <- bru_mapper_aggregate(rescale = TRUE)

  expect_equal(ibm_n(mapper), NA_integer_)
  expect_equal(
    ibm_n(mapper, state = c(1, 1)),
    2
  )
  expect_equal(
    ibm_n_output(mapper, input = list(block = c(1, 2, 2))),
    2
  )
  expect_equal(ibm_values(mapper, state = c(1, 1)), seq_len(2))

  val <- c(3, 2.5, 5)
  A <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2, 3),
    j = c(1, 4, 2, 3, 5),
    x = c(1 / 3, 2 / 3, 1 / 2, 1 / 2, 3 / 3),
    dims = c(3, 5)
  )
  expect_equal(ibm_eval(mapper, input = input, state = state), val)
  expect_equal(ibm_jacobian(mapper, input = input, state = state), A)

  set.seed(123L)
  delta <- (runif(length(state)) * 2 - 1) * 1e-6
  num_deriv <-
    (ibm_eval(mapper, input = input, state = state + delta) -
      ibm_eval(mapper, input = input, state = state - delta)) /
      (2 * sum(delta^2)^0.5)
  jac_deriv <-
    as.vector(ibm_jacobian(mapper, input = input, state = state) %*%
      (delta / sum(delta^2)^0.5))
  expect_equal(
    num_deriv,
    jac_deriv,
    tolerance = lowtol
  )
})


test_that("logsumexp mapper", {
  mapper <- bru_mapper_logsumexp(rescale = FALSE)

  expect_equal(ibm_n(mapper), NA_integer_)
  expect_equal(
    ibm_n(mapper, state = c(1, 1)),
    2
  )
  expect_equal(
    ibm_n_output(mapper, input = list(block = c(1, 2, 2))),
    2
  )
  expect_equal(ibm_values(mapper, state = c(1, 1)), seq_len(2))

  state <- c(1, 2, 3, 4, 5)

  # Test default input handling
  # Should treat NULL group and weights as all-1, but input needs to be supplied
  expect_equal(
    ibm_eval(mapper, input = NULL, state = state),
    5.451914,
    tolerance = lowtol
  )
  expect_equal(
    ibm_eval(mapper, input = NULL, state = state),
    5.451914,
    tolerance = lowtol
  )


  input <-
    list(block = c(1, 2, 2, 1, 3), weights = c(1, 1, 1, 2, 3))

  val <- c(4.717736, 3.313262, 6.098612)
  A <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2, 3),
    j = c(1, 4, 2, 3, 5),
    x = c(0.0242889, 0.9757111, 0.2689414, 0.7310586, 1),
    dims = c(3, 5)
  )
  expect_equal(
    ibm_eval(mapper, input = input, state = state),
    val,
    tolerance = lowtol
  )
  expect_equal(
    ibm_eval(mapper, input = input, state = state, log = FALSE),
    exp(val),
    tolerance = lowtol
  )
  expect_equal(
    ibm_jacobian(mapper, input = input, state = state),
    A,
    tolerance = midtol
  )


  mapper <- bru_mapper_logsumexp(rescale = TRUE)

  expect_equal(ibm_n(mapper), NA_integer_)
  expect_equal(
    ibm_n(mapper, state = c(1, 1)),
    2
  )
  expect_equal(
    ibm_n_output(mapper, input = list(block = c(1, 2, 2))),
    2
  )
  expect_equal(ibm_values(mapper, state = c(1, 1)), seq_len(2))

  val <- c(3.619124, 2.620115, 5)
  A <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2, 3),
    j = c(1, 4, 2, 3, 5),
    x = c(0.0242889, 0.9757111, 0.2689414, 0.7310586, 1),
    dims = c(3, 5)
  )
  expect_equal(
    ibm_eval(mapper, input = input, state = state),
    val,
    tolerance = lowtol
  )
  expect_equal(
    ibm_eval(mapper, input = input, state = state, log = FALSE),
    exp(val),
    tolerance = lowtol
  )
  expect_equal(
    ibm_jacobian(mapper, input = input, state = state),
    A,
    tolerance = midtol
  )

  set.seed(123L)
  delta <- (runif(length(state)) * 2 - 1) * 1e-6
  num_deriv <-
    (ibm_eval(mapper, input = input, state = state + delta) -
      ibm_eval(mapper, input = input, state = state - delta)) /
      (2 * sum(delta^2)^0.5)
  jac_deriv <-
    as.vector(ibm_jacobian(mapper, input = input, state = state) %*%
      (delta / sum(delta^2)^0.5))
  expect_equal(
    num_deriv,
    jac_deriv,
    tolerance = lowtol
  )
})


test_that("Multi-mapper bru input", {
  mapper <- bru_mapper_multi(list(
    space = bru_mapper_index(4),
    time = bru_mapper_index(3)
  ))
  expect_equal(ibm_n(mapper), 12)
  expect_equal(ibm_n(mapper, multi = 1), list(space = 4, time = 3))
  expect_equal(ibm_values(mapper), seq_len(12))
  expect_equal(
    as.data.frame(ibm_values(mapper, multi = 1)),
    expand.grid(space = seq_len(4), time = seq_len(3)),
    ignore_attr = TRUE
  )

  list_data <- list(time = 1:3, space = 2:4)
  olist_data <- list(space = 2:4, time = 1:3)
  df_data <- as.data.frame(list_data)
  matrix_data <- cbind(time = 1:3, space = 2:4)
  omatrix_data <- cbind(2:4, 1:3)
  A <- Matrix::sparseMatrix(
    i = 1:3,
    j = c(2, 7, 12),
    x = 1,
    dims = c(3, 12)
  )
  expect_equal(ibm_jacobian(mapper, list_data), A)
  expect_equal(ibm_jacobian(mapper, olist_data), A)
  expect_equal(ibm_jacobian(mapper, df_data), A)
  expect_equal(ibm_jacobian(mapper, matrix_data), A)
  expect_equal(ibm_jacobian(mapper, omatrix_data), A)

  data <- cbind(df_data, y = 1:3)
})




test_that("User defined mappers", {
  # .S3method was unavailable in R 3.6!
  skip_if_not(utils::compareVersion("4", R.Version()$major) <= 0)

  # User defined mapper objects

  ibm_jacobian.bm_test <- function(mapper, input, ...) {
    message("---- IBM_JACOBIAN from inner environment ----")
    Matrix::sparseMatrix(
      i = seq_along(input),
      j = input,
      x = rep(2, length(input)),
      dims = c(length(input), mapper$n)
    )
  }
  .S3method("ibm_jacobian", "bm_test", ibm_jacobian.bm_test)

  bm_test <- function(n, ...) {
    bru_mapper_define(
      list(n = n),
      new_class = "bm_test"
    )
  }

  m <- bm_test(n = 20)
  cmp <- y ~ -1 + indep(x, model = "iid", mapper = m)
  mydata <- data.frame(y = rnorm(15) + 2 * (1:15), x = 1:15)
  expect_message(
    object = {
      ibm_jacobian(m, mydata$x)
    },
    "---- IBM_JACOBIAN from inner environment ----",
    label = "ibm_jacobian generic call"
  )

  skip_on_cran()
  local_bru_safe_inla()

  expect_message(
    object = {
      fit <- bru(cmp, data = mydata, family = "gaussian")
    },
    "---- IBM_JACOBIAN from inner environment ----",
    label = "Non-interactive bru() call"
  )
})



test_that("Collect mapper, direct construction", {
  skip_on_cran()
  set.seed(1234L)

  mapper <- bru_mapper_collect(
    list(
      u = bru_mapper_index(4),
      v = bru_mapper_index(4)
    ),
    hidden = TRUE
  )
  expect_equal(ibm_n(mapper), 8)
  expect_equal(ibm_n(mapper, inla_f = TRUE), 4)
  expect_equal(ibm_n(mapper, multi = 1), list(u = 4, v = 4))
  expect_equal(ibm_values(mapper), seq_len(8))
  expect_equal(ibm_values(mapper, inla_f = TRUE), seq_len(4))
  expect_equal(
    as.data.frame(ibm_values(mapper, multi = 1)),
    list(u = seq_len(4), v = seq_len(4)),
    ignore_attr = TRUE
  )

  list_data <- list(u = 1:3, v = 2:4)
  A <- Matrix::bdiag(
    Matrix::sparseMatrix(
      i = 1:3,
      j = 1:3,
      x = 1,
      dims = c(3, 4)
    ),
    Matrix::sparseMatrix(
      i = 1:3,
      j = 2:4,
      x = 1,
      dims = c(3, 4)
    )
  )
  A <- as(as(as(A, "dMatrix"), "generalMatrix"), "TsparseMatrix")
  expect_equal(ibm_jacobian(mapper, list_data), A)
  expect_equal(
    as(
      ibm_jacobian(mapper, list_data[["u"]], inla_f = TRUE)[
        , ibm_inla_subset(mapper),
        drop = FALSE
      ],
      "dgTMatrix"
    ),
    as(A[
      seq_along(list_data[["u"]]),
      ibm_inla_subset(mapper),
      drop = FALSE
    ], "dgTMatrix")
  )
})


test_that("Collect mapper, automatic construction", {
  skip_on_cran()

  data <- data.frame(val = 1:3, y = 1:3)

  mapper <- bru_mapper_collect(
    list(
      u = bru_mapper_index(4),
      v = bru_mapper_index(4)
    ),
    hidden = TRUE
  )

  cmp1 <- y ~
    -1 +
    indep(val,
      model = "bym",
      mapper = mapper,
      graph = Matrix::Diagonal(4) + 1
    )
  # index mapper

  cmp2 <- y ~
    -1 +
    indep(val,
      model = "bym",
      n = 4,
      graph = Matrix::Diagonal(4) + 1
    )
  # fm_mesh_1d mapper

  lik <- like(formula = y ~ ., data = data)

  # These tests do not trigger INLA:inla.mesh.1d usage:
  cmp1 <- component_list(cmp1, lhoods = like_list(list(lik)))
  cmp2 <- component_list(cmp2, lhoods = like_list(list(lik)))

  for (inla_f in c(FALSE, TRUE)) {
    expect_identical(
      ibm_n(cmp1$indep$mapper, inla_f = inla_f),
      ibm_n(cmp2$indep$mapper, inla_f = inla_f)
    )
    expect_identical(
      ibm_values(cmp1$indep$mapper, inla_f = inla_f),
      ibm_values(cmp2$indep$mapper, inla_f = inla_f)
    )
    if (inla_f) {
      input <- list(mapper = list(
        main = data$val,
        group = rep(1, 3),
        replicate = rep(1, 3)
      ))
    } else {
      input <- list(mapper = list(
        main = list(u = data$val),
        group = rep(1, 3),
        replicate = rep(1, 3)
      ))
    }
    expect_identical(
      as.matrix(ibm_jacobian(cmp1$indep$mapper,
        input = input,
        inla_f = inla_f
      )),
      as.matrix(ibm_jacobian(cmp2$indep$mapper,
        input = input,
        inla_f = inla_f
      ))
    )
  }
})



test_that("Collect mapper works", {
  skip_on_cran()
  local_bru_safe_inla()

  graph <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 2, 3, 4),
    j = c(2, 3, 4, 1, 2, 3),
    x = rep(1, 6),
    dims = c(4, 4)
  )

  set.seed(12345L)
  data <- data.frame(
    y = 4 + rnorm(6),
    x = c(1, 2, 3, 2, 3, 4)
  )

  fit_inla <- INLA::inla(y ~ 1 + f(x, model = "bym", graph = graph),
    data = data
  )

  expect_error(
    {
      fit_bru <-
        bru(y ~ Intercept(1) + field(x, model = "bym", graph = graph),
          data = data,
          options = list(bru_initial = list(field = rep(10, 8)))
        )
    },
    NA
  )
})



test_that("Marginal mapper", {
  m1 <- bru_mapper_marginal(qexp, pexp, rate = 1 / 8)
  state0 <- -5:5
  val1 <- ibm_eval(m1, state = state0)

  state1 <- ibm_eval(m1, state = val1, reverse = TRUE)
  expect_equal(state1, state0)

  m2 <- bru_mapper_marginal(qexp, pexp, rate = 1 / 8, inverse = TRUE)
  state2 <- ibm_eval(m2, state = val1)
  expect_equal(state2, state0)

  val2 <- ibm_eval(m2, state = state2, reverse = TRUE)
  expect_equal(val2, val1)

  dqexp <- function(p, rate = 1, lower.tail = TRUE, log.p = FALSE, log = FALSE) {
    if (log.p) {
      if (lower.tail) {
        val <- log1p(-exp(p)) + log(rate)
      } else {
        val <- p + log(rate)
      }
    } else {
      if (lower.tail) {
        val <- log1p(-p) + log(rate)
      } else {
        val <- log(p) + log(rate)
      }
    }
    if (log) {
      return(val)
    }
    return(exp(val))
  }
  m1_d <- bru_mapper_marginal(qexp, pexp, dexp, rate = 1 / 8)
  m1_dq <- bru_mapper_marginal(qexp, pexp, NULL, dqexp, rate = 1 / 8)
  jac1 <- ibm_jacobian(m1, state = state0)
  jac1_d <- ibm_jacobian(m1_d, state = state0)
  jac1_dq <- ibm_jacobian(m1_dq, state = state0)
  # rbind(diag(jac1), diag(jac1_d), diag(jac1_dq))
  expect_equal(sum(abs(diag(jac1) - diag(jac1_d))), 0, tolerance = lowtol)
  expect_equal(sum(abs(diag(jac1_d) - diag(jac1_dq))), 0, tolerance = lowtol)
})


test_that("Mapper lists", {
  m1 <- bru_mapper_index(4L)
  m2 <- bru_mapper_index(3L)

  expect_s3_class(c(m1), "bm_list")
  expect_s3_class(c(m1, m2), "bm_list")
  expect_s3_class(c(c(m1, m2), c(m2, m1)), "bm_list")
})

test_that("Mesh 1d mapper", {
  m <- bru_mapper(fmesher::fm_mesh_1d(1:5, boundary = "f"), indexed = TRUE)

  loc <- seq(-2, 7, by = 0.5)
  val <- ibm_eval(m, input = loc, state = seq_len(ibm_n(m)))
  expect_length(val, length(loc))
  expect_equal(val, loc)
})

test_that("Mesh 2d mapper", {
  m <- bru_mapper(fmesher::fmexample$mesh)

  loc <- fm_pixels(fmesher::fmexample$mesh, dims = c(5, 5), mask = FALSE)
  val <- ibm_eval(m, input = loc, state = seq_len(ibm_n(m)))
  expect_length(val, 25)

  loc <- fm_pixels(fmesher::fmexample$mesh, dims = c(5, 5), mask = TRUE)
  val <- ibm_eval(m, input = loc, state = seq_len(ibm_n(m)))
  expect_length(val, 9)
})
