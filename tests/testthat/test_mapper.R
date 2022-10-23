local_bru_testthat_setup()

test_that("Multi-mapper bru input", {
  skip_on_cran()
  local_bru_safe_inla()
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

  cmp1 <- y ~ indep(list(time = time, space = space),
    model = "ar1", mapper = mapper
  ) - 1
  fit1 <- bru(cmp1, family = "gaussian", data = data)

  cmp2 <- y ~ indep(list(space, time),
    model = "ar1", mapper = mapper
  ) - 1
  fit2 <- bru(cmp2, family = "gaussian", data = data)

  cmp3 <- y ~ indep(data.frame(time = time, space = space),
    model = "ar1", mapper = mapper
  ) - 1
  fit3 <- bru(cmp3, family = "gaussian", data = data)

  cmp4 <- y ~ indep(cbind(space, time),
    model = "ar1", mapper = mapper
  ) - 1
  fit4 <- bru(cmp4, family = "gaussian", data = data)

  expect_equal(fit1$summary.random, fit2$summary.random)
  expect_equal(fit1$summary.random, fit3$summary.random)
  expect_equal(fit1$summary.random, fit4$summary.random)
})


# Define in outer environment:


test_that("User defined mappers", {
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

  bm_test <- function(n, ...) {
    bru_mapper_define(
      list(n = n),
      new_class = "bm_test",
      methods = list(
        ibm_jacobian = ibm_jacobian.bm_test
      )
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


test_that("User defined mappers 2", {
  # .S3method was unavailable in R 3.6!
  skip_if_not(utils::compareVersion("4", R.Version()$major) <= 0)

  # User defined mapper objects

  ibm_jacobian.bm_test <- function(mapper, input, ...) {
    message("---- IBM_JACOBIAN from inner environment 2 ----")
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
    "---- IBM_JACOBIAN from inner environment 2 ----",
    label = "ibm_jacobian generic call"
  )

  skip_on_cran()
  local_bru_safe_inla()

  expect_message(
    object = {
      fit <- bru(cmp, data = mydata, family = "gaussian")
    },
    "---- IBM_JACOBIAN from inner environment 2 ----",
    label = "Non-interactive bru() call"
  )
})



test_that("mapper collection direct construction consistency", {
  skip_on_cran()
  local_bru_safe_inla()
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


test_that("mapper collection automatic construction consistency", {
  local_bru_safe_inla()

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
  # inla.mesh.1d mapper

  lik <- like(formula = y ~ ., data = data)

  cmp1 <- component_list(cmp1, lhoods = list(lik))
  cmp2 <- component_list(cmp2, lhoods = list(lik))

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
      as.matrix(ibm_jacobian(cmp1$indep$mapper, input = input, inla_f = inla_f)),
      as.matrix(ibm_jacobian(cmp2$indep$mapper, input = input, inla_f = inla_f))
    )
  }
})
