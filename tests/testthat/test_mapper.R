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
  expect_equal(ibm_amatrix(mapper, list_data), A)
  expect_equal(ibm_amatrix(mapper, olist_data), A)
  expect_equal(ibm_amatrix(mapper, df_data), A)
  expect_equal(ibm_amatrix(mapper, matrix_data), A)
  expect_equal(ibm_amatrix(mapper, omatrix_data), A)

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

  ibm_amatrix.bm_test <- function(mapper, input, ...) {
    message("---- IBM_AMATRIX from inner environment ----")
    Matrix::sparseMatrix(
      i = seq_along(input),
      j = input,
      x = rep(2, length(input)),
      dims = c(length(input), mapper$n)
    )
  }

  bm_test <- function(n, ...) {
    bru_mapper(
      list(n = n),
      new_class = "bm_test",
      ibm_amatrix = ibm_amatrix.bm_test
    )
  }

  m <- bm_test(n = 20)
  cmp <- y ~ -1 + indep(x, model = "iid", mapper = m)
  mydata <- data.frame(y = rnorm(15) + 2 * (1:15), x = 1:15)

  skip_on_cran()
  local_bru_safe_inla()

  expect_message(
    {
      fit <- bru(cmp, data = mydata, family = "gaussian")
    },
    "---- IBM_AMATRIX from inner environment ----",
    label = "Non-interactive bru() call"
  )
})
