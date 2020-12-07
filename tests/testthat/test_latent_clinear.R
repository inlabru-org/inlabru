local_bru_testthat_setup()

test_that("bru: clinear component", {
  skip_on_cran()
  local_bru_safe_inla()

  # Seed influences data as well as predict()!
  set.seed(123)

  input.df <- data.frame(x = cos(1:100))
  input.df <- within(input.df, y <- 5 + 2 * x + rnorm(100, mean = 0, sd = 0.1))

  fit <- bru(
    y ~ myLin(main = x, model = "clinear", range = c(0, Inf)) +
      Intercept(1),
    family = "gaussian",
    data = input.df
  )

  expect_equal(fit$summary.random[["myLin"]][1, "mean"], 2.002517, tolerance = midtol)
  expect_equal(fit$summary.random[["myLin"]][1, "sd"], 0.01269209, tolerance = midtol)

  pr <- predict(
    fit,
    data.frame(x = c(1, 2)),
    ~ myLin + 2,
    n.samples = 5,
    seed = 1L
  )

  expect_equal(pr[, "mean"], c(4.0, 6.0), tolerance = midtol)
})
