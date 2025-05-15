test_that("bru: clinear component", {
  skip_on_cran()
  local_bru_safe_inla()

  # Seed influences data as well as predict()!
  withr::local_seed(123)

  input.df <- data.frame(x = cos(1:100))
  input.df <- within(input.df, y <- 5 + 2 * x + rnorm(100, mean = 0, sd = 0.1))

  fit <- bru(
    y ~ myLin(main = x, model = "clinear", range = c(0, Inf)) +
      Intercept(1),
    family = "gaussian",
    data = input.df,
    options = list(bru_initial = list(myLin = 2))
  )

  expect_equal(
    fit$summary.random[["myLin"]][1, "mean"],
    2.002517,
    tolerance = hitol
  )
  expect_equal(
    fit$summary.random[["myLin"]][1, "sd"],
    0.013,
    tolerance = hitol
  )

  skip_if_not_installed("sn")
  pr <- predict(
    fit,
    data.frame(x = c(1, 2)),
    ~ myLin + 2,
    n.samples = 5,
    seed = 1L
  )

  expect_equal(pr[, "mean"], c(4.0, 6.0), tolerance = midtol)
})
