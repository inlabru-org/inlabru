context("Latent models - linear (test_latent_linear.R)")

test_that("bru: linear component", {
  skip_on_cran()
  skip_if_not(bru_safe_inla())

  # Seed influences data as well as predict()!
  set.seed(123)

  input.df <- data.frame(x = cos(1:100))
  input.df <- within(input.df, y <- 5 + 2 * x + rnorm(100, mean = 0, sd = 0.1))

  fit <- bru(y ~ myLin(main = x, model = "linear") + Intercept,
    family = "gaussian", data = input.df,
    options = list(num.threads = "1:1")
  )

  expect_equal(fit$summary.fixed["myLin", "mean"], 2.02, midtol)
  expect_equal(fit$summary.fixed["myLin", "sd"], 0.0126, midtol)

  pr <- predict(
    fit,
    data.frame(x = c(1, 2)),
    ~ myLin + 1,
    n.samples = 5,
    seed = 1L
  )

  expect_equal(pr[, "mean"], c(3.005013, 5.010026), midtol)
})
