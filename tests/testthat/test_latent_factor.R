context("Latent models - factor (test_latent_factor.R)")

test_that("bru: factor component", {

  # Seed influences data as well as predict()!
  set.seed(123)

  # Factorial data
  input.df <- data.frame(x = c(rep("Alpha", 100), rep("Beta", 100)))
  input.df$y <- c(rep(1, 100), rep(-2, 100)) + rnorm(200, mean = 0, sd = 0.1)

  # Fit the model
  fit <- bru(y ~ fac(map = x, model = "factor") - Intercept,
    family = "gaussian",
    data = input.df
  )

  # Check fixed effect results
  expect_equal(fit$summary.fixed[1, "mean"], 1.009040, midtol)
  expect_equal(fit$summary.fixed[1, "sd"], 0.009721, midtol)
  expect_equal(fit$summary.fixed[2, "mean"], -2.01075, midtol)
  expect_equal(fit$summary.fixed[2, "sd"], 0.009721, midtol)

  # Check if prediction works
  pr <- predict(fit, data.frame(x = c("Alpha", "Beta")), ~ fac + 1, seed = 1)
  expect_equal(pr[, "mean"], c(2.008706, -1.010480), hitol)
})
