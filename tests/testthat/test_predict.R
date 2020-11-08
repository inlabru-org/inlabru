test_that("bru: factor component", {
  skip_on_cran()
  skip_if_not(bru_safe_inla())
  
  # Required for reproducible predict() and generate() output.
  set.seed(1234L)

  input.df <- data.frame(x = cos(1:10))
  input.df <- within(input.df, y <- 5 + 2 * cos(1:10) + rnorm(10, mean = 0, sd = 0.1))

  # Fit a model with fixed effect 'x' and intercept 'Intercept'

  fit <- bru(y ~ x, family = "gaussian", data = input.df)

  # Predict posterior statistics of 'x'

  xpost <- predict(fit, data = NULL, formula = ~x_latent,
                   n.samples = 5, seed = 12345L)

  # The statistics include mean, standard deviation, the 2.5% quantile, the median,
  # the 97.5% quantile, minimum and maximum sample drawn from the posterior as well as
  # the coefficient of variation and the variance.


  # The predict function can also be used to simulataneously estimate posteriors
  # of multiple variables:

  xipost <- generate(fit,
    data = NULL,
    formula = ~ c(
      Intercept = Intercept_latent,
      x = x_latent
    ),
    n.samples = 5,
    seed = 12345L
  )

  expect_equal(is.matrix(xipost), TRUE)
  expect_equal(rownames(xipost), c("Intercept", "x"))

  xipost <- predict(fit,
    data = NULL,
    formula = ~ c(
      Intercept = Intercept_latent,
      x = x_latent
    ),
    n.samples = 5,
    seed = 12345L
  )

  expect_equal(is.data.frame(xipost), TRUE)
  expect_equal(rownames(xipost), c("Intercept", "x"))
})
