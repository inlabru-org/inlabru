local_bru_testthat_setup()

test_that("bru: factor component", {
  skip_on_cran()
  local_bru_safe_inla()

  # Required for reproducible predict() and generate() output.
  set.seed(1234L)

  input.df <- data.frame(x = cos(1:10))
  input.df <- within(input.df, y <- 5 + 2 * cos(1:10) + rnorm(10, mean = 0, sd = 0.1))

  # Fit a model with fixed effect 'x' and intercept 'Intercept'

  fit <- bru(y ~ x + z(1:10, model="iid"), family = "gaussian", data = input.df)

  # Predict posterior statistics of 'x'

  xpost <- predict(
    fit,
    data = NULL,
    formula = ~x_latent,
    n.samples = 5,
    seed = 12345L
  )

  xpost2 <- predict(
    fit,
    data = NULL,
    formula = ~ {
      tmp <- c(
        a = x_latent,
        b = exp(x_latent)
      )
      c(tmp, a_b = sum(tmp), c = as.vector(diff(tmp)))
    },
    n.samples = 5,
    seed = 12345L
  )
  
  xpost3 <- predict(
    fit,
    data = NULL,
    formula = ~ {
      tmp <- c(
        a = x_latent,
        b = exp(x_latent)
      )
      data.frame(A = 2, B = tmp, C = c(sum(tmp), as.vector(diff(tmp))))
    },
    n.samples = 5,
    seed = 12345L
  )
  
  # The statistics include mean, standard deviation, the 2.5% quantile, the median,
  # the 97.5% quantile, minimum and maximum sample drawn from the posterior as well as
  # the coefficient of variation and the variance.

  expect_equal(is.data.frame(xpost), TRUE)
  expect_equal(nrow(xpost), 1)

  expect_equal(is.data.frame(xpost2), TRUE)
  expect_equal(nrow(xpost2), 4)
  expect_equal(rownames(xpost2), c("a", "b", "a_b", "c"))


  # The predict function can also be used to simultaneously estimate posteriors
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
})
