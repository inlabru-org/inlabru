local_bru_testthat_setup()

test_that("bru: factor component", {
  skip_on_cran()
  local_bru_safe_inla()

  # Required for reproducible predict() and generate() output.
  set.seed(1234L)

  input.df <- data.frame(x = cos(1:10), zz = rep(c(1, 10), each = 5))
  input.df <- within(input.df, y <- 5 + 2 * cos(1:10) +
    rnorm(10, mean = 0, sd = 1)[zz] +
    rnorm(10, mean = 0, sd = 0.1))

  # Fit a model with fixed effect 'x' and intercept 'Intercept'

  fit <- bru(y ~ x + z(zz, model = "iid", mapper = bru_mapper_index(10)),
    family = "gaussian", data = input.df
  )

  # Predict posterior statistics of 'x'

  xpost <- predict(
    fit,
    data = NULL,
    formula = ~ x_latent + fit$summary.random$z$mean[1],
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


  # Evaluate effect with _eval feature

  xpost4 <- generate(
    fit,
    data = NULL,
    formula = ~ c(
      z_eval(c(1, 2, 11, 12, 12)),
      z_eval(c(1, 2, 11, 12, 12))
    ),
    n.samples = 5,
    seed = 12345L
  )

  # The first four rows should equal the last 5 rows.
  expect_equal(
    xpost4[1:5, , drop = FALSE],
    xpost4[6:10, , drop = FALSE]
  )
  # The index 12 values should be equal.
  expect_equal(
    xpost4[4, , drop = FALSE],
    xpost4[5, , drop = FALSE]
  )
  # The columns should all be different
  expect_equal(
    sum(xpost4[, 1] == xpost4[, 2]),
    0
  )


  xpost5 <- predict(
    fit,
    data = NULL,
    formula = ~ z_eval(c(1, 11)) * Precision_for_z^0.5,
    n.samples = 500,
    seed = 12345L
  )

  # sd for z(11) should be close to 1
  expect_equal(mean(xpost5[2, "sd"]),
    1,
    tolerance = hitol
  )
})
