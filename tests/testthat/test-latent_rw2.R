test_that("Latent models: RW2 mapping", {
  skip_on_cran()
  local_bru_safe_inla()

  set.seed(123L)

  data1 <- data.frame(
    time = rep(c(1, 2, 4, 8, 16), times = 4),
    obs = rep(c(1, 2, 1, 4, 2), times = 4) + rnorm(20, sd = 0.5)
  )

  cmp1 <- obs ~ time(time,
    model = "rw2", values = 2^(0:4),
    constr = FALSE, scale.model = TRUE
  ) - Intercept
  fit1 <- bru(cmp1, data = data1, family = "gaussian")

  expect_equal(
    fit1$summary.random$time$mean,
    c(1.631781, 1.895681, 1.025671, 3.959789, 1.841140),
    tolerance = midtol
  )
  expect_equal(
    fit1$summary.random$time$sd,
    c(0.2350089, 0.2438991, 0.2749531, 0.2542052, 0.2473781),
    tolerance = hitol
  )
})


test_that("Latent models: RW2 mapping, data is list with different I/O sizes", {
  skip_on_cran()
  local_bru_safe_inla()

  set.seed(123L)

  data1 <- list(
    time = c(1, 2, 4, 8, 16),
    obs = rep(c(1, 2, 1, 4, 2), times = 4) + rnorm(20, sd = 0.5)
  )

  cmp1 <- obs ~ time(time,
    model = "rw2", values = 2^(0:4),
    constr = FALSE, scale.model = TRUE
  ) - Intercept
  formula <- obs ~ rep(time, times = 4)
  fit1 <- bru(cmp1,
    formula = formula, data = data1, family = "gaussian",
    allow_combine = TRUE
  )

  expect_equal(
    fit1$summary.random$time$mean,
    c(1.631781, 1.895681, 1.025671, 3.959789, 1.841140),
    tolerance = midtol
  )
  expect_equal(
    fit1$summary.random$time$sd,
    c(0.2350011, 0.2438999, 0.2749631, 0.2542013, 0.2473662),
    tolerance = hitol
  )
})
