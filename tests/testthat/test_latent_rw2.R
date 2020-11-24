local_bru_testthat_setup()

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
    c(0.2365901, 0.2463170, 0.2807306, 0.2577222, 0.2494937),
    tolerance = midtol
  )
})
