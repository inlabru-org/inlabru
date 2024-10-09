test_that("Linear predictor offset", {
  skip_on_cran()
  local_bru_safe_inla()

  dat <- data.frame(
    id = c(1:8),
    deaths = c(5, 30, 2, 4, 5, 4, 7, 10),
    pop = c(1000, 2300, 300, 400, 500, 700, 1000, 700)
  )

  # Used to fail
  m1 <- bru(
    deaths ~ 1 + myoffset(log(pop), model = "offset"),
    data = dat,
    family = "poisson"
  )
  m2 <- bru(
    ~ 1 + myoffset(log(pop), model = "offset"),
    formula = deaths ~ Intercept + myoffset,
    data = dat,
    family = "poisson"
  )

  expect_equal(
    m1$summary.fixed$mean,
    m2$summary.fixed$mean,
    ignore_attr = TRUE,
    tolerance = lowtol
  )
  expect_equal(
    m1$summary.fixed$sd,
    m2$summary.fixed$sd,
    ignore_attr = TRUE,
    tolerance = midtol
  )

  # Give an error if an 'offset' option is specified.
  # Note the since 2.11.1.9013 there is no need to make the have the full
  # predictor length for a non-linear model.
  expect_error(
    {
      m3 <- bru(
        ~ Intercept(1),
        formula = deaths ~ Intercept,
        data = dat,
        family = "poisson",
        options = list(offset = log(dat$pop))
      )
    },
    "An offset option was specified which may interfere"
  )
})
