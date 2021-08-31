local_bru_testthat_setup()

test_that("basic intercept model", {
  skip_on_cran()
  local_bru_safe_inla()
  options <- list(
    control.inla = list(h = 0.005)
  )
  mycomp <- y ~ 1
  mydata <- local_basic_intercept_testdata()
  fit <- bru(mycomp,
    family = "normal",
    data = mydata,
    options = options
  )

  expect_equal(
    fit$summary.fixed["Intercept", ]$mean,
    0.090405156,
    tolerance = lowtol
  )
})

test_that("basic intercept model, spatial data", {
  skip_on_cran()
  local_bru_safe_inla()
  options <- list(
    control.inla = list(h = 0.005)
  )
  mycomp <- y ~ 1
  mydata <- local_basic_intercept_testdata()
  mydata$coord1 <- 11
  mydata$coord2 <- 12
  coordinates(mydata) <- c("coord1", "coord2")

  fit <- bru(mycomp,
    family = "normal",
    data = mydata,
    options = options
  )

  expect_equal(
    fit$summary.fixed["Intercept", ]$mean,
    0.090405156,
    tolerance = lowtol
  )
})

test_that("basic fixed effect model", {
  skip_on_cran()
  local_bru_safe_inla()
  options <- list(
    control.inla = list(h = 0.005)
  )
  mycomp <- y ~ 1 + x1
  mydata <- local_basic_fixed_effect_testdata()
  fit <- bru(mycomp,
    family = "normal",
    data = mydata,
    options = options
  )

  expect_equal(
    fit$summary.fixed["Intercept", ]$mean,
    0.08537663,
    tolerance = lowtol
  )
})

test_that("basic fixed effect model, order relevance", {
  skip_on_cran()
  local_bru_safe_inla()
  options <- list(
    control.inla = list(h = 0.005)
  )
  mydata <- local_basic_fixed_effect_testdata()
  mycomp1 <- y ~ Intercept + x1
  fit1 <- bru(mycomp1,
    family = "normal",
    data = mydata,
    options = options
  )
  mycomp2 <- y ~ x1 + Intercept
  fit2 <- bru(mycomp2,
    family = "normal",
    data = mydata,
    options = options
  )

  expect_equal(
    fit2$summary.fixed["Intercept", ],
    fit1$summary.fixed["Intercept", ],
    tolerance = lowtol
  )
})
