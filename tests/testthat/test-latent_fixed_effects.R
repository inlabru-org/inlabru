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
    0.09140515,
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
    0.09140515,
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
    0.08637662,
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
  mycomp1 <- y ~ Intercept(1) + x1
  fit1 <- bru(mycomp1,
    family = "normal",
    data = mydata,
    options = options
  )
  mycomp2 <- y ~ x1 + Intercept(1)
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

test_that("interaction fixed effect model", {
  skip_on_cran()
  local_bru_safe_inla()
  options <- list()
  mydata <- local_basic_fixed_effect_testdata()
  set.seed(123L)
  mydata <- cbind(mydata, x2 = sample(x = factor(c("A", "B")), size = nrow(mydata), replace = TRUE))
  mycomp1 <- y ~ -1 + mix(~ -1 + x1:x2, model = "fixed")
  fit0 <- INLA::inla(y ~ -1 + x1:x2, data = mydata, family = "normal")
  fit1 <- bru(mycomp1,
    family = "normal",
    data = mydata,
    options = options
  )
  expect_equal(
    fit1$summary.random$mix$ID,
    rownames(fit0$summary.fixed)
  )
  expect_equal(
    fit1$summary.random$mix$ID,
    c("x1:x2A", "x1:x2B")
  )
  expect_equal(
    fit1$summary.random$mix$mean,
    fit0$summary.fixed$mean,
    tolerance = lowtol
  )
  expect_equal(
    fit1$summary.random$mix$sd,
    fit0$summary.fixed$sd,
    tolerance = midtol
  )
})
