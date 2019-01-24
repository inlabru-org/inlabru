context("Latent models - fixed effects (test_latent_fixed_effects.R)")

test_that("basic intercept model", {
  mycomp <- y ~ 1
  mydata <- basic_intercept_testdata()
  fit <- bru(mycomp, family = "normal",
             data = mydata)

  expect_equal(fit$summary.fixed["Intercept",]$mean, 0.090405156)
})

test_that("basic fixed effect model", {
  mycomp <- y ~ 1 + x1
  mydata <- basic_fixed_effect_testdata()
  fit <- bru(mycomp, family = "normal",
             data = mydata)

  expect_equal(fit$summary.fixed["Intercept",]$mean, 0.08537663)
})

test_that("basic fixed effect model, order relevance", {
  mydata <- basic_fixed_effect_testdata()
  mycomp1 <- y ~ Intercept + x1
  fit1 <- bru(mycomp1, family = "normal",
              data = mydata)
  mycomp2 <- y ~ x1 + Intercept
  fit2 <- bru(mycomp2, family = "normal",
              data = mydata)

  expect_equal(fit2$summary.fixed["Intercept",],
               fit1$summary.fixed["Intercept",],
               tolerance = 1e-9,
               scale = 1)
})
