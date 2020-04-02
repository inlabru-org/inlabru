context("2D LGCP fitting and prediction (test_lgcp_2d.R)")
library(INLA)

data <- gorillas_lgcp_2d_testdata()

test_that("2D LGCP fitting: Result object", {
  expect_is(data$fit, "bru")
})


test_that("2D LGCP fitting: INLA intercept", {
  expect_equal(data$fit$summary.fixed["Intercept", "mean"], 1.121929, tolerance = midtol)
  expect_equal(data$fit$summary.fixed["Intercept", "sd"], 0.5799173, tolerance = midtol)
})

test_that("2D LGCP fitting: INLA random field", {
  expect_equal(data$fit$summary.random$mySmooth$mean[c(1, 456, 789, 1058, 1479)],
    c(-2.0224597, 0.3874104, -0.4473572, 0.4019972, -1.7000660),
    tolerance = midtol
  )
  expect_equal(data$fit$summary.random$mySmooth$sd[c(1, 436, 759, 1158, 1279)],
    c(1.5924485, 0.8243210, 0.8209047, 0.7928983, 1.0671142),
    tolerance = midtol
  )
})

test_that("2D LGCP fitting: predicted random field", {
  loc <- SpatialPoints(gorillas$mesh$loc[, c(1, 2)])
  proj4string(loc) <- CRS(proj4string(gorillas$nests))
  pr <- predict(data$fit, loc, ~mySmooth, n.samples = 500, seed = 5657L)
  expect_equal(pr$mean[c(1, 255, 778, 1000)],
    c(-2.057675, -1.766163, -1.512785, -1.488362),
    tolerance = hitol
  )
  expect_equal(pr$sd[c(2, 215, 656, 1010)],
    c(0.0000000, 0.6629913, 0.9822118, 1.2876455),
    tolerance = hitol
  )
})

test_that("2D LGCP fitting: predicted intensity integral", {
  ips <- ipoints(data$gorillas$boundary, data$gorillas$mesh)
  Lambda <- predict(data$fit, ips, ~ sum(weight * exp(mySmooth + Intercept)),
    n.samples = 500, seed = 5657L
  )

  expect_equal(Lambda$mean, 647.4751, tolerance = hitol)
  expect_equal(Lambda$sd, 25.54122, tolerance = hitol)
})

test_that("Supplying integration points instead of samplers", {
  ips <- ipoints(data$gorillas$boundary, data$gorillas$mesh)
  fit_ips <- lgcp(data$cmp, data$gorillas$nests, ips = ips)

  expect_equal(fit_ips$summary.fixed["Intercept", "mean"],
    data$fit$summary.fixed["Intercept", "mean"],
    tolerance = midtol
  )
  expect_equal(fit_ips$summary.fixed["Intercept", "sd"],
    data$fit$summary.fixed["Intercept", "sd"],
    tolerance = hitol
  )
  expect_equal(fit_ips$summary.random$mySmooth$mean,
    data$fit$summary.random$mySmooth$mean,
    tolerance = midtol
  )
  expect_equal(fit_ips$summary.random$mySmooth$sd,
    data$fit$summary.random$mySmooth$sd,
    tolerance = hitol
  )
})


test_that("Supplying domain definition", {
  fit_dom <- lgcp(data$cmp,
    data$gorillas$nests,
    samplers = data$gorillas$boundary,
    domain = list(coordinates = data$gorillas$mesh),
    options = list(control.inla = list(int.strategy = "eb"))
  )

  expect_equal(fit_dom$summary.fixed["Intercept", "mean"],
    data$fit$summary.fixed["Intercept", "mean"],
    tolerance = midtol
  )
  expect_equal(fit_dom$summary.fixed["Intercept", "sd"],
    data$fit$summary.fixed["Intercept", "sd"],
    tolerance = midtol
  )
  expect_equal(fit_dom$summary.random$mySmooth$mean,
    data$fit$summary.random$mySmooth$mean,
    tolerance = midtol
  )
  expect_equal(fit_dom$summary.random$mySmooth$sd,
    data$fit$summary.random$mySmooth$sd,
    tolerance = midtol
  )
})
