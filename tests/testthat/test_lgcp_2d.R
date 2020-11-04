context("2D LGCP fitting and prediction (test_lgcp_2d.R)")

test_that("2D LGCP fitting", {
  skip_on_cran()
  skip_if_not(bru_safe_inla())
  disable_PROJ6_warnings()

  options <- list(
    control.inla = list(
      int.strategy = "eb",
      h = 0.005
    ),
    num.threads = "1:1"
  )

  data <- gorillas_lgcp_2d_testdata()

  expect_s3_class(data$fit, "bru")

  # test_that("2D LGCP fitting: INLA intercept", {
  expect_equal(data$fit$summary.fixed["Intercept", "mean"], 1.121929, tolerance = midtol)
  expect_equal(data$fit$summary.fixed["Intercept", "sd"], 0.5799173, tolerance = midtol)

  index <- c(1, 456, 789, 1058, 1479)
  # test_that("2D LGCP fitting: INLA random field", {
  expect_equal(data$fit$summary.random$mySmooth$mean[index],
    c(-2.0212259, 0.3862700, -0.4485562, 0.4033954, -1.7001208),
    tolerance = midtol
  )
  expect_equal(data$fit$summary.random$mySmooth$sd[index],
    c(1.5936344, 0.8840181, 0.8289204, 0.7830659, 1.0539057),
    tolerance = midtol
  )
  expect_equal(data$fit$summary.hyperpar["Range for mySmooth", "mean"],
    2.117453,
    tolerance = midtol
  )
  expect_equal(data$fit$summary.hyperpar["Stdev for mySmooth", "mean"],
    1.105217,
    tolerance = midtol
  )

  # test_that("2D LGCP fitting: predicted random field", {
  loc <- SpatialPoints(data$gorillas$mesh$loc[, c(1, 2)],
    proj4string = INLA::inla.sp_get_crs(data$gorillas$nests)
  )
  pr <- predict(data$fit, loc, ~mySmooth, n.samples = 500, seed = 5657L)
  expect_equal(pr$mean[c(1, 255, 778, 1000)],
    c(-2.057675, -1.766163, -1.512785, -1.488362),
    tolerance = hitol
  )
  expect_equal(pr$sd[c(2, 215, 656, 1010)],
    c(1.5174210, 0.6629913, 0.9822118, 1.2876455),
    tolerance = hitol
  )

  # test_that("2D LGCP fitting: predicted intensity integral", {
  ips <- ipoints(data$gorillas$boundary, data$gorillas$mesh)
  Lambda <- predict(data$fit, ips, ~ sum(weight * exp(mySmooth + Intercept)),
    n.samples = 500, seed = 5657L
  )

  expect_equal(Lambda$mean, 647.4751, tolerance = hitol)
  expect_equal(Lambda$sd, 25.54122, tolerance = hitol)

  # test_that("Supplying integration points instead of samplers", {
  ips <- ipoints(data$gorillas$boundary, data$gorillas$mesh)
  fit_ips <- lgcp(data$cmp, data$gorillas$nests,
    ips = ips,
    options = options
  )

  expect_equal(fit_ips$summary.fixed["Intercept", "mean"],
    data$fit$summary.fixed["Intercept", "mean"],
    tolerance = lowtol
  )
  expect_equal(fit_ips$summary.fixed["Intercept", "sd"],
    data$fit$summary.fixed["Intercept", "sd"],
    tolerance = lowtol
  )
  expect_equal(fit_ips$summary.random$mySmooth$mean,
    data$fit$summary.random$mySmooth$mean,
    tolerance = lowtol
  )
  expect_equal(fit_ips$summary.random$mySmooth$sd,
    data$fit$summary.random$mySmooth$sd,
    tolerance = lowtol
  )

  # test_that("Supplying domain definition", {
  fit_dom <- lgcp(data$cmp,
    data$gorillas$nests,
    samplers = data$gorillas$boundary,
    domain = list(coordinates = data$gorillas$mesh),
    options = options
  )

  expect_equal(fit_dom$summary.fixed["Intercept", "mean"],
    data$fit$summary.fixed["Intercept", "mean"],
    tolerance = lowtol
  )
  expect_equal(fit_dom$summary.fixed["Intercept", "sd"],
    data$fit$summary.fixed["Intercept", "sd"],
    tolerance = lowtol
  )
  expect_equal(fit_dom$summary.random$mySmooth$mean,
    data$fit$summary.random$mySmooth$mean,
    tolerance = lowtol
  )
  expect_equal(fit_dom$summary.random$mySmooth$sd,
    data$fit$summary.random$mySmooth$sd,
    tolerance = lowtol
  )
})
