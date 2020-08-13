context("2D LGCP fitting and prediction (test_lgcp_2d.R)")
library(INLA)

test_that("2D LGCP fitting", {
  
  options <- list(control.inla = list(int.strategy = "eb",
                                      h = 0.005),
                  num.threads = 1)
  
  data <- gorillas_lgcp_2d_testdata()

  expect_s3_class(data$fit, "bru")

  # test_that("2D LGCP fitting: INLA intercept", {
  expect_equal(data$fit$summary.fixed["Intercept", "mean"], 1.121929, tolerance = hitol)
  expect_equal(data$fit$summary.fixed["Intercept", "sd"], 0.5799173, tolerance = hitol)

  index <- c(1, 456, 789, 1058, 1479)
  # test_that("2D LGCP fitting: INLA random field", {
  expect_equal(data$fit$summary.random$mySmooth$mean[index],
               c(-2.0259879,  0.3891452, -0.4462794,  0.4028431, -1.7016834),
               tolerance = midtol
  )
  expect_equal(data$fit$summary.random$mySmooth$sd[index],
    c(1.5934275, 0.8846649, 0.8297741, 0.7842688, 1.0545877),
    tolerance = hitol
  )

  # test_that("2D LGCP fitting: predicted random field", {
  loc <- SpatialPoints(data$gorillas$mesh$loc[, c(1, 2)])
  proj4string(loc) <- INLA::inla.sp_get_crs(data$gorillas$nests)
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
  fit_ips <- lgcp(data$cmp, data$gorillas$nests, ips = ips,
                  options = options)

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

  # test_that("Supplying domain definition", {
  fit_dom <- lgcp(data$cmp,
    data$gorillas$nests,
    samplers = data$gorillas$boundary,
    domain = list(coordinates = data$gorillas$mesh),
    options = options
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
