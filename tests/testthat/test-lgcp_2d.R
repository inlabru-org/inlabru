test_that("2D LGCP fitting", {
  skip_on_cran()
  local_bru_safe_inla()

  set.seed(123L)

  options <- list(
    control.inla = list(
      int.strategy = "eb",
      h = 0.005
    )
  )

  data(gorillas, package = "inlabru", envir = environment())

  matern <- INLA::inla.spde2.pcmatern(gorillas$mesh,
    prior.sigma = c(0.1, 0.01),
    prior.range = c(5, 0.01)
  )
  cmp <- coordinates ~ mySmooth(main = coordinates, model = matern) +
    Intercept(1)

  fit <- lgcp(
    cmp,
    data = gorillas$nests,
    samplers = gorillas$boundary,
    domain = list(coordinates = gorillas$mesh),
    options = options
  )

  expect_s3_class(fit, "bru")

  # test_that("2D LGCP fitting: INLA intercept", {
  expect_equal(
    fit$summary.fixed["Intercept", "mean"],
    1.1109,
    tolerance = midtol
  )
  expect_equal(
    fit$summary.fixed["Intercept", "sd"],
    0.579,
    tolerance = midtol
  )

  index <- c(1, 456, 789, 1058, 1479)
  # test_that("2D LGCP fitting: INLA random field", {
  expect_snapshot_value(
    fit$summary.random$mySmooth$mean[index],
    tolerance = midtol,
    style = "serialize"
  )
  expect_snapshot_value(
    fit$summary.random$mySmooth$sd[index],
    tolerance = midtol,
    style = "serialize"
  )
  expect_equal(
    fit$summary.hyperpar["Range for mySmooth", "mean"],
    2.1137,
    tolerance = midtol
  )
  expect_equal(
    fit$summary.hyperpar["Stdev for mySmooth", "mean"],
    1.1049,
    tolerance = midtol
  )

  # test_that("2D LGCP fitting: predicted random field", {
  loc <- SpatialPoints(gorillas$mesh$loc[, c(1, 2)],
    proj4string = fm_CRS(gorillas$nests)
  )
  set.seed(123L)
  skip_if_not_installed("sn")
  pr <- predict(fit, loc, ~mySmooth,
    n.samples = 5, seed = 5657L,
    parallel.configs = FALSE
  )
  # Prediction variability includes reordering differences, so need large
  # tolerances unless n.samples is large
  expect_snapshot_value(
    pr$mean[c(1, 255, 778, 1000)],
    tolerance = 1,
    style = "serialize"
  )
  expect_snapshot_value(
    pr$sd[c(2, 215, 656, 1010)],
    tolerance = 1,
    style = "serialize"
  )

  # test_that("2D LGCP fitting: predicted intensity integral", {
  ips <- fm_int(gorillas$mesh, gorillas$boundary)
  set.seed(123L)
  Lambda <- predict(fit, ips, ~ sum(weight * exp(mySmooth + Intercept)),
    n.samples = 10, seed = 5657L
  )

  expect_equal(
    Lambda$mean,
    670,
    tolerance = max(Lambda$mean.mc_std_err) * 6 / 670
  )
  expect_equal(
    Lambda$sd,
    27,
    tolerance = max(Lambda$sd.mc_std_err) * 6 / 27
  )

  # test_that("Supplying integration points instead of samplers&domains", {
  ips <- fm_int(gorillas$mesh, gorillas$boundary)
  fit_ips <- lgcp(
    cmp,
    gorillas$nests,
    ips = ips,
    options = options
  )

  expect_equal(
    fit_ips$summary.fixed["Intercept", "mean"],
    fit$summary.fixed["Intercept", "mean"],
    tolerance = midtol
  )
  expect_equal(
    fit_ips$summary.fixed["Intercept", "sd"],
    fit$summary.fixed["Intercept", "sd"],
    tolerance = midtol
  )
  expect_equal(
    fit_ips$summary.random$mySmooth$mean,
    fit$summary.random$mySmooth$mean,
    tolerance = midtol
  )
  expect_equal(
    fit_ips$summary.random$mySmooth$sd,
    fit$summary.random$mySmooth$sd,
    tolerance = midtol
  )
})
