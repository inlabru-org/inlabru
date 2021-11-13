local_bru_testthat_setup()

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
  cmp <- coordinates ~ mySmooth(main = coordinates, model = matern) + Intercept(1)

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
    1.0912,
    tolerance = midtol
  )
  expect_equal(
    fit$summary.fixed["Intercept", "sd"],
    0.5774,
    tolerance = midtol
  )

  index <- c(1, 456, 789, 1058, 1479)
  # test_that("2D LGCP fitting: INLA random field", {
  expect_equal(
    fit$summary.random$mySmooth$mean[index],
    c(-1.7432, 0.5305, -0.2831, 0.5830, -1.4833),
    tolerance = midtol
  )
  expect_equal(
    fit$summary.random$mySmooth$sd[index],
    c(1.5914, 0.8826, 0.8273, 0.7815, 1.0522),
    tolerance = midtol
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
    proj4string = INLA::inla.sp_get_crs(gorillas$nests)
  )
  set.seed(123L)
  pr <- predict(fit, loc, ~mySmooth,
    n.samples = 5, seed = 5657L,
    parallel.configs = FALSE
  )
  # Prediction variability includes reordeing differences, so need large
  # tolerances unless n.samples is large
  expect_equal(
    pr$mean[c(1, 255, 778, 1000)],
    c(-2.7279196, -1.7653978, -1.7801908, -0.6127033),
    tolerance = 1
  )
  expect_equal(
    pr$sd[c(2, 215, 656, 1010)],
    c(0.5892231, 0.4663549, 1.0070714, 1.2699464),
    tolerance = 1
  )

  # test_that("2D LGCP fitting: predicted intensity integral", {
  ips <- ipoints(gorillas$boundary, gorillas$mesh)
  set.seed(123L)
  Lambda <- predict(fit, ips, ~ sum(weight * exp(mySmooth + Intercept)),
    n.samples = 10, seed = 5657L
  )

  expect_equal(
    Lambda$mean,
    669.5,
    tolerance = midtol
  )
  expect_equal(
    Lambda$sd,
    40,
    tolerance = hitol * 5
  )

  # test_that("Supplying integration points instead of samplers&domains", {
  ips <- ipoints(gorillas$boundary, gorillas$mesh)
  fit_ips <- lgcp(
    cmp,
    gorillas$nests,
    ips = ips,
    options = options
  )

  expect_equal(
    fit_ips$summary.fixed["Intercept", "mean"],
    fit$summary.fixed["Intercept", "mean"],
    tolerance = lowtol
  )
  expect_equal(
    fit_ips$summary.fixed["Intercept", "sd"],
    fit$summary.fixed["Intercept", "sd"],
    tolerance = lowtol
  )
  expect_equal(
    fit_ips$summary.random$mySmooth$mean,
    fit$summary.random$mySmooth$mean,
    tolerance = lowtol
  )
  expect_equal(
    fit_ips$summary.random$mySmooth$sd,
    fit$summary.random$mySmooth$sd,
    tolerance = lowtol
  )
})
