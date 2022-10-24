local_bru_testthat_setup()

test_that("2D LGCP fitting", {
  skip_on_cran()
  local_bru_safe_inla()

  # test_that("2D LGCP fitting: Factor covariate (as SpatialPixelsDataFrame)", {
  data(gorillas, package = "inlabru", envir = environment())

  # Uses the component label to pick the covariate layer to extract,
  # so doesn't need an explicit main_layer="vegetation".

  mdl <- coordinates ~ vegetation(
    main = gorillas$gcov$vegetation,
    model = "iid"
  ) - Intercept

  fit <- lgcp(
    mdl, gorillas$nests,
    samplers = gorillas$boundary,
    domain = list(coordinates = gorillas$mesh),
    options = list(
      control.inla = list(
        int.strategy = "eb",
        h = 0.005
      )
    )
  )

  expect_snapshot_value(
    fit$summary.random$veg$mean,
    tolerance = midtol,
    style = "serialize"
  )
  expect_snapshot_value(
    fit$summary.random$veg$sd,
    tolerance = midtol,
    style = "serialize"
  )

  # test_that("2D LGCP fitting: Continuous covariate (as function)", {
  elev <- gorillas$gcov$elevation
  elev$elevation <- elev$elevation - mean(elev$elevation, na.rm = TRUE)

  mdl2 <- coordinates ~ beta.elev(
    main = elev,
    main_layer = "elevation",
    model = "linear"
  ) + Intercept(1)
  fit2 <- lgcp(mdl2, gorillas$nests,
    samplers = gorillas$boundary,
    domain = list(coordinates = gorillas$mesh),
    options = list(
      control.inla = list(
        int.strategy = "eb",
        h = 0.005
      )
    )
  )

  expect_equal(fit2$summary.fixed["beta.elev", "mean"], 0.004192824, tolerance = midtol)
  expect_equal(fit2$summary.fixed["beta.elev", "sd"], 0.00249103, tolerance = midtol)
  expect_equal(fit2$summary.fixed["Intercept", "mean"], 3.069781, tolerance = midtol)
  expect_equal(fit2$summary.fixed["Intercept", "sd"], 0.05587102, tolerance = midtol)

  f.elev <- function(x, y) {
    spp <- SpatialPoints(data.frame(x = x, y = y),
      proj4string = INLA::inla.sp_get_crs(elev)
    )
    v <- over(spp, elev)
    v[is.na(v)] <- 0 # NAs are a problem! Remove them
    return(v$elevation)
  }

  mdl3 <- coordinates ~ beta.elev(main = f.elev(x, y), model = "linear") + Intercept(1)
  fit3 <- lgcp(mdl3, gorillas$nests,
    samplers = gorillas$boundary,
    domain = list(coordinates = gorillas$mesh),
    options = list(
      control.inla = list(
        int.strategy = "eb",
        h = 0.005
      )
    )
  )

  expect_equal(fit2$summary.fixed["beta.elev", "mean"], 0.004192824, tolerance = midtol)
  expect_equal(fit2$summary.fixed["beta.elev", "sd"], 0.00249103, tolerance = midtol)
  expect_equal(fit2$summary.fixed["Intercept", "mean"], 3.069781, tolerance = midtol)
  expect_equal(fit2$summary.fixed["Intercept", "sd"], 0.05587102, tolerance = midtol)
})
