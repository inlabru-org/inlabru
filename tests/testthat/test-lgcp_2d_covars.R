test_that("2D LGCP fitting", {
  skip_on_cran()
  local_bru_safe_inla()

  # test_that("2D LGCP fitting: Factor covariate (as SpatialPixelsDataFrame)", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  gorillas <- gorillas_sf
  gorillas$gcov <- gorillas_sf_gcov()

  # Uses the component label to pick the covariate layer to extract,
  # so doesn't need an explicit main_layer="vegetation".

  mdl <- geometry ~ vegetation(
    main = gorillas$gcov$vegetation,
    model = "iid"
  ) - Intercept

  fit <- lgcp(
    mdl, gorillas$nests,
    samplers = gorillas$boundary,
    domain = list(geometry = gorillas$mesh),
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
  elev <- elev - mean(terra::values(elev), na.rm = TRUE)

  mdl2 <- geometry ~ beta.elev(
    main = elev,
    main_layer = "elevation",
    model = "linear"
  ) + Intercept(1)
  fit2 <- lgcp(mdl2, gorillas$nests,
    samplers = gorillas$boundary,
    domain = list(geometry = gorillas$mesh),
    options = list(
      control.inla = list(
        int.strategy = "eb",
        h = 0.005
      )
    )
  )

  expect_equal(fit2$summary.fixed["beta.elev", "mean"], 0.004192824,
    tolerance = midtol
  )
  expect_equal(fit2$summary.fixed["beta.elev", "sd"], 0.00249103,
    tolerance = midtol
  )
  expect_equal(fit2$summary.fixed["Intercept", "mean"], 3.069781,
    tolerance = midtol
  )
  expect_equal(fit2$summary.fixed["Intercept", "sd"], 0.05587102,
    tolerance = midtol
  )

  f.elev <- function(xy) {
    eval_spatial(elev, xy)
  }

  mdl3 <- geometry ~ beta.elev(
    main = f.elev(.data.$geometry),
    model = "linear"
  ) +
    Intercept(1)
  fit3 <- lgcp(mdl3, gorillas$nests,
    samplers = gorillas$boundary,
    domain = list(geometry = gorillas$mesh),
    options = list(
      control.inla = list(
        int.strategy = "eb",
        h = 0.005
      )
    )
  )

  expect_equal(fit2$summary.fixed["beta.elev", "mean"], 0.004192824,
    tolerance = midtol
  )
  expect_equal(fit2$summary.fixed["beta.elev", "sd"], 0.00249103,
    tolerance = midtol
  )
  expect_equal(fit2$summary.fixed["Intercept", "mean"], 3.069781,
    tolerance = midtol
  )
  expect_equal(fit2$summary.fixed["Intercept", "sd"], 0.05587102,
    tolerance = midtol
  )
})

# > bru_timings(fit)
# Task Iteration       Time     System    Elapsed
# user.self Preprocess         0 0.195 secs 0.002 secs 0.222 secs
# 1         Preprocess         1 0.454 secs 0.002 secs 0.479 secs
# 2         Run inla()         1 0.866 secs 0.016 secs 0.707 secs
# > bru_timings(fit2)
# Task Iteration       Time     System    Elapsed
# user.self Preprocess         0 0.420 secs 0.006 secs 0.430 secs
# 1         Preprocess         1 0.267 secs 0.008 secs 0.272 secs
# 2         Run inla()         1 0.695 secs 0.018 secs 0.439 secs
# > bru_timings(fit3)
# Task Iteration       Time     System    Elapsed
# user.self Preprocess         0 0.067 secs 0.002 secs 0.068 secs
# 1         Preprocess         1 0.327 secs 0.002 secs 0.328 secs
# 2         Run inla()         1 0.619 secs 0.016 secs 0.407 secs

# > bru_timings(fit)
# Task Iteration       Time     System    Elapsed
# 1 Preprocess         0 0.449 secs 0.022 secs 0.514 secs
# 2 Preprocess         1 0.199 secs 0.062 secs 0.231 secs
# 3 Run inla()         1 0.962 secs 0.044 secs 0.753 secs
# > bru_timings(fit2)
# Task Iteration       Time     System    Elapsed
# 1 Preprocess         0 0.507 secs 0.028 secs 0.525 secs
# 2 Preprocess         1 0.220 secs 0.000 secs 0.219 secs
# 3 Run inla()         1 0.648 secs 0.012 secs 0.405 secs
# > bru_timings(fit3)
# Task Iteration       Time     System    Elapsed
# 1 Preprocess         0 0.137 secs 0.021 secs 0.148 secs
# 2 Preprocess         1 0.060 secs 0.000 secs 0.059 secs
# 3 Run inla()         1 0.626 secs 0.020 secs 0.455 secs
