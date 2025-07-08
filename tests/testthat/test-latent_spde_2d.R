test_that("Georeferenced data with sf", {
  skip_on_cran()
  local_bru_safe_inla()

  withr::local_seed(123)
  mydata <- expand.grid(
    Easting = seq(5, 45, by = 20),
    Northing = seq(10, 30, by = 10),
    KEEP.OUT.ATTRS = FALSE
  )
  mydata[["obs"]] <- (mydata$Easting - 20) / 10 + rnorm(NROW(mydata))
  mydata <- sf::st_as_sf(mydata, coords = c("Easting", "Northing"))

  mesh <- fm_mesh_2d_inla(
    loc = mydata,
    offset = 5,
    max.edge = 4,
    n = 16
  )

  matern <- INLA::inla.spde2.pcmatern(
    mesh,
    prior.sigma = c(10, 0.01),
    prior.range = c(4, 0.01)
  )

  # Check that mistaken empty or unnamed arguments are detected
  cmp <- obs ~ Intercept(1) + field(geometry, model = matern, )
  expect_error(
    bru_comp_list(cmp),
    "Unnamed arguments detected in component .* position\\(s\\) 3"
  )

  cmp <- obs ~ Intercept(1) + field(geometry, model = matern)

  fit <- bru(
    cmp,
    data = mydata,
    options = list(
      control.inla = list(
        int.strategy = "eb"
      )
    )
  )

  # Check Intercept
  expect_equal(
    fit$summary.fixed["Intercept", "mean"],
    0.5398535,
    tolerance = midtol
  )

  # Check SPDE
  expect_equal(
    fit$summary.random$field$mean[mesh$idx$loc[1:3]],
    c(-2.6003077, -0.2699909, 3.5188725),
    tolerance = midtol
  )


  # Check that explicit access to the data object works
  cmp <- obs ~ Intercept(1) + field(sf::st_coordinates(.data.), model = matern)

  fit <- bru(
    cmp,
    data = mydata,
    options = list(
      control.inla = list(
        int.strategy = "eb"
      )
    )
  )

  # Check Intercept
  expect_equal(
    fit$summary.fixed["Intercept", "mean"],
    0.5398535,
    tolerance = midtol
  )

  # Check SPDE
  expect_equal(
    fit$summary.random$field$mean[mesh$idx$loc[1:3]],
    c(-2.6003077, -0.2699909, 3.5188725),
    tolerance = midtol
  )

  pred_df <- fm_pixels(mesh, dims = c(8, 8), format = "sf")

  skip_if_not_installed("sn")
  pred <- predict(fit, pred_df, ~ exp(Intercept + field), n.samples = 5)
  expect_s3_class(pred, "sf")
})


test_that("Georeferenced data with sf, with groups", {
  skip_on_cran()
  local_bru_safe_inla()

  withr::local_seed(123)
  mydata <- expand.grid(
    Easting = seq(5, 45, by = 20),
    Northing = seq(10, 30, by = 10),
    KEEP.OUT.ATTRS = FALSE
  )
  mydata <- rbind(
    cbind(mydata, season = 1L),
    cbind(mydata, season = 2L)
  )
  mydata[["obs"]] <- (mydata$Easting - 20) / 10 + rnorm(NROW(mydata))
  mydata <- sf::st_as_sf(mydata, coords = c("Easting", "Northing"))

  mesh <- fm_mesh_2d_inla(
    loc = mydata,
    offset = 5,
    max.edge = 4,
    n = 16
  )

  matern <- INLA::inla.spde2.pcmatern(
    mesh,
    prior.sigma = c(10, 0.01),
    prior.range = c(4, 0.01)
  )

  cmp <- obs ~ Intercept(1) + field(geometry,
    model = matern,
    group = season,
    ngroup = 2
  )

  fit <- bru(
    cmp,
    data = mydata,
    options = list(
      control.inla = list(
        int.strategy = "eb"
      )
    )
  )

  # Check Intercept
  expect_equal(
    fit$summary.fixed["Intercept", "mean"],
    0.6343082,
    tolerance = midtol
  )

  # Check SPDE
  expect_equal(
    fit$summary.random$field$mean[mesh$idx$loc[1:3]],
    c(-2.6947679, -0.3644125, 3.4244099),
    tolerance = midtol
  )

  pred_df <- fm_pixels(mesh, dims = c(8, 8), format = "sf")
  expect_s3_class(pred_df, "sf")
  pred_df <- fm_cprod(pred_df, data.frame(season = seq_len(2)))

  skip_if_not_installed("sn")
  pred <- predict(fit, pred_df, ~ exp(Intercept + field), n.samples = 5)
  expect_s3_class(pred, "sf")
})
