test_that("fm_pixels sp vs sf", {
  skip_on_cran()
  local_bru_safe_inla()
  skip_if_not_installed("sn")

  mesh <- fm_mesh_2d_inla(cbind(0, 0),
    offset = 10, max.edge = 1,
    crs = fm_CRS("longlat_globe")
  )

  mydata <- sp::SpatialPointsDataFrame(
    mesh$loc,
    data = data.frame(y = rnorm(mesh$n) + 10),
    proj4string = fm_CRS("longlat_globe")
  )

  fit <- bru(
    ~ 0 +
      Intercept(1) +
      field(
        sp::coordinates,
        model = INLA::inla.spde2.pcmatern(
          mesh,
          prior.range = c(1, 0.01),
          prior.sigma = c(1, 0.01)
        )
      ),
    like(
      y ~ .,
      family = "gaussian",
      data = mydata
    )
  )

  system.time({
    set.seed(1234L)
    surface1 <- fm_pixels(mesh, dims = c(5, 5), mask = TRUE, format = "sp")
    density1 <- predict(fit,
      surface1,
      ~ exp(field_eval(sp::coordinates(.data.)) + Intercept),
      n.samples = 10,
      seed = 12345L
    )
  })

  system.time({
    set.seed(1234L)
    surface2 <- fm_pixels(mesh, dims = c(5, 5), mask = TRUE, format = "sf")
    density2 <- predict(fit,
      surface2,
      ~ exp(field_eval(geometry) + Intercept),
      n.samples = 10,
      seed = 12345L
    )
  })

  expect_equal(
    density1$mean,
    density2$mean
  )
})
