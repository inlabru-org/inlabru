local_bru_testthat_setup()

test_that("2D modelling on the globe", {
  skip_on_cran()
  local_bru_safe_inla()
  skip_if_not(fm_has_PROJ6())

  set.seed(123L)

  options <- list(
    control.inla = list(
      int.strategy = "eb"
    )
  )

  mesh <- INLA::inla.mesh.create(globe = 2)

  data <- data.frame(
    Long = seq(-179, 179, length.out = 100),
    Lat = seq(-89, 89, length.out = 100)
  )
  data <- within(
    data,
    {
      y <- rpois(nrow(data), exp(1 + sin(Lat / 180 * pi) * 2))
    }
  )
  coordinates(data) <- c("Long", "Lat")
  proj4string(data) <- fm_CRS("longlat_globe")
  data <- fm_spTransform(data, CRSobj = fm_CRS("sphere"))

  matern <- INLA::inla.spde2.pcmatern(
    mesh,
    prior.range = c(0.1, 0.01),
    prior.sigma = c(0.1, 0.01)
  )
  cmp <- y ~ mySmooth(main = coordinates, model = matern) - 1

  fit <- bru(
    cmp,
    data = data,
    family = "poisson",
    options = options
  )

  expect_s3_class(fit, "bru")

  expect_equal(
    fit$summary.hyperpar["Range for mySmooth", "mean"],
    2.5367148,
    tolerance = hitol
  )

  expect_equal(
    fit$summary.hyperpar["Stdev for mySmooth", "mean"],
    0.5494477,
    tolerance = hitol
  )
})


test_that("2D LGCP modelling on the globe", {
  skip_on_cran()
  local_bru_safe_inla()
  skip_if_not(fm_has_PROJ6())
  
  set.seed(123L)

  options <- list(
    control.inla = list(
      int.strategy = "eb"
    )
  )

  mesh <- INLA::inla.mesh.create(globe = 2)

  data <- data.frame(
    Long = seq(-179, 179, length.out = 100),
    Lat = seq(-89, 89, length.out = 100)
  )
  coordinates(data) <- c("Long", "Lat")
  proj4string(data) <- fm_CRS("longlat_globe")
  data <- fm_spTransform(data, CRSobj = fm_CRS("sphere"))

  matern <- INLA::inla.spde2.pcmatern(
    mesh,
    prior.range = c(0.1, 0.01),
    prior.sigma = c(0.1, 0.01)
  )
  cmp <- coordinates ~ mySmooth(main = coordinates, model = matern) - 1

  expect_equal(
    sum(ipoints(mesh)$weight),
    4 * pi
  )

  fit <- bru(
    cmp,
    data = data,
    family = "cp",
    domain = list(coordinates = mesh),
    options = options
  )

  expect_s3_class(fit, "bru")

  expect_equal(
    fit$summary.hyperpar["Range for mySmooth", "mean"],
    1.5192092,
    tolerance = hitol
  )

  expect_equal(
    fit$summary.hyperpar["Stdev for mySmooth", "mean"],
    0.8297119,
    tolerance = hitol
  )
})
