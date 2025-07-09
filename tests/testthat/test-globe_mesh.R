test_that("2D modelling on the globe", {
  skip_on_cran()
  local_bru_safe_inla()

  withr::local_seed(123L)

  options <- list(
    control.inla = list(
      int.strategy = "eb"
    )
  )

  mesh <- fm_rcdt_2d_inla(globe = 2)

  data <- data.frame(
    Long = rep(seq(-179, 179, length.out = 10), times = 10),
    Lat = rep(seq(-89, 89, length.out = 10), each = 10)
  )
  data <- within(
    data,
    expr = {
      y <- rpois(nrow(data), exp(1 + sin(Lat / 180 * pi) * 2))
    }
  )
  data <- sf::st_as_sf(data,
    coords = c("Long", "Lat"),
    crs = fm_crs("longlat_globe")
  )
  data <- fm_transform(data, crs = fm_crs("sphere"))

  matern <- INLA::inla.spde2.pcmatern(
    mesh,
    prior.range = c(0.1, 0.01),
    prior.sigma = c(0.1, 0.01)
  )
  cmp <- y ~ mySmooth(main = geometry, model = matern) - 1

  fit <- bru(
    cmp,
    data = data,
    family = "poisson",
    options = options
  )

  expect_s3_class(fit, "bru")

  expect_equal(
    fit$summary.hyperpar["Range for mySmooth", "mean"],
    2.7076690,
    tolerance = midtol
  )

  expect_equal(
    fit$summary.hyperpar["Stdev for mySmooth", "mean"],
    0.5524327,
    tolerance = midtol
  )
})


test_that("2D LGCP modelling on the globe", {
  skip_on_cran()
  local_bru_safe_inla()

  withr::local_seed(123L)

  options <- list(
    control.inla = list(
      int.strategy = "eb"
    )
  )

  mesh <- fm_rcdt_2d_inla(globe = 2, crs = fm_crs("sphere"))

  data <- data.frame(
    Long = rep(seq(0, 360 * 9 / 10, length.out = 10), times = 10),
    Lat = 180 / pi * asin(rep(
      seq(1 / 90, 89 / 90, length.out = 10)^0.5,
      each = 10
    ))
  )
  data <- sf::st_as_sf(data,
    coords = c("Long", "Lat"),
    crs = fm_crs("longlat_globe")
  )
  data <- fm_transform(data, crs = fm_crs("sphere"))

  cmp <- geometry ~
    field(main = sf::st_coordinates(geometry), model = "fixed") +
    Intercept(1)

  expect_equal(
    sum(fm_int(mesh)$weight),
    4 * pi
  )

  fit <- bru(
    cmp,
    data = data,
    family = "cp",
    domain = list(geometry = mesh),
    options = options
  )

  expect_s3_class(fit, "bru")

  expect_equal(
    fit$summary.random$field[, "mean"],
    c(0, 0, 2.840911),
    tolerance = midtol
  )
  expect_equal(
    fit$summary.random$field[, "sd"],
    c(0.2089972, 0.2089972, 0.2968185),
    tolerance = midtol
  )
})
