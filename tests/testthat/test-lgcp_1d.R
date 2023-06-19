test_data <- function() {
  data(Poisson2_1D, package = "inlabru", envir = environment())
  x <- seq(0, 55, length = 50)
  mesh1D <- INLA::inla.mesh.1d(x, boundary = "free", degree = 2)
  matern <- INLA::inla.spde2.pcmatern(
    mesh1D,
    prior.range = c(150, 0.75),
    prior.sigma = c(0.1, 0.75),
    constr = TRUE
  )
  mdl <- x ~ spde1D(main = x, model = matern) + Intercept(1)
  fit <- lgcp(
    mdl,
    pts2,
    domain = list(x = mesh1D),
    options = list(
      control.inla = list(int.strategy = "eb")
    )
  )
  list(
    mesh1D = mesh1D,
    fit = fit
  )
}


test_that("1D LGCP fitting", {
  skip_on_cran()
  local_bru_safe_inla()

  result <- test_data()
  mesh1D <- result$mesh1D
  fit <- result$fit

  # Needed for reproducible predict
  set.seed(123L)

  expect_s3_class(fit, "bru")

  expect_equal(
    fit$summary.fixed["Intercept", "mean"],
    0.60974,
    tolerance = midtol
  )
  expect_equal(
    fit$summary.fixed["Intercept", "sd"],
    0.10916,
    tolerance = midtol
  )

  expect_snapshot_value(
    fit$summary.random$spde1D$mean[c(1, 27, 50)],
    tolerance = midtol,
    style = "deparse"
  )
  expect_snapshot_value(
    fit$summary.random$spde1D$sd[c(2, 32, 29)],
    tolerance = midtol,
    style = "deparse"
  )

  skip_if_not_installed("sn")

  # predicted intensity integral
  ips <- fm_int(mesh1D, name = "x")
  Lambda <- predict(
    fit,
    ips,
    ~ sum(weight * exp(spde1D + Intercept)),
    n.samples = 100,
    seed = 4354
  )
  expect_equal(Lambda$mean, 130.25, tolerance = hitol)
  expect_equal(Lambda$sd, 12.28, tolerance = 1)
})





test_data_discrete <- function() {
  data(Poisson2_1D, package = "inlabru", envir = environment())
  xx <- ceiling(pts2$x)
  data <- data.frame(x = xx)
  x <- seq(1, 55, length = 55)
  mesh1D <- INLA::inla.mesh.1d(x, boundary = "free")
  matern <- INLA::inla.spde2.pcmatern(mesh1D,
    prior.range = c(0.01, 0.01),
    prior.sigma = c(1, 0.01),
    constr = TRUE
  )
  mdl <- x ~ spde1D(main = x, model = matern) + Intercept(1)
  fit <- lgcp(mdl,
    data = data,
    domain = list(x = x),
    options = list(
      control.inla = list(int.strategy = "eb")
    )
  )
  list(
    mesh1D = mesh1D,
    fit = fit
  )
}


test_that("1D LGCP fitting, discrete point domain", {
  skip_on_cran()
  local_bru_safe_inla()

  result <- test_data_discrete()
  mesh1D <- result$mesh1D
  fit <- result$fit

  # Needed for reproducible predict
  set.seed(123L)

  expect_s3_class(fit, "bru")

  expect_equal(
    fit$summary.fixed["Intercept", "mean"],
    0.605155,
    tolerance = hitol
  )
  expect_equal(
    fit$summary.fixed["Intercept", "sd"],
    0.11224,
    tolerance = hitol
  )

  expect_snapshot_value(
    fit$summary.random$spde1D$mean[c(1, 27, 50)],
    tolerance = hitol,
    style = "deparse"
  )
  expect_snapshot_value(
    fit$summary.random$spde1D$sd[c(2, 32, 29)],
    tolerance = hitol,
    style = "deparse"
  )

  skip_if_not_installed("sn")

  # predicted intensity integral
  ips <- data.frame(x = 1:55, weight = 1)
  Lambda <- predict(
    fit,
    ips,
    ~ sum(weight * exp(spde1D + Intercept)),
    n.samples = 100,
    seed = 4354
  )
  expect_equal(Lambda$mean, 130.127, tolerance = hitol)
  expect_equal(Lambda$sd, 11.09, tolerance = 1)
})





test_that("1D LGCP fitting, compressed format", {
  skip_on_cran()
  local_bru_safe_inla()
  local_bru_options_set(
    control.inla = list(int.strategy = "eb"),
    control.compute = list(dic = FALSE, waic = FALSE)
  )

  data(Poisson2_1D, package = "inlabru", envir = environment())
  x <- seq(0, 55, length = 50)
  mesh1D <- INLA::inla.mesh.1d(x, boundary = "free")
  matern <- INLA::inla.spde2.pcmatern(
    mesh1D,
    prior.range = c(1, 0.01),
    prior.sigma = c(0.1, 0.75),
    constr = TRUE
  )

  mdl <- ~ spde1D(main = x, model = matern) + Intercept(1)

  fit1 <- bru(
    mdl,
    like(
      formula = x ~ .,
      family = "cp",
      data = pts2,
      domain = list(x = mesh1D),
      options = list(bru_compress_cp = FALSE)
    )
  )
  fit2 <- bru(
    mdl,
    like(
      formula = x ~ .,
      family = "cp",
      data = pts2,
      domain = list(x = mesh1D),
      options = list(bru_compress_cp = TRUE)
    )
  )

  expect_equal(
    fit1$summary.hyperpar,
    fit2$summary.hyperpar,
    tolerance = hitol
  )
  expect_equal(
    fit1$summary.fixed,
    fit2$summary.fixed,
    tolerance = midtol
  )
  expect_equal(
    fit1$summary.random,
    fit2$summary.random,
    tolerance = midtol
  )

  fit3 <- bru(
    mdl,
    like(
      formula = x ~ spde1D + Intercept,
      family = "cp",
      data = pts2,
      domain = list(x = mesh1D),
      options = list(bru_compress_cp = FALSE)
    )
  )
  fit4 <- bru(
    mdl,
    like(
      formula = x ~ spde1D + Intercept,
      family = "cp",
      data = pts2,
      domain = list(x = mesh1D),
      options = list(bru_compress_cp = TRUE)
    )
  )

  expect_equal(
    fit3$summary.hyperpar,
    fit4$summary.hyperpar,
    tolerance = hitol
  )
  expect_equal(
    fit3$summary.fixed,
    fit4$summary.fixed,
    tolerance = midtol
  )
  expect_equal(
    fit3$summary.random,
    fit4$summary.random,
    tolerance = midtol
  )
})
