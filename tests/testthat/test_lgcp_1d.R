local_bru_testthat_setup()

test_data <- function() {
  data(Poisson2_1D, package = "inlabru", envir = environment())
  x <- seq(0, 55, length = 50)
  mesh1D <- INLA::inla.mesh.1d(x, boundary = "free")
  matern <- INLA::inla.spde2.pcmatern(
    mesh1D,
    prior.range = c(150, 0.75),
    prior.sigma = c(0.1, 0.75)
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
    1.054199,
    tolerance = hitol
  )
  expect_equal(
    fit$summary.fixed["Intercept", "sd"],
    0.4200817,
    tolerance = hitol
  )

  expect_equal(
    fit$summary.random$spde1D$mean[c(1, 27, 50)],
    c(-0.4146421, 0.1373841, -3.2198098),
    tolerance = hitol
  )
  expect_equal(
    fit$summary.random$spde1D$sd[c(2, 32, 29)],
    c(0.5887868, 0.4267676, 0.4288160),
    tolerance = hitol
  )

  pr <- predict(fit,
    data = data.frame(x = mesh1D$loc), formula = ~spde1D,
    n.samples = 100, seed = 84354
  )
  expect_equal(
    pr$mean,
    fit$summary.random$spde1D$mean,
    tolerance = hitol
  )
  expect_equal(pr$sd, fit$summary.random$spde1D$sd, tolerance = hitol)

  # predicted intensity integral
  ips <- ipoints(c(0, 55), 100, name = "x")
  Lambda <- predict(
    fit,
    ips,
    ~ sum(weight * exp(spde1D + Intercept)),
    n.samples = 100,
    seed = 4354
  )
  expect_equal(Lambda$mean, 131.0741, tolerance = hitol)
  expect_equal(Lambda$sd, 11.35888, tolerance = 1)
})





test_data_discrete <- function() {
  data(Poisson2_1D, package = "inlabru", envir = environment())
  xx <- ceiling(pts2$x)
  data <- data.frame(x = rep(xx, xx))
  x <- seq(1, 55, length = 55)
  mesh1D <- INLA::inla.mesh.1d(x, boundary = "free")
  matern <- INLA::inla.spde2.pcmatern(mesh1D,
    prior.range = c(0.01, 0.01),
    prior.sigma = c(1, 0.01)
  )
  mdl <- x ~ spde1D(main = x, model = matern) + Intercept(1)
  fit <- lgcp(mdl,
    data = pts2,
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


test_that("1D LGCP fitting", {
  skip_on_cran()
  local_bru_safe_inla()

  result <- test_data_discrete()
  mesh1D <- result$mesh1D
  fit <- result$fit

  # Needed for reproducible predict
  set.seed(123L)

  expect_s3_class(fit, "bru")

  expect_snapshot_value(
    fit$summary.fixed["Intercept", "mean"],
    tolerance = midtol,
    style = "serialize"
  )
  expect_snapshot_value(
    fit$summary.fixed["Intercept", "sd"],
    tolerance = midtol,
    style = "serialize"
  )

  expect_snapshot_value(
    fit$summary.random$spde1D$mean[c(1, 27, 50)],
    tolerance = midtol,
    style = "serialize"
  )
  expect_snapshot_value(
    fit$summary.random$spde1D$sd[c(2, 32, 29)],
    tolerance = midtol,
    style = "serialize"
  )

  pr <- predict(
    fit,
    data = data.frame(x = mesh1D$loc),
    formula = ~ spde1D + Intercept,
    n.samples = 100,
    seed = 84354
  )
  expect_equal(
    pr$mean,
    fit$summary.random$spde1D$mean + fit$summary.fixed$mean,
    tolerance = hitol
  )

  # predicted intensity integral
  ips <- ipoints(c(0, 55), 100, name = "x")
  Lambda <- predict(
    fit,
    ips,
    ~ sum(weight * exp(spde1D + Intercept)),
    n.samples = 100,
    seed = 4354
  )
  expect_equal(Lambda$mean, 131.5858, tolerance = hitol)
  expect_equal(Lambda$sd, 12.37687, tolerance = 1)
})
