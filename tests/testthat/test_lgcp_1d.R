context("1D LGCP fitting and prediction (test_lgcp_1d.R)")

test_data <- function() {
  data(Poisson2_1D, package = "inlabru")
  x <- seq(0, 55, length = 50)
  mesh1D <- INLA::inla.mesh.1d(x, boundary = "free")
  matern <- INLA::inla.spde2.pcmatern(mesh1D,
    prior.range = c(150, 0.75),
    prior.sigma = c(0.1, 0.75)
  )
  mdl <- x ~ spde1D(main = x, model = matern) + Intercept
  fit <- lgcp(mdl, pts2,
    ips = ipoints(c(0, 55), 50, name = "x"),
    options = list(
      control.inla = list(int.strategy = "eb"),
      num.threads = "1:1"
    )
  )
  list(
    mesh1D = mesh1D,
    fit = fit
  )
}


test_that("1D LGCP fitting", {
  skip_on_cran()
  skip_if_not(bru_safe_inla())

  result <- test_data()
  mesh1D <- result$mesh1D
  fit <- result$fit

  # Needed for reproducible predict
  set.seed(123L)

  expect_is(fit, "bru")

  expect_equal(fit$summary.fixed["Intercept", "mean"], 1.08959, tolerance = midtol)
  expect_equal(fit$summary.fixed["Intercept", "sd"], 0.4206289, tolerance = midtol)

  expect_equal(fit$summary.random$spde1D$mean[c(1, 27, 50)], c(-0.46315457, 0.09792757, -3.25164489), tolerance = midtol)
  expect_equal(fit$summary.random$spde1D$sd[c(2, 32, 29)], c(0.5887868, 0.4267676, 0.4288160), tolerance = midtol)

  pr <- predict(fit,
    data = data.frame(x = mesh1D$loc), formula = ~spde1D,
    n.samples = 100, seed = 84354
  )
  expect_equal(pr$mean, fit$summary.random$spde1D$mean, tolerance = hitol)
  expect_equal(pr$sd, fit$summary.random$spde1D$sd, tolerance = hitol)

  # predicted intensity integral
  ips <- ipoints(c(0, 55), 100, name = "x")
  Lambda <- predict(fit, ips, ~ sum(weight * exp(spde1D + Intercept)), n.samples = 100, seed = 4354)
  expect_equal(Lambda$mean, 131.5858, tolerance = hitol)
  expect_equal(Lambda$sd, 12.37687, tolerance = 1)
})





test_data_discrete <- function() {
  data(Poisson2_1D)
  xx <- ceiling(pts2$x)
  data <- data.frame(x = rep(xx, xx))
  x <- seq(0, 55, length = 56)
  mesh1D <- INLA::inla.mesh.1d(x, boundary = "free")
  matern <- INLA::inla.spde2.pcmatern(mesh1D,
    prior.range = c(150, 0.75),
    prior.sigma = c(0.1, 0.75)
  )
  mdl <- x ~ spde1D(main = x, model = matern) + Intercept
  fit <- lgcp(mdl,
    data = pts2,
    domain = list(x = mesh1D),
    options = list(
      control.inla = list(int.strategy = "eb"),
      num.threads = "1:1"
    )
  )
  fit2 <- lgcp(mdl,
    data = pts2,
    domain = list(x = (0:55)),
    options = list(
      control.inla = list(int.strategy = "eb"),
      num.threads = "1:1"
    )
  )
  list(
    mesh1D = mesh1D,
    fit = fit
  )
}


test_that("1D LGCP fitting", {
  skip_on_cran()
  skip_if_not(bru_safe_inla())

  result <- test_data_discrete()
  mesh1D <- result$mesh1D
  fit <- result$fit

  # Needed for reproducible predict
  set.seed(123L)

  expect_is(fit, "bru")

  expect_equal(fit$summary.fixed["Intercept", "mean"], 1.08959, tolerance = midtol)
  expect_equal(fit$summary.fixed["Intercept", "sd"], 0.4206289, tolerance = midtol)

  expect_equal(fit$summary.random$spde1D$mean[c(1, 27, 50)],
    c(-0.4619925, 0.2925785, -1.7602729),
    tolerance = midtol
  )
  expect_equal(fit$summary.random$spde1D$sd[c(2, 32, 29)],
    c(0.5905830, 0.4206042, 0.4219461),
    tolerance = midtol
  )

  pr <- predict(fit,
    data = data.frame(x = mesh1D$loc), formula = ~spde1D,
    n.samples = 100, seed = 84354
  )
  expect_equal(pr$mean, fit$summary.random$spde1D$mean, tolerance = hitol)
  expect_equal(pr$sd, fit$summary.random$spde1D$sd, tolerance = hitol)

  # predicted intensity integral
  ips <- ipoints(c(0, 55), 100, name = "x")
  Lambda <- predict(fit, ips, ~ sum(weight * exp(spde1D + Intercept)), n.samples = 100, seed = 4354)
  expect_equal(Lambda$mean, 131.5858, tolerance = hitol)
  expect_equal(Lambda$sd, 12.37687, tolerance = 1)
})
