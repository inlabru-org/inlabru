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

  expect_snapshot_value(
    fit$summary.fixed["Intercept", "mean"],
    tolerance = midtol,
    style = "deparse"
  )
  expect_snapshot_value(
    fit$summary.fixed["Intercept", "sd"],
    tolerance = midtol,
    style = "deparse"
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

  pr <- predict(fit,
    data = data.frame(x = mesh1D$loc),
    formula = ~ list(
      Intercept = Intercept,
      spde1D = spde1D,
      both = Intercept + spde1D
    ),
    n.samples = 100, seed = 84354
  )

  # ggplot(data.frame(idx = seq_along(pr$both$mean),
  #                   fit = fit$summary.random$spde1D$mean+fit$summary.fixed$mean,
  #                   est = pr$both$mean,
  #                   std_err = pr$both$mean.mc_std_err)) +
  #   geom_point(aes(idx, fit)) +
  #   geom_line(aes(idx, est)) +
  #   geom_ribbon(aes(idx, ymin=est-2*std_err, ymax=est+2*std_err), alpha = 0.2)
  #
  # ggplot(data.frame(idx = seq_along(pr$spde1D$mean),
  #                   fit = fit$summary.random$spde1D$mean,
  #                   est = pr$spde1D$mean,
  #                   std_err = pr$spde1D$mean.mc_std_err)) +
  #   geom_point(aes(idx, fit)) +
  #   geom_line(aes(idx, est)) +
  #   geom_ribbon(aes(idx, ymin=est-2*std_err, ymax=est+2*std_err), alpha = 0.2)
  #
  # ggplot(data.frame(idx = seq_along(pr$both$sd),
  #                   fit = fit$summary.random$spde1D$sd+fit$summary.fixed$sd,
  #                   est = pr$both$sd,
  #                   std_err = pr$both$sd.mc_std_err)) +
  #   geom_point(aes(idx, fit)) +
  #   geom_line(aes(idx, est)) +
  #   geom_ribbon(aes(idx, ymin=est-2*std_err, ymax=est+2*std_err), alpha = 0.2)
  #
  # ggplot(data.frame(idx = seq_along(pr$both$mean),
  #                   fit = fit$summary.random$spde1D$mean + fit$summary.fixed$mean,
  #                   est = pr$both$mean,
  #                   std_err = pr$both$mean.mc_std_err)) +
  #   geom_point(aes(idx, abs(est - fit) / std_err))

  expect_true(
    all(abs(pr$both$mean -
      (fit$summary.random$spde1D$mean + fit$summary.fixed$mean)) /
      pr$both$mean.mc_std_err <= 3)
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
    tolerance = hitol,
    style = "deparse"
  )
  expect_snapshot_value(
    fit$summary.fixed["Intercept", "sd"],
    tolerance = hitol,
    style = "deparse"
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
    prior.sigma = c(0.1, 0.75)
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
    tolerance = hitol
  )
  expect_equal(
    fit1$summary.random,
    fit2$summary.random,
    tolerance = hitol
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
