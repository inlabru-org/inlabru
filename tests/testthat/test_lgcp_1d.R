context("1D LGCP fitting and prediction (test_lgcp_1d.R)")
library(INLA)

test_data <- function() {
  data(Poisson2_1D)
  x <- seq(0, 55, length = 50)
  mesh1D <- inla.mesh.1d(x, boundary = "free")
  matern <- inla.spde2.pcmatern(mesh1D, prior.range = c(150, 0.75), prior.sigma = c(0.1, 0.75))
  mdl <- x ~ spde1D(main = x, model = matern) + Intercept
  fit <- lgcp(mdl, pts2,
              ips = ipoints(c(0, 55), 50, name = "x"),
              options = list(control.inla = list(int.strategy = "eb"))
  )
  list(mesh1D = mesh1D,
       fit = fit)
}


test_that("1D LGCP fitting", {
  result <- test_data()
  mesh1D <- result$mesh1D
  fit <- result$fit
  
  expect_is(fit, "bru")

  expect_equal(fit$summary.fixed["Intercept", "mean"], 1.08959, tolerance = midtol)
  expect_equal(fit$summary.fixed["Intercept", "sd"], 0.4206289, tolerance = midtol)

  expect_equal(fit$summary.random$spde1D$mean[c(1, 27, 50)], c(-0.46315457, 0.09792757, -3.25164489), tolerance = midtol)
  expect_equal(fit$summary.random$spde1D$sd[c(2, 32, 29)], c(0.5887868, 0.4267676, 0.4288160), tolerance = midtol)

  pr <- predict(fit, data = data.frame(x = mesh1D$loc), formula = ~spde1D,
                n.samples = 100, seed = 84354)
  expect_equal(pr$mean, fit$summary.random$spde1D$mean, tolerance = hitol)
  expect_equal(pr$sd, fit$summary.random$spde1D$sd, tolerance = hitol)

  # predicted intensity integral
  ips <- ipoints(c(0, 55), 100, name = "x")
  Lambda <- predict(fit, ips, ~ sum(weight * exp(spde1D + Intercept)), n.samples = 100, seed = 4354)
  expect_equal(Lambda$mean, 131.5858, tolerance = hitol)
  expect_equal(Lambda$sd, 12.37687, tolerance = 1)
})
