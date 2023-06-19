latent_spde1D_testdata <- function() {
  local_bru_safe_inla()
  data(Poisson2_1D, package = "inlabru", envir = environment())
  x <- seq(0, 55, length = 50)
  mesh1D <- INLA::inla.mesh.1d(x, boundary = "free")

  matern <- INLA::inla.spde2.pcmatern(mesh1D,
    prior.range = c(1, 0.01),
    prior.sigma = c(1, 0.01)
  )

  cmp <- count ~ field(main = x, model = matern) + Intercept(1)
  # This model is sensitive to the integration strategy; "eb" is too smooth.
  fit <- bru(cmp,
    data = countdata2, family = "poisson",
    options = list(
      E = countdata2$exposure,
      control.inla = list(h = 0.005)
    )
  )

  list(
    data = countdata2,
    matern = matern,
    cmp = cmp,
    fit = fit
  )
}

test_that("Latent models: SPDE 1D", {
  skip_on_cran()
  local_bru_safe_inla()
  data <- latent_spde1D_testdata()

  # Check Intercept + SPDE (highly negatively correlated)
  expect_snapshot_value(
    data$fit$summary.fixed["Intercept", "mean"] +
      data$fit$summary.random$field$mean[c(1, 25, 50)],
    tolerance = hitol,
    style = "serialize"
  )
})
