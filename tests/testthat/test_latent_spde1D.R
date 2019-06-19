context("Latent models - 1D SPDE (test_latent_spde1D.R)")

latent_spde1D_testdata <- function() {
  data(Poisson2_1D)
  x <- seq(0, 55, length = 50)
  mesh1D <- inla.mesh.1d(x, boundary = "free")
  
  matern <- inla.spde2.pcmatern(mesh1D,
                                  prior.range=c(1, 0.01),
                                  prior.sigma=c(10, 0.01))
  
  cmp <- count ~ field(map = x, model = matern) + Intercept
  fit <- bru(cmp, data = countdata2, family = "poisson",
             options = list(E = countdata2$exposure))
  
  list(data = countdata2,
       matern = matern,
       cmp = cmp,
       fit = fit)
}

test_that("Latent models: SPDE 1D", {
  skip_on_cran()
  library(INLA)
  data <- latent_spde1D_testdata()
  
  # Check Intercept
  expect_equal(data$fit$summary.fixed["Intercept", "mean"], 5.982867, midtol)
  
  # Check SPDE
  expect_equal(data$fit$summary.random$field$mean[c(1,25,50)],
               c(-5.195488, -4.480026, -7.150024), midtol)
  
})

