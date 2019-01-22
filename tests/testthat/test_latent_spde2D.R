context("Latent models - 2D SPDE - Group parameter")
library(INLA)

latent_spde2D_group_testdata <- function() {
  
  # Load and reduce data set
  data(mrsea)
  mrsea$points = mrsea$points[mrsea$points$season==1 | mrsea$points$season==2,]
  mrsea$samplers = mrsea$samplers[mrsea$samplers$season==1 | mrsea$samplers$season==2,]
  
  # Integration points
  ips <- ipoints(mrsea$samplers, group = "season")
  
  # Run the model
  matern <- inla.spde2.pcmatern(mrsea$mesh, 
                                prior.sigma = c(0.1, 0.01), 
                                prior.range = c(10000, 0.01))
  
  warning('Using workaround for known bug (fixed in new backend)')
  season = mrsea$points$season 
  
  cmp <- coordinates + season ~ mySmooth(map = coordinates, model = matern, group = season, ngroup = 2) + Intercept
  fit <- lgcp(cmp, mrsea$points, ips = ips, options = list(control.inla = list(int.strategy = "eb")))
  
  list(mrsea = mrsea,
       matern = matern,
       cmp = cmp,
       fit = fit)
}

data <- latent_spde2D_group_testdata()

test_that("Latent models: SPDE with group parameter (spatiotemporal)", {
  
  # Check Intercept
  expect_equal(data$fit$summary.fixed["Intercept", "mean"], -9.231591, midtol)
  
  # Check SPDE
  expect_equal(data$fit$summary.random$mySmooth$mean[c(1,250,550)], c(-1.2750692, -1.9794945, 0.7566357), midtol)
  
})

