context("2D LGCP fitting and prediction")
library(INLA)


data(gorillas, package = "inlabru")
matern <- inla.spde2.pcmatern(gorillas$mesh, prior.sigma = c(0.1, 0.01), prior.range = c(5, 0.01))
cmp <- coordinates ~ mySmooth(map = coordinates,model = matern) +Intercept
init.tutorial()
fit <- lgcp(cmp, gorillas$nests, samplers = gorillas$boundary)

test_that("2D LGCP fitting: Result object", {
  expect_is(fit, "bru")
})


test_that("2D LGCP fitting: INLA intercept", {

  expect_equal(fit$summary.fixed["Intercept","mean"], 1.121929, tolerance = midtol)
  expect_equal(fit$summary.fixed["Intercept","sd"], 0.5799173, tolerance = midtol)

})

test_that("2D LGCP fitting: INLA random field", {
  expect_equal(fit$summary.random$mySmooth$mean[c(1, 456, 789, 1058, 1479)],
               c(-2.0224597, 0.3874104, -0.4473572, 0.4019972, -1.7000660),
               tolerance = midtol)
  expect_equal(fit$summary.random$mySmooth$sd[c(1, 436, 759, 1158, 1279)],
               c(1.5924485,0.8243210,0.8209047,0.7928983,1.0671142),
               tolerance = midtol)
})

# test_that("2D LGCP fitting: predicted random field", {
#   
#   warning("This test needs to be improved by passing a seed to inla.posterior.sample()")
#   
#   loc = SpatialPoints(gorillas$mesh$loc[,c(1,2)])
#   proj4string(loc) = CRS(proj4string(gorillas$nests))
#   pr <- predict(fit, loc,  ~ mySmooth, n.samples = 500)
#   expect_equal(pr$mean[c(1,255,778,1000)], c(-2.057675,-1.766163,-1.512785,-1.488362), tolerance = hitol)
#   expect_equal(pr$sd[c(2,215,656,1010)], c(0.0000000,0.6629913,0.9822118,1.2876455), tolerance = hitol)
#   
# })

test_that("2D LGCP fitting: predicted intensity integral", {
  
  warning("This test needs to be improved by passing a seed to inla.posterior.sample()")
  ips = ipoints(gorillas$boundary, gorillas$mesh)
  Lambda <- predict(fit, ips, ~ sum(weight * exp(mySmooth + Intercept)), n.samples = 500) 
  
  expect_equal(Lambda$mean, 647.4751, tolerance = hitol)
  expect_equal(Lambda$sd, 25.54122, tolerance = hitol)
  
})

# test_that("Supplying integration points instead of samplers", {
#   ips = ipoints(gorillas$boundary, gorillas$mesh)
#   fit_ips <- lgcp(cmp, gorillas$nests, ips = ips)
#   
#   expect_equal(fit_ips$summary.fixed["Intercept","mean"], fit$summary.fixed["Intercept","mean"], tolerance = lowtol)
#   expect_equal(fit_ips$summary.fixed["Intercept","sd"], fit$summary.fixed["Intercept","sd"], tolerance = lowtol)
#   expect_equal(fit_ips$summary.random$mySmooth$mean, fit$summary.random$mySmooth$mean, tolerance = lowtol)
#   expect_equal(fit_ips$summary.random$mySmooth$sd, fit$summary.random$mySmooth$sd, tolerance = lowtol)
# })


test_that("Supplying domain definition", {
  fit_dom <- lgcp(cmp, gorillas$nests, samplers = gorillas$boundary, domain = list(coordinates = gorillas$mesh))
  
  expect_equal(fit_dom$summary.fixed["Intercept","mean"], fit$summary.fixed["Intercept","mean"], tolerance = midtol)
  expect_equal(fit_dom$summary.fixed["Intercept","sd"], fit$summary.fixed["Intercept","sd"], tolerance = midtol)
  expect_equal(fit_dom$summary.random$mySmooth$mean, fit$summary.random$mySmooth$mean, tolerance = midtol)
  expect_equal(fit_dom$summary.random$mySmooth$sd, fit$summary.random$mySmooth$sd, tolerance = midtol)
})

