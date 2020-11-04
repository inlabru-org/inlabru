context("1D LGCP fitting and prediction - nonlinear (test_lgcp_1d_nonlinear.R)")

test_that("Mexdolphin: Hazard rate detection function", {
  skip_if_not(bru_safe_inla())
  data(mexdolphin, package = "inlabru")
  
  hr <- function(distance, lsig) {
    1 - exp(-(distance / (exp(lsig)))^-1)
  }

  cmp <- ~ lsig + Intercept
  formula <- distance ~ log(hr(distance, lsig)) + Intercept
  ips <- ipoints(inla.mesh.1d(seq(0, 8, by = 0.1)), name = "distance")

  fit <- lgcp(
    components = cmp,
    mexdolphin$points,
    ips = ips,
    formula = formula,
    options = list(control.inla = list(int.strategy = "eb"))
  )

  # plot(hr(ips$distance, fit$summary.fixed["lsig", "mean"]))
  expect_equal(fit$summary.fixed["lsig", "mean"], 1.038281, tolerance = lowtol)
  expect_equal(fit$summary.fixed["lsig", "sd"], 0.5183252, tolerance = lowtol)
  expect_equal(fit$summary.fixed["Intercept", "mean"], 2.325408, tolerance = lowtol)
  expect_equal(fit$summary.fixed["Intercept", "sd"], 0.2900139, tolerance = lowtol)
})
