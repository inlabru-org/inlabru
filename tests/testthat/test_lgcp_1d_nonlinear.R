local_bru_testthat_setup()

test_that("Mexdolphin: Hazard rate detection function", {
  skip_on_cran()
  local_bru_safe_inla()
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

  #  ggplot(data.frame(distance = c(
  #    mexdolphin$points$distance,
  #    -mexdolphin$points$distance))) +
  #    geom_density(aes(distance, after_stat(count))) +
  #    geom_line(aes(distance, est),
  #              data = data.frame(distance = seq(-8,8, by = 0.01)) %>%
  #                mutate(est = hr(abs(distance),
  #                                fit$summary.fixed["lsig","mean"]) *
  #                         exp(fit$summary.fixed["Intercept","mean"])))
  #
  #  plot(ips$distance, hr(ips$distance, fit$summary.fixed["lsig", "mean"]))

  expect_equal(fit$summary.fixed["lsig", "mean"], 1.038281, tolerance = lowtol)
  expect_equal(fit$summary.fixed["lsig", "sd"], 0.5183252, tolerance = lowtol)
  expect_equal(fit$summary.fixed["Intercept", "mean"], 2.325408, tolerance = lowtol)
  expect_equal(fit$summary.fixed["Intercept", "sd"], 0.2900139, tolerance = lowtol)
})
