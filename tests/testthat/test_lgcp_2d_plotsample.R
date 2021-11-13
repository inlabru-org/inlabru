local_bru_testthat_setup()

test_that("2D LGCP fitting and prediction: Plot sampling", {
  skip_on_cran()
  local_bru_safe_inla()

  options <- list(
    control.inla = list(
      int.strategy = "eb"
    )
  )
  data(gorillas, package = "inlabru", envir = environment())

  matern <- INLA::inla.spde2.pcmatern(gorillas$mesh,
    prior.sigma = c(0.1, 0.01),
    prior.range = c(5, 0.01)
  )

  cmp <- coordinates ~ my.spde(main = coordinates, model = matern)
  fit <- lgcp(cmp,
    data = gorillas$plotsample$nests,
    samplers = gorillas$plotsample$plots,
    domain = list(coordinates = gorillas$mesh),
    options = options
  )

  expect_equal(
    sum(fit$bru_info$lhoods[[1]]$E),
    7.096605,
    tolerance = lowtol
  )
  # Test points selected with sd+sd less than 1.2, for a more stable check.
  expect_snapshot_value(
    fit$summary.fixed["Intercept", "mean"] + fit$summary.random$my.spde$mean[c(19, 100, 212)],
    tolerance = midtol,
    style = "serialize"
  )
  expect_snapshot_value(
    fit$summary.fixed["Intercept", "sd"] + fit$summary.random$my.spde$sd[c(19, 100, 212)],
    tolerance = hitol,
    style = "serialize"
  )
})
