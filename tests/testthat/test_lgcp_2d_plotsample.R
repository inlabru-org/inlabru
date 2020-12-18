local_bru_testthat_setup()

test_that("2D LGCP fitting and prediction: Plot sampling", {
  skip_on_cran()
  local_bru_safe_inla()

  options <- list(
    bru_verbose = TRUE,
    verbose = TRUE,
    control.inla = list(
      int.strategy = "eb"
    )
  )
  data(gorillas, package = "inlabru")

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
    7.095665,
    tolerance = midtol
  )
  expect_equal(
    fit$summary.fixed["Intercept", "mean"],
    1.796279,
    tolerance = midtol
  )
  expect_equal(
    fit$summary.fixed["Intercept", "sd"],
    0.5221979,
    tolerance = midtol
  )
  expect_equal(
    fit$summary.random$my.spde$mean[c(1, 100, 300)],
    c(-1.566168, 1.177564, -1.584715),
    tolerance = midtol
  )
  expect_equal(
    fit$summary.random$my.spde$sd[c(1, 100, 300)],
    c(1.3523279, 0.6516389, 1.2547961),
    tolerance = midtol
  )
})
