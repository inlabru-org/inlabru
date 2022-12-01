local_bru_testthat_setup()

test_that("1D integration points can be generated", {
  local_bru_safe_inla()

  obs <- data.frame(y = c(10, 20))
  pred <- data.frame(
    x = c(1, 2, 3, 4, 5),
    weights = c(1, 1, 1, 1, 1),
    grp = c(1, 1, 2, 2, 2)
  )

  agg <- bru_mapper_logsumexp(rescale = FALSE, n_block = nrow(obs))

  comp <- ~ Intercept(1) + x

  fit <- bru(
    comp,
    like(
      y ~ ibm_eval(
        agg,
        input = list(weights = weights, block = grp),
        state = Intercept + x
      ),
      family = "poisson",
      response_data = obs,
      data = pred
    )
  )

  expect_equal(
    fit$summary.fixed$mean,
    c(1.371, 0.124),
    tolerance = midtol
  )
})
