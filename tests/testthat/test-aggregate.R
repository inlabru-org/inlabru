test_that("Aggregated Gaussian observations", {
  local_bru_safe_inla()

  obs <- data.frame(
    x = c(10, 20, 30),
    y = c(10, 20, 30),
    z = c(10, 20, 30)
  )
  pred <- data.frame(
    x = c(1, 2, 3, 4, 5, 6),
    y = c(1, 20, 3, 40, 5, 60),
    weights = c(1, 1, 1, 1, 1, 1),
    grp = c(1, 1, 2, 2, 2, 3)
  )

  # Aggregation by average:
  agg <- bru_mapper_aggregate(rescale = TRUE, n_block = nrow(obs))

  comp <- ~ Intercept(1) + x

  fit <- bru(
    comp,
    bru_obs(
      z ~ ibm_eval(
        agg,
        input = list(weights = weights, block = grp),
        state = Intercept + x
      ),
      family = "normal",
      response_data = obs,
      data = pred,
      control.family = list(
        hyper = list(
          prec = list(
            initial = 6,
            fixed = TRUE
          )
        )
      ),
      allow_combine = TRUE
    )
  )

  expect_equal(
    fit$summary.fixed$mean,
    c(3.033, 4.426),
    tolerance = midtol
  )

  # With basic sf storage:

  obs_sf <- sf::st_as_sf(obs, coords = c("x", "y"))
  pred_sf <- sf::st_as_sf(pred, coords = c("x", "y"))

  comp_sf <- ~ Intercept(1) + x(sf::st_coordinates(pred_sf)[, "X"])

  fit_sf <- bru(
    comp_sf,
    bru_obs(
      z ~ ibm_eval(
        agg,
        input = list(weights = weights, block = grp),
        state = Intercept + x
      ),
      family = "normal",
      response_data = obs_sf,
      data = pred_sf,
      control.family = list(
        hyper = list(prec = list(initial = 6, fixed = TRUE))
      ),
      allow_combine = TRUE
    )
  )

  expect_equal(
    fit_sf$summary.fixed$mean,
    c(3.033, 4.426),
    tolerance = midtol
  )
})



test_that("Aggregated Poisson observations", {
  local_bru_safe_inla()

  obs <- data.frame(y = c(10, 20, 30))
  pred <- data.frame(
    x = c(1, 2, 3, 4, 5, 6),
    weights = c(1, 1, 1, 1, 1, 1),
    grp = c(1, 1, 2, 2, 2, 3)
  )

  # Aggregation by summation on the intensity/expectation scale
  # (log-sum-exp since the predictor is log-intensity)
  agg <- bru_mapper_logsumexp(rescale = FALSE, n_block = nrow(obs))

  comp <- ~ Intercept(1) + x

  fit <- bru(
    comp,
    bru_obs(
      y ~ ibm_eval(
        agg,
        input = list(weights = weights, block = grp),
        state = Intercept + x
      ),
      family = "poisson",
      response_data = obs,
      data = pred,
      allow_combine = TRUE
    )
  )

  expect_equal(
    fit$summary.fixed$mean,
    c(0.337, 0.470),
    tolerance = midtol
  )

  # With E specification:

  obs <- data.frame(
    y = c(10, 20, 30),
    E = c(1, 2, 3)
  )

  fit <- bru(
    comp,
    bru_obs(
      y ~ ibm_eval(
        agg,
        input = list(weights = weights, block = grp),
        state = Intercept + x
      ),
      family = "poisson",
      E = E,
      response_data = obs,
      data = pred,
      allow_combine = TRUE
    )
  )

  expect_equal(
    fit$summary.fixed$mean,
    c(0.639, 0.237),
    tolerance = midtol
  )
})
