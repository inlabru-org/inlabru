test_that("Aggregated Gaussian observations, using fm_block_eval", {
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

  comp <- ~ Intercept(1) + x

  fit <- bru(
    comp,
    bru_obs(
      z ~ fm_block_eval(
        block = grp,
        n_block = nrow(obs),
        weights = weights,
        rescale = TRUE,
        values = Intercept + x
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
      z ~ fm_block_eval(
        block = grp,
        weights = weights,
        rescale = TRUE,
        n_block = nrow(obs),
        values = Intercept + x
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



test_that("Aggregated Gaussian observations, using aggregate feature", {
  local_bru_safe_inla()

  obs <- data.frame(
    x = c(10, 20, 30),
    y = c(1000, 2000, 3000),
    z = c(10, 20, 30)
  )
  pred <- data.frame(
    x = c(1, 2, 3, 4, 5, 6),
    y = c(1, 20, 3, 40, 5, 60),
    weights = c(1, 1, 1, 1, 1, 1),
    grp = c(1, 1, 2, 2, 2, 3)
  )
  domain <- list()

  comp <- ~ Intercept(1) + x + y

  fit <- bru(
    comp,
    bru_obs(
      z ~ Intercept + x + y,
      family = "normal",
      response_data = obs,
      data = pred,
      aggregate = "average",
      aggregate_input = list(
        weights = weights,
        block = grp,
        n_block = bru_response_size(.response_data.)
      ),
      control.family = list(
        hyper = list(
          prec = list(
            initial = 6,
            fixed = TRUE
          )
        )
      )
    )
  )

  expect_equal(
    fit$summary.fixed$mean,
    c(3.636, 3.889, 0.0505),
    tolerance = midtol
  )

  expect_no_error({
    bru_obs(
      z ~ Intercept + x + y,
      family = "normal",
      response_data = obs,
      data = pred,
      aggregate = "average",
      aggregate_input = list(
        weights = weights,
        block = grp
      ),
      control.family = list(
        hyper = list(
          prec = list(
            initial = 6,
            fixed = TRUE
          )
        )
      )
    )
  })

  expect_error(
    {
      bru_obs(
        z ~ Intercept + x + y,
        family = "normal",
        response_data = obs,
        data = pred,
        aggregate = "average",
        aggregate_input = list(
          weights = weights,
          n_block = bru_response_size(.response_data.)
        ),
        control.family = list(
          hyper = list(
            prec = list(
              initial = 6,
              fixed = TRUE
            )
          )
        )
      )
    },
    paste0(
      "Aggregation requested, but `aggregate_input[['block']]` ",
      "evaluates to NULL."
    ),
    fixed = TRUE
  )
})


test_that("Aggregated Gaussian observations, using domain/samplers feature", {
  local_bru_safe_inla()

  obs <- data.frame(
    x = c(10, 20, 30),
    z = c(10, 20, 30)
  )
  pred <- data.frame(
    y = c((1 + 20) / 2, (3 + 40 + 5) / 3, 60)
  )
  domain <- list(x = 1:6)
  samplers <- list(x = list(1:2, 3:5, 6))

  comp <- ~ Intercept(1) + x + y

  fit <- bru(
    comp,
    bru_obs(
      z ~ Intercept + x + y,
      family = "normal",
      response_data = obs,
      data = pred,
      aggregate = "average",
      domain = domain,
      samplers = samplers,
      control.family = list(
        hyper = list(
          prec = list(
            initial = 6,
            fixed = TRUE
          )
        )
      )
    ),
    options = list(control.inla = list(
      int.strategy = "eb"
    ))
  )

  expect_equal(
    fit$summary.fixed$mean,
    c(3.636, 3.889, 0.0505),
    tolerance = midtol
  )

  expect_error(
    {
      bru_model(
        comp,
        bru_obs(
          z ~ Intercept + x + y,
          family = "normal",
          response_data = obs,
          data = NULL,
          aggregate = "average",
          domain = domain,
          samplers = samplers,
          control.family = list(
            hyper = list(
              prec = list(
                initial = 6,
                fixed = TRUE
              )
            )
          )
        )
      )
    },
    paste0(
      "The input evaluation 'y' for 'y' failed. ",
      "Perhaps the data object doesn't contain the needed variables?"
    ),
    fixed = TRUE
  )

  skip_if(utils::packageVersion("fmesher") >= "0.4.0.9006")
  # For fmesher < 0.4.0.9006, detect character .block info
  domain <- list(x = 1:6, y = 2:4)
  samplers <- list(x = list(1:2, 3:5, 6))
  expect_error(
    {
      bru_obs(
        z ~ Intercept + x + y,
        family = "normal",
        response_data = obs,
        data = NULL,
        aggregate = "average",
        domain = domain,
        samplers = samplers,
        control.family = list(
          hyper = list(
            prec = list(
              initial = 6,
              fixed = TRUE
            )
          )
        )
      )
    },
    "'character' aggregation block information detected.",
    fixed = TRUE
  )
})



test_that("Aggregated Poisson observations, using mapper", {
  local_bru_safe_inla()

  obs <- data.frame(y = c(10, 20, 30))
  pred <- data.frame(
    x = c(1, 2, 3, 4, 5, 6),
    weights = c(1, 1, 1, 1, 1, 1),
    grp = c(1, 1, 2, 2, 2, 3)
  )

  # Aggregation by summation on the intensity/expectation scale
  # (log-sum-exp since the predictor is log-intensity)
  agg <- bm_logsumexp(rescale = FALSE, n_block = nrow(obs))

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

test_that("Aggregated Gaussian observations, using mapper", {
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
  agg <- bm_aggregate(rescale = TRUE, n_block = nrow(obs))

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
