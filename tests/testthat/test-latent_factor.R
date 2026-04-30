test_that("bru: factor component", {
  skip_on_cran()
  local_bru_safe_inla()

  # Seed influences data as well as predict()!
  withr::local_seed(123L)

  # Factor models
  input.df <- data.frame(
    x1 = factor(rep(c("Alpha", "Beta", "Mu"), c(70, 60, 70))),
    x2 = factor(rep(
      c("Delta", "Gamma", "H1", "H2", "H3", "H4", "H5"),
      c(75, 75, 10, 10, 10, 10, 10)
    ))
  )
  input.df$y <-
    c(1, -2, 0)[as.numeric(input.df$x1)] +
    c(1, -2, 1, -2, 1, -2, 1)[as.numeric(input.df$x2)] +
    rnorm(200, mean = 0, sd = 0.1)

  # Fit the model
  fit <- bru(
    components = ~ Intercept(main = 1) +
      fac1(main = x1, model = "factor_contrast") +
      fac2(main = x2, model = "iid"),
    formula = y ~ Intercept + fac1 + fac2,
    family = "gaussian",
    data = input.df,
    options = list(
      verbose = FALSE,
      control.inla = list(
        int.strategy = "eb"
      )
    )
  )

  #  fit0$summary.hyperpar
  #  fit$summary.hyperpar
  #  rbind(
  #    cbind(ID = rownames(fit0$summary.fixed), fit0$summary.fixed),
  #    fit0$summary.random$x2
  #  )
  #  rbind(
  #    cbind(ID = rownames(fit$summary.fixed), fit$summary.fixed),
  #    fit$summary.random$fac1,
  #    fit$summary.random$fac2
  #  )

  # Check factor effect results
  expect_equal(
    fit$summary.random$fac1$mean,
    c(-3.072345, -1.092111),
    tolerance = midtol
  )
  expect_equal(
    fit$summary.random$fac1$sd,
    c(0.04367061, 0.05013299),
    tolerance = midtol
  )
  expect_equal(
    fit$summary.random$fac2$mean,
    c(
      1.213695, -1.723471,  1.310547, -1.671748,
      1.302147, -1.700772,  1.271110
    ),
    tolerance = midtol
  )
  expect_equal(
    fit$summary.random$fac2$sd,
    c(
      0.5225545, 0.5212273, 0.5216542, 0.5216542,
      0.5216542, 0.5216542, 0.5216542
    ),
    tolerance = hitol
  )

  # Check if prediction works
  skip_if_not_installed("sn")
  pr <- predict(
    fit,
    data.frame(
      x1 = factor(c("Alpha", "Beta", "Mu")),
      x2 = factor(c("Delta", "Delta", "Gamma"))
    ),
    ~ fac1 + fac2,
    n.samples = 500,
    seed = 1,
    num.threads = "1:1:1"
  )
  expect_equal(
    pr[, "mean"],
    c(1.223399, -1.849207, -2.807144),
    tolerance = midtol
  )
  expect_equal(
    pr[, "sd"],
    c(0.5064775, 0.5113702, 0.5070317),
    tolerance = hitol
  )
})


test_that("bru: indexed factor component", {
  skip_on_cran()
  local_bru_safe_inla()

  # Seed influences data as well as predict()!
  withr::local_seed(123)

  # Factor models
  input.df <- data.frame(
    x3 = rep(c(1, 2, 3, 4), times = 2)
  )
  input.df$y <-
    c(11, 12, 13, 14)[input.df$x3] +
    rnorm(8, mean = 0, sd = 0.1)

  # Fit the model
  fit <- bru(
    components = ~ Intercept(main = 1) +
      fac3(main = x3, model = "factor_contrast"),
    formula = y ~ Intercept + fac3,
    family = "gaussian",
    data = input.df,
    options = list(
      verbose = FALSE,
      control.inla = list(
        int.strategy = "eb",
        h = 0.005
      )
    )
  )


  # Check factor effect results
  expect_equal(
    fit$summary.random$fac3$mean,
    c(1.095837, 2.122574, 2.961865),
    tolerance = midtol
  )
  expect_equal(
    fit$summary.random$fac3$sd,
    c(0.07778122, 0.07778122, 0.07778122),
    tolerance = midtol
  )
})
