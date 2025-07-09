test_that("bru: linear component", {
  skip_on_cran()
  local_bru_safe_inla()

  # Seed influences data as well as predict()!
  withr::local_seed(123)

  input.df <- data.frame(x = cos(1:100))
  input.df <- within(input.df, y <- 5 + 2 * x + rnorm(100, mean = 0, sd = 0.1))

  fit <- bru(
    y ~ myLin(main = x, model = "linear") + Intercept(1),
    family = "gaussian",
    data = input.df
  )

  expect_equal(fit$summary.fixed["myLin", "mean"], 2.002273, tolerance = midtol)
  expect_equal(fit$summary.fixed["myLin", "sd"], 0.01323361, tolerance = hitol)

  skip_if_not_installed("sn")
  pr <- predict(
    fit,
    data.frame(x = c(1, 2)),
    ~ myLin + 2,
    n.samples = 5,
    seed = 1L
  )

  expect_equal(pr[, "mean"], c(4.005013, 6.010026), tolerance = midtol)
})


test_that("bru: linear predictor detection", {
  skip_on_cran()
  local_bru_safe_inla()

  # Seed influences data as well as predict()!
  withr::local_seed(123)

  input.df <- data.frame(x = cos(1:100), z = sin(1:100))
  input.df <- within(input.df, y <- 5 + 2 * x + rnorm(100, mean = 0, sd = 0.1))

  fit <- bru(
    ~ x(x) + z(z) + Intercept(1),
    formula = y ~ .,
    family = "gaussian",
    data = input.df,
    options = list(bru_run = FALSE)
  )
  expect_equal(
    as.character(bru_obs_expr(
      fit$bru_info$lhoods[[1]],
      fit$bru_info$effects
    )),
    "x + z + Intercept",
    info = "Should construct correct expression"
  )
  expect_true(
    fit$bru_info$lhoods[[1]][["is_additive"]],
    info = "Should detect additive predictor"
  )
  expect_true(
    fit$bru_info$lhoods[[1]][["linear"]],
    info = "Should detect linear predictor"
  )
  expect_equal(
    bru_used(fit),
    bru_used(effect = c("x", "z", "Intercept"), latent = character(0))
  )

  fit <- bru(
    ~ x(x) + z(z) + Intercept(1),
    formula = y ~ x + Intercept,
    family = "gaussian",
    data = input.df,
    options = list(bru_run = FALSE)
  )
  expect_equal(
    as.character(bru_obs_expr(
      fit$bru_info$lhoods[[1]],
      fit$bru_info$effects
    )),
    "x + Intercept",
    info = "Should construct correct expression"
  )
  expect_true(
    fit$bru_info$lhoods[[1]][["is_additive"]],
    info = "Should detect additive predictor"
  )
  expect_true(
    fit$bru_info$lhoods[[1]][["linear"]],
    info = "Should detect linear predictor"
  )
  expect_equal(
    bru_used(fit),
    bru_used(effect = c("x", "Intercept"), latent = character(0))
  )

  # Predictor with use of _latent should be detected as non-additive
  fit <- bru(
    ~ x(x) + z(z) + Intercept(1),
    formula = y ~ x + z_latent + Intercept,
    family = "gaussian",
    data = input.df,
    options = list(bru_run = FALSE)
  )
  expect_equal(
    as.character(bru_obs_expr(
      fit$bru_info$lhoods[[1]],
      fit$bru_info$effects
    )),
    "x + Intercept + z_latent",
    info = "Should construct correct expression"
  )
  expect_false(
    fit$bru_info$lhoods[[1]][["is_additive"]],
    info = "Should detect non-additive predictor"
  )
  expect_false(
    fit$bru_info$lhoods[[1]][["linear"]],
    info = "Should detect non-linear predictor, due to non-pure additivity"
  )
  expect_equal(
    bru_used(fit),
    bru_used(effect = c("x", "Intercept"), latent = "z"),
    info = "Should detect two components, 'x' and 'Intercept', and latent 'z'"
  )

  # Predictor with missing component/variable should be detected as
  # non-additive, but only give errors when attempting to evaluate the model
  fit <- bru(
    ~ x(x) + z(z) + Intercept(1),
    formula = y ~ x + Intercept + something,
    family = "gaussian",
    data = input.df,
    options = list(bru_run = FALSE)
  )
  expect_equal(
    as.character(bru_obs_expr(
      fit$bru_info$lhoods[[1]],
      fit$bru_info$effects
    )),
    "x + Intercept + something",
    info = "Should construct correct expression"
  )
  expect_false(
    fit$bru_info$lhoods[[1]][["is_additive"]],
    info = "Should detect non-additive predictor"
  )
  expect_false(
    fit$bru_info$lhoods[[1]][["linear"]],
    info = "Should detect non-linear predictor, due to missing variable"
  )
  expect_equal(
    bru_used(fit),
    bru_used(effect = c("x", "Intercept"), latent = character(0)),
    info = "Should detect two components, 'x' and 'Intercept'"
  )
})
