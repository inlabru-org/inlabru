test_that("tidy.bru returns correct columns for fixed effects", {
  local_bru_safe_inla()

  fit <- bru(
    components = y ~ Intercept(1) + x1(main = x),
    family = "gaussian",
    data = data.frame(x = 1:10, y = rnorm(10)),
    options = list(bru_run = FALSE)
  )
  fe <- data.frame(
    mean = c(0.5, 1.2),
    sd = c(0.1, 0.2),
    `0.025quant` = c(0.3, 0.8),
    `0.975quant` = c(0.7, 1.6),
    check.names = FALSE,
    row.names = c("Intercept", "x1")
  )
  fake_fit <- fit
  fake_fit$summary.fixed <- fe

  result <- tidy(fake_fit)
  expect_s3_class(result, "tbl_df")
  expect_named(
    result,
    c("term", "estimate", "std.error", "conf.low", "conf.high")
  )
  expect_equal(result$term, c("Intercept", "x1"))
  expect_equal(result$estimate, c(0.5, 1.2))
})

test_that("tidy.bru errors on unknown effects argument", {
  local_bru_safe_inla()

  fit <- bru(
    components = y ~ Intercept(1) + x1(main = x),
    family = "gaussian",
    data = data.frame(x = 1:10, y = rnorm(10)),
    options = list(bru_run = FALSE)
  )
  fake_fit <- fit

  expect_error(tidy(fake_fit, effects = "latent"), "must be")
})

test_that("glance.bru returns one-row tibble with expected columns", {
  local_bru_safe_inla()

  fit <- bru(
    components = y ~ Intercept(1) + x1(main = x),
    family = "gaussian",
    data = data.frame(x = 1:10, y = rnorm(10)),
    options = list(bru_run = FALSE)
  )
  fake_fit <- fit
  fake_fit$dic <- list(dic = 123.4)
  fake_fit$waic <- list(waic = 125.0)
  fake_fit$mlik <- matrix(c(-61.0, 0), nrow = 1)
  fake_fit$bru_timings <- data.frame(
    Elapsed = as.difftime(c(1.0, 2.0), units = "secs")
  )

  result <- glance(fake_fit)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1L)
  expect_equal(result$dic, 123.4)
  expect_equal(result$waic, 125.0)
  expect_equal(result$nobs, 10L)
  expect_equal(result$elapsed, 3.0)
})

test_that("glance.bru returns NA nobs for any point process fit", {
  local_bru_safe_inla()

  # nobs is ill-defined for "cp" likelihoods.
  fit <- bru(
    components = y ~ Intercept(1) + x1(main = y),
    family = "cp",
    data = data.frame(y = runif(10)),
    domain = list(y = fmesher::fm_mesh_1d(0:1)),
    options = list(bru_run = FALSE)
  )
  expect_true(is.na(glance(fit)$nobs))
})

test_that("glance.bru returns NA for missing fields gracefully", {
  local_bru_safe_inla()

  fit <- bru(
    components = y ~ Intercept(1) + x1(main = x),
    family = "gaussian",
    data = data.frame(x = 1:10, y = rnorm(10)),
    options = list(bru_run = FALSE)
  )
  fake_fit <- fit

  result <- glance(fake_fit)
  expect_true(is.na(result$dic))
  expect_true(is.na(result$waic))
})
