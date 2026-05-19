test_that("tidy.bru returns correct columns for fixed effects", {
  fe <- data.frame(
    mean = c(0.5, 1.2),
    sd = c(0.1, 0.2),
    `0.025quant` = c(0.3, 0.8),
    `0.975quant` = c(0.7, 1.6),
    check.names = FALSE,
    row.names = c("Intercept", "x1")
  )
  fake_fit <- structure(
    list(summary.fixed = fe, bru_info = list(lhoods = list())),
    class = c("bru", "inla")
  )

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
  fake_fit <- structure(list(summary.fixed = NULL), class = c("bru", "inla"))
  expect_error(tidy(fake_fit, effects = "latent"), "must be")
})

test_that("glance.bru returns one-row tibble with expected columns", {
  fake_fit <- structure(
    list(
      dic = list(dic = 123.4),
      waic = list(waic = 125.0),
      mlik = matrix(c(-61.0, 0), nrow = 1),
      bru_timings = data.frame(
        Elapsed = as.difftime(c(1.0, 2.0), units = "secs")
      ),
      bru_info = list(
        lhoods = list(
          list(family = "gaussian", data = data.frame(x = 1:15))
        )
      )
    ),
    class = c("bru", "inla")
  )

  result <- glance(fake_fit)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1L)
  expect_equal(result$dic, 123.4)
  expect_equal(result$waic, 125.0)
  expect_equal(result$nobs, 15L)
  expect_equal(result$elapsed, 3.0)
})

test_that("glance.bru returns NA nobs for any multi-likelihood fit", {
  # nobs is ill-defined for joint models regardless of whether families
  # match — possibly-shared observations and different supports across
  # likelihoods make the row sum meaningless.
  mixed <- structure(
    list(
      bru_info = list(
        lhoods = list(
          list(family = "gaussian", data = data.frame(x = 1:10)),
          list(family = "poisson", data = data.frame(x = 1:5))
        )
      )
    ),
    class = c("bru", "inla")
  )
  same <- structure(
    list(
      bru_info = list(
        lhoods = list(
          list(family = "gaussian", data = data.frame(x = 1:10)),
          list(family = "gaussian", data = data.frame(x = 1:5))
        )
      )
    ),
    class = c("bru", "inla")
  )
  expect_true(is.na(glance(mixed)$nobs))
  expect_true(is.na(glance(same)$nobs))
})

test_that("glance.bru returns NA for missing fields gracefully", {
  fake_fit <- structure(
    list(bru_info = list(lhoods = list())),
    class = c("bru", "inla")
  )
  result <- glance(fake_fit)
  expect_true(is.na(result$dic))
  expect_true(is.na(result$waic))
})
