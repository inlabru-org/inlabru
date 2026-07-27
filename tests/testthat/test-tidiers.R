# TODO: Replace the brittle mockups with actual bru/bru_info/bru_model/bru_obs/etc
# objects so that accessor methods can work properly when the data structure
# templates change.

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
    list(summary.fixed = fe, bru_info = list(model = list(lhoods = list()))),
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
      bru_info = structure(
        list(
          inlabru_version = as.character(utils::packageVersion("inlabru")),
          model = structure(
            list(
              lhoods = structure(
                list(
                  structure(
                    list(
                      response_data = data.frame(x = 1:15),
                      response = "x",
                      family = "normal"
                    ),
                    class = "bru_obs"
                  )
                ),
                class = "bru_obs_list"
              )
            ),
            class = "bru_model"
          )
        ),
        class = "bru_info"
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

test_that("glance.bru returns NA nobs for any point process fit", {
  # nobs is ill-defined for "cp" likelihoods.
  bru_ver <- as.character(utils::packageVersion("inlabru"))
  fit <- structure(
    list(
      bru_info = structure(
        list(
          inlabru_version = bru_ver,
          model = structure(
            list(
              lhoods = structure(
                list(
                  structure(
                    list(
                      response_data = data.frame(x = 1:10),
                      response = "x",
                      family = "cp",
                      inla.family = "xpoisson"
                    ),
                    class = "bru_obs"
                  )
                ),
                class = "bru_obs_list"
              )
            ),
            class = "bru_model"
          )
        ),
        class = "bru_info"
      )
    ),
    class = c("bru", "inla")
  )
  expect_true(is.na(glance(fit)$nobs))
})

test_that("glance.bru returns NA for missing fields gracefully", {
  fake_fit <- structure(
    list(
      bru_info = structure(
        list(
          inlabru_version = as.character(utils::packageVersion("inlabru")),
          model = structure(
            list(
              lhoods = list()
            )
          ),
          class = "bru_model"
        ),
        class = "bru_info"
      )
    ),
    class = c("bru", "inla")
  )
  result <- glance(fake_fit)
  expect_true(is.na(result$dic))
  expect_true(is.na(result$waic))
})
