test_that("INLA loading", {
  skip_on_cran()
  # On non-CRAN systems, INLA should be installed when running the tests
  expect_error(local_bru_safe_inla(), NA, class = "skip")
})
