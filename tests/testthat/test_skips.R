test_that("INLA loading", {
  skip_on_cran()
  # Note:
  # On non-CRAN systems, INLA should be installed when running the tests
  # but all inla-requiring tests should be protected specifically for missing
  # INLA.  The check below should NOT be done, as it will fail when testing
  # without Suggests packages, which includes INLA.
  # expect_error(local_bru_safe_inla(), NA, class = "skip")
})
