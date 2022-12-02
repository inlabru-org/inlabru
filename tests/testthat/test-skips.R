# Note:
# On non-CRAN systems, INLA should be installed when running the tests
# but all inla-requiring tests should be protected specifically for missing
# INLA.  The check below should NOT be done, as it will fail when testing
# without Suggests packages, which includes INLA.
# The test is now disabled since it's empty
#
# Note 2: This test was likely dubious, as testthat might run it in the same
#         thread as some previous tests.

# test_that("INLA loading", {
#   skip_on_cran()
#   # expect_error(local_bru_safe_inla(), NA, class = "skip")
# })
