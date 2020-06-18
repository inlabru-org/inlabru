
# Do not run tests if INLA is not available (e.g. on CRAN)

if (require("INLA", quietly = TRUE)) {
  library(testthat)
  library(INLA)
  library(rgdal)

  test_check("inlabru")
}
