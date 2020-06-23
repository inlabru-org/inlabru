
# Do not run tests if INLA is not available (e.g. on CRAN)

if (require("INLA", quietly = TRUE)) {
  library(testthat)
  library(INLA)
  library(PROJ6INLA200618)
  library(rgdal)

  test_check("inlabru")
}
