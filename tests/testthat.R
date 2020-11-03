
# Do not run tests if INLA is not available (e.g. on CRAN)

if (require("INLA", quietly = TRUE)) {
  library(testthat)
  library(INLA)
  options("rgdal_show_exportToProj4_warnings" = "none")
  library(rgdal)

  test_check("inlabru")
}
