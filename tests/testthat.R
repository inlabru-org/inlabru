# Do not run tests if INLA is not available (e.g. on CRAN)
if (bru_safe_inla()) {
  library(testthat)
  options("rgdal_show_exportToProj4_warnings" = "none")

  test_check("inlabru")
}
