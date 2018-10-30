
# Do not run tests if INLA is not available (e.g. on CRAN)

if (require("INLA", quietly = TRUE)) {

  library(testthat)
  library(INLA)
  
  lowtol = 1E-5
  midtol = 1E-3
  hitol = 1E-1
  
  
  test_check("inlabru")
}