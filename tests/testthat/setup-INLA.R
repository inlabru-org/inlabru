options("rgdal_show_exportToProj4_warnings" = "none")
#library(inlabru)

if (requireNamespace("INLA", quietly = TRUE)) {
  # Save the num.threads option so it can be restored in teardown-INLA.R
  testthat_inla_num_threads <- INLA::inla.getOption("num.threads")
  # Tests should set num.threads = "1:1" to ensure within-system repeatability
  # by calling skip_if_not(bru_safe_inla())
}