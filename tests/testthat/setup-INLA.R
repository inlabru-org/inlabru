if (requireNamespace("INLA", quietly = TRUE)) {
  # Save the num.threads option so it can be restored in teardown-INLA.R
  inlabru_testthat_inla_num_threads <- INLA::inla.getOption("num.threads")
  # Tests should set num.threads = "1:1" to ensure within-system repeatability
  # by calling skip_if_not(bru_safe_inla())
  withr::defer(
    INLA::inla.setOption(num.threads = inlabru_testthat_inla_num_threads),
    teardown_env()
  )
}
