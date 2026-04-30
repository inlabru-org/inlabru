local_bru_testthat_setup(envir = testthat::teardown_env())
bru_options_set_local(
  "bru_compat_pre_2_14_enable" = FALSE,
  .envir = testthat::teardown_env()
)
