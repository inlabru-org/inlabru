test_that("sp data detected", {
  skip_on_cran()
  local_bru_safe_inla()
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  local_bru_options_set(bru_run = FALSE)
  withr::local_options(lifecycle_verbosity = "error")
  expect_error(
    {
      fit_sp <- bru(
        as.integer(date) ~ field(sp::coordinates, model = "linear"),
        data = gorillas_sp()$nests
      )
    },
    paste0(
      "`data` argument has deprecated support for `Spatial` input; ",
      "was deprecated in inlabru 2.12.0.9023."
    )
  )
})
