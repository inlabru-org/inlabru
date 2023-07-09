test_that("bru_log", {
  local_bru_options_set(
    # Show messages up to and including level 2 (default 0)
    bru_verbose = 2,
    # Store messages to an including level 3 (default Inf, storing all)
    bru_verbose_store = 3
  )

  bru_log_bookmark("bookmark 1")
  bru_log_message("Test message 1", verbosity = 1)
  bru_log_message("Test message 2", verbosity = 2)
  bru_log_bookmark("bookmark 2")
  bru_log_message("Test message 3", verbosity = 3)
  bru_log_message("Test message 4", verbosity = 4)

  expect_equal(length(bru_log()["bookmark 1"]), 2L)
  expect_equal(length(bru_log()["bookmark 2"]), 1L)
})
