test_that("multiplication works", {
  local_bru_options_set(
    # Show messages up to and including level 2 (default 0)
    bru_verbose = 2,
    # Store messages to an including level 3 (default Inf, storing all)
    bru_verbose_store = 3
  )

  bru_log_bookmark("test bookmark 1")
  expect_message(
    bru_log_message("Test message 1", verbosity = 1),
    regexp = ".*"
  )
  expect_message(
    bru_log_message("Test message 2", verbosity = 2),
    regexp = ".*"
  )
  bru_log_bookmark("test bookmark 2")
  expect_message(
    bru_log_message("Test message 3", verbosity = 3),
    regexp = NA
  )
  expect_message(
    bru_log_message("Test message 4", verbosity = 4),
    regexp = NA
  )
  expect_equal(length(bru_log_get(bookmark = "test bookmark 1")), 3)
  expect_equal(length(bru_log_get(bookmark = "test bookmark 2")), 1)
})
