test_that("Missing data infilling", {
  skip_on_cran()
  local_bru_safe_inla()
  skip_if_not(bru_safe_sp())
  skip_if_not_installed("terra")

  points <-
    sp::SpatialPointsDataFrame(
      matrix(1:6, 3, 2),
      data = data.frame(val = c(NA, NA, NA))
    )
  input_coord <- expand.grid(x = 0:6, y = 0:7)
  input_data <- data.frame(val = as.vector(input_coord$y))
  input <-
    sp::SpatialPixelsDataFrame(
      input_coord,
      data = input_data
    )

  val0 <- bru_fill_missing(input, points, points$val)

  input <-
    sp::SpatialPointsDataFrame(
      input_coord,
      data = input_data
    )
  val2 <- bru_fill_missing(input, points, points$val)

  expect_equal(val2, val0)

  input <- sp::SpatialGrid(sp::GridTopology(c(0, 0), c(1, 1), c(7, 8)))
  input_data <- data.frame(val = as.vector(sp::coordinates(input)[, 2]))
  input <-
    sp::SpatialGridDataFrame(
      input,
      data = input_data
    )
  val2 <- bru_fill_missing(input, points, points$val)

  expect_equal(val2, val0)
})



test_that("Laplace distribution", {
  q <- -5:5
  rate <- 2
  p <- plaplace(q, rate = rate)

  expect_equal(plaplace(q, rate = rate, lower.tail = TRUE), p)
  expect_equal(plaplace(-q, rate = rate, lower.tail = FALSE), p)
  expect_equal(
    plaplace(q, rate = rate, lower.tail = TRUE, log.p = TRUE),
    log(p)
  )
  expect_equal(
    plaplace(-q, rate = rate, lower.tail = FALSE, log.p = TRUE),
    log(p)
  )
  expect_equal(qlaplace(p, rate = rate, lower.tail = TRUE), q)
  expect_equal(qlaplace(p, rate = rate, lower.tail = FALSE), -q)
  expect_equal(
    qlaplace(log(p), rate = rate, lower.tail = TRUE, log.p = TRUE),
    q
  )
  expect_equal(
    qlaplace(log(p), rate = rate, lower.tail = FALSE, log.p = TRUE),
    -q
  )
})
