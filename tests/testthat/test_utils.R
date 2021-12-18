local_bru_testthat_setup()

test_that("Missing data infilling", {
  skip_on_cran()
  local_bru_safe_inla()

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

  input <- SpatialGrid(GridTopology(c(0, 0), c(1, 1), c(7, 8)))
  input_data <- data.frame(val = as.vector(coordinates(input)[, 2]))
  input <-
    sp::SpatialGridDataFrame(
      input,
      data = input_data
    )
  val2 <- bru_fill_missing(input, points, points$val)

  expect_equal(val2, val0)
})
