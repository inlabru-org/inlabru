test_that("eval_spatial.SpatRaster", {
  skip_if_not_installed("terra")

  # Load the Gorilla data
  data(gorillas, package = "inlabru", envir = environment())

  covs <- c(
    terra::rast(gorillas$gcov$elevation),
    terra::rast(gorillas$gcov$elevation)
  )
  names(covs) <- c("A", "B")
  covs$B <- covs$B + 1

  # Multirow
  where <- gorillas$nests[seq_len(5), , drop = FALSE]
  A_val <- eval_spatial(covs,
    where = where,
    layer = rep("A", nrow(where))
  )
  B_val <- eval_spatial(covs,
    where = where,
    layer = rep("B", nrow(where))
  )
  AB_val <- eval_spatial(covs,
    where = where,
    layer = rep(
      c("A", "B"),
      c(2, nrow(where) - 2)
    )
  )

  expect_equal(A_val + 1, B_val)
  expect_equal(c(A_val[1:2], B_val[3:5]), AB_val)

  # Multirow
  where <- gorillas$nests[seq_len(5), , drop = FALSE]
  A1_val <- eval_spatial(covs,
    where = where,
    layer = "A"
  )
  B1_val <- eval_spatial(covs,
    where = where,
    layer = "B"
  )

  expect_equal(A1_val, A_val)
  expect_equal(B1_val, B_val)

  # Single row
  where <- gorillas$nests[5, , drop = FALSE]
  A2_val <- eval_spatial(covs,
    where = where,
    layer = "A"
  )
  B2_val <- eval_spatial(covs,
    where = where,
    layer = "B"
  )

  expect_equal(A2_val, A_val[5])
  expect_equal(B2_val, B_val[5])
})


test_that("eval_spatial.Spatial*", {
  # Load the Gorilla data
  data(gorillas, package = "inlabru", envir = environment())

  covs <- gorillas$gcov$elevation
  names(covs) <- "A"
  covs$B <- covs$A + 1

  # Multirow
  where <- gorillas$nests[seq_len(5), , drop = FALSE]
  A_val <- eval_spatial(covs,
    where = where,
    layer = rep("A", nrow(where))
  )
  B_val <- eval_spatial(covs,
    where = where,
    layer = rep("B", nrow(where))
  )
  AB_val <- eval_spatial(covs,
    where = where,
    layer = rep(
      c("A", "B"),
      c(2, nrow(where) - 2)
    )
  )

  expect_equal(A_val + 1, B_val)
  expect_equal(c(A_val[1:2], B_val[3:5]), AB_val)

  # Multirow
  where <- gorillas$nests[seq_len(5), , drop = FALSE]
  A1_val <- eval_spatial(covs,
    where = where,
    layer = "A"
  )
  B1_val <- eval_spatial(covs,
    where = where,
    layer = "B"
  )

  expect_equal(A1_val, A_val)
  expect_equal(B1_val, B_val)

  # Single row
  where <- gorillas$nests[5, , drop = FALSE]
  A2_val <- eval_spatial(covs,
    where = where,
    layer = "A"
  )
  B2_val <- eval_spatial(covs,
    where = where,
    layer = "B"
  )

  expect_equal(A2_val, A_val[5])
  expect_equal(B2_val, B_val[5])
})





test_that("eval_spatial.sf", {
  # Load the Gorilla data
  data(gorillas, package = "inlabru", envir = environment())

  nests <- sf::st_as_sf(gorillas$nests)
  plots <- sf::st_as_sf(gorillas$plotsample$plots)
  plots$something_num <- seq_len(nrow(plots))
  nests$something_num <- eval_spatial(plots, nests, layer = "something_num")

  expect_equal(
    sort(unique(nests$something), na.last = TRUE),
    c(7, 12, 13, 14, 15, 17, 18, 19, 20, 23, 24, NA)
  )

  plots$something_char <- as.character(seq_len(nrow(plots)))
  nests$something_char <- eval_spatial(plots, nests, layer = "something_char")

  expect_equal(nests$something_char, as.character(nests$something_num))
})



test_that("eval_spatial.stars", {
  skip_if_not_installed("stars")

  # Load the Gorilla data
  data(gorillas, package = "inlabru", envir = environment())

  covs <- c(
    stars::st_as_stars(gorillas$gcov$elevation),
    stars::st_as_stars(gorillas$gcov$elevation)
  )
  names(covs) <- c("A", "B")
  covs$B <- covs$B + 10000L

  # Multirow
  where <- gorillas$nests[seq_len(5), , drop = FALSE]
  A_val <- eval_spatial(covs,
    where = where,
    layer = rep("A", nrow(where))
  )
  B_val <- eval_spatial(covs,
    where = where,
    layer = rep("B", nrow(where))
  )
  AB_val <- eval_spatial(covs,
    where = where,
    layer = rep(
      c("A", "B"),
      c(2, nrow(where) - 2)
    )
  )

  expect_equal(A_val + 10000L, B_val)
  expect_equal(c(A_val[1:2], B_val[3:5]), AB_val)

  # Multirow
  where <- gorillas$nests[seq_len(5), , drop = FALSE]
  A1_val <- eval_spatial(covs,
    where = where,
    layer = "A"
  )
  B1_val <- eval_spatial(covs,
    where = where,
    layer = "B"
  )

  expect_equal(A1_val, A_val)
  expect_equal(B1_val, B_val)

  # Single row
  where <- gorillas$nests[5, , drop = FALSE]
  A2_val <- eval_spatial(covs,
    where = where,
    layer = "A"
  )
  B2_val <- eval_spatial(covs,
    where = where,
    layer = "B"
  )

  expect_equal(A2_val, A_val[5], ignore_attr = "class")
  expect_equal(B2_val, B_val[5], ignore_attr = "class")
})
