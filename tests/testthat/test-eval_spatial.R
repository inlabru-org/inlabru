test_that("eval_spatial.SpatRaster", {
  skip_if_not_installed("terra")

  # Load the Gorilla data
  data(gorillas_sf, package = "inlabru", envir = environment())
  gcov <- gorillas_sf_gcov()

  covs <- c(
    gcov$elevation,
    gcov$elevation
  )
  names(covs) <- c("A", "B")
  covs$B <- covs$B + 1

  # Multirow
  where <- gorillas_sf$nests[seq_len(5), , drop = FALSE]
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
  where <- gorillas_sf$nests[seq_len(5), , drop = FALSE]
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
  where <- gorillas_sf$nests[5, , drop = FALSE]
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
  skip_if_not(bru_safe_sp())

  # Load the Gorilla data
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  gorillas <- gorillas_sp()

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
  data(gorillas_sf, package = "inlabru", envir = environment())

  nests <- gorillas_sf$nests
  plots <- gorillas_sf$plotsample$plots
  plots$something_num <- seq_len(nrow(plots))
  nests$something_num <- eval_spatial(plots, nests, layer = "something_num")
  something_num <- eval_spatial(plots, plots, layer = "something_num")

  expect_equal(
    sort(unique(nests$something_num), na.last = TRUE),
    c(3, 6, 8, 9, 16, 17, 19, 20, 21, 22, 25, NA)
  )

  plots$something_char <- as.character(seq_len(nrow(plots)))
  nests$something_char <- eval_spatial(plots, nests, layer = "something_char")

  expect_equal(nests$something_char, as.character(nests$something_num))
})



test_that("eval_spatial.stars", {
  skip_if_not_installed("terra")
  skip_if_not_installed("stars")

  # Load the Gorilla data
  data(gorillas_sf, package = "inlabru", envir = environment())
  gcov <- gorillas_sf_gcov()

  covs <- c(
    gcov$elevation,
    gcov$elevation
  )
  names(covs) <- c("A", "B")
  covs$B <- covs$B + 10000L

  # Multirow
  where <- gorillas_sf$nests[seq_len(5), , drop = FALSE]
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
  where <- gorillas_sf$nests[seq_len(5), , drop = FALSE]
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
  where <- gorillas_sf$nests[5, , drop = FALSE]
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
