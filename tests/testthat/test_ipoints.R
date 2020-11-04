context("Integration point construction (test_ipoints.R)")

test_that("1D integration points can be generated", {
  disable_PROJ6_warnings()
  ips <- ipoints(c(0, 10), 3, name = "myDim")

  expect_is(ips, "data.frame")
  expect_equal(nrow(ips), 3)
  expect_equal(ncol(ips), 2)
  expect_equal(names(ips), c("myDim", "weight"))
  expect_equal(as.numeric(ips[1, ]), c(0, 2.5))
  expect_equal(as.numeric(ips[2, ]), c(5, 5))
})


test_that("conversion of 1D mesh to integration points", {
  disable_PROJ6_warnings()
  mesh <- inla.mesh.1d(seq(0, 10, by = 1))
  ips <- ipoints(mesh, name = "time")

  expect_is(ips, "data.frame")
  expect_equal(nrow(ips), 11)
  expect_equal(ncol(ips), 2)
  expect_equal(names(ips), c("time", "weight"))
  expect_equal(as.numeric(ips[1, ]), c(0, 0.5))
  expect_equal(as.numeric(ips[5, ]), c(4, 1))
  expect_equal(as.numeric(ips[11, ]), c(10, 0.5))
})

test_that("conversion of SpatialPolygon to integration points", {
  disable_PROJ6_warnings()
  data(gorillas, package = "inlabru")
  ips <- ipoints(gorillas$boundary)

  expect_is(ips, "SpatialPointsDataFrame")
  expect_equal(names(ips), "weight")
  expect_equal(colnames(data.frame(ips)), c("weight", "x", "y", "optional"))
  expect_equal(sum(ips$weight), 19.87366, tolerance = lowtol)
})

test_that("conversion of SpatialPolygon to integration points when domain is defined via a mesh", {
  disable_PROJ6_warnings()
  data(gorillas, package = "inlabru")
  ips <- ipoints(gorillas$boundary, gorillas$mesh)

  expect_is(ips, "SpatialPointsDataFrame")
  expect_equal(colnames(data.frame(ips)), c("weight", "x", "y", "optional"))
  expect_equal(sum(ips$weight), sum(ipoints(gorillas$boundary)$weight), tolerance = hitol)
})

test_that("conversion of 2D mesh to integration points", {
  disable_PROJ6_warnings()
  data(gorillas, package = "inlabru")
  ips <- ipoints(gorillas$mesh)

  expect_is(ips, "SpatialPointsDataFrame")
  expect_equal(colnames(data.frame(ips)), c("vertex", "weight", "x", "y", "optional"))
  expect_equal(sum(ips$weight), 27.65967, tolerance = lowtol)
})

test_that("SpatialLinesDataFrame to integration points using grouping parameter", {
  disable_PROJ6_warnings()
  data(mrsea, package = "inlabru")
  mrsea <- mrsea_rebuild_CRS(mrsea, use_km = FALSE)
  ips <- ipoints(mrsea$samplers, mrsea$mesh, group = "season")
  
  expect_is(ips, "SpatialPointsDataFrame")
  expect_equal(colnames(data.frame(ips)),
               c("weight", "vertex", "season", "x", "y", "coordinateZ", "optional"))
  expect_equal(sum(ips$weight) / 2288791, 1, tolerance = midtol)

  data(mrsea, package = "inlabru")
  mrsea <- mrsea_rebuild_CRS(mrsea, use_km = TRUE)
  ips <- ipoints(mrsea$samplers, mrsea$mesh, group = "season")
  
  expect_is(ips, "SpatialPointsDataFrame")
  expect_equal(colnames(data.frame(ips)),
               c("weight", "vertex", "season", "x", "y", "z", "optional"))
  expect_equal(sum(ips$weight) / 2288791, 1, tolerance = midtol)
})
