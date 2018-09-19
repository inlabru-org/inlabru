
test_that("1D integration points can be generated", {
  ips = ipoints(c(0,10), 3, name = "myDim")
  
  expect_match(class(ips), "data.frame")
  expect_equal(nrow(ips), 3)
  expect_equal(ncol(ips), 2)
  expect_equal(names(ips), c("myDim", "weight"))
  expect_equal(as.numeric(ips[1,]), c(0, 2.5))
  expect_equal(as.numeric(ips[2,]), c(5, 5))
  
})


test_that("conversion of 1D mesh to integration points", {
  mesh = inla.mesh.1d(seq(0,10,by = 1))
  ips = ipoints(mesh, name = "time")
  
  expect_match(class(ips), "data.frame")
  expect_equal(nrow(ips), 11)
  expect_equal(ncol(ips), 2)
  expect_equal(names(ips), c("time", "weight"))
  expect_equal(as.numeric(ips[1,]), c(0, 0.5))
  expect_equal(as.numeric(ips[5,]), c(4, 1))
  expect_equal(as.numeric(ips[11,]), c(10, 0.5))
  
})

test_that("conversion of SpatialPolygon to integration points", {
  data(gorillas, package = "inlabru")
  ips = ipoints(gorillas$boundary)
  
  expect_match(class(ips), "SpatialPointsDataFrame")
  expect_equal(names(ips), "weight")
  expect_equal(colnames(data.frame(ips)), c("weight", "x", "y", "optional"))
  expect_equal(sum(ips$weight), 19.87366, tolerance = 1E-5)
})

test_that("conversion of SpatialPolygon to integration points when domain is defined via a mesh", {
  data(gorillas, package = "inlabru")
  ips = ipoints(gorillas$boundary, gorillas$mesh)
  
  expect_match(class(ips), "SpatialPointsDataFrame")
  expect_equal(colnames(data.frame(ips)), c("weight", "x", "y", "optional"))
  expect_equal(sum(ips$weight), 19.87366, tolerance = 1E-5)
  
})

test_that("conversion of 2D mesh to integration points", {
  data(gorillas, package = "inlabru")
  ips = ipoints(gorillas$mesh)
  
  expect_match(class(ips), "SpatialPointsDataFrame")
  expect_equal(colnames(data.frame(ips)), c("vertex", "weight", "x", "y", "optional"))
  expect_equal(sum(ips$weight), 27.65967, tolerance = 1E-5)
  
})

test_that("SpatialLinesDataFrame to integration points using grouping parameter", {
  data(mrsea, package = "inlabru")
  ips = ipoints(mrsea$samplers, mrsea$mesh, group = "season")
  
  expect_match(class(ips), "SpatialPointsDataFrame")
  expect_equal(colnames(data.frame(ips)), c("weight", "vertex", "season", "x", "y", "optional"))
  expect_equal(floor(sum(ips$weight)), 2288791)
})

