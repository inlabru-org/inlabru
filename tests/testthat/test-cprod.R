local_bru_testthat_setup()

# sp sf test
test_that("cprod(..., na.rm = TRUE) sf output can be generated", {
  local_bru_safe_inla()
  sf_obj1 <- sf::st_as_sf(data.frame(x=1:3, y=3:5), coords=c("x","y"))
  sf_obj2 <- sf::st_as_sf(data.frame(x=3:6, y=5:8), coords=c("x","y"))
  ips <- cprod(sf_obj1, sf_obj2, na.rm = TRUE)

  expect_s3_class(ips, "sf")
  expect_equal(nrow(ips), 1)
  expect_equal(ncol(ips), 2)
  expect_equal(names(ips), c("geometry", "weight"))
  expect_equal(as.numeric(unlist(sf::st_geometry(ips))), c(3, 5))
  expect_equal(as.numeric(ips[[2]]), 1)
})

test_that("cprod(..., na.rm = FALSE) sf output can be generated", {
  local_bru_safe_inla()
  sf_obj1 <- sf::st_as_sf(data.frame(x=1:3, y=3:5), coords=c("x","y"))
  sf_obj2 <- sf::st_as_sf(data.frame(x=3:6, y=5:8), coords=c("x","y"))
  ips <- cprod(sf_obj1, sf_obj2, na.rm = FALSE)

  expect_s3_class(ips, "sf")
  expect_equal(nrow(ips), 6)
  expect_equal(ncol(ips), 2)
  expect_equal(names(ips), c("geometry", "weight"))
  expect_equal(as.numeric(unlist(sf::st_geometry(ips))),
               c(1, 3, 2, 4, 3, 5, 4, 6, 5, 7, 6, 8))
  expect_equal(as.numeric(ips[[2]]), rep(1,6))
})

test_that("cprod(na.rm = TRUE) sp output can be generated", {
  local_bru_safe_inla()
  sf_obj1 <- sf::st_as_sf(data.frame(x=1:3, y=3:5), coords=c("x","y"))
  sf_obj2 <- sf::st_as_sf(data.frame(x=3:6, y=5:8), coords=c("x","y"))
  if(require(sp, quietly = TRUE)){
    sp_obj1 <- as(sf_obj1, 'Spatial')
    sp_obj2 <- as(sf_obj2, 'Spatial')
  }
  ips <- cprod(sp_obj1, sp_obj2, na.rm = TRUE)

  expect_s4_class(ips, "Spatial") # or precisely SpatialPointsDataFrame
  expect_equal(nrow(ips), 1)
  expect_equal(ncol(ips), 1)
  expect_equal(names(ips), "weight")
  expect_equal(as.numeric(ips@coords), c(3, 5))
  expect_equal(as.numeric(ips[[1]]), 1)
})

test_that("cprod(na.rm = FALSE) sp output can be generated", {
  local_bru_safe_inla()
  sf_obj1 <- sf::st_as_sf(data.frame(x=1:3, y=3:5), coords=c("x","y"))
  sf_obj2 <- sf::st_as_sf(data.frame(x=3:6, y=5:8), coords=c("x","y"))
  if(require(sp, quietly = TRUE)){
    sp_obj1 <- as(sf_obj1, 'Spatial')
    sp_obj2 <- as(sf_obj2, 'Spatial')
  }
  ips <- cprod(sp_obj1, sp_obj2, na.rm = FALSE)

  expect_s4_class(ips, "Spatial") # or precisely SpatialPointsDataFrame
  expect_equal(nrow(ips), 6)
  expect_equal(ncol(ips), 1)
  expect_equal(names(ips), "weight")
  expect_equal(as.numeric(ips@coords),
               c(1, 3, 2, 4, 3, 5, 4, 6, 5, 7, 6, 8))
  expect_equal(as.numeric(ips[[1]]), rep(1,6))
})

# ips <- ipoints(c(0, 10), 50, name = "myDim")
# ips <- ipoints(matrix(c(0, 3, 5, 10), nrow = 2, byrow = TRUE), 50)
# mesh <- inla.mesh.1d(seq(0, 10, by = 1))
# ips <- ipoints(mesh, name = "time")
#
# ips <- cprod(
#        ipoints(c(0, 10), 10, name = "x"),
#        ipoints(c(1, 5), 10, name = "y")
#        )
#
# ips1 <- ipoints(rbind(c(0, 3), c(3, 8)), 17, name = "myDim")
# ips2 <- ipoints(domain = c(1, 2, 4), name = "myDiscreteDim")
# cprod(ips1, ips2)
#
#
# data(gorillas, package = "inlabru")
# ips1 <- ipoints(gorillas$boundary)
# cprod(ips1,ips1)
# ips2 <- ipoints(gorillas$boundary, domain = gorillas$mesh)
