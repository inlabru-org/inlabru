local_bru_testthat_setup()

test_that("1D integration points can be generated", {
  local_bru_safe_inla()
  ips <- ipoints(c(0, 10), 3, name = "myDim")

  expect_s3_class(ips, "data.frame")
  expect_equal(nrow(ips), 3)
  expect_equal(ncol(ips), 2)
  expect_equal(names(ips), c("myDim", "weight"))
  expect_equal(as.numeric(ips[1, ]), c(5 / 3, 10 / 3))
  expect_equal(as.numeric(ips[2, ]), c(15 / 3, 10 / 3))
  expect_equal(as.numeric(ips[3, ]), c(25 / 3, 10 / 3))
})


test_that("conversion of 1D mesh to integration points", {
  local_bru_safe_inla()
  mesh <- INLA::inla.mesh.1d(seq(0, 10, by = 1))
  ips <- ipoints(mesh, name = "time")

  expect_s3_class(ips, "data.frame")
  expect_equal(nrow(ips), 11)
  expect_equal(ncol(ips), 2)
  expect_equal(names(ips), c("time", "weight"))
  expect_equal(as.numeric(ips[1, ]), c(0, 0.5))
  expect_equal(as.numeric(ips[5, ]), c(4, 1))
  expect_equal(as.numeric(ips[11, ]), c(10, 0.5))
})

test_that("conversion of SpatialPolygon to integration points", {
  local_bru_safe_inla()
  data(gorillas, package = "inlabru", envir = environment())
  expect_warning(
    ips <- ipoints(gorillas$boundary),
    "Computing integration points from polygon"
  )

  expect_s4_class(ips, "SpatialPointsDataFrame")
  expect_equal(names(ips), "weight")
  expect_equal(colnames(data.frame(ips)), c("weight", "x", "y", "optional"))
  expect_equal(sum(ips$weight), 19.87366, tolerance = lowtol)
})

test_that("conversion of SpatialPolygon to integration points when domain is defined via a mesh", {
  local_bru_safe_inla()
  data(gorillas, package = "inlabru", envir = environment())
  ips <- ipoints(gorillas$boundary, gorillas$mesh)
  expect_warning(
    ips_nodomain <- ipoints(gorillas$boundary),
    "Computing integration points from polygon"
  )

  expect_s4_class(ips, "SpatialPointsDataFrame")
  expect_equal(colnames(data.frame(ips)), c("weight", "x", "y", "optional"))
  expect_equal(sum(ips$weight),
    sum(ips_nodomain$weight),
    tolerance = midtol
  )
})

test_that("conversion of 2D mesh to integration points", {
  local_bru_safe_inla()
  data(gorillas, package = "inlabru", envir = environment())
  ips <- ipoints(gorillas$mesh)

  expect_s4_class(ips, "SpatialPointsDataFrame")
  expect_equal(colnames(data.frame(ips)), c("vertex", "weight", "x", "y", "optional"))
  expect_equal(sum(ips$weight), 27.64229, tolerance = lowtol)
})

test_that("SLDF in metres to integration points using grouping parameter", {
  local_bru_safe_inla()
  skip_if_not(fm_has_PROJ6())

  data(mrsea, package = "inlabru", envir = environment())
  mrsea <- local_mrsea_convert(mrsea, use_km = FALSE)
  suppressWarnings(
    ips <- ipoints(mrsea$samplers, mrsea$mesh, group = "season")
  )
  # # For rgdal 1.5-23, we had this version
  #  expect_warning(
  #    ...,
  #    "export to PROJ failed: generic error of unknown origin"
  #  )

  expect_s4_class(ips, "SpatialPointsDataFrame")
  expect_equal(
    colnames(data.frame(ips)),
    c("weight", "vertex", "season", "x", "y", "coordinateZ", "optional")
  )
  # Should be a factor 1000 relative to the kilometre scale, since both schemes us
  # CRS information to convert to km, but the weight information is in metres here
  expect_equal(sum(ips$weight) / 2293712, 1, tolerance = midtol)
})

test_that("SLDF in kilometres to integration points using grouping parameter", {
  local_bru_safe_inla()
  skip_if_not(fm_has_PROJ6())

  data(mrsea, package = "inlabru", envir = environment())
  mrsea <- local_mrsea_convert(mrsea, use_km = TRUE)
  suppressWarnings(
    ips <- ipoints(mrsea$samplers, mrsea$mesh, group = "season")
  )

  expect_s4_class(ips, "SpatialPointsDataFrame")
  expect_equal(
    colnames(data.frame(ips)),
    c("weight", "vertex", "season", "x", "y", "coordinateZ", "optional")
  )
  expect_equal(sum(ips$weight) / 2293.712, 1, tolerance = midtol)
})


test_that("Polygon integration with holes", {
  local_bru_safe_inla()
  set.seed(123L)

  plyA <- sp::SpatialPolygons(list(
    sp::Polygons(
      list(
        sp::Polygon(matrix(c(0, 3, 3, 0, 0, 0, 3, 3) + runif(8) * 0.01, 4, 2),
                    hole = FALSE),
        sp::Polygon(matrix(c(1, 2, 2, 1, 1, 1, 2, 2) + runif(8) * 0.01, 4, 2),
                    hole = TRUE)
      ),
      ID = "A"
    )
  ))

  bndA <- INLA::inla.sp2segment(plyA)
  m <- INLA::inla.mesh.2d(
    loc.domain = bndA$loc,
    max.edge = 1
  )
  ipA1 <- ipoints(plyA, m, int.args = list(poly_method = "legacy",
                                           method = "direct",
                                           nsub2 = 1))
  ipA2 <- ipoints(plyA, m, int.args = list(poly_method = "legacy",
                                           method = "stable",
                                           nsub2 = 1))
  ipA3 <- ipoints(plyA, m, int.args = list(method = "direct",
                                           nsub2 = 1))
  ipA4 <- ipoints(plyA, m, int.args = list(method = "stable",
                                           nsub2 = 1))
  ipA1$test <- "A1"
  ipA2$test <- "A2"
  ipA3$test <- "A3"
  ipA4$test <- "A4"

  if (FALSE) {
    pl <- ggplot2::ggplot() +
      gg(m) +
      gg(plyA)
    pl

    pl +
      gg(ipA1, mapping = aes(col = weight, size = weight)) +
      gg(ipA2, mapping = aes(col = weight, size = weight)) +
      gg(ipA3, mapping = aes(col = weight, size = weight)) +
      gg(ipA4, mapping = aes(col = weight, size = weight)) +
      ggplot2::facet_wrap(vars(test))
  }

  expect_equal(sum(ipA1$weight), 8.779508, tolerance = midtol)
  expect_equal(sum(ipA2$weight), 8.779508, tolerance = midtol)

  expect_equal(sum(ipA3$weight), 7.780959, tolerance = midtol)
  expect_equal(sum(ipA4$weight), 7.780959, tolerance = midtol)
})


test_that("Integration line splitting", {
  local_bru_safe_inla()

  mesh <- INLA::inla.mesh.2d(
    loc.domain = cbind(0, 0),
    offset = 2,
    max.edge = 0.5
  )

  expect_error(
    object = {
      sl <- split_lines(
        mesh,
        sp = rbind(c(-1, 0), c(-1, 1)),
        ep = rbind(c(1, 0), c(1, 1))
      )
    },
    NA
  )

  # Check issue #63 (problem for single line input), fixed
  expect_error(
    object = {
      sl <- split_lines(
        mesh,
        sp = cbind(-1, 0),
        ep = cbind(1, 0)
      )
    },
    NA
  )

  # Check if empty input is ok
  expect_error(
    object = {
      sl <- split_lines(
        mesh,
        sp = matrix(0, 0, 2),
        ep = matrix(0, 0, 2)
      )
    },
    NA
  )
})
