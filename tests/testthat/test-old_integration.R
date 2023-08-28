test_that("Continuous integration", {
  withr::local_options(lifecycle_verbosity = "quiet")

  domain <- fm_mesh_1d(2:5)

  samplers <- c(3, 5)
  ips_ <- data.frame(
    x = seq(3, 5, by = 0.5),
    weight = c(1 / 6, 2 / 3, 1 / 3, 2 / 3, 1 / 6),
    .block = 1L
  )
  ips_ <- ips_[order(ips_$x), ]
  rownames(ips_) <- NULL

  ips <- ipoints(domain = domain, samplers = samplers)
  ips <- ips[order(ips$x), ]
  rownames(ips) <- NULL
  expect_identical(ips, ips_, tolerance = lowtol, )

  # Check blockwise integration
  samplers <- rbind(c(3, 5), c(2, 4.5))
  ips <- ipoints(domain = domain, samplers = samplers)
  sums <- as.vector(by(ips$weight, ips$.block, sum))
  expect_equal(sums, c(2, 2.5))

  # degree = 2
  domain <- fm_mesh_1d(2:5, degree = 2)

  samplers <- c(3, 5)
  ips_ <- data.frame(
    x = c(3:5, 3.5, 4.5),
    weight = c(1 / 6, 1 / 3, 1 / 6, 2 / 3, 2 / 3),
    .block = 1L
  )
  ips_ <- ips_[order(ips_$.block, ips_$x), ]

  ips <- ipoints(domain = domain, samplers = samplers)
  ips <- ips[order(ips$.block, ips$x), ]
  expect_identical(ips, ips_)

  # Check blockwise integration
  samplers <- rbind(c(3, 5), c(2, 4.5))
  ips <- fm_int(domain, samplers = samplers)
  sums <- as.vector(by(ips$weight, ips$.block, sum))
  expect_equal(sums, c(2, 2.5))
})





test_that("conversion of SpatialPolygon to integration points when domain is defined via a mesh", {
  skip_if_not(bru_safe_sp())
  withr::local_options(lifecycle_verbosity = "quiet")

  ips <- ipoints(domain = fmexample$mesh, samplers = fmexample$boundary_sp[[1]])

  expect_s4_class(ips, "SpatialPointsDataFrame")
  expect_equal(
    sort(colnames(as.data.frame(ips))),
    sort(c("weight", ".block", "x", "y", "z"))
  )
  expect_equal(sum(ips$weight), 18.33349, tolerance = lowtol)
})

test_that("conversion of whole 2D mesh to integration points", {
  skip_if_not(bru_safe_sp())
  withr::local_options(lifecycle_verbosity = "quiet")

  ips <- sf::st_as_sf(ipoints(fmexample$mesh))

  expect_s3_class(ips, "sf")
  expect_setequal(colnames(ips), c("weight", ".vertex", "geometry"))
  expect_equal(sum(ips$weight), 64.58135, tolerance = lowtol)

  ips <- ipoints(fmexample$mesh)

  expect_s4_class(ips, "SpatialPointsDataFrame")
  expect_setequal(colnames(as.data.frame(ips)), c("weight", ".vertex", "x", "y", "z"))
  expect_equal(sum(ips$weight), 64.58135, tolerance = lowtol)
})


test_that("Polygon integration with holes", {
  withr::local_options(lifecycle_verbosity = "quiet")

  set.seed(123L)

  plyA <- sp::SpatialPolygons(list(
    sp::Polygons(
      list(
        sp::Polygon(matrix(c(0, 3, 3, 0, 0, 0, 3, 3) + runif(8) * 0.01, 4, 2),
          hole = FALSE
        ),
        sp::Polygon(matrix(c(1, 2, 2, 1, 1, 1, 2, 2) + runif(8) * 0.01, 4, 2),
          hole = TRUE
        )
      ),
      ID = "A"
    )
  ))

  bndA <- fm_as_segm(plyA)
  m <- fm_mesh_2d_inla(
    loc.domain = bndA$loc,
    max.edge = 1
  )
  ipA1 <- ipoints(domain = m, samplers = plyA, int.args = list(
    method = "direct",
    nsub2 = 1
  ))
  ipA2 <- ipoints(domain = m, samplers = plyA, int.args = list(
    method = "stable",
    nsub2 = 1
  ))
  ipA3 <- ipoints(domain = m, samplers = plyA, int.args = list(
    method = "direct"
  ))
  ipA4 <- ipoints(domain = m, samplers = plyA, int.args = list(
    method = "stable"
  ))
  ipA5 <- fm_int(domain = m, samplers = plyA, int.args = list(
    method = "direct"
  ))
  ipA6 <- fm_int(domain = m, samplers = plyA, int.args = list(
    method = "stable"
  ))
  ipA1$type <- "low"
  ipA2$type <- "low"
  ipA3$type <- "ipoints"
  ipA4$type <- "ipoints"
  ipA5$type <- "fm_int"
  ipA6$type <- "fm_int"
  ipA1$method <- "direct"
  ipA2$method <- "stable"
  ipA3$method <- "direct"
  ipA4$method <- "stable"
  ipA5$method <- "direct"
  ipA6$method <- "stable"

  # if (FALSE) {
  #   require("ggplot2")
  #   pl <- ggplot2::ggplot() +
  #     geom_fm(data = m, alpha = 0.2) +
  #     gg(plyA)
  #   pl
  #
  #   pl +
  #     gg(ipA1, mapping = aes(col = weight, size = weight)) +
  #     gg(ipA3, mapping = aes(col = weight, size = weight)) +
  #     gg(ipA5, mapping = aes(col = weight, size = weight)) +
  #     gg(ipA2, mapping = aes(col = weight, size = weight)) +
  #     gg(ipA4, mapping = aes(col = weight, size = weight)) +
  #     gg(ipA6, mapping = aes(col = weight, size = weight)) +
  #     ggplot2::facet_grid(rows = vars(method), cols = vars(type)) +
  #     scale_size_area()
  # }

  # > sf::st_area(sf::st_as_sf(plyA))
  # [1] 8.006112

  expect_equal(sum(ipA1$weight), 7.781, tolerance = midtol)
  expect_equal(sum(ipA2$weight), 7.781, tolerance = midtol)
  expect_equal(sum(ipA3$weight), 8.004, tolerance = midtol)
  expect_equal(sum(ipA4$weight), 8.004, tolerance = midtol)
  expect_equal(sum(ipA5$weight), 8.004, tolerance = midtol)
  expect_equal(sum(ipA6$weight), 8.004, tolerance = midtol)
})



test_that("SpatialLines integration", {
  skip_if_not(bru_safe_sp())
  withr::local_options(lifecycle_verbosity = "quiet")

  bnd <- fm_segm(cbind(c(0, 1, 0.9, 0), c(0, 0, 1, 1)), is.bnd = TRUE)
  mesh <- fm_mesh_2d_inla(boundary = bnd)

  theline <-
    sf::as_Spatial(sf::st_zm(fm_as_sfc(fm_segm(rbind(
      c(0.8, 0.2), c(0.2, 0.8)
    ), is.bnd = FALSE))))

  ips <- ipoints(domain = mesh, samplers = theline)
  ips_fmesher <- fm_int(domain = mesh, samplers = theline)

  expect_equal(as.data.frame(ips), as.data.frame(ips_fmesher))

  theline <- SpatialLinesDataFrame(
    theline,
    data = data.frame(theblock = 4),
    match.ID = FALSE
  )

  ips <- ipoints(domain = mesh, samplers = theline, group = "theblock")
  ips_fmesher <- fm_int(
    domain = list(
      coordinates = mesh,
      theblock = sort(unique(theline[["theblock"]]))
    ),
    samplers = theline
  )

  expect_equal(as.data.frame(ips), as.data.frame(ips_fmesher))
})
