test_that("Discrete integration", {
  domain <- 2:5
  samplers <- 3:7
  ips_ <- data.frame(x = 3:5, weight = rep(1, 3), .block = 1L:3L)

  ips <- fm_int(domain, samplers = samplers)
  expect_identical(ips, ips_)

  domain <- as.character(domain)
  samplers <- as.character(samplers)
  ips <- fm_int(domain, samplers)
  ips_$x <- as.character(ips_$x)
  expect_identical(ips, ips_)

  domain <- factor(domain, levels = domain)
  samplers <- factor(samplers, levels = samplers)
  ips <- fm_int(domain, samplers)
  ips_$x <- factor(ips_$x, levels = domain)
  expect_identical(ips, ips_)
})


test_that("Continuous integration", {
  local_bru_safe_inla()

  domain <- INLA::inla.mesh.1d(2:5)

  samplers <- c(3, 7)
  ips_ <- data.frame(
    x = c(3:5, 3.5, 4.5),
    weight = c(1 / 6, 1 / 3, 1 / 6, 2 / 3, 2 / 3),
    .block = 1L
  )

  ips <- fm_int(domain, samplers = samplers)
  expect_identical(ips, ips_)

  # Check blockwise integration
  samplers <- rbind(c(3, 7), c(0, 10))
  ips <- fm_int(domain, samplers = samplers)
  sums <- as.vector(by(ips$weight, ips$.block, sum))
  expect_equal(sums, c(2, 3))
})


test_that("Tensor space integration", {
  local_bru_safe_inla()

  mesh_time <- INLA::inla.mesh.1d(1:5)
  mesh_space <- INLA::inla.mesh.1d(c(0, 5, 10))
  domain <- list(space = mesh_space, time = mesh_time)
  samplers1 <- tibble::tibble(
    time = rbind(c(1, 3), c(2, 4), c(3, 5)),
    space = rbind(c(0, 10), c(0, 5), c(5, 10)),
    weight = c(1, 10, 100)
  )

  ips1 <- fm_int(domain, samplers1)

  samplers2 <- tibble::tibble(
    space = samplers1$space,
    time = samplers1$time,
    weight = samplers1$weight
  )

  ips2 <- fm_int(domain, samplers2)

  expect_equal(sort(names(ips1)), sort(names(ips2)))
  expect_equal(
    dplyr::arrange(ips1, .block, time, space),
    dplyr::arrange(ips2[names(ips1)], .block, time, space)
  )
})



# From old ipoints tests

test_that("conversion of SpatialPolygon to integration points when domain is defined via a mesh", {
  local_bru_safe_inla()
  data(gorillas, package = "inlabru", envir = environment())
  ips <- fm_int(gorillas$mesh, gorillas$boundary)

  expect_s4_class(ips, "SpatialPointsDataFrame")
  expect_equal(colnames(as.data.frame(ips)), c("weight", ".block", "x", "y", "z"))
  expect_equal(sum(ips$weight),
    gorillas$boundary@polygons[[1]]@area,
    tolerance = lowtol
  )
})

test_that("conversion of whole 2D mesh to integration points", {
  local_bru_safe_inla()
  data(gorillas, package = "inlabru", envir = environment())

  ips <- fm_int(gorillas$mesh, format = "sp")

  expect_s4_class(ips, "SpatialPointsDataFrame")
  expect_equal(colnames(as.data.frame(ips)), c("weight", ".block", "x", "y", "z"))
  expect_equal(sum(ips$weight), 27.64229, tolerance = lowtol)

  ips <- fm_int(gorillas$mesh, format = "sf")

  expect_s3_class(ips, "sf")
  expect_equal(colnames(as.data.frame(ips)), c("weight", ".block", "geometry"))
  expect_equal(sum(ips$weight), 27.64229, tolerance = lowtol)
})

test_that("SLDF in metres to integration points using grouping parameter", {
  local_bru_safe_inla()

  data(mrsea, package = "inlabru", envir = environment())
  mrsea <- local_mrsea_convert(mrsea, use_km = FALSE)
  ips <- fm_int(
    list(coordinates = mrsea$mesh, season = 1:4),
    mrsea$samplers
  )

  expect_s4_class(ips, "SpatialPointsDataFrame")
  expect_true("weight" %in% colnames(as.data.frame(ips)))
  expect_true("season" %in% colnames(as.data.frame(ips)))

  # Should be a factor 1000 relative to the kilometre scale, since both schemes us
  # CRS information to convert to km, but the weight information is in metres here
  expect_equal(sum(ips$weight) / 2293712, 1, tolerance = midtol)
})

test_that("SLDF in kilometres to integration points using grouping parameter", {
  local_bru_safe_inla()

  data(mrsea, package = "inlabru", envir = environment())
  mrsea <- local_mrsea_convert(mrsea, use_km = TRUE)
  ips <- fm_int(
    list(coordinates = mrsea$mesh, season = 1:4),
    mrsea$samplers
  )

  expect_s4_class(ips, "SpatialPointsDataFrame")
  expect_true("weight" %in% colnames(as.data.frame(ips)))
  expect_true("season" %in% colnames(as.data.frame(ips)))

  expect_equal(sum(ips$weight) / 2293.712, 1, tolerance = midtol)
})


test_that("Polygon integration with holes", {
  local_bru_safe_inla()
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

  bndA <- INLA::inla.sp2segment(plyA)
  m <- INLA::inla.mesh.2d(
    loc.domain = bndA$loc,
    max.edge = 1
  )
  ipA3 <- fm_int(m, plyA, int.args = list(
    method = "direct",
    nsub2 = 1
  ))
  ipA4 <- fm_int(m, plyA, int.args = list(
    method = "stable",
    nsub2 = 1
  ))
  ipA3$test <- "A3"
  ipA4$test <- "A4"

  if (FALSE) {
    pl <- ggplot2::ggplot() +
      gg(m) +
      gg(plyA)
    pl

    pl +
      gg(ipA3, mapping = aes(col = weight, size = weight)) +
      gg(ipA4, mapping = aes(col = weight, size = weight)) +
      ggplot2::facet_wrap(vars(test))
  }

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




# Additional mesh integration tests

test_that("flat mesh integration", {
  skip_on_cran()
  local_bru_safe_inla()

  mesh <- INLA::inla.mesh.2d(cbind(0, 0), offset = 1, max.edge = 2)

  ips0 <- fm_int(mesh, int.args = list(nsub2 = 0))
  ips9 <- fm_int(mesh, int.args = list(nsub2 = 9))

  expect_equal(sum(ips0$weight), sum(ips9$weight))
})

test_that("sphere and globe mesh integration", {
  skip_on_cran()
  local_bru_safe_inla()

  mesh <- INLA::inla.mesh.create(globe = 1)

  ips0 <- fm_int(mesh, int.args = list(nsub2 = 0))
  ips9 <- fm_int(mesh, int.args = list(nsub2 = 9))

  expect_equal(sum(ips0$weight), 4 * pi)
  expect_equal(sum(ips9$weight), 4 * pi)

  mesh_ <- mesh
  mesh_$loc <- mesh_$loc * 1000

  ips0_ <- fm_int(mesh_, int.args = list(nsub2 = 0))
  ips9_ <- fm_int(mesh_, int.args = list(nsub2 = 9))

  expect_equal(sum(ips0_$weight), 4 * pi * 1e6)
  expect_equal(sum(ips9_$weight), 4 * pi * 1e6)

  suppressWarnings(
    mesh_2 <- INLA::inla.mesh.create(globe = 1, crs = fm_CRS("globe"))
  )

  ips0_2 <- fm_int(mesh_2, int.args = list(nsub2 = 0))
  ips9_2 <- fm_int(mesh_2, int.args = list(nsub2 = 9))

  expect_equal(sum(ips0_2$weight), 4 * pi * 6370.997^2)
  expect_equal(sum(ips9_2$weight), 4 * pi * 6370.997^2)

  # Corrected spherical/globe code updated in INLA fmesher on 22.11.27
  skip_if_not_installed("INLA", "22.11.28")

  ips0_3 <- fm_int(mesh_2, int.args = list(nsub2 = 0))
  ips9_3 <- fm_int(mesh_2, int.args = list(nsub2 = 9))

  expect_equal(sum(ips0_3$weight), 4 * pi * 6370.997^2)
  expect_equal(sum(ips9_3$weight), 4 * pi * 6370.997^2)
})

test_that("flat SpatialPolygons integration", {
  skip_on_cran()
  local_bru_safe_inla()

  mesh <- INLA::inla.mesh.2d(cbind(0, 0), offset = 2, max.edge = 4)

  poly <- SpatialPolygons(list(Polygons(list(Polygon(rbind(
    c(-1, -1), c(-1, 1), c(1, 1), c(1, -1)
  ))), ID = "A")))
  poly <- sf::st_as_sf(poly)

  ips0 <- fm_int(mesh, samplers = poly, int.args = list(nsub2 = 0, method = "direct"))
  ips1 <- fm_int(mesh, samplers = poly, int.args = list(nsub2 = 1, method = "direct"))
  ips9 <- fm_int(mesh, samplers = poly, int.args = list(nsub2 = 9, method = "direct"))
  ips19 <- fm_int(mesh, samplers = poly, int.args = list(nsub2 = 19, method = "direct"))

  #  ggplot() + gg(mesh) +
  #    geom_point(aes(x=x,y=y,size=weight,colour="0"),data=ips0) +
  #    geom_point(aes(x=x,y=y,size=weight,colour="1"),data=ips1) +
  #    geom_point(aes(x=x,y=y,size=weight,colour="9"),data=ips9) +
  #    geom_point(aes(x=x,y=y,size=weight,colour="19"),data=ips19)

  expect_equal(sum(ips0$weight), 6.627417, tolerance = midtol)
  expect_equal(sum(ips1$weight), 4.970563, tolerance = midtol)
  expect_equal(sum(ips9$weight), 4.108999, tolerance = midtol)
  expect_equal(sum(ips19$weight), 4.009587, tolerance = midtol)
})

test_that("globe polygon integration", {
  skip_on_cran()
  local_bru_safe_inla()

  suppressWarnings(
    mesh <- INLA::inla.mesh.create(globe = 1, crs = fm_CRS("globe"))
  )

  poly <- sp::SpatialPolygons(
    list(sp::Polygons(list(sp::Polygon(rbind(
      c(-45, -45), c(-45, 45), c(45, 45), c(45, -45)
    ))), ID = "A")),
    proj4string = fm_CRS("longlat_globe")
  )

  ips1 <- fm_int(mesh, samplers = poly, int.args = list(nsub = 2))

  expect_equal(
    nrow(ips1),
    9
  )
})
