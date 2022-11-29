local_bru_testthat_setup()

test_that("flat mesh integration", {
  skip_on_cran()
  local_bru_safe_inla()

  mesh <- INLA::inla.mesh.2d(cbind(0, 0), offset = 1, max.edge = 2)

  ips0 <- mesh_triangle_integration(mesh, nsub = 0)
  ips9 <- mesh_triangle_integration(mesh, nsub = 9)

  expect_equal(sum(ips0$weight), sum(ips9$weight))
})

test_that("sphere and globe mesh integration", {
  skip_on_cran()
  local_bru_safe_inla()

  mesh <- INLA::inla.mesh.create(globe = 1)

  ips0 <- mesh_triangle_integration(mesh, nsub = 0)
  ips9 <- mesh_triangle_integration(mesh, nsub = 9)

  expect_equal(sum(ips0$weight), 4 * pi)
  expect_equal(sum(ips9$weight), 4 * pi)

  mesh_ <- mesh
  mesh_$loc <- mesh_$loc * 1000

  ips0_ <- mesh_triangle_integration(mesh_, nsub = 0)
  ips9_ <- mesh_triangle_integration(mesh_, nsub = 9)

  expect_equal(sum(ips0_$weight), 4 * pi * 1e6)
  expect_equal(sum(ips9_$weight), 4 * pi * 1e6)

  suppressWarnings(
    mesh_2 <- INLA::inla.mesh.create(globe = 1, crs = fm_CRS("globe"))
  )

  ips0_2 <- mesh_triangle_integration(mesh_2, nsub = 0)
  ips9_2 <- mesh_triangle_integration(mesh_2, nsub = 9)

  expect_equal(sum(ips0_2$weight), 4 * pi * 6370.997^2)
  expect_equal(sum(ips9_2$weight), 4 * pi * 6370.997^2)

  # Corrected spherical/globe code updated in INLA fmesher on 22.11.27
  skip_if_not_installed("INLA", "22.11.28")

  ips0_3 <- bru_int_polygon(mesh_2, nsub = 0)
  ips9_3 <- bru_int_polygon(mesh_2, nsub = 9)

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

  ips0 <- bru_int_polygon(mesh, samplers = poly, nsub = 0, method = "direct")
  ips1 <- bru_int_polygon(mesh, samplers = poly, nsub = 1, method = "direct")
  ips9 <- bru_int_polygon(mesh, samplers = poly, nsub = 9, method = "direct")
  ips19 <- bru_int_polygon(mesh, samplers = poly, nsub = 19, method = "direct")

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

  poly <- SpatialPolygons(
    list(Polygons(list(Polygon(rbind(
      c(-45, -45), c(-45, 45), c(45, 45), c(45, -45)
    ))), ID = "A")),
    proj4string = fm_CRS("longlat_globe")
  )

  ips0 <- bru_int_polygon(mesh, samplers = poly, nsub = 2)

  expect_equal(
    nrow(ips0),
    9
  )
})
