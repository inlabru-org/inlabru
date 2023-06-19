test_that("sf gorillas lgcp vignette", {
  ##  skip("Feature not yet implemented")

  skip_on_cran()
  local_bru_safe_inla()
  skip_if_not_installed("sf")

  # Code adapted from the lgcp vignette with an additional mesh construction
  # step instead of using the mesh provided in the gorillas data.

  # All sp objects in the example are replaced with sf equivalents.

  data(gorillas, package = "inlabru", envir = environment())
  gorillas_sf <- list()
  gorillas_sf$nests <- sf::st_as_sf(gorillas$nests)

  # Bug in st_as_sf making check_ring_dir=TRUE have no effect, as st_as_sfc.SpatialPolygons
  # ignores it. To get around it, need to convert to sfc_POLYGON first, and then do a separate
  # st_sfc call, which does use check_ring_dir.
  # No effect:
  b1 <- sf::st_as_sf(gorillas$boundary, check_ring_dir = TRUE)
  # st_sfc gives the proper effect:
  b2 <- b1
  b2$geometry <- sf::st_sfc(sf::st_geometry(b2$geometry), check_ring_dir = TRUE)

  gorillas_sf$boundary <- b2

  ## Build boundary information:
  ## inla.mesh.2d doesn't support sf yet.
  boundary <- list(
    fm_as_inla_mesh_segment(gorillas_sf$boundary),
    NULL
  )

  ## Build the mesh:
  ## Suppress PROJ support warnings, until the inla.mesh functions are updated
  suppressWarnings(
    mesh_sf <- INLA::inla.mesh.2d(
      boundary = boundary,
      max.edge = c(0.54, 0.97),
      min.angle = c(30, 21),
      ## Safeguard against large meshes.
      max.n = c(48000, 16000),
      ## Don't build a huge mesh!
      max.n.strict = c(128000, 128000),
      ## Filter away adjacent points.
      cutoff = 0.01,
      ## Offset for extra boundaries, if needed.
      offset = c(0.73, 1.55),
      ## Build mesh in this crs:
      crs = fm_CRS(gorillas$nests)
    )
  )

  # library(ggplot2)
  # ggplot() +
  #   gg(mesh_sf) +
  #   geom_sf(data = gorillas_sf$boundary, alpha = 0.2, fill = "blue")

  matern <- INLA::inla.spde2.pcmatern(mesh_sf,
    prior.sigma = c(0.1, 0.01),
    prior.range = c(5, 0.01)
  )

  cmp <- geometry ~ mySmooth(geometry, model = matern) +
    Intercept(1)

  # Check integration construction
  ips_sp <- fm_int(mesh_sf, gorillas$boundary)
  ips_sf <- fm_int(mesh_sf, gorillas_sf$boundary)

  expect_equal(
    ips_sp$weight,
    ips_sf$weight,
    tolerance = 1e-3
  )

  fit <- lgcp(
    cmp,
    data = gorillas_sf$nests,
    samplers = gorillas_sf$boundary,
    domain = list(geometry = mesh_sf)
  )
})
