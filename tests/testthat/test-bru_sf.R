test_that("sf gorillas lgcp vignette", {
  ##  skip("Feature not yet implemented")

  skip_on_cran()
  local_bru_safe_inla()
  skip_if_not_installed("sf")

  # Code adapted from the lgcp vignette with an additional mesh construction
  # step instead of using the mesh provided in the gorillas data.

  # All sp objects in the example are replaced with sf equivalents.

  data(gorillas_sf, package = "inlabru", envir = environment())

  # ## Not needed anymore for this test, but kept for future reference
  # # Bug in st_as_sf making check_ring_dir=TRUE have no effect, as
  # # st_as_sfc.SpatialPolygons ignores it. To get around it, need to convert to
  # # sfc_POLYGON first, and then do a separate st_sfc call, which does use
  # # check_ring_dir.
  # # No effect:
  # b1 <- sf::st_as_sf(gorillas$boundary, check_ring_dir = TRUE)
  # # st_sfc gives the proper effect:
  # b2 <- b1
  # b2$geometry <- sf::st_sfc(sf::st_geometry(b2$geometry),
  #   check_ring_dir = TRUE)
  #
  # gorillas_sf$boundary <- b2

  ## Build boundary information:
  boundary <- list(
    gorillas_sf$boundary,
    NULL
  )

  ## Build the mesh:
  mesh_sf <- fm_mesh_2d_inla(
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
    crs = fm_crs(gorillas_sf$nests)
  )

  # library(ggplot2)
  # ggplot() +
  #   gg(mesh_sf) +
  #   geom_sf(data = gorillas_sf$boundary, alpha = 0.2, fill = "blue")

  matern <- INLA::inla.spde2.pcmatern(mesh_sf,
    prior.sigma = c(0.1, 0.01),
    prior.range = c(0.1, 0.01)
  )

  cmp <- geometry ~ mySmooth(geometry, model = matern) +
    Intercept(1)

  expect_error(
    fit <- lgcp(
      cmp,
      data = gorillas_sf$nests,
      samplers = gorillas_sf$boundary,
      domain = list(geometry = mesh_sf)
    ),
    NA
  )
})
