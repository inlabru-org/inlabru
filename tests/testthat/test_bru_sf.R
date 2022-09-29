test_that("sf gorillas lgcp vignette", {

  # Code adapted from the lgcp vignette with an additional mesh construction
  # step instead of using the mesh provided in the gorillas data.

  # All sp objects in the example are replaced with sf equivalents.

  data(gorillas, package = "inlabru")
  gorillas_sf = list()
  gorillas_sf$nests = sf::st_as_sf(gorillas$nests)
  gorillas_sf$boundary = sf::st_as_sf(gorillas$boundary)

  ## Build boundary information:
  ## (fmesher supports SpatialPolygons, but this app is not (yet) intelligent enough for that.)
  boundary <- list(
    fm_as_inla_mesh_segment(gorillas_sf$boundary),
    NULL)

  ## Build the mesh:
  mesh_sf <- INLA::inla.mesh.2d(
    boundary = boundary,
    max.edge = c(0.24, 0.97),
    min.angle = c(30, 21),
    max.n = c(48000, 16000),
    ## Safeguard against large meshes.
    max.n.strict = c(128000, 128000),
    ## Don't build a huge mesh!
    cutoff = 0.01,
    ## Filter away adjacent points.
    offset = c(0.73, 1.55)
  ) ## Offset for extra boundaries, if needed.

  matern <- INLA::inla.spde2.pcmatern(gorillas$mesh,
                                      prior.sigma = c(0.1, 0.01),
                                      prior.range = c(5, 0.01))

  cmp <- coordinates ~ mySmooth(main = coordinates, model = matern) + Intercept(1)

  fit <- lgcp(
    cmp,
    data = gorillas_sf$nests,
    samplers = gorillas_sf$boundary,
    domain = list(coordinates = gorillas$mesh)
  )

})
