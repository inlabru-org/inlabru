local_bru_testthat_setup()

test_that("sf standards compliance: basic polygons", {
  out <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0), ncol = 2, byrow = TRUE)
  hole1 <- matrix(c(1, 1, 1, 2, 2, 2, 2, 1, 1, 1), ncol = 2, byrow = TRUE)
  hole2 <- matrix(c(5, 5, 5, 6, 6, 6, 6, 5, 5, 5), ncol = 2, byrow = TRUE)

  goodp1 <- sf::st_polygon(list(out))
  expect_true(st_check_polygon(goodp1))

  goodp2 <- sf::st_polygon(list(out, hole1, hole2))
  expect_true(st_check_polygon(goodp2))

  # A "hole" outside the main ring
  badp1 <- sf::st_polygon(list(out, hole1, hole2, out + 15))
  expect_false(st_check_polygon(badp1))

  # A "hole" overlapping another hole
  badp2 <- sf::st_polygon(list(out, hole1, hole1 + 1))
  expect_false(st_check_polygon(badp2))
})

# FL note: see the inla.mesh.segment documentation for the idx and is.bnd arguments
# sf stores polygons with repeated coordinates for the first and last point,
# whereas inla.mesh.segment requires the index into loc to do that job;
# For is.bnd=TRUE, the idx default closes the polygon, but not for is.bnd=FALSE,
# which is used for "linestring" type information.
# Further note: inla.mesh.segment has two ways of specifying the index;
# as a sequence, or as a two-column matrix. But it is always stored as a two column matrix,
# with no general guarantee that one line connects to the one in the next row.
# This makes conversion to sp and sf polygons more difficult, which is why there
# isn't a general fm_as_sp.inla.mesh.segment method. There is some code in various places,
# including in 'excursions' that could be used as a starting point to doing it properly
# for sf conversion.
#
# The fm_as_inla_mesh_segment.SpatialPoints method had a bug, w.r.t is.bnd
# handling, and has now been fixed.

test_that("Conversion from sfc_POINT to inla.mesh.segment", {

  local_bru_safe_inla()

  ## sfc_POINT ##

  # compare inla.mesh.segment with matrix input
  # to fm_as_inla_mesh_segment with sf point input

  # matrix version
  loc.bnd <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
  segm.bnd <- INLA::inla.mesh.segment(loc.bnd,
                                is.bnd = TRUE
  )

  # sf version
  loc.sf <- sf::st_as_sf(as.data.frame(loc.bnd),
                         coords = c(1, 2)
  )

  segm.bnd.sf <- fm_as_inla_mesh_segment(loc.sf, is.bnd = TRUE)

  expect_identical(segm.bnd.sf, segm.bnd)
  str(segm.bnd)
  str(segm.bnd.sf)

  crs <- sf::st_crs(sf::st_geometry(loc.sf))

  # check warning message for xyz (there shouldn't be one)
  loc.sf.xyz <- sf::st_as_sf(as.data.frame(cbind(loc.bnd, 1)),
                             coords = c(1, 2, 3)
  )
  class(loc.sf.xyz)
  class(sf::st_geometry(loc.sf.xyz))
  class(sf::st_geometry(loc.sf.xyz)[[1]])

  segm.bnd.sf.xym <- fm_as_inla_mesh_segment(loc.sf.xyz)

})


test_that("Conversion from sfc_LINESTRING to inla.mesh.segment", {

  local_bru_safe_inla()

  ## scf_LINESTRING ##

  pts1 <- rbind(c(0, 3), c(0, 4), c(1, 5), c(2, 5))
  pts2 <- rbind(c(1, 1), c(0, 0), c(0, -1), c(-2, -2))
  seg1 <- INLA::inla.mesh.segment(
    loc = pts1,
    idx = seq_len(nrow(pts1)),
    is.bnd = FALSE
  )

  seg2 <- INLA::inla.mesh.segment(
    loc = pts2,
    idx = seq_len(nrow(pts2)),
    is.bnd = FALSE
  )

  seg <- fm_internal_sp2segment_join(list(seg1, seg2),
                                     grp = seq_len(2))
  expect_identical(seg$grp, as.matrix(rep(1:2, each = 3)))

  line_str1 <- sf::st_linestring(pts1)
  line_str2 <- sf::st_linestring(pts2)
  line_sfc <- sf::st_as_sfc(list(line_str1, line_str2))
  line_sf <- sf::st_sf(geometry = line_sfc)
  seg_sf <- fm_as_inla_mesh_segment(line_sf)

  expect_identical(seg_sf, seg)

  #  str(seg)
  #  str(seg_sf)

})


test_that("Conversion from sfc_POLYGON to inla.mesh.segment", {

  local_bru_safe_inla()

  ## scf_POLYGON ##

  pts1 <- rbind(c(0, 3), c(0, 4), c(1, 5), c(2, 5), c(0, 3))
  pts2 <- rbind(c(1, 2), c(0, 0), c(0, -1), c(-2, -2), c(1, 2))
  seg1 <- INLA::inla.mesh.segment(
    loc = pts1[1:4, , drop = FALSE],
    is.bnd = TRUE
  )

  seg2 <- INLA::inla.mesh.segment(
    loc = pts2[1:4, , drop = FALSE],
    is.bnd = TRUE
  )

  seg <- fm_internal_sp2segment_join(list(seg1, seg2),
                                     grp = seq_len(2))
  expect_identical(seg$grp, as.matrix(rep(1:2, each = 4)))

  line_str1 <- sf::st_polygon(list(pts1))
  line_str2 <- sf::st_polygon(list(pts2))
  line_sfc <- sf::st_as_sfc(list(line_str1, line_str2))
  line_sf <- sf::st_sf(geometry = line_sfc)
  seg_sf <- fm_as_inla_mesh_segment(line_sf)

  expect_identical(seg_sf, seg)

  #  str(seg)
  #  str(seg_sf)

})

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
    as.inla.mesh.segment(gorillas_sf$boundary),
    NULL)

  ## Build the mesh:
  mesh_sf <- inla.mesh.2d(
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


