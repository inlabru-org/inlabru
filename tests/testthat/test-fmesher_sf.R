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
#
# When relying on ring orientation, need to ensure that stored sf objects
# are in correct sf standard CCW orientation; sf doesn't take into account
# that geos has CW as canonical orientation, so output from st_* function
# do not always follow the sf standard:
#   sfc <- sf::st_sfc(x, check_ring_dir = TRUE)
# See https://github.com/r-spatial/sf/issues/2096


test_that("Conversion from sfc_POINT to inla.mesh.segment", {
  local_bru_safe_inla()

  ## sfc_POINT ##

  # compare inla.mesh.segment with matrix input
  # to fm_as_inla_mesh_segment with sf point input

  # matrix version
  loc.bnd <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
  segm.bnd <- INLA::inla.mesh.segment(loc.bnd,
    is.bnd = TRUE,
    crs = fm_CRS()
  )

  # sf version
  loc.sf <- sf::st_as_sf(as.data.frame(loc.bnd),
    coords = c(1, 2)
  )

  segm.bnd.sf <- fm_as_inla_mesh_segment(loc.sf, is.bnd = TRUE)

  expect_identical(segm.bnd.sf, segm.bnd)
  #  str(segm.bnd)
  #  str(segm.bnd.sf)

  crs <- sf::st_crs(sf::st_geometry(loc.sf))

  # check warning message for xyz (there shouldn't be one)
  loc.sf.xyz <- sf::st_as_sf(as.data.frame(cbind(loc.bnd, 1)),
    coords = c(1, 2, 3)
  )
  class(loc.sf.xyz)
  class(sf::st_geometry(loc.sf.xyz))
  class(sf::st_geometry(loc.sf.xyz)[[1]])

  expect_no_warning(fm_as_inla_mesh_segment(loc.sf.xyz))
})


test_that("Conversion from sfc_LINESTRING to inla.mesh.segment", {
  local_bru_safe_inla()

  ## scf_LINESTRING ##

  pts1 <- rbind(c(0, 3), c(0, 4), c(1, 5), c(2, 5))
  pts2 <- rbind(c(1, 1), c(0, 0), c(0, -1), c(-2, -2))
  seg1 <- INLA::inla.mesh.segment(
    loc = pts1,
    idx = seq_len(nrow(pts1)),
    is.bnd = FALSE,
    crs = fm_CRS()
  )

  seg2 <- INLA::inla.mesh.segment(
    loc = pts2,
    idx = seq_len(nrow(pts2)),
    is.bnd = FALSE,
    crs = fm_CRS()
  )

  seg <- fm_internal_sp2segment_join(list(seg1, seg2),
    grp = seq_len(2)
  )
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

  pts0 <- rbind(c(-7, -7), c(7, -7), c(7, 7), c(-7, 7), c(-7, -7)) # covering (CCW)
  pts0b <- pts0 + 10
  pts1 <- rbind(c(0, 3), c(0, 4), c(1, 5), c(2, 5), c(0, 3)) # hole (CW)
  pts2 <- rbind(c(1, 2), c(0, 0), c(0, -1), c(-2, -2), c(1, 2)) # hole (CW)
  seg0 <- INLA::inla.mesh.segment(
    loc = pts0[1:4, , drop = FALSE],
    is.bnd = TRUE,
    crs = fm_CRS()
  )
  seg0b <- INLA::inla.mesh.segment(
    loc = pts0b[1:4, , drop = FALSE],
    is.bnd = TRUE,
    crs = fm_CRS()
  )

  seg1 <- INLA::inla.mesh.segment(
    loc = pts1[1:4, , drop = FALSE],
    is.bnd = TRUE,
    crs = fm_CRS()
  )

  seg2 <- INLA::inla.mesh.segment(
    loc = pts2[1:4, , drop = FALSE],
    is.bnd = TRUE,
    crs = fm_CRS()
  )

  seg <- fm_internal_sp2segment_join(list(seg0, seg1, seg2, seg0b),
    grp = c(1, 1, 1, 2)
  )
  expect_identical(seg$grp, as.matrix(rep(c(1L, 1L, 1L, 2L), each = 4)))

  line_str1 <- sf::st_polygon(list(pts0, pts1, pts2))
  line_str2 <- sf::st_polygon(list(pts0b))
  line_sfc <- sf::st_as_sfc(list(line_str1, line_str2))
  line_sf <- sf::st_sf(geometry = line_sfc)
  seg_sf <- fm_as_inla_mesh_segment(line_sf)

  expect_identical(seg_sf, seg)

  #  str(seg)
  #  str(seg_sf)
})



test_that("Conversion from sfc_MULTIPOLYGON to inla.mesh.segment", {
  local_bru_safe_inla()

  ## scf_MULTIPOLYGON ##

  pts0 <- rbind(c(-7, -7), c(7, -7), c(7, 7), c(-7, 7), c(-7, -7)) # covering (CCW)
  pts0b <- pts0 + 15
  pts1 <- rbind(c(0, 3), c(0, 4), c(1, 5), c(2, 5), c(0, 3)) # hole (CW)
  pts2 <- rbind(c(1, 2), c(0, 0), c(0, -1), c(-2, -2), c(1, 2)) # hole (CW)
  seg0 <- INLA::inla.mesh.segment(
    loc = pts0[1:4, , drop = FALSE],
    is.bnd = TRUE,
    crs = fm_CRS()
  )
  seg0b <- INLA::inla.mesh.segment(
    loc = pts0b[1:4, , drop = FALSE],
    is.bnd = TRUE,
    crs = fm_CRS()
  )

  seg1 <- INLA::inla.mesh.segment(
    loc = pts1[1:4, , drop = FALSE],
    is.bnd = TRUE,
    crs = fm_CRS()
  )

  seg2 <- INLA::inla.mesh.segment(
    loc = pts2[1:4, , drop = FALSE],
    is.bnd = TRUE,
    crs = fm_CRS()
  )

  seg_1 <- fm_internal_sp2segment_join(list(seg0, seg1, seg2, seg0b),
    grp = c(1, 1, 1, 1)
  )
  expect_identical(seg_1$grp, as.matrix(rep(c(1L), each = 16)))
  seg_2 <- fm_internal_sp2segment_join(
    list(
      seg0, seg1, seg2, seg0b,
      seg0, seg1, seg2, seg0b
    ),
    grp = c(1, 1, 1, 1, 2, 2, 2, 2)
  )
  expect_identical(seg_2$grp, as.matrix(rep(c(1L, 2L), each = 16)))

  line_str1 <- sf::st_polygon(list(pts0, pts1, pts2))
  line_str2 <- sf::st_polygon(list(pts0b))

  line_sfc <- sf::st_sfc(list(line_str1, line_str2), check_ring_dir = TRUE)

  line_sfc_c <- sf::st_combine(line_sfc)
  line_sfc_c <- c(line_sfc_c, line_sfc_c)
  seg_sf_c <- fm_as_inla_mesh_segment(line_sfc_c)

  expect_identical(seg_sf_c, seg_2)

  # sf::st_union might return CW order since that's the canonical geos ordering,
  # and sf doesn't account for that, despite the sf standard being CCW ordering.
  # The conversion functions start by ensuring appropriate ordering.
  line_sfc_u <- sf::st_union(line_sfc)
  seg_sf_u <- fm_as_inla_mesh_segment(line_sfc_u)

  minimal_shift <- function(a, b) {
    minimal_shift <- NA
    minimum <- Inf
    n <- nrow(b)
    for (shift in seq_len(n) - 1L) {
      val <- max(abs(a - b[(shift + seq_len(n) - 1L) %% n + 1L, ]))
      if (val < minimum) {
        minimal_shift <- shift
        minimum <- val
      }
    }
    b[(minimal_shift + seq_len(n) - 1L) %% n + 1L, ]
  }
  expect_equal(nrow(seg_sf_u$loc), 4 * 4)
  B <- minimal_shift(seg_sf_u$loc[1:4, ], seg_1$loc[1:4, ])
  expect_identical(seg_sf_u$loc[1:4, ], B)
  B <- minimal_shift(seg_sf_u$loc[4 + 1:4, ], seg_1$loc[4 + 1:4, ])
  expect_identical(seg_sf_u$loc[4 + 1:4, ], B)
  B <- minimal_shift(seg_sf_u$loc[8 + 1:4, ], seg_1$loc[8 + 1:4, ])
  expect_identical(seg_sf_u$loc[8 + 1:4, ], B)
  B <- minimal_shift(seg_sf_u$loc[12 + 1:4, ], seg_1$loc[12 + 1:4, ])
  expect_identical(seg_sf_u$loc[12 + 1:4, ], B)
})




test_that("Conversion from sfc_GEOMETRY to inla.mesh.segment", {
  local_bru_safe_inla()

  ## scf_GEOMETRY ##

  pts0 <- rbind(c(-7, -7), c(7, -7), c(7, 7), c(-7, 7), c(-7, -7)) # covering (CCW)
  pts0b <- pts0 + 15
  pts1 <- rbind(c(0, 3), c(0, 4), c(1, 5), c(2, 5), c(0, 3)) # hole (CW)
  pts2 <- rbind(c(1, 2), c(0, 0), c(0, -1), c(-2, -2), c(1, 2)) # hole (CW)
  seg0 <- INLA::inla.mesh.segment(
    loc = pts0[1:4, , drop = FALSE],
    is.bnd = TRUE,
    crs = fm_CRS()
  )
  seg0b <- INLA::inla.mesh.segment(
    loc = pts0b[1:4, , drop = FALSE],
    is.bnd = TRUE,
    crs = fm_CRS()
  )

  seg1 <- INLA::inla.mesh.segment(
    loc = pts1[1:4, , drop = FALSE],
    is.bnd = TRUE,
    crs = fm_CRS()
  )

  seg2 <- INLA::inla.mesh.segment(
    loc = pts2[1:4, , drop = FALSE],
    is.bnd = TRUE,
    crs = fm_CRS()
  )

  seg_1 <- fm_internal_sp2segment_join(list(seg0, seg1, seg2, seg0b, seg0b),
    grp = c(1, 1, 1, 1, 2)
  )
  expect_identical(seg_1$grp, as.matrix(rep(c(1L, 2L), times = c(16, 4))))

  line_str1 <- sf::st_polygon(list(pts0, pts1, pts2))
  line_str2 <- sf::st_polygon(list(pts0b))

  line_sfc1 <- sf::st_combine(sf::st_sfc(list(line_str1, line_str2), check_ring_dir = TRUE))
  line_sfc2 <- sf::st_sfc(list(line_str2), check_ring_dir = TRUE)
  line_sfc <- c(line_sfc1, line_sfc2)

  expect_s3_class(line_sfc, "sfc_GEOMETRY")
  seg_sf <- fm_as_inla_mesh_segment(line_sfc)

  expect_identical(seg_sf, seg_1)
})
