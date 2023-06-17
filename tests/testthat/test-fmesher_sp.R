# FL note: see the inla.mesh.segment documentation for the idx and is.bnd arguments
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

test_that("Conversion from matrix to inla.mesh.segment", {
  local_bru_safe_inla()

  ## sfc_POINT ##

  # compare inla.mesh.segment with matrix input
  # to fm_as_inla_mesh_segment with sf point input

  # matrix version
  loc.bnd <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
  segm.bnd <- INLA::inla.mesh.segment(loc.bnd,
    is.bnd = TRUE
  )

  segm.bnd.sp <- fm_as_inla_mesh_segment(loc.bnd, is.bnd = TRUE, closed = TRUE)

  expect_identical(segm.bnd.sp, segm.bnd)
})


test_that("Conversion from Lines to inla.mesh.segment", {
  local_bru_safe_inla()

  ## sp::Lines ##

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
    grp = seq_len(2)
  )
  expect_identical(seg$grp, as.matrix(rep(1:2, each = 3)))

  seg_sp <- fm_as_inla_mesh_segment(
    sp::Lines(list(sp::Line(pts1), sp::Line(pts2)), ID = "A"),
    grp = 1:2
  )

  expect_identical(seg_sp, seg)
})


test_that("Conversion from Polygons to inla.mesh.segment", {
  extract_sequences <- function(seg) {
    sequences <- list()
    unused_edges <- rep(TRUE, nrow(seg$idx))
    i <- integer(0)
    while (any(unused_edges)) {
      edge <- min(which(unused_edges))
      i <- seg$idx[edge, 1]
      while (length(edge) > 0) {
        edge <- min(edge)
        i <- c(i, seg$idx[edge, 2])
        unused_edges[edge] <- FALSE
        edge <- which(unused_edges & (seg$idx[, 1] == i[length(i)]))
      }
      sequences[[length(sequences) + 1]] <- i
      i <- integer(0)
    }
    sequences
  }

  local_bru_safe_inla()

  ## Polygon ##
  pts1 <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)) # solid
  pts2 <- rbind(c(0, 0), c(0, 1), c(1, 1), c(1, 0), c(0, 0)) # hole
  seg1 <- INLA::inla.mesh.segment(
    loc = pts1[1:4, , drop = FALSE],
    is.bnd = TRUE
  )

  seg2 <- INLA::inla.mesh.segment(
    loc = pts2[1:4, , drop = FALSE],
    is.bnd = TRUE
  )

  poly1 <- sp::Polygon(pts1[5:1, ], hole = FALSE)
  poly2 <- sp::Polygon(pts2[5:1, ], hole = TRUE)
  seg1_sp <- fm_as_inla_mesh_segment(poly1)
  seg2_sp <- fm_as_inla_mesh_segment(poly2)

  expect_identical(seg1_sp$loc[seg1_sp$idx[, 1], ], seg1$loc[seg1$idx[, 1], ])
  expect_identical(seg2_sp$loc[seg2_sp$idx[, 1], ], seg2$loc[seg2$idx[, 1], ])

  seq_seg1 <- extract_sequences(seg1)
  seq_seg1_sp <- extract_sequences(seg1_sp)
  expect_identical(
    seg1_sp$loc[seq_seg1_sp[[1]], ],
    seg1$loc[seq_seg1[[1]], ]
  )
  seq_seg2 <- extract_sequences(seg2)
  seq_seg2_sp <- extract_sequences(seg2_sp)
  expect_identical(
    seg2_sp$loc[seq_seg2_sp[[1]], ],
    seg2$loc[seq_seg2[[1]], ]
  )

  ## Polygons ##

  # Winding order and hold status is messy for sp.
  # Should focus on the sf conversions instead.
  if (FALSE) {
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
      grp = seq_len(2)
    )
    expect_identical(seg$grp, as.matrix(rep(1:2, each = 4)))

    poly_sp <- sp::Polygons(list(
      sp::Polygon(pts1, hole = TRUE),
      sp::Polygon(pts2, hole = FALSE)
    ), ID = "A")
    seg_sp <- fm_as_inla_mesh_segment(
      poly_sp,
      grp = 1:2
    )

    seq_seg <- extract_sequences(seg)
    seq_seg_sp <- extract_sequences(seg_sp)
    # Matches:
    expect_identical(
      seg_sp$loc[seq_seg_sp[[1]], ],
      seg$loc[seq_seg[[1]], ]
    )
    # Doesn't match:
    expect_identical(
      seg_sp$loc[seq_seg_sp[[2]], ],
      seg$loc[seq_seg[[2]], ]
    )
  }
})



# if (FALSE) {
#  # Testing for fmesher functions using sf objects.
#  # For testing while editing code.  Eventually these will be the
#  # basis for proper package tests.  For now should work
#  # standalone by sourcing the relevant fmesher_*.R scripts directly
#
#  library(inlabru)
#  library(INLA)
#  library(sf)
#  source(here::here("R", "fmesher_crs.R"))
#  source(here::here("R", "sf_utils.R"))
#  gorillas_sf <- readRDS(here::here("sf", "Data", "gorillas_sf.rds"))
#
# }
