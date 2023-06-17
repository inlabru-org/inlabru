# GENERICS
# refine = function(...) {UseMethod("refine")}
# tsplit = function(...) {UseMethod("tsplit")}


#' Query if a point is inside the mesh boundary
#'
#'
#' @aliases is.inside
#' @export
#' @keywords internal
#' @param mesh an inla.mesh object.
#' @param loc Points in space stored either as data.frame, a two-column matrix of x and y coordinates or a SpatialPoints object.
#' @param mesh.coords Coordinate names of the mesh. Use only if loc is a data.frame with respective column names.
#' @return Single column matrix of Boolean values indicating if a point is inside the mesh.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'
#' @examples
#' \dontrun{
#' if (bru_safe_inla(quietly = TRUE)) {
#'   # Load Gorilla data
#'
#'   data("gorillas", package = "inlabru")
#'
#'   # Check if all Gorilla nests are inside the mesh
#'
#'   all(is.inside(gorillas$mesh, gorillas$nests))
#'
#'   # Also works for locations not stored as SpatialPoints object
#'
#'   loc <- coordinates(gorillas$nests)
#'   all(is.inside(gorillas$mesh, loc))
#' }
#' }
#'
is.inside <- function(mesh, loc, mesh.coords = NULL) {
  lifecycle::deprecate_soft(
    "2.8.0",
    "is.inside()",
    "fm_is_within()",
    details = "`is.inside(mesh, loc)` becomes `fm_is_within(loc, mesh)`"
  )

  if (!is.null(mesh.coords)) {
    loc <- as.matrix(loc[, mesh.coords])
  }
  fm_is_within(loc, mesh)
}





#' @describeIn inlabru-deprecated Extract vertex locations from an `inla.mesh`.
#' Converts the vertices of an `inla.mesh` object into a `SpatialPointsDataFrame`.
#' Deprecated in favour of [fm_vertices()]
#'
#' @export
vertices.inla.mesh <- function(...) {
  lifecycle::deprecate_warn(
    "2.8.0",
    "vertices.inla.mesh()",
    "fm_vertices()",
    details =
      c(
        "!" = "fm_vertices() by default returns 'sf' instead of SPDF.",
        "!" = "fm_vertices() includes a '.vertex' column instead of a 'vertex' column."
      )
  )

  vrt <- fm_vertices(..., format = "sp") %>%
    dplyr::rename(vertex = .data$.vertex)

  vrt
}


#' @title Generate `SpatialPixels` covering an `inla.mesh`
#'
#' @description Generate `SpatialPixels` covering an `inla.mesh`
#'
#' @export
#' @seealso [fm_pixels()]
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'
#' @param mesh An `inla.mesh` object
#' @param nx Number of pixels in x direction
#' @param ny Number of pixels in y direction
#' @param mask If logical and TRUE, remove pixels that are outside the mesh.
#' If `mask` is a `Spatial` object, only return pixels covered by this object.
#' @return `SpatialPixelsDataFrame` covering the mesh
#'
#' @examples
#' \donttest{
#' if (require(ggplot2, quietly = TRUE)) {
#'   data("mrsea", package = "inlabru")
#'   pxl <- fm_pixels(mrsea$mesh,
#'     nx = 50, ny = 50, mask = mrsea$boundary,
#'     format = "sp"
#'   )
#'   ggplot() +
#'     gg(pxl, fill = "grey", alpha = 0.5) +
#'     gg(mrsea$mesh)
#' }
#' }
#'
pixels <- function(mesh, nx = 150, ny = 150, mask = TRUE) {
  lifecycle::deprecate_soft(
    "2.8.0",
    "pixels()",
    "fm_pixels(format = 'sp')",
    details = "The fm_pixels() function can generate sf, terra, and sp output."
  )
  fm_pixels(mesh, nx = nx, ny = ny, mask = mask, format = "sp")
}





#' Refine an inla.mesh object
#'
#'
#' @aliases refine.inla.mesh
#' @keywords internal
#'
#' @param mesh an inla.mesh object
#' @param refine A list of refinement options passed on to
#' `INLA::inla.mesh.create`
#' @return mesh A refined inla.mesh object
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}

refine.inla.mesh <- function(mesh, refine = list(max.edge = 1)) {
  lifecycle::deprecate_warn(
    "2.7.0",
    "refine.inla.mesh()",
    details = c(
      "!" = "This function is experimental and will be replaced by a new method"
    )
  )

  rmesh <- INLA::inla.mesh.create(
    loc = mesh$loc,
    interior = INLA::inla.mesh.interior(mesh),
    boundary = INLA::inla.mesh.boundary(mesh),
    refine = refine
  )
  return(rmesh)
}

#' Split triangles of a mesh into four triangles
#'
#' * Warning: The original triangle edges are not always preserved.
#' * Warning: Works in euclidean coordinates. Not suitable for sphere.
#' * Warning: Experimental; will be replaced by a new method
#'
#' @aliases tsplit.inla.mesh
#' @keywords internal
#' @param mesh an inla.mesh object
#' @param n number of splitting recursions
#' @return mesh A refined inla.mesh object
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}

tsplit.inla.mesh <- function(mesh, n = 1) {
  lifecycle::deprecate_warn(
    "2.7.0",
    "tsplit.inla.mesh()",
    details = c(
      "!" = "This function is experimental and will be replaced by a new method."
    )
  )

  p1 <- mesh$loc[mesh$graph$tv[, 1], ]
  p2 <- mesh$loc[mesh$graph$tv[, 2], ]
  p3 <- mesh$loc[mesh$graph$tv[, 3], ]

  m1 <- (p2 + p3) / 2
  m2 <- (p1 + p3) / 2
  m3 <- (p1 + p2) / 2
  tri.loc <- rbind(mesh$loc, m1, m2, m3)

  split.edges <- function(segm) {
    if (is.null(segm) || (nrow(segm$idx) == 0)) {
      return(segm)
    }
    n.loc <- nrow(segm$loc)
    n.idx <- nrow(segm$idx)
    loc <- rbind(
      segm$loc,
      (segm$loc[segm$idx[, 1], ] +
        segm$loc[segm$idx[, 2], ]) / 2
    )
    idx <- rbind(
      cbind(segm$idx[, 1], n.loc + seq_len(n.idx)),
      cbind(n.loc + seq_len(n.idx), segm$idx[, 2])
    )

    segm2 <-
      INLA::inla.mesh.segment(
        loc = loc,
        idx = idx,
        grp = c(segm$grp, segm$grp),
        is.bnd = segm$is.bnd,
        crs = fm_CRS(segm)
      )

    segm2
  }

  interior2 <- split.edges(INLA::inla.mesh.interior(mesh)[[1]])
  boundary2 <- split.edges(INLA::inla.mesh.boundary(mesh)[[1]])

  mesh2 <- INLA::inla.mesh.create(
    loc = tri.loc,
    interior = interior2,
    boundary = boundary2,
    crs = fm_CRS(mesh)
  )

  if (n <= 1) {
    return(mesh2)
  } else {
    return(tsplit.inla.mesh(mesh2, n - 1))
  }
}
