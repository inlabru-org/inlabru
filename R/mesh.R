# GENERICS
# refine = function(...) {UseMethod("refine")}
# tsplit = function(...) {UseMethod("tsplit")}


#' @describeIn inlabru-deprecated
#' Find out which points are inside a mesh.
#' `r lifecycle::badge("deprecated")` in favour of [fm_is_within()].
#' Replace `is.inside(mesh, loc)` with `fm_is_within(loc, mesh)`.
#'
#' @seealso [fm_is_within()]
#' @export
#' @param mesh an inla.mesh object.
#' @param loc Points in space stored either as data.frame, a two-column matrix
#'   of x and y coordinates or a SpatialPoints object.
#' @param mesh.coords Coordinate names of the mesh. Use only if loc is a
#'   data.frame with respective column names.
#' @return `is.inside()`: Single column matrix of Boolean values indicating if a point is
#'   inside the mesh.
#'
is.inside <- function(mesh, loc, mesh.coords = NULL) {
  lifecycle::deprecate_warn(
    "2.8.0",
    "is.inside()",
    "fmesher::fm_is_within()",
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
#' @export vertices.inla.mesh
vertices.inla.mesh <- function(...) {
  lifecycle::deprecate_warn(
    "2.8.0",
    "vertices.inla.mesh()",
    "fmesher::fm_vertices()",
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
#' @description
#' `r lifecycle::badge("deprecated")` in favour of [fmesher::fm_pixels()]
#'
#' Generate `SpatialPixels` covering an `inla.mesh`.
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
#'   pxl <- fm_pixels(
#'     mrsea$mesh,
#'     dims = c(50, 50),
#'     mask = mrsea$boundary,
#'     format = "sp",
#'     minimal = TRUE
#'   )
#'   ggplot() +
#'     gg(pxl, fill = "blue", alpha = 0.75) +
#'     gg(mrsea$mesh)
#'
#'   pxl <- fm_pixels(
#'     mrsea$mesh,
#'     dims = c(50, 50),
#'     mask = mrsea$boundary,
#'     format = "sf",
#'     minimal = TRUE
#'   )
#'   ggplot() +
#'     gg(pxl, geom = "tile", fill = "blue", alpha = 0.75) +
#'     gg(mrsea$mesh)
#' }
#' }
#'
pixels <- function(mesh, nx = 150, ny = 150, mask = TRUE) {
  lifecycle::deprecate_warn(
    "2.8.0",
    "pixels()",
    "fmesher::fm_pixels(format = 'sp')",
    details = "The `fm_pixels()` function can generate `sf`, `terra`, and `sp` output."
  )
  fm_pixels(mesh, dims = c(nx, ny), mask = mask, format = "sp", minimal = FALSE)
}
