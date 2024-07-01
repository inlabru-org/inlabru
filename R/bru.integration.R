#' @describeIn inlabru-deprecated
#' `r lifecycle::badge("deprecated")` in favour of [fmesher::fm_int()]
#'
#' @param samplers Description of the integration region boundary.
#' In 1D, a length 2 vector or two-column matrix where each row describes an interval,
#' or `NULL`
#' In 2D either a `SpatialPolygon` or a `SpatialLinesDataFrame` with a `weight` column
#' defining the width of the a transect line, and optionally further columns used by the
#' `group` argument, or `NULL`.  When `domain` is `NULL`, `samplers` may also
#' be an `inla.mesh.1d` or `inla.mesh` object, that is then treated as a `domain`
#' argument instead.
#' @param domain Either
#' * when `samplers` is a 1D interval(s) definition only, `domain` can be
#'   a single integer for the number of integration points to place in each 1D
#'   interval, overriding `int.args[["nsub1"]]`, and otherwise
#' * when `samplers` is `NULL`, `domain` can be a numeric vector of points,
#'   each given integration weight 1 (and no additional points are added
#'   in between),
#' * an `inla.mesh.1d` object for continuous 1D integration, or
#' * an `inla.mesh.2d` object for continuous 2D integration.
#' @param name Character array stating the name of the domains dimension(s).
#' If `NULL`, the names are taken from coordinate names from `samplers` for
#' `Spatial*` objects, otherwise "x", "y", "z" for 2D regions and
#' `"x"` for 1D regions
#' @param group Column names of the `samplers` object (if applicable) for which
#' the integration points are calculated independently and not merged when
#' aggregating to mesh nodes.
#' @param int.args List of arguments passed to `bru_int_polygon`.
#' * `method`: "stable" (to aggregate integration weights onto mesh nodes)
#'   or "direct" (to construct a within triangle/segment integration scheme
#'   without aggregating onto mesh nodes)
#' * `nsub1`, `nsub2`: integers controlling the number of internal integration
#'   points before aggregation. Points per triangle: `(nsub2+1)^2`.
#'   Points per knot segment: `nsub1`
#' * `poly_method`: if set to "legacy", selects an old polygon integration method
#'   that doesn't handle holes. No longer supported, and will generate an error.
#' @param project `r lifecycle::badge("deprecated")` Deprecated in favour of `int.args(method=...)`.
#' If TRUE, aggregate the integration points to mesh vertices. Default:
#' `project = (identical(int.args$method, "stable"))`
#'
#' @return `ipoints()`: A `data.frame`, `tibble`, `sf`, or `SpatialPointsDataFrame` of 1D and
#' 2D integration points, including a `weight` column and `.block` column.
#'
#' @importFrom sp coordnames coordinates
ipoints <- function(samplers = NULL, domain = NULL, name = NULL, group = NULL,
                    int.args = NULL,
                    project = deprecated()) {
  lifecycle::deprecate_stop(
    "2.8.0.9004",
    "ipoints()",
    "fmesher::fm_int()",
    details = c(
      "`ipoints(samplers, domain)` has been replaced by more versatile `fm_int(domain, samplers, ...)` methods."
    )
  )
}




#' @describeIn inlabru-deprecated
#' (Blockwise) cross product of integration points.
#'
#' Calculates the groupwise cross product of integration points in different
#' dimensions and multiplies their weights accordingly.
#' If the object defining points in a particular dimension has no
#' weights attached to it all weights are assumed to be 1.
#'
#' Legacy wrapper for [fm_cprod()]
#'
#' @seealso [fm_cprod()]
#' @keywords internal
#'
#'
#' @param na.rm logical; if `TRUE`, the rows with weight `NA` from the
#' non-overlapping full_join will be removed; if `FALSE`, set the undefined weights to `NA`.
#' If `NULL` (default), act as `TRUE`, but warn if any elements needed removing.
#' @param .blockwise logical; if `FALSE`, computes full tensor product integration.
#' If `TRUE`, computes within-block tensor product integration (used internally
#' by [fm_int()]).
#' Default `FALSE`
#' @return A `data.frame`, `sf`, or `SpatialPointsDataFrame` of multidimensional
#' integration points and their weights
cprod <- function(..., na.rm = NULL, .blockwise = FALSE) {
  lifecycle::deprecate_warn(
    "2.8.0",
    "cprod()",
    "fmesher::fm_cprod()"
  )
  fm_cprod(..., na.rm = na.rm, .blockwise = .blockwise)
}
