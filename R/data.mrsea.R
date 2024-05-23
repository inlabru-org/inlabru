#' @name mrsea
#' @title Marine renewables strategic environmental assessment
#' @docType data
#' @description Data imported from package MRSea, see <https://www.creem.st-andrews.ac.uk/software/>
#'
#' @format A list of objects:
#'  \describe{
#'    \item{`points`:}{ A `SpatialPointsDataFrame` object containing the locations of
#'    XXXXX.}
#'    \item{`samplers`:}{ A `SpatialLinesDataFrame` object containing the transect lines
#'    that were surveyed.}
#'    \item{`mesh`:}{ An `fm_mesh_2d` object containing a Delaunay triangulation
#'    mesh (a type of discretization of continuous space) covering the survey region.}
#'    \item{`boundary`:}{ An `SpatialPolygonsDataFrame` object defining the boundary of the
#'    survey region.}
#'    \item{`covar`:}{ An `SpatialPointsDataFrame` containing sea depth estimates.}
#'  }
#' @source
#' Library `MRSea`.
#'
#' @references
#'
#' NONE YET
#'
#' @examples
#' if (bru_safe_inla() &&
#'   require(ggplot2, quietly = TRUE)) {
#'   ggplot() +
#'     gg(mrsea$mesh) +
#'     gg(mrsea$samplers) +
#'     gg(mrsea$points) +
#'     gg(mrsea$boundary)
#' }
"mrsea"
