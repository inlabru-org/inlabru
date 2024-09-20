#' @name mrsea
#' @title Marine renewables strategic environmental assessment
#' @docType data
#' @description Data imported from package MRSea, see
#'   <https://www.creem.st-andrews.ac.uk/software/>
#'
#' @format A list of objects:
#'  \describe{
#'    \item{`points`}{ A `sf` object containing the locations of
#'    XXXXX.}
#'    \item{`samplers`}{ A `sf` object containing the transect lines
#'    that were surveyed.}
#'    \item{`mesh`}{ An `fm_mesh_2d` object containing a Delaunay triangulation
#'    mesh (a type of discretization of continuous space) covering the survey
#'      region.}
#'    \item{`boundary`}{ An `sf` object defining the boundary polygon of the
#'    survey region.}
#'    \item{`covar`}{ An `sf` containing sea depth estimates.}
#'  }
#' @source
#' Library `MRSea`.
#'
#' @references
#'
#' NONE YET
#'
#' @examples
#' if (require(ggplot2, quietly = TRUE)) {
#'   ggplot() +
#'     geom_fm(data = mrsea$mesh) +
#'     gg(mrsea$samplers) +
#'     gg(mrsea$points) +
#'     gg(mrsea$boundary)
#' }
"mrsea"
