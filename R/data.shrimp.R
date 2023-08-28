#' @name shrimp
#' @title Blue and red shrimp in the Western Mediterranean Sea
#' @docType data
#' @description Blue and red shrimp in the Western Mediterranean Sea.
#'
#' @usage data(shrimp)
#'
#' @format A list of objects:
#'  \describe{
#'    \item{`haul`:}{ A `SpatialPointsDataFrame` object containing haul locations}
#'    \item{`mesh`:}{ An `inla.mesh` object containing a Delaunay triangulation
#'    mesh (a type of discretization of continuous space) covering the haul locations.}
#'    \describe{
#'       \item{`catch`}{Catch in Kg.}
#'       \item{`landing`}{Landing in Kg.}
#'       \item{`depth`}{Mean depth of the fishery haul.}
#'     }
#'  }
#' @source
#' Pennino, Maria Grazia. Personal communication.
#'
#' @references
#' Pennino, M. G., Paradinas, I., Munoz, F., Illian, J.,Quilez-Lopez, A., Bellido, J.M., Conesa,
#' D. Accounting for preferential sampling in species distribution models. Ecology and Evolution,  In Press.
#'
#' @examples
#' \donttest{
#' if (require(ggplot2, quietly = TRUE)) {
#'   data(shrimp, package = "inlabru")
#'   ggplot() +
#'     gg(shrimp$mesh) +
#'     gg(shrimp$hauls) +
#'     coord_equal()
#' }
#' }
"shrimp"
