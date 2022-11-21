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
NULL

#' Shrimp data import
#'
#'
#' Load `shrimp` data stored as file `gamba.Rdata` and construct spatial object
#'
#' @aliases import.shrimp
#' @keywords internal
#' @return The [shrimp] data set
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'

import.shrimp <- function() {
  # Load the raw data
  load(file = paste0(system.file("inst/extdata", package = "inlabru"), "/gamba.Rdata"))

  # Remove ambiguous coordinates
  gamba <- gamba[, c(1, 2, 3, 6, 7)]

  # Rename columns
  colnames(gamba) <- c("catch", "landing", "depth", "lat", "lon")

  # Turn into spatial object
  coordinates(gamba) <- c("lon", "lat")
  proj4string(gamba) <- CRS("+proj=longlat")

  # Make a mesh
  mesh <- INLA::inla.mesh.2d(loc.domain = gamba, max.edge = 0.15, offse = 0.15)

  # Final shrimp object
  shrimp <- list(hauls = gamba, mesh = mesh)

  return(shrimp)
}

# shrimp <- import.shrimp()
# use_data(shrimp)
