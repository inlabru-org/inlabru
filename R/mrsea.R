#' @name mrsea
#' @title Marine renewables strategic environmental assessment
#' @docType data
#' @description Data imported from package MRSea, see http://creem2.st-andrews.ac.uk/software/
#' 
#' @usage mrsea
#' 
#' @format A list of objects:
#'  \describe{
#'    \item{\code{points}:}{ A \code{SpatialPointsDataFrame} object containing the locations of 
#'    XXXXX.}
#'    \item{\code{samplers}:}{ A \code{SpatialLinesDataFrame} object containing the transect lines
#'    that were surveyed.}
#'    \item{\code{mesh}:}{ An \code{inla.mesh} object containing a Delaunay triangulation 
#'    mesh (a type of discretization of continuous space) covering the survey region.}
#'    \item{\code{boundary}:}{ An \code{SpatialPolygonsDataFrame} object defining the boundary of the 
#'    survey region.}
#'    \item{\code{covar}:}{ An \code{SpatialPointsDataFrame} containing sea depth estimates.}
#'  }
#' @source 
#' Library \code{MRSea}.
#' 
#' @references 
#' 
#' NONE YET
#' 
#' @examples
#'  data(mrsea)
#'  ggplot() + gg(mrsea$mesh) + gg(mrsea$samplers) + gg(mrsea$points) + gg(mrsea$boundary)
NULL
