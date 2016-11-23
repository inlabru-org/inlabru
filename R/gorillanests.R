#' @name gorillanests
#' @title Gorilla Nesting Sites.
#' @docType data
#' @description This a subset of the \code{gorillas} dataset from the package \code{spatstat}, reformatted
#' as point process data for use with \code{inlabru}. 
#' 
#' @usage data(gorillanests)
#' 
#' @format The data contain these objects:
#'  \describe{
#'    \item{\code{gnests}:}{ A \code{SpatialPointsDataFrame} object containing the locations of 
#'    the gorilla nests.}
#'    \item{\code{gnestboundary}:}{ An \code{SpatialPolygonsDataFrame} object defining the boundary
#'    of the regoin that was searched for the nests.}
#'    \item{\code{gnestmesh}:}{ An \code{inla.mesh} object containing a mesh that can be used
#'    with function \code{lgcp} to fit a LGCP to the nest data.}
#'  }
#' @source 
#' Library \code{spatstat}.
#' 
#' @references 
#' Funwi-Gabga, N. (2008) A pastoralist survey and fire impact assessment in the Kagwene Gorilla 
#' Sanctuary, Cameroon. M.Sc. thesis, Geology and Environmental Science, University of Buea, 
#' Cameroon.
#' 
#' Funwi-Gabga, N. and Mateu, J. (2012) Understanding the nesting spatial behaviour of gorillas 
#' in the Kagwene Sanctuary, Cameroon. Stochastic Environmental Research and Risk Assessment 
#' 26 (6), 793–811.
#' 
#' @examples
#'  data(gorillanests)
#'  ggplot()+gg(gnestmesh,lwd=0.1)+gg(gnestboundary)+gg(gnests,pch="+",cex=2)
#'  
NULL