#' Convert SpatialPoints and boundary polygon to spatstat ppp object
#' 
#' Spatstat point pattern objects consist of points and an observation windows. This
#' function uses a SpatialPoints object and a SpatialPolygon object to generate the points
#' and the window. Lastly, the ppp() function is called to create the \code{ppp} object.
#'
#' @aliases spatial.to.ppp
#' @export
#' @param points A SpatialPoints[DataFrame] object describing the point pattern.
#' @param samplers A SpatialPolygons[DataFrame] object describing the observation window.
#' @return A spatstat \code{spatstat} \code{ppp} object
#' 
#' @examples 
#' 
#' \donttest{
#' # Load Gorilla data
#' 
#' data("gorillas", package = "inlabru")
#' 
#' # Use nest locations and survey boundary to create a spatstat ppp object
#' 
#' gp <- spatial.to.ppp(gorillas$nests, gorillas$boundary)
#' class(gp)
#' 
#' # Plot it
#' 
#' plot(gp)
#' }


spatial.to.ppp = function(points, samplers) {

  bnd = samplers@polygons[[1]]@Polygons[[1]]@coords
  bnd = bnd[1:(nrow(bnd)-1),]
  gp = spatstat::ppp(x = coordinates(points)[,1], 
           y = coordinates(points)[,2], 
           window = spatstat::owin(poly = list(x=rev(bnd[,1]), y = rev(bnd[,2]))) )
  
}