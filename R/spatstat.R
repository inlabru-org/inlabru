#' Convert SpatialPoints and boundary polygon to  spatstat ppp object
#'
#' @aliases spatial.to.ppp
#' @export
#' @param points A SpatialPoints[DataFrame] describing the point pattern
#' @param samplers A SpatialPolygons[DataFrame] describing the observation window
#' @param ... arguments passed on to \link{predict}
#' @return A spatstat \code{spatstat} \code{ppp} object
#' 
#' @examples 
#' 
#' # Load Gorilla data
#' 
#' data("gorillas")
#' 
#' # Use nest locations and survey boundary to create a spatstat ppp object
#' 
#' gp <- spatial.to.ppp(gorillas$nests, gorillas$boundary)
#' class(gp)
#' 
#' # Plot it
#' 
#' plot(gp)


spatial.to.ppp = function(points, samplers) {

  bnd = samplers@polygons[[1]]@Polygons[[1]]@coords
  bnd = bnd[1:(nrow(bnd)-1),]
  gp = spatstat::ppp(x = coordinates(points)[,1], 
           y = coordinates(points)[,2], 
           window = spatstat::owin(poly = list(x=rev(bnd[,1]), y = rev(bnd[,2]))) )
  
}