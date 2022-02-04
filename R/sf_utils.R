# helper functions for working with sf objects

#'Calculate signed area for polygon
#'
#'@aliases st_signed_area
#'@export
#'@param sfg A POLYGON sfg object
#'@return Returns the signed area.  Negative values indicate
#'anti-clockwise winding direction.
#'@author Andrew Seaton \email{Andrew.Seaton.2@@glasgow.ac.uk}
#'@author Finn Lindgren \email{finn.lindgren@@gmail.com}

st_signed_area = function(sfg){

  if (!inherits(sfg, c("POLYGON", "sfg"))){
    stop("Signed area only implemented for POLYGON sfg objects")
  }

  coords <- as.matrix(sfg)
  i <- seq_len(nrow(coords)-1)
  edges <- cbind(coords[i,,drop=FALSE], coords[i+1,,drop=FALSE])
  area <- sum((edges[,3] - edges[,1]) * (edges[,2] + edges[,4]) / 2)
  return(area)

}

#' Check for "XYZ", "XYM" and "XYZM" sfg classes
#'
#' @aliases sf_check_dim
#' @param sfc An sfc object
#' @return LOGICAL indicating if any sfg element of the sfc object has class "XYZ", "XYM" or "XYZM". Internal function used to check for 3 and 4 dimensional objects.

sf_check_dim = function(sfc){

  check = sapply(sfc,
                 FUN = function(x) inherits(x, c("XYZ", "XYM", "XYZM")))

  sum(check)
}

