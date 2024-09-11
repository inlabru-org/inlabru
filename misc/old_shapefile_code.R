# Construct a covariate (function) returning region identifiers of locations
#
# @export
# @param shapefile Either a character array identifying the shapefile or a shapefile
# @param coords Coordinates to to be extracted from the location data provided to the constructed covariate
# @examples \dontrun{ }
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#

shapefile.to.covariate <- function(shapefile, coords = c("x", "y")) {
  if (is.character(shapefile)) {
    shapefile <- maptools::readShapeSpatial(shapefile)
  }

  fun <- function(loc) {
    sp::over(sp::SpatialPoints(as.data.frame(loc[, coords])), shapefile, fn = NULL)[, 1]
  }
}
