# Convert data frame to SpatialLinesDataFrame
#   
#
# @aliases sfill
# @export
# @param data A SpatialGridDataFrame or SpatialPixelDataFrame
# @param where Spatial*DataFrame of locations for which to fill in values from \code{data}. If NULL, use \code{data} to determine the locations.
# @return Spatial object

sfill = function(data, where = NULL) {
  
  if ( is.null(where) ) { where = data }
  vallist = list()
  for (k in 1:ncol(data@data)) {
    
    dpoints = SpatialPoints(data)
    vals = data@data[,k]
    dpoints = dpoints[!is.na(vals),]
    vals = vals[!is.na(vals)]
    
    data.ow = spatstat::owin(range(coordinates(dpoints)[,1]), range(coordinates(dpoints)[,2]))
    data.ppp = spatstat::as.ppp(coordinates(dpoints), data.ow)
    where.ow = spatstat::owin(range(coordinates(where)[,1]), range(coordinates(where)[,2]))
    where.ppp = spatstat::as.ppp(coordinates(where), where.ow)
    
    nn = spatstat::nncross(where.ppp, data.ppp)[,"which"]
    vallist[[k]] = vals[nn]
  }
  ret = data.frame(do.call(data.frame, vallist))
  colnames(ret) = colnames(data@data)
  ret = sp::SpatialPixelsDataFrame(where, data = ret)


}


#' Convert data frame to SpatialLinesDataFrame
#'
#' A line in 2D space is defined by a start and an and point, each associated with 2D coordinates. 
#' This function takes a /code{data.frame} as input and assumes that each row defines a line in space.
#' In order to do so, the data frame must have at least four columns and the \code{start.cols} and
#' \code{end.cols} parameters must be used to point out the names of the columns that define
#' the start and end coordinates of the line. The data is then converted to a 
#' \code{SpatialLinesDataFrame} \code{DF}. If a coordinate reference system \code{crs} is provided
#' it is attached to \code{DF}. If also \code{to.crs} is provided, the coordinate system of \code{DF}
#' is transfromed accordingly. Additional columns of the input data, e.g. covariates,
#' are retained and attached to \code{DF}.
#' 
#'
#' @aliases sline
#' @export
#' @param data A data.frame
#' @param start.cols Character array poitning out the columns of \code{data} that hold the start points of the lines
#' @param end.cols Character array poitning out the columns of \code{data} that hold the end points of the lines
#' @param crs Coordinate reference system of the original \code{data}
#' @param to.crs Coordinate reference system for the SpatialLines ouput. 
#' @return SpatialLinesDataFrame 
#' 
#' @examples 
#'
#' \donttest{
#' # Create a data frame defining three lines 
#' lns = data.frame(xs = c(1,2,3), ys = c(1,1,1), # start points
#'                  xe = c(2,3,4), ye = c(2,2,2)) # end points
#' 
#' 
#' # Conversion to SpatialLinesDataFrame without CRS
#' spl = sline(lns, start.cols = c("xs","ys"), 
#'             end.cols = c("xe","ye"))
#' 
#' # Plot the lines
#' ggplot() + gg(spl)
#' }

sline = function(data, start.cols, end.cols, crs = CRS(as.character(NA)), to.crs = NULL) {
  
  sp = as.data.frame(data[,start.cols])
  ep = as.data.frame(data[,end.cols])
  
  colnames(sp) = c("x","y")
  colnames(ep) = c("x","y")
  
  lilist = lapply(1:nrow(sp), function(k) { Lines(list(Line(rbind(sp[k,],ep[k,]))), ID = k) } )
  spl = SpatialLines(lilist, proj4string = crs)
  
  df = data[,setdiff(names(data), c(start.cols, end.cols))]
  rownames(df) = 1:nrow(df)
  
  slines = SpatialLinesDataFrame(spl, data = df)
  
  # If requested, change CRS
  if ( !is.null(to.crs) ) slines = spTransform(slines, to.crs)
  
  slines
}


#' Convert a data.frame of boundary points into a SpatialPolgonsDataFrame
#' 
#' A polygon can be described as a sequence of points defining the polygon's boundary.
#' When given such a sequence (anti clockwise!) this function creates a 
#' SpatialPolygonsDataFrame holding the polygon decribed. By default, the
#' first two columns of \code{data} are assumed to define the x and y coordinates
#' of the points. This behavior can ba changed using the \code{cols} parameter, which
#' points out the names of the columns holding the coordinates. The coordinate
#' reference system of the resulting spatial polygon can be set via the \code{crs}
#' paraemter. Posterior conversion to a different CRS is supported using the
#' \code{to.crs} parameter.
#'
#' @aliases spoly
#' @export
#' @param data A data.frame of points describing the boundary of the polygon
#' @param cols Column names of the x and y coordinates within the data
#' @param crs Coordinate reference system of the points
#' @param to.crs Coordinate reference system for the SpatialLines ouput. 
#' @return SpatialPolygonsDataFrame 
#' 
#' @examples 
#' \donttest{
#' # Create data frame of boundary points (anti clockwise!)
#' pts = data.frame(x = c(1,2,1.7,1.3),
#'                  y = c(1,1,2,2))
#' 
#' # Convert to SpatialPolygonsDataFrame
#' pol = spoly(pts, crs = CRS(as.character(NA)))
#' 
#' # Plot it!
#' ggplot() + gg(pol)
#' }

spoly = function(data, cols = colnames(data)[1:2], crs = CRS(as.character(NA)), to.crs = NULL) {
  
  po = Polygon(data[, cols], hole=FALSE)
  pos = Polygons(list(po), ID = "tmp")
  predpoly = SpatialPolygons(list(pos), proj4string = crs)
  df = data.frame(weight = 1)
  rownames(df) = "tmp"
  spoly = SpatialPolygonsDataFrame(predpoly, data = df)
  
  # If requested, change CRS
  if ( !is.null(to.crs) ) spoly = spTransform(spoly, to.crs)
  spoly
}


#' Coordinate transformation for spatial objects
#' 
#' This is a wrapper for the \link{spTransform} function provided by the \code{sp} package.
#' Given a spatial object (or a list thereof) it will transform the coordinate system according
#' to the parameter \code{crs}. In addition to the usual spatial objects this function is
#' also capables of transforming \link[INLA]{inla.mesh} objects that are equipped with a coordinate
#' system.#'
#' @aliases stransform
#' @export
#' @param splist list of Spatial* objects
#' @param crs Coordinate reference system to change to
#' @return List of Spatial* objects
#' 
#' @examples
#' \donttest{
#' # Load Gorilla data
#' data("gorillas", package = "inlabru")
#' 
#' # Take the mesh and transform it to latitude/longitude
#' tmesh = stransform(gorillas$mesh, crs = CRS("+proj=longlat"))
#'
#' # Compare original and transformed mesh
#' 
#' multiplot(ggplot() + gg(gorillas$mesh) + ggtitle("Original mesh"),
#'           ggplot() + gg(tmesh) + ggtitle("Transformed mesh"))
#' }         

stransform = function(splist, crs) {
  if (!is.null(crs)) {
    if ( class(splist)[[1]] == "list" ) {
      for (k in 1:length(splist)) {
        if (inherits(splist[[k]], "Spatial")) {
          # cn = coordnames(splist[[k]])
          splist[[k]] = sp::spTransform(splist[[k]], crs)
          #coordnames(splist[[k]]) = cn
        } else if (inherits(splist[[k]], "inla.mesh")) {
          splist[[k]] = INLA::inla.spTransform(splist[[k]], CRSobj = crs)
        }
      }
    } else { 
      splist = stransform(list(splist), crs = crs)[[1]]
    } 
    splist
  } else {splist}
}
