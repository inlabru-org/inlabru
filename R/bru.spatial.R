#' Convert data frame to SpatialLinesDataFrame
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

sline = function(data, start.cols, end.cols, crs, to.crs = NULL) {
  
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


#' Create a SpatialPolygonsDataFrame from a boundary description
#'
#' @aliases spoly
#' @export
#' @param points A matrix or data.frame of points describing the boundary of the polygon
#' @param crs Coordinate reference system of the points
#' @param to.crs Coordinate reference system for the SpatialLines ouput. 
#' @return SpatialPolygonsDataFrame 

spoly = function(points, crs, to.crs = NULL) {
  
  po = Polygon(points, hole=FALSE)
  pos = Polygons(list(po), ID = "tmp")
  predpoly = SpatialPolygons(list(pos), proj4string = crs)
  df = data.frame(weight = 1)
  rownames(df) = "tmp"
  spoly = SpatialPolygonsDataFrame(predpoly, data = df)
  
  # If requested, change CRS
  if ( !is.null(to.crs) ) spoly = spTransform(spoly, to.crs)
  spoly
}


#' Apply spTransform to all Spatial* elements of a list 
#'
#' @aliases stransform
#' @export
#' @param splist list of Spatial* objects
#' @param crs Coordinate reference system to change to
#' @return List of Spatial* objects

stransform = function(splist, crs) {
  for (k in 1:length(splist)) {
    if (inherits(splist[[k]], "Spatial")) {
      splist[[k]] = spTransform(splist[[k]], crs)
    } else if (inherits(splist[[k]], "inla.mesh")) {
      splist[[k]] = inla.spTransform(splist[[k]], CRSobj = crs)
    }
  }
  splist
}
