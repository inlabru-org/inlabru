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
