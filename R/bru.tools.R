#' Project data to mesh vertices under the assumption of lineariity
#'   
#'
#' @aliases vertex.projections
#' @export
#' @param mesh An inla.mesh object
#' @param points A SpatialPointsDataFrame object
#' @param columns A character array of the points columns which whall be projected
#' @return SpatialPointsDataFrame of mesh vertices with projected data attached

vertex.projection = function(mesh, points, columns = names(points)){
  
  res = inla.fmesher.smorg(mesh$loc, mesh$graph$tv, points2mesh = coordinates(points))
  tri = res$p2m.t 
  
  data = list()
  for (k in 1:length(columns)){
    cn = columns[k]
    nw = points@data[,columns] * res$p2m.b
    w.by = by(as.vector(nw), as.vector(mesh$graph$tv[tri,]), sum, simplify = TRUE)
    data[[cn]] = as.vector(w.by) 
  }
  
  data = data.frame(data)
  coords = mesh$loc[as.numeric(names(w.by)),c(1,2)]
  data$vertex = as.numeric(names(w.by))
  ret = SpatialPointsDataFrame(coords, proj4string = CRS(proj4string(points)), data = data)
  return(ret)
}
