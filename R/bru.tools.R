# Project data to mesh vertices under the assumption of lineariity
#   
#
# @aliases vertex.projection
# @export
# @param points A SpatialPointsDataFrame object
# @param mesh An inla.mesh object
# @param columns A character array of the points columns which whall be projected
# @param group Character array identifying columns in \code{points}. These coloumns are interpreted as factors and the projection is performed independently for eah combination of factor levels.
# @return SpatialPointsDataFrame of mesh vertices with projected data attached

vertex.projection = function(points, mesh, columns = names(points), group = NULL){
  
  if ( is.null(group) ) {
    
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
    coordnames(ret) = coordnames(points)
    
  } else {
    fn = function(X) {
      coordinates(X) = coordnames(points)
      ret = vertex.projection(X, mesh, columns = columns)
      for (g in group) { ret[[g]] = X[[g]][1] }
      ret
    }
    idx = as.list(data.frame(points)[,group,drop=FALSE])
    ret = by(points, idx, fn)
    ret = do.call(rbind, ret)
    proj4string(ret) = proj4string(points)
  }
  ret
}