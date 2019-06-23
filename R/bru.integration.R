#' @title Generate integration points 
#'
#' @description 
#' This function generates points in one or two dimensions with a weight attached to each point.
#' The weighted sum of a function evaluated at these points is the integral of that function approximated
#' by linear basis functions. The parameter \code{region} describes the area(s) integrated over. 
#' 
#' In case of a single dimension \code{region} is supposed to be a two-column \code{matrix} where
#' each row describes the start and end point of the interval to integrate over. In the two-dimensional
#' case \code{region} can be either a \code{SpatialPolygon}, an \code{inla.mesh} or a 
#' \code{SpatialLinesDataFrame} describing the area to integrate over. If a \code{SpatialLineDataFrame}
#' is provided it has to have a column called 'weight' in order to indicate the width of the line.
#' 
#' The domain parameter is an \code{inla.mesh.1d} or \code{inla.mesh} object that can be employed to 
#' project the integration points to the vertices of the mesh. This reduces the final number of
#' integration points and reduces the computational cost of the integration. The projection can also 
#' prevent numerical issues in spatial LGCP models where each observed point is ideally surrounded
#' by three integration point sitting at the coresponding mesh vertices. For convenience, the
#' \code{domain} parameter can also be a single integer setting the number of equally spaced integration
#' points in the one-dimensional case. 
#' 
#' @aliases ipoints
#' @export
#' 
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 
#' @param region Description of the integration region boundary. 
#' In 1D either a vector of two numerics or a two-column matrix where each row describes and interval. 
#' In 2D either a \code{SpatialPolygon} or a \code{SpatialLinesDataFrame} with a weight column defining the width of the line.
#' @param domain In 1D a single numeric setting the numer of integration points or an \code{inla.mesh.1d} 
#' defining the locations to project the integration points to. In 2D \code{domain} has to be an
#' \code{inla.mesh} object describing the projection and granularity of the integration.
#' @param name Character array stating the name of the domains dimension(s)
#' @param group Column names of the \code{region} object (if applicable) for which the integration points are calculated independently and not merged by the projection.
#' @param project If TRUE, project the integration points to mesh vertices
#' 
#' @return A \code{data.frame} or \code{SpatialPointsDataFrame} of 1D and 2D integration points, respectively.
#' 
#' @examples
#' \donttest{
#' if (require("INLA", quietly = TRUE)) {
#' 
#' # Create 50 integration points covering the dimension 'myDim' between 0 and 10. 
#' 
#' ips = ipoints(c(0,10), 50, name = "myDim")
#' plot(ips)
#' 
#' # Create integration points for the two intervals [0,3] and [5,10]
#' 
#' ips = ipoints(matrix(c(0,3, 5,10), nrow = 2, byrow = TRUE), 50)
#' plot(ips)
#' 
#' # Convert a 1D mesh into integration points
#' mesh = inla.mesh.1d(seq(0,10,by = 1))
#' ips = ipoints(mesh, name = "time")
#' plot(ips)
#'
#' 
#' # Obtain 2D integration points from a SpatialPolygon
#' 
#' data(gorillas, package = "inlabru")
#' ips = ipoints(gorillas$boundary)
#' ggplot() + gg(gorillas$boundary) + gg(ips, aes(size = weight))
#' 
#' 
#' #' Project integration points to mesh vertices
#' 
#' ips = ipoints(gorillas$boundary, domain = gorillas$mesh)
#' ggplot() + gg(gorillas$mesh) +  gg(gorillas$boundary) + gg(ips, aes(size = weight))
#' 
#' 
#' # Turn a 2D mesh into integration points
#' 
#' ips = ipoints(gorillas$mesh)
#' ggplot() + gg(gorillas$boundary) + gg(ips, aes(size = weight))
#' }
#' }

ipoints = function(region = NULL, domain = NULL, name = "x", group = NULL, project) {
  
  pregroup = NULL
  
  # If region is null treat domain as the region definition
  if ( is.null(region) ) {
    if ( is.null(domain) ) { stop("regio and domain can not be NULL at the same time.") }
    else { region = domain ; domain = NULL }
  }
  
  
  if ( is.data.frame(region) ) {
    if (!("weight" %in% names(region))) { region$weight = 1 }
    ips = region
  }
  
  else if (is.integer(region)){
    
    ips = data.frame(weight = rep(1,length(region)))
    ips[name] = region
  }
  
  else if (is.numeric(region)) {
    
    if ( is.null(dim(region)) ){ region = matrix(region, nrow = 1) }
    
    if ( ncol(region) == 1) {
      
      ips = data.frame(x = region[,1], weight = 1)
      colnames(ips) = c(name, "weight")
      
    } else {
      
      ips = list()
      
      for (j in 1:nrow(region) ) {
        
        subregion = region[j,]
        
        # If domain is NULL set domain to a 1D mesh with 30 equally spaced vertices and boundary according to region
        # If domain is a single numeric set domain to a 1D mesh with n=domain vertices and boundary according to region
        if ( is.null(domain) ) { subdomain = INLA::inla.mesh.1d(seq(min(subregion), max(subregion), length.out = 30)) }
        else if ( is.numeric(domain)) { subdomain = INLA::inla.mesh.1d(seq(min(subregion), max(subregion), length.out = domain)) }
        else { subdomain = stop("1D weight projection not yet implemented") }
        
        fem = INLA::inla.mesh.1d.fem(subdomain)
        ips[[j]] = data.frame(weight = Matrix::diag(fem$c0))
        ips[[j]][name] = subdomain$loc
        ips[[j]] = ips[[j]][,c(2,1)] # make weights second column
      }
      
      ips = do.call(rbind, ips)
    }
    
  } else if ( inherits(region, "inla.mesh") ){
    
    # If domain is provided: break
    if ( !is.null(domain) ) stop("Integration region provided as 2D and domain is not NULL.")
    
    # transform to equal area projection
    if ( !is.null(region$crs) && !(is.na(region$crs@projargs))) {
      crs = region$crs
      region = stransform(region, crs = CRS("+proj=cea +units=km"))
    }
    
    ips = vertices(region)
    ips$weight = INLA::inla.mesh.fem(region, order = 1)$va
    
    # backtransform
    if ( !is.null(region$crs) && !(is.na(region$crs@projargs))) { ips = stransform(ips, crs = crs) }
    
  } else if ( inherits(region, "inla.mesh.1d") ){
    
    ips = data.frame(x = region$loc)
    colnames(ips) = name
    ips$weight = Matrix::diag(INLA::inla.mesh.fem(region)$c0)
    
  } else if ( class(region) == "SpatialPoints" ){
    
    ips = region
    ips$weight = 1
    
  } else if ( class(region) == "SpatialPointsDataFrame" ){
    
      if (!("weight" %in% names(region))) { 
        warning("The integration points provided have no weight column. Setting weights to 1.")
        region$weight = 1
      }
      
      ips = region
  
  } else if ( inherits(region, "SpatialLines") || inherits(region, "SpatialLinesDataFrame") ){
    
    # If SpatialLines are provided convert into SpatialLinesDataFrame and attach weight = 1
    if ( class(region)[1] == "SpatialLines" ) { 
      region = SpatialLinesDataFrame(region, data = data.frame(weight = rep(1, length(region)))) 
    }
    
    # Set weights to 1 if not provided
    if (!("weight" %in% names(region))) { 
      warning("The integration points provided have no weight column. Setting weights to 1.")
      region$weight = 1
    }
    
    ips = int.slines(region, domain, group = group)
  
  } else if (inherits(region,"SpatialPolygons")){
    
    # If SpatialPolygons are provided convert into SpatialPolygonsDataFrame and attach weight = 1
    if ( class(region)[1] == "SpatialPolygons" ) { 
      region = SpatialPolygonsDataFrame(region,
                                        data = data.frame(weight = rep(1, length(region))),
                                        match.ID = FALSE)
    }
    
    cnames = coordnames(region)
    p4s = proj4string(region)
    
    # Convert region and domain to equal area CRS
    if ( !is.null(domain$crs) && !is.na(domain$crs@projargs)){
      region = stransform(region, crs = CRS("+proj=cea +units=km"))
    }
    
    polyloc = do.call(rbind, lapply(1:length(region), 
                                    function(k) cbind(
                                      x = rev(coordinates(region@polygons[[k]]@Polygons[[1]])[,1]),
                                      y = rev(coordinates(region@polygons[[k]]@Polygons[[1]])[,2]),
                                      group = k)))
    
    # If domain is NULL, make a mesh with the polygons as boundary
    if ( is.null(domain) ) {
      max.edge = max(diff(range(polyloc[,1])), diff(range(polyloc[,2])))/20
      domain = INLA::inla.mesh.2d(boundary = region, max.edge = max.edge)
      domain$crs = CRS(proj4string(region))
    } else {
      if ( !is.null(domain$crs) && !is.na(domain$crs@projargs)) 
        domain = stransform(domain, crs = CRS("+proj=cea +units=km"))
    }
    
    ips = int.polygon(domain, loc = polyloc[,1:2], group = polyloc[,3])
    df = data.frame(region@data[ips$group, pregroup, drop = FALSE],
                    weight = ips[,"weight"])
    ips = SpatialPointsDataFrame(ips[,c("x","y")], data = df, match.ID = FALSE)
    proj4string(ips) = proj4string(region)
    
    if ( !is.na(p4s) ) {
      ips = stransform(ips, crs = CRS(p4s))  
    }
    
  }
  
  ips
  
}

#' @title Cross product of integration points
#'
#' @description 
#' Calculates the cross product of integration points in different dimensions
#' and multiplies their weights accordingly. If the object defining points in a particular
#' dimension has no weights attached to it all weights are assumend to be 1.
#' 
#' @aliases cprod
#' @export
#' 
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 
#' @param ... \code{data.frame} or \code{SpatialPointsDataFrame} objects, each one usually obtained by a call to the \link{ipoints} function.
#' @return A \code{data.frame} or \code{SpatialPointsDataFrame} of multidimensional integration points and their weights
#' 
#' @examples
#' \donttest{
#' # ipoints needs INLA
#' if (require("INLA", quietly = TRUE)) {
#' # Create integration points in dimension 'myDim' and 'myDiscreteDim' 
#' ips1 = ipoints(c(0,8), name = "myDim")
#' ips2 = ipoints(as.integer(c(1,2,3)), name = "myDiscreteDim")
#' 
#' # Calculate the cross product
#' ips = cprod(ips1, ips2)
#' 
#' # Plot the integration points
#' plot(ips$myDim, ips$myDiscreteDim, cex = 10*ips$weight)
#' }
#' }

cprod = function(...) {
  ipl = list(...)
  ipl = ipl[!vapply(ipl, is.null, TRUE)]
  if ( length(ipl) == 0 ) return(NULL)
  
  if ( length(ipl) == 1 ) {
    ips = ipl[[1]]
  } else {
    ips1 = ipl[[1]]
    ips2 = do.call(cprod, ipl[2:length(ipl)])
    if (! "weight" %in% names(ips1) ) { ips1$weight = 1 }
    if (! "weight" %in% names(ips2) ) { ips2$weight = 1 }
    
    loc1 = ips1[,setdiff(names(ipl[[1]]),"weight"), drop = FALSE]
    w1 = data.frame(weight = ips1$weight)
    loc2 = ips2[,setdiff(names(ips2),"weight"), drop = FALSE]
    w2 = data.frame(weight2 = ips2[,"weight"])
    
    # Merge the locations. In case of Spatial objects we need to use the sp:merge
    # function. Unfortunately sp::merge replicates entries in a different order than
    # base merge so we need to reverse the order of merging the weights
    
    if ( inherits(loc1, "Spatial") ) {
      ips = sp::merge(loc1, loc2, duplicateGeoms = TRUE) 
      weight = merge(w2, w1)
    } else if ( inherits(loc2, "Spatial") ){
      ips = sp::merge(loc2, loc1, duplicateGeoms = TRUE)
      weight = merge(w2, w1)
    } else {
      ips = merge(loc1, loc2)
      weight = merge(w1, w2)
    }
    ips$weight = weight$weight * weight$weight2
    
  }
  ips
}

# Integration points for log Gaussian Cox process models using INLA
# 
# prerequisits:
# 
# - List of integration dimension names, extend and quadrature 
# - Samplers: These may live in a subset of the dimensions, usually space and time 
#             ("Where and wehen did a have a look at the point process")
# - Actually this is a simplified view. Samplers should have start and end time !
# 
# Procedure:
# - Select integration strategy by type of samplers: 
#       1) SpatialPointsDataFrame: Assume these are already integration points
#       2) SpatialLinesDataFrame: Use simplified integration along line with (width provided by samplers)
#       3) SpatialPolygonDataFrame: Use full integration over polygons
#
# - Create integration points from samplers. Do NOT perform simplification projection here!
# - Simplify integration points. 
#   1) Group by non-mesh dimensions, e.g. time, weather
#   2) For each group simplify with respect to mesh-dimensions, e.g. space
#   3) Merge
#   
#   Dependencies (iDistance):
#   int.points(), int.polygon(), int.1d(), int.expand(), recurse.rbind()
#
# @aliases ipoints
# @export
# @param samplers A Spatial[Points/Lines/Polygons]DataFrame objects
# @param points A SpatialPoints[DataFrame] object
# @param config An integration configuration. See \link{iconfig}
# @return Integration points


ipmaker = function(samplers, domain, dnames, model = NULL, data = NULL) {
  
  # Fill missing domain definitions using meshes from effects where map equals the domain name
  meshes = list()
  for (e in effect(model)) {meshes[[paste0(as.character(e$map), collapse ="")]] = e$mesh}
  for ( nm in dnames) {
    if ( is.null(domain[[nm]]) ) { domain[[nm]] = meshes[[nm]] }
  }
  
  # Fill missing domain definitions with data ranges
  for ( nm in dnames) {
    if ( !(nm %in% names(domain)) & !is.null(data) & !(nm %in% names(samplers))){
      if ( nm == "coordinates" ) {
        domain[["coordinates"]] = INLA::inla.mesh.2d(loc.domain = coordinates(data), max.edge = diff(range(coordinates(data)[,1]))/10)
        domain[["coordinates"]]$crs = INLA::inla.CRS(proj4string(data))
      } else {
        domain[[nm]] = range(data[[nm]])
      }
    }
  }
  
  
  if ( "coordinates" %in% dnames ) { spatial = TRUE } else { spatial = FALSE }
  
  # Dimensions provided via samplers (except "coordinates")
  samp.dim = intersect(names(samplers), dnames)
  
  # Dimensions provided via domain but not via samplers
  nosamp.dim = setdiff(names(domain), c(samp.dim, "coordinates"))
  
  # Check if a domain definition is missing
  missing.dims = setdiff(dnames, c(names(domain), samp.dim))
  if ( length(missing.dims > 0) ) stop(paste0("Domain definitions missing for dimensions: ", paste0(missing.dims, collapse = ", ")))
  
  if ( spatial ) {
    ips = ipoints(samplers, domain$coordinates, project = TRUE, group = samp.dim)
  } else { 
    ips = NULL
  }
  
  
  lips = lapply(nosamp.dim, function(nm) ipoints(NULL, domain[[nm]], name = nm))
  ips = do.call(cprod, c(list(ips), lips))
  
}




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

vertex.projection = function(points, mesh, columns = names(points), group = NULL, fill = NULL){
  
  if ( is.null(group) | (length(group) == 0) ) {
    
    res = INLA::inla.fmesher.smorg(mesh$loc, mesh$graph$tv, points2mesh = coordinates(points))
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
    
    ret = SpatialPointsDataFrame(coords,
                                 proj4string = CRS(proj4string(points)),
                                 data = data,
                                 match.ID = FALSE)
    coordnames(ret) = coordnames(points)
    
    # If null is not not NULL, add vertices to which no data was projected
    # and set their projected data according to `fill`
    
    if ( !is.null(fill) ) {
      vrt = vertices(mesh)
      vrt = vrt[setdiff(vrt$vertex, data$vertex),]
      if ( nrow(vrt) > 0 ){
        for (nm in setdiff(names(data), "vertex")) vrt[[nm]] = fill
        ret = rbind(ret, vrt)
      }
      ret = ret[match(1:mesh$n, ret$vertex),]
    } 
    
    
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
# @example
#
# pts = data.frame(x = 50 * runif(10), weight = abs(rnorm(100)))
# msh = inla.mesh.1d(seq(0,50,by=1))
# pts$year = c(rep(1,5), rep(2,5))
# ip =  vertex.projection.1d(pts, msh)
# ggplot(ip) + geom_point(aes(x=x, y=weight))
#
# ip =  vertex.projection.1d(pts, msh, group = "year", fill = 0, column = "weight")
# head(ip)
# ggplot(ip) + geom_point(aes(x=x, y=weight, color = year))

vertex.projection.1d = function(points, mesh, group = NULL, column = "weight", simplify = TRUE, fill = NULL) {

  dname = setdiff(names(points),c(column, group))
  if ( length(dname)>1 ) { dname = dname[1] }
  
  xx = points[, dname]
  ww = points[, column]
  iv = findInterval(xx, mesh$loc)
  
  # Left and right vertex location
  left = mesh$loc[iv]  
  right = mesh$loc[iv+1]  
  
  # Relative location within the two neighboring vertices
  w.right = (xx-left)/(right-left)
  w.left = 1 - w.right
  
  # Projected integration points
  ips = rbind(data.frame(x = left, vertex = iv),
              data.frame(x = right, vertex = iv+1))
  ips[column] = c(ww * w.left, ww * w.right)
  
  
  # Simplify 
  if ( simplify ) {
    bygroup = list(vertex = ips$vertex)
    if ( !is.null(group) ) { bygroup = c(bygroup, as.list(rbind(points[,group,drop=FALSE], points[,group,drop=FALSE]))) }
    ips = aggregate(ips[,column, drop = FALSE], by = bygroup, FUN = sum)
  }
  
  # Add x-coordinate
  ips[dname] = mesh$loc[ips$vertex]
  
  # Fill
  if ( !is.null(fill) ) { 
    miss = setdiff(1:length(mesh$loc), ips$vertex)
    mips = data.frame(vertex = miss, x = mesh$loc[miss])
    mips[,column] = fill
    ips = rbind(ips, merge(mips, ips[,group, drop=FALSE]))
  }
  
  ips
}


#' Weighted summation (integration) of data frame subsets
#'
#' A typical task in statistical inference to integrate a (multivariate) function along one or
#' more dimensions of its domain. For this purpose, the function is evaluated at some points
#' in the domain and the values are summed up using weights that depend on the area being 
#' integrated over. This function performs the weighting and summation conditional for each level
#' of the dimensions that are not integrated over. The parameter \code{dims} states the the 
#' dimensions to integrate over. The set of dimensions that are held fixed is the set difference
#' of all column names in \code{data} and the dimensions stated by \code{dims}.
#' 
#' @aliases int
#' @export
#' @param data A \code{data.frame} or \code{Spatial} object. Has to have a \code{weight} column with numeric values.
#' @param values Numerical values to be summed up, usually the result of function evaluations.
#' @param dims Column names (dimension names) of the \code{data} object to integrate over.
#' @return A \code{data.frame} of integrals, one for each level of the cross product of all dimensions not being integrated over.
#' 
#' @examples 
#' \donttest{
#' # ipoints needs INLA
#' if (require("INLA", quietly = TRUE)) {
#' # Create integration points in two dimensions, x and y
#'
#' ips = cprod(ipoints(c(0,10), 10, name = "x"),
#'             ipoints(c(1,5), 10, name = "y"))
#'
#' # The sizes of the domains are 10 and 4 for x and y, respectively.
#' # Integrating f(x,y) = 1 along x and y should result in the total
#' # domain size 40
#'
#' int(ips, rep(1, nrow(ips)), c("x","y"))
#' }
#' }



int = function(data, values, dims = NULL) {
  keep = setdiff(names(data), c(dims, "weight"))
  if (length(keep) > 0 & !is.null(dims)) { 
    agg = aggregate(values * data$weight, by = as.list(data[,keep,drop=FALSE]), FUN = sum)
    names(agg)[ncol(agg)] = "integral" # paste0("integral_{",dims,"}(",deparse(values),")")
  } else { 
    agg = sum(data$weight * values) 
  }
  
  agg
}

