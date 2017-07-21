#' @title Generate integration points 
#'
#' @description Generate integration points based on a desciption of the integration region.
#' @aliases ipoints
#' @export
#' 
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 
#' @param region Description of the integration region boundary
#' @param domain An object describing a discretization of the domain
#' @param name Character array stating the name of the domains dimension(s)
#' @param group Column names of the \code{region} object (if applicable) for which the integration points are calculated independently and not merged.
#' @param project If TRUE, project the integration points to mesh vertices
#' @return A \code{data.frame} or \code{SpatialPointsDataFrame}
#' 
#' @examples
#' 
#' ips = ipoints(c(0,10), name = "myDim")
#' ips
#' 

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
        if ( is.null(domain) ) { subdomain = inla.mesh.1d(seq(min(subregion), max(subregion), length.out = 30)) }
        else if ( is.numeric(domain)) { subdomain = inla.mesh.1d(seq(min(subregion), max(subregion), length.out = domain)) }
        else { subdomain = stop("1D weight projection not yet implemented") }
        
        fem = inla.mesh.1d.fem(subdomain)
        ips[[j]] = data.frame(weight = diag(as.matrix(fem$c0)))
        ips[[j]][name] = subdomain$loc
        ips[[j]] = ips[[j]][,c(2,1)] # make weights second column
      }
      
      ips = do.call(rbind, ips)
    }
    
  } else if ( inherits(region, "inla.mesh") ){
    
    # If domain is provided: break
    if ( !is.null(domain) ) stop("Integration region provided as 2D and domain is not NULL.")
    
    # transform to equal area projection
    if ( !is.null(region$crs) ) {
      crs = region$crs
      region = stransform(region, crs = CRS("+proj=cea +units=km"))
    }
    
    ips = vertices(region)
    ips$weight = diag(as.matrix(inla.mesh.fem(region)$c0))
    
    # backtransform
    if ( !is.null(region$crs) ) { ips = stransform(ips, crs = crs) }
    
  } else if ( inherits(region, "inla.mesh.1d") ){
    
    ips = data.frame(x = region$loc)
    colnames(ips) = name
    ips$weight = diag(as.matrix(inla.mesh.fem(region)$c0))
    
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
      region = SpatialPolygonsDataFrame(region, data = data.frame(weight = rep(1, length(region)))) 
    }
    
    cnames = coordnames(region)
    p4s = proj4string(region)
    
    # Convert region and domain to equal area CRS
    if ( !is.null(domain$crs) ){
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
      domain = inla.mesh.2d(boundary = region, max.edge = max.edge)
      domain$crs = CRS(proj4string(region))
    } else {
      if ( !is.null(domain$crs) ) domain = stransform(domain, crs = CRS("+proj=cea +units=km"))
    }
    
    ips = int.polygon(domain, loc = polyloc[,1:2], group = polyloc[,3])
    df = data.frame(region@data[ips$group, pregroup, drop = FALSE], weight = ips[,"weight"])
    ips = SpatialPointsDataFrame(ips[,c("x","y")],data = df)
    proj4string(ips) = proj4string(region)
    ips = stransform(ips, crs = CRS(p4s))  
  }
  
  ips
  
}

#' @title Cross product of integration points
#'
#' @description Calculates the dimensional cross product of integration points and multiply their weights accordingly.
#' @aliases cprod
#' @export
#' 
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 
#' @param ... \code{data.frame} or \code{SpatialPointsDataFrame} objects
#' @return A \code{data.frame} or \code{SpatialPointsDataFrame} object
#' 
#' @examples
#' 
#' ips1 = ipoints(c(0,8), name = "myDim")
#' ips2 = ipoints(as.integer(c(0,2,5)), name = "myDiscreteDim")
#' ips = cprod(ips1, ips2)
#' plot(ips$myDim) ; points(ips$myDiscreteDim, col = "red")

cprod = function(...) {
  ipl = list(...)
  ipl = ipl[!sapply(ipl, is.null)]
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
#' 
#' @aliases int
#' @export
#' @param data A \code{data.frame} or \code{Spatial} object. Has to have a weight column with numeric values.
#' @param values Numerical values to be summed up.
#' @param dims Columns of the \code{data} obect to integrate over
#' @return A \code{data.frame} of integrated values


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

