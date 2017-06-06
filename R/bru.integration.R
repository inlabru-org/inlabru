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

ipoints = function(region, domain = NULL, name = "x", group = NULL, project) {
  
  pregroup = NULL
  
  if ( is.data.frame(region) ) {
    if (!("weight" %in% names(region))) { region$weight = 1 }
    ips = region
  }
  
  else if (is.integer(region)){
    
    ips = data.frame(weight = rep(1,length(region)))
    ips[name] = region
  }
  
  else if (is.numeric(region)) {
    # If domain is NULL set domain to a 1D mesh with 30 equally spaced vertices and boundary according to region
    # If domain is a single numeric set domain to a 1D mesh with n=domain vertices and boundary according to region
    if ( is.null(domain) ) { domain = inla.mesh.1d(seq(region[1], region[2], length.out = 30)) }
    else if ( is.numeric(domain)) { domain = inla.mesh.1d(seq(region[1], region[2], length.out = domain)) }
    
    fem = inla.mesh.1d.fem(domain)
    ips = data.frame(weight = diag(as.matrix(fem$c0)))
    ips[name] = domain$loc
    ips = ips[,c(2,1)] # make weights second column
  
  } else if ( inherits(region, "inla.mesh") ){
    
    ips = vertices(region)
    ips$weight = diag(as.matrix(inla.mesh.fem(region)$c0))
    
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
      domain = stransform(domain, crs = CRS("+proj=cea +units=km"))
    }
      
    polyloc = do.call(rbind, lapply(1:length(region), 
                                    function(k) cbind(
                                      x = rev(coordinates(region@polygons[[k]]@Polygons[[1]])[,1]),
                                      y = rev(coordinates(region@polygons[[k]]@Polygons[[1]])[,2]),
                                      group = k)))
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

ipmaker = function(samplers, config) {
  
  # First step: Figure out which marks are already included in the samplers (e.g. time)
  # For each of the levels of these marks the integration points will be computed independently
  pregroup = intersect(setdiff(names(samplers), "weight"), names(config))
  postgroup = setdiff(names(config), c(pregroup, "coordinates"))
  
  # By default we assume we are handling a Spatial* object
  spatial = TRUE
  
  #
  # Inital integration points
  # 
  
  if ( !is.null(samplers) ) {
    ips = ipoints(samplers, config[["coordinates"]]$mesh, group = pregroup, project = FALSE)
    
    # Reduce number of integration points by projection onto mesh vertices
    all.groups = intersect(names(ips), c(postgroup,pregroup))
    if ( length(all.groups) == 0 ) all.groups = NULL
    if ("coordinates" %in% names(config)) {
      ips = vertex.projection(ips, config[["coordinates"]]$mesh, columns = "weight", group = all.groups)
    }
    
  } else {
    ips = NULL
  }
  
  #
  # Expand the integration points over postgroup dimensions
  #
  
  for ( k in seq_len(length(postgroup)) ){
    gname = postgroup[[k]]
    cfg = config[[gname]]
    nips = ipoints(c(cfg$sp, cfg$ep), cfg$n, name = gname)
    if ( is.null(ips) ) { ips = nips } else { ips = cprod(ips, nips) }
  }
  
  ips
}


#' Integration point configuration generator
#'   
#'
#' @aliases iconfig
#' @export
#' @param samplers A Spatial[Points/Lines/Polygons]DataFrame objects
#' @param points A SpatialPoints[DataFrame] object
#' @param model A \link{model}
#' @param dim.names Dimension names (character array)
#' @param mesh default spatial mesh used for integration
#' @return An integration configuration

iconfig = function(samplers, points, model, dim.names = NULL, mesh = NULL, domain = NULL) {
  
  # Obtain dimensions to integrate over. 
  # These are provided as the left hand side of model$formula
  if ( is.null(dim.names) ) { dim.names = model$dim.names }
  
  # Default mesh for spatial integration:
  meshes = lapply(effect(model), function(e) e$mesh)
  names(meshes) = elabels(model)
  meshes = meshes[!unlist(lapply(meshes, is.null))]
  issp = as.vector(unlist(lapply(meshes, function(m) m$manifold == "R2")))
  if (!is.null(issp) & any(issp)) { 
    spmesh = meshes[[which(issp)[[1]]]] } 
  else { spmesh = mesh }
  
  # Function that maps each dimension name to a setup
  ret = make.setup = function(nm) {
    ret = list()
    
    # Set the name of the coordinate system
    ret$name = nm
    ret$mesh = NULL
    ret$project = FALSE
    
    # Create a function to fetch coordinates. 
    # "coordinates" is a special case for which we extract an integration mesh from the model
    if ( nm == "coordinates" ) {
      ret$get.coord = get0(nm)
      ret$n.coord = ncol(ret$get.coord(points))
      # Use the first spatial mesh that we find
      ret$mesh = spmesh
      ret$class = "matrix"  
      ret$project = TRUE
      ret$p4s = proj4string(points)
    } else {
      # Spatial* object or data.frame?
      if ( is.data.frame(points) ) { ret$get.coord = function(sp) { sp[,nm, drop = FALSE] }} 
      else { ret$get.coord = function(sp) { sp@data[,nm, drop = FALSE] } }
      ret$n.coord = 1
      ret$class = class(ret$get.coord(points)[,nm])
    }
    
    # If there is a mesh in the model that has the same name as the dimension
    # use the mesh knots/vertices as integration points. Otherwise use some
    # standard setting depending on the class of the dimension.
    
    if ( nm %in% names(meshes) ) {
      if ( meshes[[nm]]$manifold == "R1" ) {
        # The mesh is 1-dimensional. 
        # By default a two-point Gaussian quadrature is used for each interval between the mesh knots.
        # This intergration is exact for poylnomials of second degree. 
        ret$mesh = meshes[[nm]]
        ret$min = min(ret$mesh$loc)
        ret$max = max(ret$mesh$loc)
        ret$cnames = colnames(ret$get.coord(points))
        ret$scheme = "gaussian"
        ret$sp = ret$mesh$loc[1:(length(ret$mesh$loc)-1)]
        ret$ep = ret$mesh$loc[2:length(ret$mesh$loc)] 
        ret$n.points = 2
      } else { 
        stop( "not implemented: 2d mesh auto-integration") 
      }
    } else {
      # If the dimension is defined via the formula extract the respective information.
      # If not, extract information from the data
      if ( nm %in% names(domain) ) {
        ret$min = min(domain[[nm]][1])
        ret$max = max(domain[[nm]][2])
        ret$n.points = 30
      } else {
        ret$min = apply(ret$get.coord(points), MARGIN = 2, min)
        ret$max = apply(ret$get.coord(points), MARGIN = 2, max)
      }
      
      ret$cnames = colnames(ret$get.coord(points))
      
      if ( ret$class == "integer" ) {
        ret$scheme = "fixed"
        ret$sp = ret$min:ret$max
        ret$ep = NULL
        if (is.null(ret$n.points)) { ret$n.points = length(ret$sp) }
      } else if ( ret$class == "factor" ) {
        ret$scheme = "fixed"
        stop("Not implemented.")
      } else { 
        ret$scheme = "trapezoid"
        ret$sp = ret$min
        ret$ep = ret$max
        if (is.null(ret$n.points)) { ret$n.points = 20 }
      }      
    }
    ret
  }

  # Create configurations
  ret = lapply(dim.names, make.setup)
  names(ret) = dim.names
  ret
}


# Integration point configuration generator
#   
# @keywords internal
# @return An integration configuration

iconfig2 = function(data, dims, domain, samplers) {
  icfg = list()
  
  for (dname in dims){
    ret = list()
    if (dname == "coordinates") {
      ret$name = "coordinates"
      ret$domain = domain[[dname]]
      ret$class = "inla.mesh"
    } else {
      ret$name = dname
      # If no domain specification is given use data range for this dimension to construct a 1D mesh
      if ( !(dname %in% names(domain)) ) {
        domain[[dname]] = inla.mesh.1d(seq(min(data.frame(data)[,dname]), max(data.frame(data)[,dname]),length.out = 25))
      }
      ips = NULL
      if ( class(domain[[dname]]) == "inla.mesh.1d" ) { 
        dom = domain[[dname]] 
        ips = data.frame(dom$loc, weight = diag(as.matrix(inla.mesh.1d.fem(dom)$c0)))
        colnames(ips)[1] = dname
        }
      else if ( is.integer(domain[[dname]]) ) {
        dom = inla.mesh.1d(domain[[dname]])
        ips = data.frame(dom$loc, weight = 1)
        colnames(ips)[1] = dname
        }
      else if ( is.numeric(domain[[dname]]) ) { 
        dom = inla.mesh.1d(domain[[dname]])
        ips = data.frame(dom$loc, weight = diag(as.matrix(inla.mesh.1d.fem(dom)$c0)))
        colnames(ips)[1] = dname
        }
      else if ( is.factor(domain[[dname]]) ) { dom = domain[[dname]] }
      else { stop(sprintf("Unsupported domain class '%s' for domain '%s'"), class(domain[[dname]]), dname) }
      ret$domain = dom
      ret$ips = ips
      ret$class = class(domain[[dname]])
    }
    icfg[[dname]] = ret
  }
  
  if ( inherits(samplers, "Spatial") | is.data.frame(samplers) ){
    samplerdf = samplers
    for ( dname in setdiff(dims,"coordinates") ) {
      if ( dname %in% names(samplerdf) ) { 
        sdf = data.frame(samplerdf[[dname]], weight = 1)
        colnames(sdf)[1] = dname
        icfg[[dname]]$tips = sdf
      } else if ( paste0("min.",dname) %in% names(samplerdf) ) {
        

        
        ips = mint(as.data.frame(samplerdf)[,paste0(c("min.","max."),dname)], icfg[[dname]]$domain$loc, dname)
        
        
        stop("s")
      }
      
    }
  }
  
  icfg
 
}

mint = function(df,knots,dname) {
  ips = list()
  for (k in 1:nrow(df)) {
    mi = df[k,1]
    ma = df[k,2]
    msh = inla.mesh.1d(sort(unique(c(mi,ma,knots[knots<=ma & knots>=mi]))))
    w = diag(as.matrix(inla.mesh.1d.fem(msh)$c0))
    ips[[k]] = data.frame(msh$loc, weight = w)
  }
  ips = do.call(rbind, ips)
  colnames(ips) = c(dname, "weight")
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
