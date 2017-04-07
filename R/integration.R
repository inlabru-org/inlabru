# Split lines at mesh edges
#
# @aliases split.lines
# @export
# @param mesh An inla.mesh object
# @param sp Start points of lines
# @param ep End points of lines
# @param filter.zero.length Filter out segments with zero length? (Bool)
# @param ... argments to int.quadrature
# @return List of start and end points resulting from splitting the given lines
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

split.lines = function(mesh, sp, ep, filter.zero.length = TRUE) {
  
  # locations for splitting
  loc = as.matrix(rbind(sp,ep))
  idx = 1:dim(sp)[1]
  
  # Filter out segments not on the mesh
  t1 = inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(sp,z=0)))$p2m.t
  t2 = inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(ep,z=0)))$p2m.t
  # if (any(t1==0) | any(t2==0)) { warning("points outside boundary! filtering...")}
  sp = sp[!((t1==0) | (t2==0)),]
  ep = ep[!((t1==0) | (t2==0)),]
  idx = idx[!((t1==0) | (t2==0))]
  loc = as.matrix(rbind(sp,ep))
  
  # Split them segments into parts
  if ( dim(loc)[2] == 2 ) {loc = cbind(loc,rep(0,dim(loc)[1]))}
  np = dim(sp)[1]
  sp.idx = t(rbind(1:np,np+1:np))
  splt = inla.fmesher.smorg(mesh$loc,mesh$graph$tv, splitlines=list(loc=loc, idx=sp.idx))
  #plot(data$mesh)
  #points(loc)
  #points(splt$split.loc,col="blue)
  
  sp = splt$split.loc[splt$split.idx[,1],1:dim(sp)[2]] # Start point of new segments
  ep = splt$split.loc[splt$split.idx[,2],1:dim(ep)[2]] # End points of new segments
  idx = idx[splt$split.idx[,1]]
  origin = splt$split.origin
  
  # Filter out zero length segments
  if ( filter.zero.length ) {
    sl = apply((ep-sp)^2,MARGIN=1,sum)
    sp = sp[!(sl==0),]
    ep = ep[!(sl==0),]
    origin = origin[!(sl==0)]
    idx = idx[!(sl==0)]
  }
  
  return(list(sp=sp,ep=ep,split.origin=origin,idx=idx,split.loc=splt$split.loc))
  
}

# Replicate lines with different distances to base line
#
# @aliases replicate.lines
# @export
# @param sp Start points of lines
# @param ep End points of lines
# @param truncation distance at which we truncate sightings
# @param ... argments to int.quadrature
# @return List of 1) Start and end points of replicated lines 2) their distances to the original line
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

replicate.lines = function(sp,ep,truncation,scheme="equidistant",n=3,fake.distance=FALSE,geometry="euc") {
  
  if (fake.distance) {
    quad = int.quadrature(0,truncation,scheme,n)
  } else {
    quad = int.quadrature(-truncation,truncation,scheme,n)
  }
  
  dst = quad$ips
  idx = 1:dim(sp)[1]
  
  if (geometry == "euc"){
    # Normal vectors
    dvec = cbind(-(ep[,2]-sp[,2]),ep[,1]-sp[,1])
    nrm = apply(dvec^2,MARGIN=1,function(x) {return(sqrt(sum(x)))})
    dvec = dvec/rep(nrm,dim(dvec)[2])
    
    tmp.sp = list()
    tmp.ep = list()
    distance = list()
    for (k in 1:length(dst)){
      if (fake.distance){
        tmp.sp[[k]] = sp
        tmp.ep[[k]] = ep
      }
      else {
        tmp.sp[[k]] = sp + dvec*dst[k]
        tmp.ep[[k]] = ep + dvec*dst[k]
      }
      distance[[k]] = rep(dst[k],dim(sp)[1])
    }
    sp  = do.call(rbind,tmp.sp)
    ep  = do.call(rbind,tmp.ep)
    distance = abs(do.call(c,distance))
    idx = rep(idx,length(dst))
    return(list(sp=sp,ep=ep,distance=distance,idx=idx))
  
  } else if (geometry == "geo"){
    stop("not supported: geometry geo")
  }
}

# Gaussian quadrature and other integration point constructors
# 
# Contruct integration points for each of lines defined by the start and end points provided.
# The following schemes are available: 
# "equidistant" : Equidistant integration points without boundary. All weights are identical and sum uf to the length of a line.
# "gaussian": Points and weight according to the Gaussian quadrature rule. Currently only n=1 and n=2 are supported (Exact integration for linear and quadratic functions).
# "twosided-gaussian": Experimental
#    
# @aliases int.quadrature
# @export
# @param sp Start points of lines
# @param ep End points of lines
# @param scheme Integration scheme (gaussian or equdistant)
# @param n Number of integration points
# @return List with integration poins (ips), weights (w) and weights including line length (wl)
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

int.quadrature = function(sp=NULL,ep=NULL,scheme="gaussian",n.points=1,geometry="euc",coords = NULL){

  if ( is.null(colnames(sp)) & !is.null(coords)) { sp = data.frame(sp) ; colnames(sp) = coords }
  if ( is.null(colnames(ep)) & !is.null(coords) & !(length(ep)==0) ) { ep = data.frame(ep) ; colnames(ep) = coords  }
  
  if (is.vector(sp)) { 
    n.lines = 1
    len = function(a) {abs(a)}
  }
  else { 
    n.lines = dim(sp)[1] 
    len = function(x) {apply(x^2,MARGIN=1,function(y){return(sqrt(sum(y)))})}
  }
  if (scheme == "gaussian"){
    if (n.points==1){
      # points
      if (geometry == "euc"){
        ips = sp + (ep-sp)/2
        # weights
        w = rep(1,n.lines)
        wl = w*len(ep-sp)
      }  
      else if (geometry == "geo"){
        colnames(sp) = coords
        sp.euc = geo.to.euc(data.frame(sp))
        colnames(ep) = coords
        ep.euc = geo.to.euc(data.frame(ep))
        ips = sp.euc + (ep.euc-sp.euc)/2
        ips = euc.to.geo(data.frame(ips))[,coords,drop=FALSE]
        # weights
        w = rep(1,n.lines)
        wl = w*dist.geo(data.frame(sp),data.frame(ep))
      }

      # PD vector
      # pdv = cbind(-(ep[,2]-sp[,2]),ep[,1]-sp[,1])
      # Index of line a point comes from
      line.idx = 1:n.lines
      
    } else if (n.points==2){
      if (geometry == "geo") {stop("Geometry geo not supported")}
      # points
      ips1 = sp + (-0.5*sqrt(1/3)+1/2) * (ep-sp)
      ips2 = sp + (0.5*sqrt(1/3)+1/2) * (ep-sp)
      ips = rbind(ips1,ips2)
      
      # weights
      w = rep(1,dim(ips)[1])/2
      wl = w*len(ep-sp)
      
      # PD vector
      #dvec = cbind(-(ep[,2]-sp[,2]),ep[,1]-sp[,1])
      #dvec = rbind(dvec,dvec)
      
      # Index of line a point comes from
      line.idx = c(1:n.lines,1:n.lines)
      
    } else if (n.points==3){
      if (geometry == "geo") {stop("Geometry geo not supported")}
      # points
      ips1 = sp + (-0.5*sqrt(3/5)+1/2) * (ep-sp)
      ips2 = sp + 0.5*(ep-sp)
      ips3 = sp + (0.5*sqrt(3/5)+1/2) * (ep-sp)
      ips = rbind(ips1,ips2,ips3)
      
      # weights
      w = c(rep(5/9,n.lines),rep(8/9,n.lines),rep(5/9),n.lines)/2
      wl = w*len(ep-sp)
      
      # Index of line a point comes from
      line.idx = c(1:n.lines,1:n.lines)
    }
    
    else { stop("Gaussian quadrature with n>3 not implemented") }
  }
  else if (scheme == "twosided-gaussian"){
    if (geometry == "geo") {stop("Geometry geo not supported")}
    ips1 = int.quadrature(sp,sp+0.5*(ep-sp),scheme="gaussian",n.points=n.points, geometry, coords)
    ips2 = int.quadrature(sp+0.5*(ep-sp),ep,scheme="gaussian",n.points=n.points, geometry, coords)
    ips = rbind(ips1$ips,ips2$ips)
    w = 0.5*rbind(ips1$w,ips2$w)
    wl = 0.5*rbind(ips1$wl,ips2$wl)
    line.idx = rbind(ips1$line.idx,ips2$line.idx)
  }
  else if (scheme == "equidistant"){
    if (geometry == "geo") {stop("Geometry geo not supported")}
    # points
    mult = seq(0-1/(2*n.points),1+1/(2*n.points),length.out=n.points+2)[2:(n.points+1)]
    nips = list()
    for (k in 1:n.points){
      nips[[k]] = sp + mult[k]*(ep-sp)
    }
    ips = do.call(rbind,nips)
    # weights
    w = rep(1,dim(ips)[1])/n.points
    wl = w*(len(ep-sp)[rep(1:n.lines,n.points)])
    # Index of line a point comes from
    line.idx = rep(1:n.lines,n.points)
  }
  else if (scheme == "trapezoid"){
    if (geometry == "geo") {stop("Geometry geo not supported")}
    # points
    mult = seq(0,1,length.out=n.points)
    nips = list()
    for (k in 1:n.points){
      nips[[k]] = sp + mult[k]*(ep-sp)
    }
    ips = do.call(rbind,nips)
    # weights
    w = rep(1,dim(ips)[1])/(n.points-1)
    w[1] = 0.5 * w[1]
    w[length(w)] = 0.5 * w[length(w)] 
    wl = w * (len(ep-sp)[rep(1:n.lines,n.points)])
    # Index of line a point comes from
    line.idx = rep(1:n.lines,n.points)
  }
  else if (scheme == "fixed"){ 
    ips = data.frame(tmp = sp)
    #colnames(ips) = coords
    w = rep(1,length(sp))
    wl = rep(1,length(sp))
    line.idx = rep(NaN,length(sp))
  }
  return(list(ips=ips,w=w,wl=wl,line.idx=line.idx))
}


# Integration points in one dimension
# 
# This function is a wrapper for \link{int.quadrature}. It returns a data frame instead of a list.
# 
# @aliases int.1d
# @param ... see \link{int.quadrature}
# @export 
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

int.1d = function(...){
  quad = int.quadrature(...)
  ips = data.frame(quad$ips, weight = quad$wl)
  return(ips)
}

# Line transect integration
# 
# Creates a set of integration points with weights from a \link{dsdata} structure. 
# The integration points can be based on the transect lines or the line segments of the survey (parameter "on").
# By default, the integration points are arranged on a grid and their number in direction of the transect and perpendicular
# to the transect are given by the parameters n.length and n.distance. The boundary of the integration in perpendicular
# direction is given by distance.truncation. Since transect lines can be long compared to the mesh that is used to model
# the data it can be useful to set the parameter line.split to TRUE. This means that transect lines / segments are split
# into parts at the edges of the mesh before the integration points are constructed (for each of these parts). The
# mesh that is used for this procedure is (by default) the mesh of the data set (data$mesh). However, it can by replaced
# by a mesh provided as a parameter. Additionaly, the given mesh can automatically be refined (mesh.refine) if a denser set 
# of integration points is required. The latter also holds for integration points constructed using the option 
# "projection = TRUE". Hereby, after their contruction, integration points are projected onto points at the mesh vertices.
#
# @aliases int.points
# @export
# @param data Either a dsdata/etpdata data set (e.g. whales) or a data.frame describing effort data
# @param on Either "transect" or "segment". This determines on which of these the integration is based on. Alternatively a two column index matrix, first column: column index of transect start points in effort data, second column: column index of transect end points in effort data
# @param line.split TRUE or FALSE, determines if lines that cross mesh edged should be splitted
# @param mesh Mesh used to construct the integration points. By default the mesh of the given data set.
# @param mesh.split Split mesh triangles into four sub-triangles for refined integration.
# @param mesh.coords Character description of the mesh coordinates, e.g. c("lon","lat")
# @param geometry Either "geo" (geographic)  or "euc" (Euclidean)
# @param length.scheme Integration scheme along the line (transect/segment)
# @param n.length Number of integration points along the line
# @param distance.scheme Integration scheme along perpendicular distance
# @param n.distance Number of integration points perpendicular distance
# @param distance.truncation Truncation for perpendicular distance (i.e. integration limit)
# @param fake.distance Wether or not integration points stay on the transect line and distances are faked.
# @param projection Type of projection. Currently only "linear" works.
# @param group Create independent integration points for sub-groups of the data, e.g. group.by=c("year","month")
# @return List with integration poins (ips), weights (w) and weights including line length (wl)
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>


int.points = function(data, 
                  on = "segment",
                  line.split = FALSE,
                  mesh = data$mesh, 
                  mesh.split = FALSE,
                  mesh.coords = data$mesh.coords,
                  geometry = data$geometry, 
                  length.scheme = "gaussian", 
                  n.length = 1, 
                  distance.scheme = "equidistant", 
                  n.distance = 5, 
                  distance.truncation = 1, 
                  fake.distance = FALSE,
                  projection = NULL,
                  group = NULL,
                  filter.zero.length = TRUE
){
  
  if (mesh.split) { mesh = tsplit.inla.mesh(mesh) }
  
  if (is.data.frame(data)){
    covariates = data
    # Segment/transect start and end points 
    sp = strip.coords(data[, paste0("start.",mesh.coords)])
    ep = strip.coords(data[, paste0("end.",mesh.coords)])
    idx = data.frame(1:nrow(sp), 1:nrow(sp))
  } else if ( class(data) == "SpatialLinesDataFrame" ) {
    qq = coordinates(data)
    sp = do.call(rbind, lapply(qq, function(k) do.call(rbind, lapply(k, function(x) x[1:(nrow(x)-1),]))))
    ep = do.call(rbind, lapply(qq, function(k) do.call(rbind, lapply(k, function(x) x[2:(nrow(x)),]))))
    
    idx = do.call(rbind, lapply(1:length(qq), function(k) do.call(cbind, lapply(qq[[k]], function(x) rep(k, nrow(x)-1) ))))
    idx = cbind(idx, idx)

  } else {
    # Get segments/transect indices
    if (on == "transect" ){ idx = as.transect(data) }
    else if (on == "segment" ) { idx = as.segment(data) }
    else ( idx = on )
    
    # Covariates
    covariates = data$effort
    # Segment/transect start and end points 
    sp = startpoint(idx,data)[,mesh.coords]
    ep = endpoint(idx,data)[,mesh.coords]
  }
  
  # Filter out bogus stuff
  loc = as.matrix(rbind(sp,ep))
  t1 = inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(sp,z=0)))$p2m.t
  t2 = inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(ep,z=0)))$p2m.t
  if (any(t1==0) | any(t2==0)) { warning("filtering...")}
  sp = sp[!((t1==0) | (t2==0)),]
  ep = ep[!((t1==0) | (t2==0)),]
  idx = idx[!((t1==0) | (t2==0)),]
  
  
  # Replicate segments for different distances (out:sp,ep,distance)
  line.rep = replicate.lines(sp,ep,distance.truncation,scheme=distance.scheme,n=n.distance,fake.distance=fake.distance)
  sp = line.rep$sp
  ep = line.rep$ep
  distance = line.rep$distance
  #split.origin=splt$split.origin,idx=idx
  idx = idx[line.rep$idx,]
  
  # Split lines
  if ( line.split ) {
    line.spl = split.lines(mesh, sp, ep, filter.zero.length)
    sp = line.spl$sp
    ep = line.spl$ep
    distance = distance[line.spl$split.origin]
    idx = idx[line.spl$split.origin,]
  }
  
  # Determine integration points along lines
  quad = int.quadrature(sp, ep, scheme = length.scheme, n.points = n.length, geometry=geometry,coords=mesh.coords)
  ips = quad$ips 
  w = 2*distance.truncation*quad$wl/n.distance
  distance = distance[quad$line.idx]
  idx = idx[quad$line.idx,]
  
  # Wrap everything up and perform projection according to distance and given group argument
  ips = data.frame(ips)
  colnames(ips) = mesh.coords
  if ( is.null(group) ) {
    ips= cbind(ips, distance = distance, weight = w)  
  } else {
    if (class(data)[[1]] == "dsdata") {
      ips= cbind(ips, distance = distance, weight = w, data$effort[idx[,1],group,drop=FALSE])  
    } else {
      ips= cbind(ips, distance = distance, weight = w, as.data.frame(data)[idx[,1],group,drop=FALSE])  
    }
    
  }
  
  if ( !is.null(projection)){
    if (!is.null(data$mesh.coords)) { coords = data$mesh.coords } else { coords = colnames(ips)[1:2] } # Workaround for SpatialPoints

    fn = function(x) { 
      pr =  project.weights(x, mesh, mesh.coords = data$mesh.coords)
      cbind(pr, x[rep(1, nrow(pr)),c("distance", group), drop = FALSE])
    }
    ips = recurse.rbind(fun = fn, ips, cols = c("distance", group))
    
  }
  if (length(idx[,1]) == dim(ips)[1]) { ips$idx = idx[,1] }
  return(ips)
}

# Prjection integration done FAST. Not implemented completely, has to deal with distances.

project.weights = function(ips, mesh, mesh.coords = NULL){
  
  w = ips[,"weight"]
  
  if ( is.null(mesh.coords) ) {
    ips.loc = ips
  } else {
    ips.loc = as.matrix(ips[,mesh.coords])
  }
  
  res = inla.fmesher.smorg(mesh$loc, mesh$graph$tv, points2mesh = ips.loc)
  tri = res$p2m.t 
  
  # Filter integration points outside mesh
  if ( any(tri==0) ) {
    ips = ips[!tri==0,]
    w = w[!tri==0]
    res$p2m.b = res$p2m.b[!tri==0,]
    warning("Integration points outside mesh. Omitting.")
  }
  
  nw = w * res$p2m.b
  w.by = by(as.vector(nw), as.vector(mesh$graph$tv[tri,]), sum, simplify = TRUE)
  nw = as.vector(w.by)
  
  if ( is.null(mesh.coords)) {
    ret = data.frame(mesh$loc[as.numeric(names(w.by))], w = nw)
  } else {
    nips = mesh$loc[as.numeric(names(w.by)),1:length(mesh.coords)]
    colnames(nips) = mesh.coords
    ret = cbind(nips, weight = nw)
  }
  
  return(ret)
  
}


# Expand integration points over additional dimensions
# 
# See int.quadrature on how to formulate the further arguments.
#
# @aliases int.expand
# @export
# @param ips Integration points to expand
# @param ... lists with arguments for multiple calls of int.quadrature
# @return Integration points
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

int.expand = function(ips, ...) {
  dots = list(...)
  x = list()
  weight = list()
  for (k in 1: length(dots)){
    if (is.numeric(dots[[k]])) {
      weight[[k]] = rep(1,length(dots[[k]]))
    }
    else {
      quad = do.call(int.quadrature,dots[[k]])
      x[[k]] = quad$ips
      weight[[k]] = quad$wl
      names(x)[[k]] = names(dots)[[k]]
    }
  }
  if (length(dots)>1){ 
    grd = do.call(expand.grid, x) 
    w = apply(do.call(expand.grid, weight), MARGIN = 1, prod)
  }
  else { 
    grd = data.frame(x[[1]])
    names(grd) = names(x)[[1]]
    w = weight[[1]]
  }
  rips = ips[as.vector(t(matrix(1:nrow(ips),ncol = nrow(grd),nrow = nrow(ips)))) , ]
  rgrd = grd[rep(seq_len(nrow(grd)), nrow(ips)), , drop = FALSE]
  rips$weight = rips$weight * w[rep(seq_len(nrow(grd)), nrow(ips))]
  return(cbind(rips[,!(colnames(rips)=="weight"),drop=FALSE],rgrd,weight=rips$weight))
}

