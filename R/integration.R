#' Split lines at mesh edges
#'
#' @aliases split_lines
#' @export
#' @param mesh An inla.mesh object
#' @param sp Start points of lines
#' @param ep End points of lines
#' @param filter.zero.length Filter out segments with zero length? (Bool)
#' @param ... argments to int.quadrature
#' @return List of start and end points resulting from splitting the given lines
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' @keywords internal
split_lines = function(mesh, sp, ep, filter.zero.length = TRUE) {

  # locations for splitting
  loc = as.matrix(rbind(sp,ep))
  idx = 1:dim(sp)[1]

  # Filter out segments not on the mesh
  t1 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(sp,z=0)))$p2m.t
  t2 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(ep,z=0)))$p2m.t
  # if (any(t1==0) | any(t2==0)) { warning("points outside boundary! filtering...")}
  sp = sp[!((t1==0) | (t2==0)),]
  ep = ep[!((t1==0) | (t2==0)),]
  idx = idx[!((t1==0) | (t2==0))]
  loc = as.matrix(rbind(sp,ep))

  # Split them segments into parts
  if ( dim(loc)[2] == 2 ) {loc = cbind(loc,rep(0,dim(loc)[1]))}
  np = dim(sp)[1]
  sp.idx = t(rbind(1:np,np+1:np))
  splt = INLA::inla.fmesher.smorg(mesh$loc,mesh$graph$tv, splitlines=list(loc=loc, idx=sp.idx))
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
    len = function(x) {abs(x)}
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
      else if (geometry == "geo"){ stop("Geometry geo not supported") }

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




int.slines = function(data, mesh, group = NULL, project = TRUE) {
  
  # Extract start and end coordinates
  qq = coordinates(data)
  sp = do.call(rbind, lapply(qq, function(k) do.call(rbind, lapply(k, function(x) x[1:(nrow(x)-1),]))))
  ep = do.call(rbind, lapply(qq, function(k) do.call(rbind, lapply(k, function(x) x[2:(nrow(x)),]))))
  
  idx = do.call(rbind, lapply(1:length(qq), function(k) do.call(cbind, lapply(qq[[k]], function(x) rep(k, nrow(x)-1) ))))
  idx = cbind(idx, idx)
  
  if ( !is.null(mesh) ) {
    # Filter out points outside the mesh...
    loc = as.matrix(rbind(sp,ep))
    t1 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(sp,z=0)))$p2m.t
    t2 = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=as.matrix(data.frame(ep,z=0)))$p2m.t
    if (any(t1==0) | any(t2==0)) { warning("Found spatial lines with start or end point ouside of the mesh. Omitting.")}
    sp = sp[!((t1==0) | (t2==0)),]
    ep = ep[!((t1==0) | (t2==0)),]
    idx = idx[!((t1==0) | (t2==0)),]
  
    # Split at mesh edges
    line.spl = split_lines(mesh, sp, ep, TRUE)
    sp = line.spl$sp
    ep = line.spl$ep
    idx = idx[line.spl$split.origin,]
    
  }
  
  # Determine integration points along lines
  
  sp3d = within(data.frame(sp), Z <- 0)
  colnames(sp3d) = c("X1","X2","Z")
  coordinates(sp3d) = c("X1","X2","Z")
  proj4string(sp3d) = proj4string(data)
  ep3d = within(data.frame(ep), Z <- 0)
  colnames(ep3d) = c("X1","X2","Z")
  coordinates(ep3d) = c("X1","X2","Z")
  proj4string(ep3d) = proj4string(data)
  
  sp3d = spTransform(sp3d, CRSobj = CRS("+proj=geocent +ellps=sphere +R=1.00"))
  ep3d = spTransform(ep3d, CRSobj = CRS("+proj=geocent +ellps=sphere +R=1.00"))
  mp3d = SpatialPoints((coordinates(sp3d) + coordinates(ep3d))/2, proj4string = CRS("+proj=geocent +ellps=sphere +R=1.00"))
  
  ips = coordinates(spTransform(mp3d, CRS(proj4string(data))))[,c(1,2)]
  w = spDists(coordinates(spTransform(sp3d, CRSobj = CRS("+proj=longlat")))[,c(1,2)],
              coordinates(spTransform(ep3d, CRSobj = CRS("+proj=longlat")))[,c(1,2)],
              diagonal = TRUE, longlat = TRUE)
  
  # Wrap everything up and perform projection according to distance and given group argument
  ips = data.frame(ips)
  colnames(ips) = c("x","y")
  
  # Weights
  ips = cbind(ips, weight = w)  
  if ( "weight" %in% names(data) ) { ips$weight = ips$weight * data$weight[idx[,1]] }
  
  
  if ( !is.null(group) ) { ips = cbind(ips, as.data.frame(data)[idx[,1],group,drop=FALSE]) }
  
  coordinates(ips) = c("x","y")
  if (!is.null(coordnames(data))) coordnames(ips) = coordnames(data)
  proj4string(ips) = proj4string(data)
  
  # Project to mesh vertices
  if ( project & !is.null(mesh) ) ips = vertex.projection(ips, mesh, columns = "weight", group = group)
  
  ips
}



join_segm <- function(...) {
  segm_list <- list(...)
  loc <- matrix(0, 0, 3)
  idx <- matrix(0, 0, 2)
  for (k in seq_along(segm_list)) {
    idx <- rbind(idx, segm_list[[k]]$idx + nrow(loc))
    loc <- rbind(loc, segm_list[[k]]$loc)
  }
  
  # Collapse duplicate points
  new_loc <- loc
  new_idx <- seq_len(nrow(loc))
  prev_idx <- 0
  for (k in seq_len(nrow(loc))) {
    if (any(is.na(new_loc[k, ]))) {
      new_idx[k] <- NA
    } else {
      if (prev_idx == 0) {
        prev_dist <- 1
      } else {
        prev_dist <- ((new_loc[seq_len(prev_idx), 1] - new_loc[k, 1])^2 +
                        (new_loc[seq_len(prev_idx), 2] - new_loc[k, 2])^2 +
                        (new_loc[seq_len(prev_idx), 3] - new_loc[k, 3])^2)^0.5
      }
      if (all(prev_dist > 0)) {
        prev_idx <- prev_idx + 1
        new_idx[k] <- prev_idx
        new_loc[prev_idx, ] <- new_loc[k, ]
      } else {
        new_idx[k] <- which.min(prev_dist)
      }
    }
  }
  idx <- matrix(new_idx[idx], nrow(idx), 2)
  # Remove NA and atomic lines
  ok <-
    !is.na(idx[, 1]) &
    !is.na(idx[, 2]) &
    idx[, 1] != idx[, 2]
  idx <- idx[ok, , drop=FALSE]
  # Set locations
  loc <- new_loc[seq_len(prev_idx), , drop = FALSE]
  
  INLA::inla.mesh.segment(loc = loc,
                          idx = idx,
                          is.bnd = FALSE)
}



#' Construct the intersection mesh of a mesh and a polygon
#'
#' @param mesh \code{inla.mesh} object to be intersected
#' @param poly \code{inla.mesh.segment} object with a closed polygon
#'   to intersect with the mesh 
#' @author Finn Lindgren <\email{finn.lindgren@@gmail.com}>
#' @keywords internal
intersection_mesh <- function(mesh, poly) {
  if (ncol(poly$loc) < 3) {
    poly$loc <- cbind(poly$loc, 0)
  }
  
  all_edges <- INLA::inla.mesh.segment(loc = mesh$loc,
                                 idx = cbind(as.vector(t(mesh$graph$tv)),
                                             as.vector(t(mesh$graph$tv[, c(2, 3, 1), drop=FALSE]))),
                                 is.bnd = FALSE)
  
  mesh_cover <- INLA::inla.mesh.create(loc = rbind(mesh$loc, poly$loc),
                                 interior = c(list(all_edges)))
  
  split <- INLA::inla.fmesher.smorg(mesh_cover$loc,
                              mesh_cover$graph$tv,
                              splitlines = list(loc = poly$loc,
                                                idx = poly$idx))
  split_segm <- INLA::inla.mesh.segment(loc = split$split.loc,
                                        idx = split$split.idx,
                                        is.bnd = FALSE)
  
  joint_segm <- join_segm(split_segm, all_edges)
  
  mesh_joint_cover <- INLA::inla.mesh.create(interior = list(joint_segm),
                                             extend = TRUE)
  
  mesh_poly <- INLA::inla.mesh.create(boundary = poly)
  
  loc_tri <-
    (mesh_joint_cover$loc[mesh_joint_cover$graph$tv[, 1], , drop=FALSE] +
       mesh_joint_cover$loc[mesh_joint_cover$graph$tv[, 2], , drop=FALSE] +
       mesh_joint_cover$loc[mesh_joint_cover$graph$tv[, 3], , drop=FALSE]) / 3
  ok_tri <-
    INLA::inla.mesh.projector(mesh, loc = loc_tri)$proj$ok &
    INLA::inla.mesh.projector(mesh_poly, loc = loc_tri)$proj$ok
  if (any(ok_tri)) {
    loc_subset <- unique(sort(as.vector(mesh_joint_cover$graph$tv[ok_tri, , drop=FALSE])))
    new_idx <- integer(mesh$n)
    new_idx[loc_subset] = seq_along(loc_subset)
    tv_subset <- matrix(new_idx[mesh_joint_cover$graph$tv[ok_tri, , drop=FALSE]],
                        ncol = 3)
    loc_subset <- mesh_joint_cover$loc[loc_subset, , drop=FALSE]
    mesh_subset <- INLA::inla.mesh.create(loc = loc_subset,
                                          tv = tv_subset,
                                          extend = FALSE)
  } else {
    mesh_subset <- NULL
  }
  
  mesh_subset
}

#' Safe(?) way to construct integration weights for mesh/polygon intersections
#'
#' @param mesh Mesh on which to integrate
#' @param inter \code{list} of \code{loc}, integration points,
#'   and \code{weight}, integration weights
#' @author Finn Lindgren <\email{finn.lindgren@@gmail.com}>
#' @keywords internal
integration_weight_construction <- function(mesh, inter) {
  # Safe mesh for point lookups
  all_edges <- INLA::inla.mesh.segment(loc = mesh$loc,
                                       idx = cbind(as.vector(t(mesh$graph$tv)),
                                                   as.vector(t(mesh$graph$tv[, c(2, 3, 1), drop=FALSE]))),
                                       is.bnd = FALSE)
  mesh_cover <- INLA::inla.mesh.create(loc = mesh$loc,
                                       interior = list(all_edges),
                                       extend = TRUE)
  
  proj_cover <- INLA::inla.mesh.projector(mesh_cover,
                                          loc = inter$loc)
  # Remove points not in the original 'mesh'
  A <- proj_cover$proj$A[, mesh_cover$idx$loc, drop = FALSE]
  # Renormalise
  A <- A / Matrix::rowSums(A)

  # Convert integration weights from intersection mesh points to original mesh points
  weight <- as.vector(as.vector(inter$weight) %*% A)
  
  list(loc = mesh$loc, weight = weight)
}


#' Basic robust integration weights for mesh/polygon intersections
#'
#' @param mesh Mesh on which to integrate
#' @param bnd \code{inla.mesh.segment} defining the integration domain
#' @param nsub number of subdivision points along each triangle edge, giving
#'    \code{(nsub + 1)^2} proto-integration points used to compute
#'   the vertex weights
#'   (default \code{9}, giving 100 proto-integration points)
#' @return \code{list} with elements \code{loc} and \code{weight} with
#'   one integration point for each mesh vertex of triangles overlapping
#'   the integration domain
#' @author Finn Lindgren <\email{finn.lindgren@@gmail.com}>
#' @keywords internal
make_stable_integration_points <- function(mesh, bnd, nsub = NULL) {
  # Construct a barycentric grid of subdivision triangle midpoints
  if (is.null(nsub)) {
    nsub <- 9
  }
  stopifnot(nsub >= 0)
  nB <- (nsub + 1)^2

  b <- seq(1/3, 1/3 + nsub, length = nsub + 1) / (nsub + 1)
  bb <- as.matrix(expand.grid(b, b))
  # Points above the diagonal should be reflected into the lower triangle:
  refl <- rowSums(bb) > 1
  if (any(refl)) {
    bb[refl, ] = cbind(1 - bb[refl, 2], 1 - bb[refl, 1])
  }
  # Construct complete barycentric coordinates:
  barycentric_grid <- cbind(1 - rowSums(bb), bb)

  # Construct integration points
  nT <- nrow(mesh$graph$tv)
  loc <- matrix(0.0, nT*nB, 3)
  idx_end <- 0
  for (tri in seq_len(nT)) {
    idx_start <- idx_end + 1
    idx_end <- idx_start + nB - 1
    loc[seq(idx_start, idx_end, length = nB), ] <-
      as.matrix(barycentric_grid %*%
                  mesh$loc[mesh$graph$tv[tri, ], , drop = FALSE])
  }
  
  # Construct integration weights
  weight <- rep(INLA::inla.mesh.fem(mesh, order = 1)$ta / nB, each = nB)
  
  # Filter away point outside integration domain boundary:
  mesh_bnd <- INLA::inla.mesh.create(boundary = bnd)
  ok <- INLA::inla.mesh.projector(mesh_bnd, loc = loc)$proj$ok
  
  list(loc = loc[ok, , drop = FALSE],
       weight = weight[ok])
}

#' Integration points for polygons inside an inla.mesh
#' 
#' @aliases int.polygon
#' @export
#' @param mesh An inla.mesh object
#' @param loc Locations defining the polygons
#' @param group If loc defines multiple polygons then this is the ID of the group for each location in loc
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}> and Finn Lindgren <\email{finn.lindgren@@gmail.com}>
#' @keywords internal

int.polygon = function(mesh, loc, group = NULL){
  
  if ( is.null(group) ) { group = rep(1, nrow(loc)) }
  ipsl <- list()
  # print(paste0("Number of polygons to integrate over: ", length(unique(group)) ))
  for ( g in unique(group) ) {
    gloc <- loc[group == g, , drop=FALSE]

    # Combine polygon with mesh boundary to get mesh covering the intersection.
    bnd <- INLA::inla.mesh.segment(loc = gloc, is.bnd = TRUE)

    method <- "stable"
    inter <- list(loc = matrix(0, 0, 3), weight = numeric(0))
    if (method == "stable") {
      inter <- make_stable_integration_points(mesh, bnd)
    } else {
      # Create intersection mesh
      mesh_inter <- intersection_mesh(mesh, bnd)
      if (!is.null(mesh_inter)) {
        inter <- list(loc = mesh_inter$loc,
                      weight = INLA::inla.mesh.fem(mesh_inter, order = 1)$va)
      }
    }
    
    # Compute integration locations and weights
    loc_weight <- integration_weight_construction(mesh, inter)

    ok <-
      INLA::inla.mesh.project(mesh, loc_weight$loc)$ok &
      (loc_weight$weight > 0)

    ips <- data.frame(loc_weight$loc[ok, 1:2, drop=FALSE])
    colnames(ips) <- c("x","y")
    ips$weight <- loc_weight$weight[ok]

    ips$group <- g
    ipsl <- c(ipsl, list(ips))
  }
  
  do.call(rbind, ipsl)
}



