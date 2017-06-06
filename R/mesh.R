# GENERICS
refine = function(...){UseMethod("refine")}
tsplit = function(...){UseMethod("tsplit")}

# Plot an inla.mesh using ggplot
#
# @aliases ggp.mesh
# @export
# @param mesh an inla.mesh object
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

ggp.mesh = function(mesh, col = NULL, nx = 400, add = NULL, mcol = "black") {

  xlim = range(mesh$loc[,1])
  ylim = range(mesh$loc[,2])
  aspect = diff(ylim)/diff(xlim)
  ny = round(nx * aspect)
  proj <- inla.mesh.projector(mesh, dims = c(nx,ny))
  x = seq(xlim[1], xlim[2], length.out = nx)
  y = seq(ylim[1], ylim[2], length.out = ny)
  grid = expand.grid(x=x,y=y)
  
  
  if ( !is.null(col) ) {
    mcol = "green"
    loc = data.frame(cbind(grid$x,grid$y))
    col = INLA::inla.spde.make.A(mesh, loc = as.matrix(loc)) %*% as.vector(col)
    msk = is.inside(mesh,as.matrix(loc))
    # col[msk] = NA
    df = data.frame(grid, col=as.vector(col), alpha = msk)
    gg = ggplot(df, aes(x=x,y=y) )
    gg = gg + geom_raster(aes(fill = col, alpha = alpha), hjust=0.5, vjust=0.5, interpolate = TRUE)
    gg = gg + scale_alpha_discrete(guide = 'none')
    gg = gg + theme(legend.title=element_blank())
    
  } else {
    df = data.frame(grid)
    if (!is.null(add)) { gg = add } else { gg = ggplot(df, aes(x=x,y=y) )}
  }
  
  # Plot mesh lines
  gg = gg + geom_segment(data = data.frame(a=mesh$loc[mesh$graph$tv[,1],c(1,2)],b=mesh$loc[mesh$graph$tv[,2],c(1,2)]), 
                         aes_string(x="a.1",y="a.2",xend="b.1",yend="b.2"), color = mcol)
  
  gg = gg + geom_segment(data = data.frame(a=mesh$loc[mesh$graph$tv[,2],c(1,2)],b=mesh$loc[mesh$graph$tv[,3],c(1,2)]),
                         aes_string(x="a.1",y="a.2",xend="b.1",yend="b.2"), color = mcol)
  
  gg = gg + geom_segment(data = data.frame(a=mesh$loc[mesh$graph$tv[,1],c(1,2)],b=mesh$loc[mesh$graph$tv[,3],c(1,2)]),
                         aes_string(x="a.1",y="a.2",xend="b.1",yend="b.2"), color = mcol)
  
  gg = gg + coord_fixed()

  return(gg)
}




#' Query if a point is inside the mesh boundary
#'
#'
#' @aliases is.inside
#' @export
#' @param mesh an inla.mesh object
#' @param loc points to query
#' @param mesh.coords Coordinate names of the mesh. Use only if loc is a data.frame with respective column names.
#' @return inside Boolean, TRUE if inside
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

is.inside = function(mesh, loc, mesh.coords = NULL) {
  if (!is.null(mesh.coords) & is.data.frame(loc)) { loc = as.matrix(loc[,mesh.coords,drop=FALSE])}
  p2m = inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=loc)
  return(!(p2m$p2m.t == 0))
}

# Query if a point is inside a polygon AND inside the mesh;
#
#
# @aliases is.in.polygon
# @export
# @param mesh an inla.mesh object
# @param ploc Points defining a polygon
# @param loc Points to quer
# @param mask.mesh Mask points outside mesh, default: TRUE
# @param mesh.coords Coordinate names of the mesh. Use only if loc is a data.frame with respective column names.
# @return inside Boolean, TRUE if inside polygon
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

is.inside.polygon = function(mesh, ploc, loc, mesh.coords = NULL, mask.mesh = TRUE ) {
  if (!is.null(mesh.coords) & is.data.frame(loc)) { loc = as.matrix(loc[,mesh.coords,drop=FALSE])}
  if (!is.null(mesh.coords) & is.data.frame(ploc)) { ploc = as.matrix(ploc[,mesh.coords,drop=FALSE])}
  
  mask <- sp::point.in.polygon(loc[,1], loc[,2],ploc[,1], ploc[,2]) > 0
  if (mask.mesh){
    mask2 = is.inside(mesh, loc)
    return(mask & mask2)
  } else {
    return(mask)
  }
}


#' Plot an inla.mesh using ggplot
#'
#' @export
vertices = function(...){UseMethod("vertices")}



#' @title Extract vertex locations from an \code{inla.mesh}
#'
#' @description Converts the vertices of an \code{inla.mesh} object into a \code{SpatialPointsDataFrame}.
#' 
#' @aliases vertices.inla.mesh
#' @export vertices.inla.mesh
#' 
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 
#' @param mesh An \code{inla.mesh} object
#' @return A \code{SpatialPointsDataFrame} of vertex locations
#' 
#' @examples
#' 
#' data("mrsea")
#' vrt = vertices(mrsea$mesh)
#' ggplot() + gg(mrsea$mesh) + gg(vrt, color = "red")
#' 

vertices.inla.mesh = function(mesh) {
  if (is.null(mesh$crs)) {
    mesh$loc
  } else {
    if (any(!(mesh$loc[,3]==0))) { vrt = mesh$loc } else { vrt = mesh$loc[,c(1,2)] }
    SpatialPointsDataFrame(vrt, proj4string = mesh$crs, data = data.frame(vertex = 1:nrow(mesh$loc)))
  }
  
}


#' @title Generate \code{SpatialPixels} covering an \code{inla.mesh}
#'
#' @description Generate \code{SpatialPixels} covering an \code{inla.mesh}
#' 
#' @aliases pixels
#' @export
#' 
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 
#' @param mesh An \code{inla.mesh} object
#' @param nx Number of pixels in x direction
#' @param ny Number of pixels in y direction
#' @param mask If logical and TRUE, remove pixels that are outside the mesh. If \code{mask} is a \code{Spatial} object, only return pixels covered by this object. 
#' @return \code{SpatialPixels} covering the mesh
#' 
#' @examples
#' 
#' data("mrsea")
#' pxl = pixels(mrsea$mesh, nx = 50, ny = 50)
#' ggplot() + gg(pxl) + gg(mrsea$mesh)
#' 

pixels = function(mesh, nx = 150, ny = 150, mask = TRUE) {
  lattice <- inla.mesh.lattice(x=seq(min(mesh$loc[,1]), max(mesh$loc[,1]),length = nx),
                               y=seq(min(mesh$loc[,2]), max(mesh$loc[,2]),length = ny))
  pixels <- data.frame(x=lattice$loc[,1], y=lattice$loc[,2])

  coordinates(pixels) <- c("x","y")
  
  if (!is.null(mesh$crs)) { pixels = SpatialPixels(pixels, proj4string = mesh$crs) }
  else { pixels = SpatialPixels(pixels) }
  
  if ( is.logical(mask) && (mask == TRUE) ){ 
    pixels = pixels[is.inside(mesh, coordinates(pixels))] 
  } else if ( !is.null(mask) ) {
    pixels = pixels[as.vector(!is.na(over(pixels, mask))), ]
  }
}



# Triangle indices of points given a mesh
#
# @aliases triangle
# @export
# @param mesh A inla.mesh
# @param loc Locations using the coordinate system of the mesh 
# @return tri Triangle indices
# @examples \\dontrun{}
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

triangle = function(mesh,loc){
  mcross = function(a,b) { return( (a[,1])*(b[2]) - (a[,2])*(b[1]) )}
  tri = numeric(length=dim(loc)[1])
  tv = mesh$graph$tv  
  for (j in 1:nrow(tv)) {
    v = mesh$loc[tv[j,],c(1,2)]
    
    a = v[1,]
    b = v[2,]
    c = v[3,]
    
    ap =  loc - cbind(rep(a[1],nrow(loc)),rep(a[2],nrow(loc)))
    bp =  loc - cbind(rep(b[1],nrow(loc)),rep(b[2],nrow(loc)))
    cp =  loc - cbind(rep(c[1],nrow(loc)),rep(c[2],nrow(loc)))
    
    ab = a-b
    bc = b-c
    ca = c-a
    
    c1 = mcross(ap,ab)
    c2 = mcross(bp,bc)
    c3 = mcross(cp,ca)
    
    # AP x AB, BP x BC, and CP x CA must have same sign
    inside = which( ((sign(c1) == -1) & (sign(c2)==-1) & (sign(c3)==-1)) | ((sign(c1) == 1) & (sign(c2)==1) & (sign(c3)==1)))
    
    tri[inside] = j
  }
  return(tri)
}




#' Refine an inla.mesh object
#'
#'
#' @aliases refine.inla.mesh
#' @export
#' @param mesh an inla.mesh object
#' @param refine A list of refinement options passed on to \link{inla.mesh.create}
#' @return mesh A refined inla.mesh object
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

refine.inla.mesh = function(mesh, refine = list(max.edge=1)){
  rmesh = inla.mesh.create(loc=mesh$loc,interior=inla.mesh.interior(mesh),boundary=inla.mesh.boundary(mesh),refine=refine)
  return(rmesh)
}

#' Split triangles of a mesh into four triangles
#'
#' Warning: does not reconstruct interior boundary
#' Warning2: Works in euclidean coordinates. Not suitable for sphere.
#'
#' @aliases tsplit.inla.mesh
#' @export
#' @param mesh an inla.mesh object
#' @param n number of splitting recursions
#' @return mesh A refined inla.mesh object
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

tsplit.inla.mesh = function(mesh, n = 1){
  
  n = 1
  
  p1 = mesh$loc[mesh$graph$tv[,1],]
  p2 = mesh$loc[mesh$graph$tv[,2],]
  p3 = mesh$loc[mesh$graph$tv[,3],]
  
  m1 = p1 + 0.5*(p2-p1)
  m2 = p1 + 0.5*(p3-p1)
  m3 = p2 + 0.5*(p3-p2)
  all.loc = rbind(mesh$loc,m1,m2,m3)
  
  bnd.mid = mesh$loc[mesh$segm$bnd$idx[,1],] + 0.5 * ( mesh$loc[mesh$segm$bnd$idx[,2],] - mesh$loc[mesh$segm$bnd$idx[,1],]  )
  all.bnd = rbind(mesh$segm$bnd$loc,bnd.mid)
  
  
  mesh2 = inla.mesh.create(loc = all.loc, boundary = all.bnd )
  
  if (n == 1) { return(mesh2) }
  else { return(tsplit.inla.mesh(mesh2,n-1))}
}





# Generate a simple default mesh
#
# @aliases default.mesh
# @export
# @param spObject A Spatial* object
# @param max.edge A parameter passed on to \link{inla.mesh.2d} which controls the granularity of the mesh. If NULL, 1/20 of the domain size is used.
# @return An \code{inla.mesh} object

default.mesh = function(spObject, max.edge = NULL, convex = -0.15){
  if (inherits(spObject, "SpatialPoints")) {
    x = c(bbox(spObject)[1,1], bbox(spObject)[1,2], bbox(spObject)[1,2], bbox(spObject)[1,1])
    y = c(bbox(spObject)[2,1], bbox(spObject)[2,1], bbox(spObject)[2,2], bbox(spObject)[2,2])
    # bnd = inla.mesh.segment(loc = cbind(x,y))
    # mesh = inla.mesh.2d(interior = bnd, max.edge = diff(bbox(spObject)[1,])/10)
    if ( is.null(max.edge) ) { max.edge = max.edge = diff(bbox(spObject)[1,])/20 }
    hull = inla.nonconvex.hull(points = coordinates(spObject), convex = convex)
    mesh = inla.mesh.2d(boundary = hull, max.edge = max.edge)
  } else {
    NULL
  }
}
