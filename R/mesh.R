# Needed for registerin S4 methods, e.g. vertices()
setClass("inla.mesh")


# GENERICS
# refine = function(...){UseMethod("refine")}
# tsplit = function(...){UseMethod("tsplit")}


#' Query if a point is inside the mesh boundary
#'
#'
#' @aliases is.inside
#' @export
#' @keywords internal
#' @param mesh an inla.mesh object.
#' @param loc Points in space stored either as data.frame, a two-column matrix of x and y coordinates or a SpatialPoints object.
#' @param mesh.coords Coordinate names of the mesh. Use only if loc is a data.frame with respective column names.
#' @return Single column matrix of Boolean values indicating if a point is inside the mesh.
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 
#' @examples 
#' \donttest{
#' if (require("INLA", quietly = TRUE)) {
#' # Load Gorilla data
#' 
#' data("gorillas", package = "inlabru")
#' 
#' # Check if all Gorilla nests are inside the mesh
#' 
#' all(is.inside(gorillas$mesh, gorillas$nests))
#' 
#' # Also works for locations not stored as SpatialPoints object
#' 
#' loc = coordinates(gorillas$nests)
#' all(is.inside(gorillas$mesh, loc))
#' }
#' }

is.inside = function(mesh, loc, mesh.coords = NULL) {
  if ( inherits(loc, "Spatial") ) { loc = coordinates(loc) }
  if (!is.null(mesh.coords) & is.data.frame(loc)) { loc = as.matrix(loc[,mesh.coords,drop=FALSE])}
  p2m = INLA::inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=loc)
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

#' Vertices
#' 
#' This is a generic function. The outcome depends on the \code{object} provided
#' 
#' @name vertices
#' @exportMethod vertices
#' @param object An object for which to call the particular vertices method.
#' @return The form of the value returned by vertices() depends on the class of its argument. See the documentation of the particular methods for details of what is produced by that method.

setGeneric("vertices", valueClass = "SpatialPointsDataFrame", function(object) {
  standardGeneric("vertices")
})

#' Vertices
#' 
#' @rdname vertices
setMethod("vertices", signature("inla.mesh"), function(object) vertices.inla.mesh(object))


#' @title Extract vertex locations from an \code{inla.mesh}
#'
#' @description Converts the vertices of an \code{inla.mesh} object into a \code{SpatialPointsDataFrame}.
#' 
#' @aliases vertices.inla.mesh
#' @export
#' @param object An \code{inla.mesh} object.
#' @return A SpatialPointsDataFrame of mesh vertex locations. The \code{vrt} column indicates the internal vertex id.
#' 
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#' 
#' @examples
#' \donttest{
#' data("mrsea")
#' vrt = vertices(mrsea$mesh)
#' ggplot() + gg(mrsea$mesh) + gg(vrt, color = "red")
#' }

vertices.inla.mesh = function(object) {
  
  if ( is.null(object$crs) ) { object$crs = CRS("")  }

  vrt = data.frame(object$loc)
  if (! is.null(colnames(vrt))) { colnames(vrt) = c("x","y","z") }
  if ( all(vrt[, 3]==0) ) { vrt = vrt[,1:2] }
  coordinates(vrt) = colnames(vrt)
  vrt = SpatialPointsDataFrame(vrt, data = data.frame(vertex = 1:nrow(object$loc)))
  proj4string(vrt) = object$crs
  
  vrt  # return
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
#' \donttest{
#' data("mrsea")
#' pxl = pixels(mrsea$mesh, nx = 50, ny = 50)
#' ggplot() + gg(pxl) + gg(mrsea$mesh)
#' }

pixels = function(mesh, nx = 150, ny = 150, mask = TRUE) {
  if ( length(nx)==1 ){
    x = seq(min(mesh$loc[,1]), max(mesh$loc[,1]),length = nx)
  } else { x = nx }
  if ( length(ny)==1 ){
    y = seq(min(mesh$loc[,2]), max(mesh$loc[,2]),length = ny)
  } else { y = ny }
  
  lattice <- INLA::inla.mesh.lattice(x=x,y=y)
  pixels <- data.frame(x=lattice$loc[,1], y=lattice$loc[,2])

  coordinates(pixels) <- c("x","y")
  
  if (!is.null(mesh$crs)) { pixels = SpatialPixels(pixels, proj4string = mesh$crs) }
  else { pixels = SpatialPixels(pixels) }
  
  if ( is.logical(mask) ) {
    if ( mask ) { pixels = pixels[is.inside(mesh, coordinates(pixels))] }
  } else {
    pixels = pixels[!is.na(sp::over(pixels, mask)),]
  }
  pixels
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
#' @keywords internal
#' 
#' @param mesh an inla.mesh object
#' @param refine A list of refinement options passed on to \link[INLA]{inla.mesh.create}
#' @return mesh A refined inla.mesh object
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

refine.inla.mesh = function(mesh, refine = list(max.edge=1)){
  rmesh = INLA::inla.mesh.create(loc=mesh$loc,interior=INLA::inla.mesh.interior(mesh),boundary=INLA::inla.mesh.boundary(mesh),refine=refine)
  return(rmesh)
}

#' Split triangles of a mesh into four triangles
#'
#' Warning: does not reconstruct interior boundary
#' Warning2: Works in euclidean coordinates. Not suitable for sphere.
#'
#' @aliases tsplit.inla.mesh
#' @keywords internal
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


  mesh2 = INLA::inla.mesh.create(loc = all.loc, boundary = all.bnd )

  if (n == 1) { return(mesh2) }
  else { return(tsplit.inla.mesh(mesh2,n-1))}
}


