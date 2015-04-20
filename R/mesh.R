#' Query if a point is inside the mesh boundary
#'
#'
#' @aliases is.inside
#' @export
#' @param mesh an inla.mesh object
#' @param points points to query
#' @param mesh.coords Coordinate names of the mesh. Use only if loc is a data.frame with respective column names.
#' @return inside Boolean, TRUE if inside
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

is.inside = function(mesh, loc, mesh.coords = NULL) {
  if (!is.null(mesh.coords) & is.data.frame(loc)) { loc = as.matrix(loc[,mesh.coords,drop=FALSE])}
  p2m = inla.fmesher.smorg(loc=mesh$loc,tv=mesh$graph$tv,points2mesh=loc)
  return(!(p2m$p2m.t == 0))
}

#' Query if a point is inside a polygon AND inside the mesh;
#'
#'
#' @aliases is.in.polygon
#' @export
#' @param mesh an inla.mesh object
#' @param ploc Points defining a polygon
#' @param loc Points to quer
#' @param mask.mesh Mask points outside mesh, default: TRUE
#' @param mesh.coords Coordinate names of the mesh. Use only if loc is a data.frame with respective column names.
#' @return inside Boolean, TRUE if inside polygon
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

is.inside.polygon = function(mesh, ploc, loc, mesh.coords = NULL, mask.mesh = TRUE ) {
  if (!is.null(mesh.coords) & is.data.frame(loc)) { loc = as.matrix(loc[,mesh.coords,drop=FALSE])}
  if (!is.null(mesh.coords) & is.data.frame(ploc)) { ploc = as.matrix(ploc[,mesh.coords,drop=FALSE])}
  require(sp)
  mask <- point.in.polygon(loc[,1], loc[,2],ploc[,1], ploc[,2]) > 0
  if (mask.mesh){
    mask2 = is.inside(mesh, loc)
    return(mask & mask2)
  } else {
    return(mask)
  }
}