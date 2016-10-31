#' Plot an inla.mesh using ggplot
#'
#' @aliases ggp.mesh
#' @export
#' @param mesh an inla.mesh object
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

ggp.mesh = function(mesh, col = NULL, nx = 400, add = NULL, mcol = rgb(0,0,0,0.3)) {
  if ( !require(ggplot2) ) { stop("This function requires the ggplot2 package to run.") }
  xlim = range(mesh$loc[,1])
  ylim = range(mesh$loc[,2])
  aspect = diff(ylim)/diff(xlim)
  ny = round(nx * aspect)
  proj <- inla.mesh.projector(mesh, dims = c(nx,ny))
  x = seq(xlim[1], xlim[2], length.out = nx)
  y = seq(ylim[1], ylim[2], length.out = ny)
  grid = expand.grid(x=x,y=y)
  
  
  if ( !is.null(col) ) {
    mcol = rgb(0,0,0,0.1)
    loc = data.frame(cbind(grid$x,grid$y))
    col = inla.spde.make.A(mesh, loc = as.matrix(loc)) %*% as.vector(col)
    msk = is.inside(mesh,as.matrix(loc))
    # col[msk] = NA
    df = data.frame(grid, col=as.vector(col), alpha = msk)
    gg = ggplot(df, aes(x=x,y=y) )
    gg = gg + geom_raster(aes(fill = col, alpha = alpha), hjust=0.5, vjust=0.5, interpolate = TRUE)
    gg = gg + scale_alpha_discrete(guide = 'none')
    gg = gg + scale_fill_gradientn(colours = topo.colors(100) ) + theme(legend.title=element_blank())
    
  } else {
    df = data.frame(grid)
    if (!is.null(add)) { gg = add } else { gg = ggplot(df, aes(x=x,y=y) )}
  }
  
  # Plot mesh lines
  gg = gg + geom_segment(data = data.frame(a=mesh$loc[mesh$graph$tv[,1],c(1,2)],b=mesh$loc[mesh$graph$tv[,2],c(1,2)]), 
                         aes(x=a.1,y=a.2,xend=b.1,yend=b.2), color = mcol)
  
  gg = gg + geom_segment(data = data.frame(a=mesh$loc[mesh$graph$tv[,2],c(1,2)],b=mesh$loc[mesh$graph$tv[,3],c(1,2)]),
                         aes(x=a.1,y=a.2,xend=b.1,yend=b.2), color = mcol)
  
  gg = gg + geom_segment(data = data.frame(a=mesh$loc[mesh$graph$tv[,1],c(1,2)],b=mesh$loc[mesh$graph$tv[,3],c(1,2)]),
                         aes(x=a.1,y=a.2,xend=b.1,yend=b.2), color = mcol)
  
  gg = gg + coord_fixed()

  return(gg)
}




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