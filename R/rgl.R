#' Plot a globe using rgl
#' 
#' Creates a textured sphere and lon/lat coordinate annotations.
#' 
#' This funciton requires the `rgl` and `sphereplot` packages.
#' 
#' @aliases globe
#' @name globe
#' @export
#' @param R Radius of the globe
#' @param R.grid Radius of the annotation sphere.
#' @param specular Light color of specular effect.
#' @param axes If TRUE, plot x, y and z axes.
#' @param box If TRUE, plot a box around the globe.
#' @param xlab,ylab,zlab Axes labels
#' 
#' @return No value, used for plotting side effect.
#' 
#' @family inlabru RGL tools
#' 
#' @example inst/examples/rgl.R

globe = function(R = 1, 
                 R.grid = 1.05,
                 specular = "black", 
                 axes = FALSE, 
                 box = FALSE, 
                 xlab = "", ylab= "", zlab = ""){
  
  
  # coordinates for texture
  n.smp = 50
  lat <- matrix(-asin(seq(-1,1, len=n.smp)), n.smp, n.smp, byrow=TRUE)
  long <- matrix(seq(-180, 180, len=n.smp)*pi/180, n.smp, n.smp)
  x <- R*cos(lat)*cos(long)
  y <- R*cos(lat)*sin(long)
  z <- R*sin(lat)
  
  # globe and texture
  requireNamespace("rgl")
  rgl::persp3d(x, y, z, col="white", 
          texture=system.file("misc/Lambert_ocean.png",package="inlabru"), 
          specular = "black", axes = axes, box = box, xlab=xlab, ylab=ylab, zlab=zlab,
          normal_x=x, normal_y=y, normal_z=z)
  
  # spheric grid
  requireNamespace("sphereplot")
  sphereplot::rgl.sphgrid(longtype = "D",add=TRUE,radius=R.grid)
  
}

#' Render Spatial* and inla.mesh objects using RGL
#' 
#' glplot is a generic function for renders various kinds of spatial objects, i.e. Spatial* data 
#' and inla.mesh objects. The function invokes particular methods which depend on the class of the first argument.
#' 
#' @aliases glplot
#' @name glplot
#' @export
#' @param object an object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' 
#' @family inlabru RGL tools
#' 
#' @example inst/examples/rgl.R

glplot = function(object, ...) { UseMethod("glplot") }

#' Visualize SpatialPoints using RGL
#' 
#' This function will calculate the cartesian coordinates of the points provided 
#' and use rgl.points() in order to render them.
#' 
#' @export
#' @name glplot.SpatialPoints
#' 
#' @param object a SpatialPoints or SpatialPointsDataFrame object.
#' @param add If TRUE, add the points to an existing plot. If FALSE, create new plot.
#' @param color vector of R color characters. See rgl.material() for details.
#' @param ... Parameters passed on to rgl.points()
#' 
#' @family inlabru RGL tools
#' 
#' @example inst/examples/rgl.R


glplot.SpatialPoints = function(object, add = TRUE, color = "red", ...) {
  
  if ( length(coordnames(object))<3 ) {
    ll = data.frame(object)
    ll$TMP.ZCOORD = 0
    coordinates(ll) = c(coordnames(object), "TMP.ZCOORD")
    proj4string(ll) = CRS(proj4string(object))
    object = ll
  }
  
  object = spTransform(object, CRSobj = CRS("+proj=geocent +ellps=sphere +R=1.00"))
  cc = coordinates(object)
  requireNamespace("rgl")
  rgl::rgl.points(x=cc[,1], y = cc[,2], z = cc[,3], add = add, color = color, ...)
  
}

#' Visualize SpatialLines using RGL
#' 
#' This function will calculate a cartesian representation of the lines provided 
#' and use rgl.linestrip() in order to render them.
#' 
#' 
#' @export
#' @name glplot.SpatialLines
#' 
#' @param object a SpatialLines or SpatialLinesDataFrame object.
#' @param add If TRUE, add the lines to an existing plot. If FALSE, create new plot.
#' @param ... Parameters passed on to rgl.linestrips().
#' 
#' @family inlabru RGL tools
#' 
#' @example inst/examples/rgl.R

glplot.SpatialLines = function(object, add = TRUE,  ...) {
  
  qq = coordinates(object)
  sp = do.call(rbind, lapply(qq, function(k) do.call(rbind, lapply(k, function(x) x[1:(nrow(x)-1),]))))
  ep = do.call(rbind, lapply(qq, function(k) do.call(rbind, lapply(k, function(x) x[2:(nrow(x)),]))))
  sp = data.frame(x = sp[,1], y = sp[,2], z = 0)
  ep = data.frame(x = ep[,1], y = ep[,2], z = 0)
  
  coordinates(sp) = c("x","y","z")
  coordinates(ep) = c("x","y","z")
  proj4string(sp) = CRS(proj4string(object))
  proj4string(ep) = CRS(proj4string(object))
  
  sp = spTransform(sp, CRSobj = CRS("+proj=geocent +ellps=sphere +R=1.00"))
  ep = spTransform(ep, CRSobj = CRS("+proj=geocent +ellps=sphere +R=1.00"))
  
  cs = coordinates(sp)
  ce = coordinates(ep)
  na = matrix(NA, ncol = 3, nrow = nrow(cs))

  mm = matrix(t(cbind(cs,ce,na)), ncol = 3, nrow = 3*nrow(ce), byrow=TRUE)

  requireNamespace("rgl")
  
  rgl::rgl.linestrips(mm, add = add, ...)
  
}


#' Visualize SpatialPoints using RGL
#' 
#' This function transforms the mesh to 3D cartesian coordinates and uses 
#' inla.plot.mesh() with \code{rgl=TRUE} to plot the result.
#' 
#' @export
#' @name glplot.inla.mesh
#' 
#' @param object an inla.mesh object.
#' @param add If TRUE, add the lines to an existing plot. If FALSE, create new plot.
#' @param col Color specification. A single named color, a vector of scalar values, or a matrix of RGB values.
#' @param ... Parameters passed on to plot.inla.mesh()
#' 
#' @family inlabru RGL tools
#'
#' @example inst/examples/rgl.R

glplot.inla.mesh = function(object, add = TRUE, col = NULL,...){
  if ( object$manifold  == "S2" ) {
    # mesh$loc = mesh$loc
  } else {
    ll = data.frame(object$loc)
    colnames(ll) = c("x","y","z")
    coordinates(ll) = c("x","y","z")
    proj4string(ll) = object$crs
    ll = spTransform(ll, CRSobj = CRS("+proj=geocent +ellps=sphere +R=1.00"))
    object$loc = coordinates(ll)
  }
  
  if (is.null(col)) { plot(object,rgl=TRUE,add=add,...) }
  else{ plot(object,rgl=TRUE,add=add,col=col,...) }
}



# Play (animate) spatial field
#
# Animates a spatial field using RGL. 
# 
# @aliases play.spatial
# @export
# @param group Example: group = list(year = c(1,2)) animates the field for years 1 and 2
# @param ... Parameters passed on to \link{plot.spatial}
#

play.spatial = function(group = list(), rgl, ...){
  if ( !rgl ) { rgl = TRUE}
  globe()
  sargs = list(...)
  myanim = function(time, ...) {
    rgl::par3d(skipRedraw = TRUE)
    grp = list()
    grp[[names(group)[[1]]]] = group[[1]][(floor(time) %% 2)+1]
    do.call(pixelplot.mesh, c(sargs, list(group = grp, add = TRUE, rgl = TRUE)))
    rgl::par3d(skipRedraw = FALSE)
    return("")
  }
  
  rgl::play3d(myanim, duration = 10, startTime = 0, fps = 1)
  
}
 