#' Plot a globe using rgl
#' 
#' @aliases globe
#' @name globe
#' @export

globe = function(R = 1, 
                 R.grid = 1.05,
                 specular = "black", 
                 axes = FALSE, 
                 box = FALSE, 
                 xlab = "", ylab= "", zlab = ""){
  
  require(sphereplot)
  
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
          specular = specular, axes = axes, box = box, xlab=xlab, ylab=ylab, zlab=zlab,
          normal_x=x, normal_y=y, normal_z=z)
  
  # spheric grid
  requireNamespace("sphereplot")
  sphereplot::rgl.sphgrid(longtype = "D",add=TRUE,radius=R.grid)
  
}

#' Plot sp objects and and meshes using RGL
#' 
#' 
#' @aliases rgl
#' @name rgl
#' @export

rgl = function(...){UseMethod("rgl")}



rgl.SpatialPoints = function(data, add = TRUE, color = "red", ...) {
  
  if ( length(coordnames(data))<3 ) {
    ll = data.frame(data)
    ll$TMP.ZCOORD = 0
    coordinates(ll) = c(coordnames(data), "TMP.ZCOORD")
    proj4string(ll) = CRS(proj4string(data))
    data = ll
  }
  
  data = spTransform(data, CRSobj = CRS("+proj=geocent +ellps=sphere +R=1.00"))
  cc = coordinates(data)
  requireNamespace("rgl")
  rgl::rgl.points(x=cc[,1], y = cc[,2], z = cc[,3], add = add, color = color, ...)
  
}

rgl.SpatialLines = function(data, add = TRUE,  ...) {
  
  qq = coordinates(data)
  sp = do.call(rbind, lapply(qq, function(k) do.call(rbind, lapply(k, function(x) x[1:(nrow(x)-1),]))))
  ep = do.call(rbind, lapply(qq, function(k) do.call(rbind, lapply(k, function(x) x[2:(nrow(x)),]))))
  sp = data.frame(x = sp[,1], y = sp[,2], z = 0)
  ep = data.frame(x = ep[,1], y = ep[,2], z = 0)
  
  coordinates(sp) = c("x","y","z")
  coordinates(ep) = c("x","y","z")
  proj4string(sp) = CRS(proj4string(data))
  proj4string(ep) = CRS(proj4string(data))
  
  sp = spTransform(sp, CRSobj = CRS("+proj=geocent +ellps=sphere +R=1.00"))
  ep = spTransform(ep, CRSobj = CRS("+proj=geocent +ellps=sphere +R=1.00"))
  
  cs = coordinates(sp)
  ce = coordinates(ep)
  na = matrix(NA, ncol = 3, nrow = nrow(cs))

  mm = matrix(t(cbind(cs,ce,na)), ncol = 3, nrow = 3*nrow(ce), byrow=TRUE)

  requireNamespace("rgl")
  
  rgl::rgl.linestrips(mm, add = add, ...)
  
}

rgl.inla.mesh = function(mesh, add = TRUE, col = NULL,...){
  if ( mesh$manifold  == "S2" ) {
    # mesh$loc = mesh$loc
  } else {
    ll = data.frame(mesh$loc)
    colnames(ll) = c("x","y","z")
    coordinates(ll) = c("x","y","z")
    proj4string(ll) = mesh$crs
    ll = spTransform(ll, CRSobj = CRS("+proj=geocent +ellps=sphere +R=1.00"))
    mesh$loc = coordinates(ll)
  }
  
  if (is.null(col)) { plot(mesh,rgl=TRUE,add=add,...) }
  else{ plot(mesh,rgl=TRUE,add=add,col=col,...) }
}

 