#' Plot linestrips on sphere using rgl
#'
#' @aliases rgl.sphlinestrips
#' @export
#' @examples \\dontrun{}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
#' 
rgl.sphlinestrips = function(lat,lon,radius=1,...){
  R=radius
  for (k in 1:(length(lat)-1)){
    x.start = geo.to.euc(data.frame(lat=lat[k],lon=lon[k]),R=1)
    x.end = geo.to.euc(data.frame(lat=lat[k+1],lon=lon[k+1]),R=1)
    n.pts = min(10,max(3,round(1000*dist.euc(x.start,x.end)/R)))
    xseq = t(sapply(seq(0,1,length.out=n.pts),FUN=function(a){t(x.start+a*(x.end-x.start))}))
    xseq = radius*normalize.euc(xseq)
    rgl.linestrips(xseq[,1],xseq[,2],xseq[,3],...)
  }
}

#' Plot lines on sphere using rgl
#'
#' @aliases rgl.sphlines
#' @export
#' @examples \\dontrun{}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
#' 
rgl.sphlines = function(lat,lon,radius=1,smooth=FALSE,...){
  if (smooth) {
    for (k in seq(1,(length(lat)-1),by=2)){
      x.start = geo.to.euc(data.frame(lat=lat[k],lon=lon[k]),R=radius)
      x.end = geo.to.euc(data.frame(lat=lat[k+1],lon=lon[k+1]),R=radius)
      n.pts = 2# min(10,max(3,round(1000*dist.euc(x.start,x.end)/R)))
      xseq = t(sapply(seq(0,1,length.out=n.pts),FUN=function(a){t(x.start+a*(x.end-x.start))}))
      xseq = radius*normalize.euc(xseq)
      rgl.linestrips(xseq[,1],xseq[,2],xseq[,3],...)
    }
  }
  else {
    x = geo.to.euc(data.frame(lat=lat,lon=lon),R=radius)
    rgl.lines(x[,1],x[,2],x[,3],...)
  }
}

#' Plot lines on sphere using rgl
#' 
#' This version of \link{rgl.rgl.sphlines2} uses a parametrization by the
#' start end end points of the lines
#'
#' @aliases rgl.sphlines2
#' @export
#' @examples \\dontrun{}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
#' 
rgl.sphlines2 = function(start,end,col="black",...){
  # And here comes the worst way to do this...
  slat = start[,"lat"]
  slon = start[,"lon"]
  elat = end[,"lat"]
  elon = end[,"lon"]
  lat = rep(0,2*length(slat))
  lat[seq(1,length(lat),by=2)] = slat
  lat[seq(2,length(lat),by=2)] = elat
  lon = rep(0,2*length(slon))
  lon[seq(1,length(lon),by=2)] = slon
  lon[seq(2,length(lon),by=2)] = elon
  rgl.sphlines(lat,lon,col=as.vector(rbind(col,col)),...)
}


rgl.earth = function(R=1,R.grid=1.05){

  require(sphereplot)
  
  # coordinates for texture
  n.smp = 50
  lat <- matrix(-asin(seq(-1,1, len=n.smp)), n.smp, n.smp, byrow=TRUE)
  long <- matrix(seq(-180, 180, len=n.smp)*pi/180, n.smp, n.smp)
  x <- R*cos(lat)*cos(long)
  y <- R*cos(lat)*sin(long)
  z <- R*sin(lat)
  
  # globe and texture
  persp3d(x, y, z, col="white", 
          texture=system.file("misc/Lambert_ocean.png",package="iDistance"), 
          specular="black", axes=FALSE, box=FALSE, xlab="", ylab="", zlab="",
          normal_x=x, normal_y=y, normal_z=z)
  
  # spheric grid
  rgl.sphgrid(longtype = "D",add=TRUE,radius=R.grid)

}

#' Plot an inla.mesh with geographic coordinates onto a sphere
#' 
#' @aliases rgl.sphmesh
#' @export
#' @examples \\dontrun{}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
#' 
rgl.sphmesh = function(mesh,radius=1,add=NULL,col=NULL,...){
  mesh$loc = geo.to.euc(data.frame(lat=mesh$loc[,2],lon=mesh$loc[,1]),R=radius)
  if (is.null(col)) { plot(mesh,rgl=TRUE,add=add,...) }
  else{ plot(mesh,rgl=TRUE,add=add,col=col,...) }
}

 