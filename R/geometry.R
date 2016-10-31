##############
# GENERICS
ssum = function(...){UseMethod("ssum")}
# distance = function(...){UseMethod("distance")}
normalize = function(...){UseMethod("normalize")}

### Unknown space

sdistance = function(sp, ep, geometry){
  class(sp) = c(geometry,"data.frame")
  class(ep) = c(geometry,"data.frame")
  if (geometry == "euc") { return(dist.euc(sp,ep)) }
  else if (geometry == "geo") { return(dist.geo(sp,ep)) }
} 

### A little tool to remove coordinate annotations
strip.coords = function(dframe) { names(dframe) = gsub("^.*?\\.","",names(dframe)); return(dframe)}


#################################################################################################
## TIME DOMAIN 

#' Time sum (SKELETON)
#'
#' @aliases sum.tim
#' @export
#' @return sum Sum of to times
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'
tsum = function(tim1,tim2) {return(rep(0,dim(tim1)[1]))}

#' Time distance (SKELETON)
#'
#' @aliases dist.tim
#' @export
#' @return dst
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'
tdist = function(tim1,tim2) {return(rep(0,dim(tim1)[1]))}



#################################################################################################
## EUCLIDIAN SPACE

#' Euclidian distance (CURRENTLY ONLY 2D)
#'
#' @aliases dist.euc
#' @export
#' @return distance Euclidian distance
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
dist.euc = function(euc1,euc2=NA) {
  if (!("x" %in% names(euc2))) {
    cols = which(names(euc1) %in% c("start.x","start.y","end.x","end.y"))
    euc2 = euc1[,cols[3:4]]
    euc1 = euc1[,cols[1:2]]
  }
  else {
    #euc1 = euc1[,c("x","y","z")]
    #euc2 = euc2[,c("x","y","z")]
  }
  dist = sqrt(apply((euc1-euc2)^2,1,sum))
  return(dist)
}

#' Euclidian sum (CURRENTLY ONLY 2D)
#'
#' @aliases sum.euc
#' @export
#' @return esum Sum of two points (interpreted as vectors)
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
ssum.euc = function(euc1,euc2=NA) {
  if (!("x" %in% names(euc2))) {
    cols = which(names(euc1) %in% c("start.x","start.y","end.x","end.y"))
    euc2 = euc1[,cols[3:4]]
    euc1 = euc1[,cols[1:2]]
  }
  else {
    euc1 = euc1[,c("x","y")]
    euc2 = euc2[,c("x","y")]
  }
  esum = euc1+euc2
  names(esum) = c("x","y")
  return(esum)
}



#' Normalize vector
#'
#' @aliases normalize.euc
#' @export
#' @return esum Sum of two points (interpreted as vectors)
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
normalize.euc = function(euc) {
  neuc = data.frame(t(apply(euc,MARGIN=1,function(x) {x/sqrt(sum(as.vector(x)^2))})))
  return(neuc)
}


#' Pseudo-cast data into Euclidian geometry
#'
#' Convert a data.frame with (lon/lat) or (c,z) entries to (x/y) entries. 
#' This is NOT a projection from the geographic space.
#'
#' @aliases pseudo.euc
#' @export
#' @return euc pseudo euclidian point
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
pseudo.euc = function(euc) {
  idx = which(names(euc) %in% c("start.lat","start.lon","end.lat","end.lon"))
  if (length(idx)>0){
    names(euc)[idx] = c("start.x","start.y","end.x","end.y")
  } 
  else if (any(names(euc) %in% c("lat","lon"))){
    names(euc)[which(names(euc) %in% c("lat","lon"))] = c("x","y")
  } else {
    names(euc)[which(names(euc) %in% c("c","z"))] = c("x","y")
  }
  return(euc)
}

#' Convert euclidean (x,y,z) to geographic (lat,lon) coordinates given the radius of the earth
#'
#' @aliases euc.to.geo
#' @export
#' @return cyl Cyldindric coordinates (c,z)
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
euc.to.geo = function(euc,R=6371) {
  lat = 180*asin(euc$z / R)/pi
  lon = 180*atan2(euc$y, euc$x)/pi
  return(data.frame(lat=lat,lon=lon))
}


#################################################################################################
## GEOGRAPHIC SPACE

#' Plot geographic coordinates
#'
#' @aliases plot.geo
#' @export
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
plot.geo = function(geo){plot(x=geo$lat,y=geo$lon)}

#' Distance in geographic geometry (geodesic distance) using Havesine formula.
#'
#'
#' @aliases dist.geo
#' @export
#' @return distance Distance in km
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
dist.geo = function(geo1,geo2=NA,R=6371) {
  if (!("lat" %in% names(geo2))){
    cols = which(names(geo1) %in% c("start.lat","start.lon","end.lat","end.lon"))
    geo2 = geo1[,cols[3:4]]; names(geo2) = c("lat","lon")
    geo1 = geo1[,cols[1:2]]; names(geo1) = c("lat","lon")
  }
  d2r = function(x) {x*pi/180} # degree to radians
  p = sin((d2r(geo2$lat) - d2r(geo1$lat))/2)^2 + cos(d2r(geo1$lat)) * cos(d2r(geo2$lat)) * sin((d2r(geo2$lon) - d2r(geo1$lon))/2)^2
  dist = R*2*asin(unlist(lapply(p,function(x) min(1,sqrt(x)))))
  return(dist)
}

dist2.geo <- function(geo1,geo2,R=6371) {
  d2r = function(x) {x*pi/180} # degree to radians
  distance = R * acos(sin(d2r(geo1$lat))*sin(d2r(geo2$lat)) + cos(d2r(geo1$lat))*cos(d2r(geo2$lat)) * cos(d2r(geo2$lon)-d2r(geo1$lon))) 
  return(distance)
}

#' Convert geographic (lon,lat) to cylindrical (c,z) coordinates given the radius of the earth
#'
#' @aliases geo.to.cyl
#' @export
#' @return cyl Cyldindrical coordinates (c,z) (range: c in [-pi,pi], z in [-1,1])
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
geo.to.cyl = function(geo,R=6371,B=NA) {
   if (is.numeric(B)){
    euc = geo.to.euc(geo)
    nms = names(euc)
    euc =  data.frame(as.matrix(euc)%*%B)
    names(euc) = nms
    geo = euc.to.geo(euc) 
   }
    lon.offset = 0
  #c = (pi*R/2)*geo$lat/90
  #z = (pi*R/4)*sin(2*pi*(geo$lon+lon.offset)/360))
  c = pi*geo$lon/180
  z = sin(pi*geo$lat/180)
  cyl = data.frame(c=c,z=z)
  class(cyl) = c("cyl","data.frame")
  return(cyl)
}

#' Convert geographic (lon,lat) to cylindric (c,z) coordinates given the radius of the earth
#'
#' @aliases geo.to.euc
#' @export
#' @return cyl Cyldindric coordinates (c,z)
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
geo.to.euc = function(geo,R=6371) {  
  euc = data.frame(x = R*cos(0.5*pi*geo$lat/90)*cos(0.5*pi*geo$lon/90), 
             y = R*cos(0.5*pi*geo$lat/90)*sin(0.5*pi*geo$lon/90),
             z = R*sin(0.5*pi*geo$lat/90) )
  return(euc)
}


#' Construct an orthonormal basis from two points in geographic space
#'
#' @aliases basis.geo
#' @export
#' @return B Orthonormal basis
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
basis.geo = function(geo1,geo2) {  
  euc1 = geo.to.euc(geo1)
  euc2 = geo.to.euc(geo2)
  z = cross(as.vector(t(euc1)),as.vector(t(euc2)))
  z = z/normvec(z)
  y = cross(cross(as.vector(t(euc1)),as.vector(t(euc2))),as.vector(t(euc1)))
  y = y/normvec(y)
  B = t(as.matrix(rbind(euc1/normvec(euc1),y,z)))
  return(B)
}

#################################################################################################
## CYLINDRIC SPACE

#' Plot cylindric coordinates
#'
#' @aliases plot.geo
#' @export
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
plot.cyl = function(cyl,...){plot(x=cyl$c,y=cyl$z,...)}

#' Distance in cylindric geometry
#'
#' @aliases dist.cyl
#' @export
#' @return distance
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
dist.cyl = function(cyl1,cyl2,R=6371) { return(sqrt((cyl1$c-cyl2$c)^2+(cyl1$z-cyl2$z)^2)) }

#' Spatial sum in cylindric geometry
#'
#' @aliases dist.cyl
#' @export
#' @return distance
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
ssum.cyl = function(cyl1,cyl2) { 
  cyl = cyl1+cyl2
  class(cyl) = c("cyl","data.frame")
  return(cyl) 
}

#' Convert cylindric (c,z) to geographic (lon,lat) coordinates
#'
#' @aliases cyl.to.geo
#' @export
#' @return (z,c)
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
cyl.to.geo = function(cz,B=NA) {
  lon = 180*cz$c/pi
  lat = 180*asin(cz$z)/pi
  geo = cbind(data.frame(lat=lat,lon=lon),cz)
  if (is.numeric(B)){
    euc = geo.to.euc(geo)
    nms = names(euc)
    euc =  data.frame(as.matrix(euc)%*%t(B))
    names(euc) = nms
    geo = euc.to.geo(euc) 
  }
  return(geo)
}

