##############
# GENERICS
ssum = function(...){UseMethod("ssum")}
dist = function(...){UseMethod("dist")}
normalize = function(...){UseMethod("normalize")}

### Unknown space

dist.data.frame = function(sp,ep, geometry){
  class(sp) = c(geometry,"data.frame")
  class(ep) = c(geometry,"data.frame")
  colnames(sp) = c("x","y")
  colnames(ep) = c("x","y")
  return(dist(sp,ep))
} 


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





#################################################################################
################################### JOYCE CODE ################################## 
#################################################################################

library(pracma)

otvec = function(v1, v2) {sum(v1*v2)}
normvec = function(x) sqrt(sum(x^2))


##########
## step1: transform lon&lat into sphere coordinates (x, y, z) system 
##########
# ## need to make sure long&lat are degrees 
# # Convert degrees to radians
#  deg2rad <- function(deg) return(deg*pi/180)
## for the cross product of vectors
tranformLongLatToXYZ = function(long, lat, R=6371)
{ ## the latitude and longitude values in ETP datasets are in degrees, not radians.  
  ## Positive values indicate N and E hemispheres, negative values S and W hemispheres
  ## need to transform long&lat from radiant to radiant
  long = long/180*pi 
  lat = lat/180*pi 
  ### how to deal with the negative long&lat 
  x = R*cos(lat)*cos(long)
  y = R*cos(lat)*sin(long)
  z = R*sin(lat)
  return(c(x,y,z))
}

##########
##  normailze a row vector, calcuate the norm; fun in step 2 to get the new v1 v2 v3 system
##########
normvec <- function(x) sqrt(sum(x^2))


#########
## step3: fun to tranform from (x,y,z) sphere into (v1,v2,v3) tangent plane system 
##      note the difference between tranforming points and vectors 
##      tranforming points is slightly more complicated than transforming vectors 
#########
### you only need the rotation matrix to tranform any vector in the new v1v2v3
### but this is different for transforming points, for any point p=(px,py,pz)
### --the transformed point in v1v2v3 is 
##  dot(P-A,v1)*v1 + dot(P-A,v2)*v2 + dot(P-A,v3)*v3 ] = p%*%Rmat - porig%*%Rmat
## where A is the new origianl point in the v1v2v3 system
transformPoint.SphereToTangent =function(p, porig, v1, v2, v3)
{#-- pv is the vector you want to transform, pv = p - porig 
  # porig (coordin in R3) is the new orgin used in the tangent plan system
  Rmat = cbind(v1, v2, v3)## this is the rotation matrix 
  pT = p%*%Rmat - porig%*%Rmat
  # pT=(p-porig)%*%Rmat
  return(pT)
}
### does not porig to tranform vectors 
transformVector.SphereToTangent = function(v, v1, v2, v3)
{ ## transform vector is different from tranform a point given new coordinates v1v2v3
  ## all you need is the rotaion matrix, simpler than tranforming points 
  Rmat = cbind(v1, v2, v3) 
  ## you can check the rotation matrix by  t(Rmat)%*%Rmat?=Identity matrix , 
  vT = v%*%Rmat
  return(vT)  
}

######## p0 is the corrodinates in R3

##########
## transform back from tangent plane (v1, v2, v3) to (xyz) sphere system 
## need to tranform sphere back to lon&lat system 
######### 
transformPoint.TangentToSphere.orginEarthCenter =function(pT, p0, v1, v2, v3)
{#-- pv is the vector you want to transform, pv = p - porig 
  #   porig here should the the coordinates of c(0,0,0) in the new tangent plane system
  ## 1st need to find the tangent system coordinates of earch center c(0, 0, 0)
  ####### you need p0 is not as the new orig, but to get c(0,0,0) in the v1c2v3 system
  pc.inv1v2v3 = transformPoint.SphereToTangent(p=c(0,0,0), porig=p0, v1, v2, v3)
  ## then pc.inv1v2v3 is used as the origin in the R3 when transforming back from tangent to R3 
  Rmat = cbind(v1, v2, v3)## this is the rotation matrix 
  psphere = pT%*%t(Rmat) - pc.inv1v2v3%*%t(Rmat)
  # solve(Rmat) == t(Rmat) ## solve(Rmat) gives the inverse of Rmat
  return(psphere)
}


#############################
## then need to tranform from sphere R3 to lon&lat 
#############################
convert.fromXYZtoLongLat = function(pxyz,R=6371)
{
  x=pxyz[1]; y=pxyz[2]; z=pxyz[3]
  long = atan2(y,x)
  lat = asin(z/R)
  return(c(lon=long*180/pi, lat=lat*180/pi))
}
# Using formulas with atan2() is more convenient. You don't have to add/subtract pi/2 or care about sign issues in different quadrants or division by zero.
# 
# lat will be >0 in the northern hemisphere
# lat will be <0 in the southern hemisphere
# lng will be >0 in the eastern hemisphere
# lng will be <0 in the western hemisphere


## p0 is the corrodinates of p0 in the R3 shpere system
############
## for the integration grid (intgrid) defined in the tangent plane, we need to tranform it back to long and lat 
## Given intgrid.tangent(an integration grid designed on tangent plane), 
##  -- the function to tranform the v1v2v3 points back to lon&lat
###########
transform.intgrid.fromTagentTolonlat = function(intgridOnTangent, v1, v2, v3, R, p0)
{
  #step1 is to transform back to R3 
  ngrid  = nrow(intgridOnTangent)
  transfromToR3 = t(apply(matrix(c(1:ngrid), ncol=1), 1, 
                          function(X){
                            transformPoint.TangentToSphere.orginEarthCenter(
                              pT=as.numeric(intgridOnTangent[X,]),p0, 
                              v1, v2, v3) }))
  tranfromFromR3Tolonlat = t(apply(matrix(c(1:ngrid), ncol=1), 1, 
                                   function(X){convert.fromXYZtoLongLat(
                                     pxyz=transfromToR3[X,],R=6371)}))
  return(tranfromFromR3Tolonlat)
}

#############
## final fun for obtaining the integration grid on lonlat by designing it on tangent sphere
#############
designintegration.oneLTsegement = function(p0lonlat, p1lonlat, w, dLT, dPD)
{
  ## w (km) is the trancation distance 
  ## dPD is the integration cell size along PD direction 
  ## dLT is the integration cell size along the transect line 
  ## calculate the distance of the LT segment, given start and end points p0 and p1 
  library(fields)
  L=rdist.earth(x1=t(as.numeric(p0lonlat)), 
                x2=t(as.numeric(p1lonlat)), miles = FALSE)
  ############
  ##step1 tranform sphere into R3 
  ## the axis of the new plane system 
  p0 = tranformLongLatToXYZ(long=p0lonlat$lon, lat=p0lonlat$lat, R=6371)
  p1 = tranformLongLatToXYZ(long=p1lonlat$lon, lat=p1lonlat$lat, R=6371)
  ############
  ## step2 work out the tangent plane coordinate system v1, v2, v3 (orthogonal vectors for the three axes and later the rotation matrix to tranform from sphere to tangent plane system)
  v1 =  p0 ; v1 = v1/normvec(v1)
  v3 = cross(v1, p1-p0); v3 = v3/normvec(v3)
  v2 = cross(v3, v1)
  # v1;v2;v3
  ############
  ## step3 tranform p0 and p1 to p0T and p1T (related fun can be found in rfun file)
  #-- tranform the points 
  p0T = transformPoint.SphereToTangent(p=p0, porig=p0,v1,v2,v3 )
  p1T = transformPoint.SphereToTangent(p=p1, porig=p0,v1,v2,v3 )
  #-- tranform the vectors 
  v1T = transformVector.SphereToTangent(v1, v1, v2, v3)
  v2T = transformVector.SphereToTangent(v2, v1, v2, v3)
  v3T = transformVector.SphereToTangent(v3, v1, v2, v3)
  #-- output
  out.step3orig  = rbind(p0T = p0T, p1T = p1T, v1T=v1T, v2T=v2T, v3T=v3T)
  rownames(out.step3orig)=c("p0T","p1T","v1","v2","v3" )
  out.step3round = rbind(p0T=round(p0T), p1T=round(p1T), v1T=round(v1T),
                         v2T=round(v2T), v3T=round(v3T))
  rownames(out.step3round)=c("p0T","p1T","v1","v2","v3" )
  ###########
  ## step4 then design the integration grid on the tangent plane along the   direction v2(LT) and v3(PD) 
  ## If use p0 as the origin, then (0, 0, 10) describes a point 10 km PD.
  # dLT = 1 # cell size on v2 axes 
  # dPD = 0.5 # cell size on v3 axes
  nPD = w/dPD  ## no. of cells on v3 direction (dPD=dv3)
  cellPD =c(.5*dPD + c(0:(nPD-1))*dPD, -.5*dPD - c(0:(nPD-1))*dPD) ## both sides of the LT 
  ########
  ### note that LT is not an integer, if the decimal part >0.5, then treated as another cell 
  ### if the decimal part is less or equal to 0.5, then combine it with another dLT to make one grid cell 
  if(L>dLT)
  {
    if(round(L/dLT) < (L/dLT)) {nLTminus = floor(L/dLT)-1}  ## no. of cells on v2 direction (dLT=dv2) -1
    else {nLTminus = ceiling(L/dLT)-1} 
    ## given the cell size on both directions, work out the v1v2v3 coordinates for the integration points 
    cellLT.lastcellsize = L-dLT*nLTminus
    cellcenter.lastcell = L -0.5*cellLT.lastcellsize
    cellLT =c(.5*dLT + c(0:(nLTminus-1))*dLT, cellcenter.lastcell) ## the final intPoint is not about dLT 
    intgridv1v2v3 = expand.grid(v2 = cellLT, v3 = cellPD)
    ## all the integration points on the tangent plane, coordinates for v1 = 0 
    intgridv1v2v3$v1 = 0 
    intgridv1v2v3$weight = dPD*dLT 
  }
  else
  { ## if length(LT) < dLT the cell size specified, then weight has nothing to do with dLT  
    cellLT = .5*L
    intgridv1v2v3 = expand.grid(v2 = cellLT, v3 = cellPD)
    ## all the integration points on the tangent plane, coordinates for v1 = 0 
    intgridv1v2v3$v1 = 0 
    intgridv1v2v3$weight = rep(dPD*L/2, nrow(intgridv1v2v3))
  }
  ## use expand.grid() to get all the points on integration grid 
  
  # head(intgridv1v2v3)
  # intsquare = rbind(p0N =   c(0, 0, w), p1N =   c(0, L, w),
  #                  p1S =   c(0, L, -w), p0S =   c(0, 0, -w))
  # p0N; p0S; p1N; p1S
  ### plot the integration grid on tangent plane with bound 
  #     plotpoly = rbind(intsquare[,-1],intsquare[1,-1] )
  #     quartz()
  #     pdf(paste(plotdir, "/intgridonTangentPlane.pdf", sep=""))# width=7,height=5)
  #     plot(plotpoly[,1], plotpoly[,2], type = "n", main = "integration grid on tangent plane ")
  #     abline(0,0)
  #     polygon(plotpoly[,1], plotpoly[,2])
  #     points(as.numeric(intgridv1v2v3$v2), as.numeric(intgridv1v2v3$v3))
  #     dev.off()
  ###################
  ## step 5 tranform the integration grid on the tangent plane back to R3, and then from R3 back to lon&lat, done in one fun: transform.intgrid.fromTagentTolonlat 
  #####################
  #   ## pc is short for the point of the sphere center, which is the orgin of the R3 system 
  #   ## in tangent system, p0 is the new orgin, when transfom back to R3, you need to use the  coordinates of pc=c(0, 0, 0) in the new tangent system, which is the origin when transform back to R3 
  #   p1check = transformPoint.TangentToSphere.orginEarthCenter(pT=p1T, v1, v2, v3)
  #   p1-p1check
  #   ### the above is to check the function tranform between tangent and sphere,
  #   ### p0-> p0T-->p0.check, the difference between p0.check and p0 is 
  #   # p1-p1check
  #   #      [,1]         [,2]          [,3]
  #   # [1,]    0 9.094947e-13 -9.094947e-13
  #   ## the difference between p1 and p1check(obtained from transforming) is essenstiallly zero. which means the fun works 
  ## for any point pT in the tangent plane with p0 as the origin, we transform it back to R3 with earth center being the origin 
  intgridlonlat=transform.intgrid.fromTagentTolonlat(intgridOnTangent=intgridv1v2v3[,c("v1", "v2","v3")],
                                                     v1, v2, v3, R=6371, p0)
  # ## you got to make sure the row of intgrid enters the function in the right sequence as (v1,v2,v3) !!intgrid[1,c("v1", "v2","v3")])!!
  # quartz()
  # pdf(paste(plotdir, "/intgridonSpherelonlat.pdf", sep=""))# width=7,height=5)  
  # plot(intgridlonlat[,"lon"],intgridlonlat[,"lat"], xlab="lon", ylab="lat", main="integration points on sphere")
  # arrows(p0lonlat$lon, p0lonlat$lat, p1lonlat$lon, p1lonlat$lat)
  # points(p0lonlat$lon, p0lonlat$lat,col="blue")
  # points(p1lonlat$lon, p1lonlat$lat,col="red")
  # dev.off()
  
  # p0lonlat.check = convert.fromXYZtoLongLat(pxyz=p0,R=6371)
  # p0lonlat.check
  # p0lonlat
  # ##### p0lonlat.check is consistent with p0lonlat original from the effdatset 
  # p1lonlat.check = convert.fromXYZtoLongLat(pxyz=p1,R=6371)
  # p1lonlat.check
  # p1lonlat
  
  ### final output should be a matrix with lon&lat for all the integration points projected from sphere to the sphere
  finalv1v2v3 = cbind(intgridv1v2v3, dPD=dPD, dLT =dLT, w=w, L=L)
  finallonlat = cbind(intgridlonlat, weight=intgridv1v2v3$weight)
  final=list(intgridv1v2v3 = finalv1v2v3 , intgridlonlat=finallonlat )
  return(final)
}

## note 
## 1) p0 (coordinates in R3) is used as a new origin in v1v2v3 system for any new transect lines. 
##    
