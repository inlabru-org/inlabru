#' @name mrsea
#' @title Marine renewables strategic environmental assessment
#' @docType data
#' @description Data imported from package MRSea, see http://creem2.st-andrews.ac.uk/software/
#' 
#' 
#' @format A list of objects:
#'  \describe{
#'    \item{\code{points}:}{ A \code{SpatialPointsDataFrame} object containing the locations of 
#'    XXXXX.}
#'    \item{\code{samplers}:}{ A \code{SpatialLinesDataFrame} object containing the transect lines
#'    that were surveyed.}
#'    \item{\code{mesh}:}{ An \code{inla.mesh} object containing a Delaunay triangulation 
#'    mesh (a type of discretization of continuous space) covering the survey region.}
#'    \item{\code{boundary}:}{ An \code{SpatialPolygonsDataFrame} object defining the boundary of the 
#'    survey region.}
#'    \item{\code{covar}:}{ An \code{SpatialPointsDataFrame} containing sea depth estimates.}
#'  }
#' @source 
#' Library \code{MRSea}.
#' 
#' @references 
#' 
#' NONE YET
#' 
#' @examples
#'  data(mrsea)
#'  ggplot() + gg(mrsea$mesh) + gg(mrsea$samplers) + gg(mrsea$points) + gg(mrsea$boundary)
NULL

#' MRSea data import
#' 
#' 
#' Load \link{mrsea} survey data from MRSea package and convert to spatial formats defined by the \link{sp} package.
#'
#' @aliases import.mrsea
#' @keywords internal
#' @return The \link{mrsea} data set
#' @examples \dontrun{mrsea = import.mrsea();}
#' @author Lindesay Scott-Hayward <\email{lass@@st-andrews.ac.uk}>
#'

import.mrsea = function() { 
  
  # library(MRSea)
  predict.data.re = NULL
  dis.data.re = NULL
  data("dis.data.re", package = "MRSea", envir = environment())
  data("predict.data.re", package = "MRSea", envir = environment())
  impact = dis.data.re$impact
  segment.id = dis.data.re$segment.id
  
  preddata_gpseas<-dplyr::group_by(predict.data.re, impact, segment.id )
  
  # Some housekeeping to change the name labels to the right format and correct the effort unit to match
  # the coordinate information.
  names(dis.data.re)[1:2]<-c('Transect.Label', 'Transect.label')
  names(dis.data.re)[8:9]<-c('x','y')
  names(dis.data.re)[5]<-c('Sample.Label')

  # change Effort column to same units as coordinates
  dis.data.re$Effort<-dis.data.re$length*1000
  
  segdata<-dis.data.re[,c("Transect.Label", "Transect.label" ,"season", "impact", "depth", "Sample.Label",
                          "segment.label" , "length", "Effort", 'x', 'y')]
  segdata<- dplyr::distinct(segdata, segdata$Sample.Label, .keep_all = TRUE)
  
  # effort, object and distance.
  # Not taken x and y as these are segement mid points not detection locations
  distdata<-dis.data.re[,c("object" ,"distance" , "Effort")]
  distdata<-na.omit(distdata)
  distdata$size<-rep(1, nrow(distdata))
  
  # obsdata
  obsdata<-na.omit(dis.data.re)
  obsdata<-obsdata[,c("object" , "Sample.Label", "distance" , "Effort")]
  obsdata$size<-rep(1, nrow(obsdata))

  preddata = predict.data.re
  colnames(preddata)[c(2,3)] = c("x","y")
  
  dsmdata = list(obsdata = obsdata, distdata = distdata, segdata = segdata, preddata = preddata)
  dset = import.dsmdata(dsmdata, covar.col = 5)
  
  # Depth data to the data set
  
  depth = predict.data.re[(predict.data.re$season==1) & (predict.data.re$impact==0), c("x.pos","y.pos","depth")]
  colnames(depth) = c("x","y","depth")
  

  ############ NEW FORMAT USING sp objects ##############
  
  # Transects lines
  lns = subset(dset$effort, is.na(det)) ; class(lns) = "data.frame"
  lns = sline(lns, c("start.x", "start.y"), c("end.x", "end.y"), crs = CRS("+proj=utm +zone=32"))
  lns$weight = 500
  
  # Detections
  pts = subset(dset$effort, !is.na(det))[c("x","y","season","distance")]; class(pts) = "data.frame"
  coordinates(pts) = c("x","y")
  proj4string(pts) = CRS("+proj=utm +zone=32")
  
  # Mesh
  mesh = dset$mesh
  mesh$crs = INLA::inla.CRS("+proj=utm +zone=32")
  
  
  # Boundary
  boundary = spoly(dset$mesh$loc[dset$mesh$segm$int$idx[,1], 1:2], CRS("+proj=utm +zone=32"))
  
  # Covariates
  covar = SpatialPointsDataFrame(depth[,1:2], data = depth[,3,drop=FALSE], proj4string = CRS("+proj=utm +zone=32"))
  
  # Remove `distance` column from transects
  lns$distance = NULL
  
  mrsea = list(points = pts,
               samplers = lns,
               boundary = boundary,
               mesh = mesh,
               covar = covar)
  
}


# INTERNAL DATA STORAGE
io_mrsea.getDataDir = function() {return(system.file("data",package="inlabru"))}


# Regenerate \link{mrsea} data and store it to \code{mrsea.RData}
# 
# Uses \code{\link{io_mrsea.pkgdata.load}} to load the data and stores
# the result to mrsea.RData. Thereby the data that is distributed with 
# our package is generated.
#
# @aliases io_mrsea.pkgdata.save
# @export
# @return NULL
# @examples \\dontrun{io_mrsea.pkgdata.save();}
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#

io_mrsea.pkgdata.save = function(){
  ## save the data we will include in the R package
  mrsea = import.mrsea()
  save(mrsea,file=paste0(io_mrsea.getDataDir(),"/mrsea.RData"))
}

mrseanames2dsmnames<-function(data){
  nam<-names(data)
  cols2change<-c('transect.id', 'transect.label', 'x.pos', 'y.pos', 'segment.id')
  id<-NULL
  for(i in 1:length(cols2change)){
    id<-c(id, grep(cols2change[i], nam))
  }
  names(data)[id]<-c('Transect.Label', 'Transect.label', 'x','y', 'Sample.Label')
  # make sure effort column is same unit as coordinate system
  data$Effort<-data$length
  return(data)
}
# ---------------------------------------------------------------------
prednamechange<-function(data){
  id<-c(grep('x.pos', colnames(data)), grep('y.pos', colnames(data)))
  names(data)[id]<-c('x', 'y')
  return(data)
}
# ---------------------------------------------------------------------
makedistdata<-function(data){
  # effort, object and distance.
  # Not taken x and y as these are segement mid points not detection locations
  distdata<-data[,c("object" ,"distance" , "Effort")]
  if(is.null(data$size)){
    distdata$size<-rep(1, nrow(distdata))  
  }else{
    distdata$size<-data$size
  }
  distdata<-na.omit(distdata)
  return(distdata)
}
# ---------------------------------------------------------------------
makeobsdata<-function(data){
  obsdata<-na.omit(data)
  if(is.null(data$size)){
    obsdata$size<-rep(1, nrow(obsdata))  
  }else{
    obsdata$size<-obsdata$size
  }
  obsdata<-obsdata[,c("object" , "Sample.Label", "distance" , "Effort", 'size')]
  return(obsdata)
}
# ---------------------------------------------------------------------