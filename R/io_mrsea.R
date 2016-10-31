#' Marine renewables strategic environmental assessment 
#' 
#' Data imported from package MRSea, see http://creem2.st-andrews.ac.uk/software/
#'  
#' @examples \\donttest{ data(mrsea) ; plot(mrsea)}
#' @name mrsea
NULL

# INTERNAL DATA STORAGE
io_mrsea.getDataDir = function() {return(system.file("data",package="iDistance"))}


#' Load \link{mrsea} survey data from raw data sets
#'
#' @aliases io_mrsea.pkgdata.load
#' @export
#' @return \code{mrsea} the \link{mrsea} data set
#' @examples \\dontrun{mrsea = io_mrsea.pkgdata.load();}
#' @author Lindesay Scott-Hayward <\email{lass@st-andrews.ac.uk}>
#'

io_mrsea.pkgdata.load = function() { 
  
  library(MRSea)
  data("dis.data.re")
  data("predict.data.re")
  
  library(dplyr)
  preddata_gpseas<-group_by(predict.data.re, impact, segment.id )
  
  # Some housekeeping to change the name labels to the right format and correct the effort unit to match
  # the coordinate information.
  names(dis.data.re)[1:2]<-c('Transect.Label', 'Transect.label')
  names(dis.data.re)[8:9]<-c('x','y')
  names(dis.data.re)[5]<-c('Sample.Label')

  # change Effort column to same units as coordinates
  dis.data.re$Effort<-dis.data.re$length*1000
  
  segdata<-dis.data.re[,c("Transect.Label", "Transect.label" ,"season", "impact", "depth", "Sample.Label",
                          "segment.label" , "length", "Effort", 'x', 'y')]
  segdata<- distinct(segdata, Sample.Label)
  
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
  
  # Attach depth data to the data set
  
  depth = predict.data.re[(predict.data.re$season==1) & (predict.data.re$impact==0), c("x.pos","y.pos","depth")]
  colnames(depth) = c("x","y","depth")
  
  # Attach CRS projection strings
  dset$p4s = "+proj=utm +zone=32"
  dset$mesh.p4s = "+proj=utm +zone=32"
  
  return(dset)
  
}


#' Regenerate \link{mrsea} data and store it to \code{mrsea.RData}
#' 
#' Uses \code{\link{io_mrsea.pkgdata.load}} to load the data and stores
#' the result to mrsea.RData. Thereby the data that is distributed with 
#' our package is generated.
#'
#' @aliases io_mrsea.pkgdata.save
#' @export
#' @return NULL
#' @examples \\dontrun{io_mrsea.pkgdata.save();}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

io_mrsea.pkgdata.save = function(){
  ## save the data we will include in the R package
  mrsea = io_mrsea.pkgdata.load()
  save(mrsea,file=paste0(io_mrsea.getDataDir(),"/mrsea.RData"))
}

#########################################################################################################
#'
#' Tools written by Lindesay 
#'

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