#' Generic distance sampling data format
#' 
#' Fields:
#' - effort = 
#' - geometry = "euc"
#' - mesh.coords = c("lon","lat")
#' 
#' @name dsdata
#' 
NULL

#' Plot generic distance sampling data
#' 
#' @aliases plot.dsdata 
#' @export
#' @param dsdata Distance sampling data set
#' @name blah
#' @examples \dontrun{ toy=toy1() ; plot(toy) }
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 

plot.dsdata = function(data,rgl=FALSE,det.col="red",add=FALSE,R=1,col=NULL,...){
  if (rgl==FALSE) {
    c1 = data$mesh.coords[1]
    c2 = data$mesh.coords[2]
    xlim=range(data$mesh$loc[,1])
    ylim=range(data$mesh$loc[,2])
    
    # Plot mesh
    plot(data$mesh,col=col,main="")
    
    # Plot transect lines
    spoint = startpoint(as.transect(data),data)
    epoint = endpoint(as.transect(data),data)
    lines(x=as.vector(t(cbind(spoint[,c1],epoint[,c1],rep(NA,dim(spoint)[1])))),
          y=as.vector(t(cbind(spoint[,c2],epoint[,c2],rep(NA,dim(spoint)[1])))),
          xlim=xlim,ylim=ylim,
          lwd=2)
    
    # Plot detections
    par(new=TRUE)
    plot(detdata(data)[,data$mesh.coords],xlim=xlim,ylim=ylim,col=det.col)
  }
  else {
    require(rgl)
    if (data$geometry == "euc"){
      # plot mesh
      plot(data$mesh,rgl=TRUE,col=col)
      
      # plot detections
      detections = cbind(detdata(data)[,c("x","y")],z=rep(0.01,nrow(detdata(data))))
      rgl.points(x=detections, z=0.001,size=5,col="red")
      
      # plot transect lines
      sp = startpoint(as.transect(data),data)
      ep = endpoint(as.transect(data),data)
      rgl.lines(cbind(matrix(cbind(sp,ep),ncol=2,byrow=TRUE),z=0.01),lwd=5,col="black")
    }
    else if (data$geometry == "geo"){
      # Elevation for plotting on top of earth surface
      R.delta=0.001
      
      # earth
      rgl.earth()
      rgl.viewpoint(0,-90)

      # detections
      detections = detdata(data)[,c("lat","lon")]
      rgl.sphpoints(long=detdata(data)[,"lon"]+360,lat=detdata(data)[,"lat"],radius=R+2*R.delta,col="red",size=5)
      if (is.null(col)){
        rgl.sphmesh(data$mesh,add=TRUE,radius=R+R.delta,lwd = 1,edge.color = rgb(0,0,0),alpha=0.3,draw.segments = TRUE)
      } else {
        rgl.sphmesh(data$mesh,add=TRUE,radius=R+R.delta,lwd = 1,edge.color = rgb(0,0,0),col=col)
      }
      
      do.call(par3d,data$par3d.args)
      
      # transects
      plot(as.transect(data$effort),data=data,rgl=TRUE,add=TRUE,radius=R+R.delta,col=rgb(0.1,0.1,0.7))
      
    }
  }
}


#' Summarize generic distance sampling data
#' 
#' @aliases summary.dsdata 
#' @export
#' @param dsdata Distance sampling data set
#' @examples \dontrun{ data(toy1) ; summary(toy1) }
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

summary.dsdata = function(dsdata){
  cat(paste0("Geometry: "),dsdata$geometry,"\n")
  cat(paste0("Coordinate names: "),dsdata$mesh.coords,"\n")
  cat(paste0("Effort columns names: "),"\n")
  cat(paste0(" ",names(dsdata$effort)), "\n")
  cat(paste0("Number of transects: "), dim(as.transect(dsdata))[1], "\n")
  cat(paste0("Total length of transects: "), sum(len.transect(as.transect(dsdata),dsdata)), "units","\n")
  
  # Mesh area
  A.mesh = sum(Matrix::diag(inla.mesh.fem(dsdata$mesh, order=1)$c0))
  cat(paste0("Total area of mesh: "), A.mesh, " square mesh units","\n")
  
  # Survey area
  if ("inner.boundary" %in% names(dsdata)) {
    msk =  is.inside.polygon(whales$mesh, whales$inner.boundary,sst$mesh$loc, whales$mesh.coords)
    A = sum(Matrix::diag(inla.mesh.fem(dsdata$mesh, order=1)$c0) * msk)
    cat(paste0("Total survey area: A ="), A, " square mesh units","\n")
  }
  
  # Number of transects
  if (!class(dsdata$effort)[1]=="etpeffort"){ # Note yet implemented by etpdata class
    cat(paste0("Number of segements: "), dim(as.segment(dsdata))[1], "\n")
  }
  
  # Number of detections
  n.det =  dim(detdata(dsdata))[1]
  cat(paste0("Number of detections : N = "),n.det, "\n")
  
  # Number of detections per unit mesh area inside survey area
  cat(paste0("Detections per survey area: N/A = "), n.det/A.mesh, "\n")
  
}

as.transect = function(...){UseMethod("as.transect")}
as.segment = function(...){UseMethod("as.segment")}
as.detection = function(...){UseMethod("as.detection")}
startpoint = function(...){UseMethod("startpoint")}
endpoint = function(...){UseMethod("endpoint")}
detdata = function(...){UseMethod("detdata")}
trdata = function(...){UseMethod("trdata")}

join = function(...){UseMethod("join")}
id = function(...){UseMethod("id")}

as.transect.dsdata = function(dsdata,...) {return(as.transect(dsdata$effort,...))}
as.segment.dsdata = function(dsdata,...) {return(as.segment(dsdata$effort,...))}
as.detection.dsdata = function(dsdata,...) {return(as.detection(dsdata$effort,...))}
startpoint.dsdata = function(dsdata,obj,...) { return(startpoint(obj,dsdata,...))}
endpoint.dsdata = function(dsdata,obj,...) { return(endpoint(obj,dsdata,...))}


detdata.dsdata = function(data,detection=NULL,...){ 
  if (is.null(detection)) {
    return(data$effort[as.detection(data)[,"start"],])
  } 
  else {
    return(data$effort[detection[,"start"],])
  }
}

trdata.dsdata = function(data,tr=NULL,...){ 
  if (is.null(tr)) {
    return(data$effort[as.transect(data)[,"start"],])
  } 
  else {
    return(data$effort[tr[,"start"],])
  }
}

#' Extract transect pointers from effort data
#' 
#' @aliases as.transect.effort 
#' @export
#' @param effort Effort data set
#' @return transect Transects
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

as.transect.effort = function(effort){
  end.idx = findInterval(1:length(levels(effort$tr)),as.numeric(effort$tr))
  if (length(end.idx)==1){ start.idx=1 } 
  else {start.idx = c(1,end.idx[1:(length(end.idx)-1)]+1)}
  tr = data.frame(start=start.idx,end=end.idx)
  class(tr) = c("transect","data.frame")
  return(tr)
}

#' Extract segment pointers from effort data
#' 
#' @aliases as.effort.effort 
#' @export
#' @param effort Effort data set
#' @return segmet Segments
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

as.segment.effort = function(effort){
  end.idx = findInterval(1:length(levels(effort$seg)),as.numeric(effort$seg))
  if (length(end.idx)==1){ start.idx=1 } 
  else {start.idx = c(1,end.idx[1:(length(end.idx)-1)]+1)}
  seg = data.frame(start=start.idx,end=end.idx)
  class(seg) = c("segment","data.frame")
  return(seg)
}

#' Extract detection pointers from effort data
#' 
#' @aliases as.detection.effort 
#' @export
#' @param effort Effort data set
#' @return sighting sightings
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

as.detection.effort = function(effort){
  idx = which(!(is.na(effort$det)))
  det = data.frame(start=idx,end=idx)
  class(det) = c("detection","data.frame")
  return(det)
}


#' Join an effort data set with a sighting data set
#' 
#' @aliases join.effort 
#' @export
#' @param effort Effort data set
#' @param detection Detections
#' @return effort Effort data set including provided sightings
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

join.effort = function(effort,detection){
  all = rbind(effort[,c("seg","det")],detection[,c("seg","det")])
  all = merge(x=all,y=detection,by=c("det","seg"),all.x=TRUE)
  final = merge(x=all,y=effort,by=c("seg"),all.x=TRUE,suffixes=c("",".EFFORT"))
  final = final[,!(names(final)=="det.EFFORT")]
  final$det[is.na(final$det)] = "0.0.0.0"
  final = final[order(final$det),]
  final = final[order(final$seg),]
  final$det[final$det=="0.0.0.0"] = NA
  class(final) = c("effort","data.frame")
  return(final)
}


#' Start points coordinates of transect(s)
#' 
#' @aliases startpoint.transect
#' @export
#' @param transect
#' @param data Data set with $effort
#' @return spoint
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

startpoint.transect = function(tr,data,keep=FALSE) {
  cols = names(data$effort) %in% c("start.x","start.y")
  if (keep){ return(data$effort[tr$start,]) }
  else {
    pts = data$effort[tr$start,cols]
    names(pts) = c("x","y")
    class(pts) = c(data$geometry,"data.frame")
    return(pts)
  }
}

#' End points coordinates of transect(s)
#' 
#' @aliases endpoint.transect
#' @export
#' @param transect
#' @param data Data set with $effort
#' @return epoint
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

endpoint.transect = function(tr,data,keep=FALSE) {
  cols = names(data$effort) %in% c("end.x","end.y")
  if (keep){ return(data$effort[tr$end,]) } 
  else {
    pts = data$effort[tr$end,cols]
    names(pts) = c("x","y")
    class(pts) = c(data$geometry,"data.frame")
    return(pts)
  } 
}

#' Number of transects
#' 
#' @aliases numel.transect
#' @export
#' @param transects
#' @param data Data set with $effort
#' @return n
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

numel.transect = function(tr,data=data) { return(dim(tr)[1]) }


#' Linedata
#' 
#' @aliases linedata linedata.transect
#' @export
#' @param transects
#' @param data Data set with $effort
#' @return data
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

linedata = function(...){
  #stop("This method is deprecated")
  UseMethod("linedata")
}
linedata.transect = function(tr,data,fields){
  return(data$effort[tr$start,fields])
}
linedata.segment = function(seg,data,fields){
  return(data$effort[seg$start,fields])
}


#' Start points coordinates of segment(s)
#' 
#' @aliases startpoint.segment
#' @export
#' @param segment
#' @param data Data set with $effort
#' @return spoint
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

startpoint.segment = startpoint.transect

#' End points coordinates of segment(s)
#' 
#' @aliases endpoint.segment
#' @export
#' @param segment
#' @param data Data set with $effort
#' @return epoint
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

endpoint.segment = endpoint.transect

#' Transect name according to convention, e.g. "Trans.Name" = 1.2
#' 
#' @aliases id.transect
#' @export
#' @param transect
#' @param data Data set with $effort
#' @return id
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

id.transect = function(tr,data){ return(data$effort[tr$start,"tr"]) }

#' Segment name according to convention, e.g. "Seg.Name" = 1.2.2
#' 
#' @aliases id.segment
#' @export
#' @param segment
#' @param data Data set with $effort
#' @return id
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

id.segment = function(seg,data){ return(data$effort[seg$start,"seg"]) }

