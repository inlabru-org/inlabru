#' iDistance data format for distance sampling surveys
#' 
#' A \code{dsdata} object is a \code{list} with the following fields:
#' \itemize{
#'  \item{effort: }{A data.frame with distance sampling data.}
#'  \item{mesh: }{An \code{inla.mesh} modeling the survey area.}
#'  \item{geometry: }{String determining the geometry of the data. 
#'  Choices are \code{"euc"} for two-dimensional data sets in the Euclidean plane and \code{"geo"} for a geographic geometry on the sphere}
#'  \item{mesh.coords: }{Column names of the effort table that describe coordinates of detections and transects. 
#'  These are arbitrary but for Euclidean data you might want to use \code{c("x","y")} and \code{c("lon","lat"}) for geographic coordinates.}
#'  }
#' Data sets included in iDistance:
#' \itemize{
#'  \item{\link{toy1}: }{A toy data set with a single transect line}
#'  \item{\link{whales}: }{Blue whales (line transect)}
#'  \item{\link{strdolphin}: }{Striped dolphins (line transect)}
#'  \item{\link{sbdolphin}: }{Short beaked dolphins (line transect)}
#' }
#' 
#' @name dsdata
#' 
NULL


#' Make dsdata 
#' 
#' @aliases make.dsdata 
#' @export
#' @examples \dontrun{  }
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 

make.dsdata = function(effort = NULL, geometry = NULL, mesh = NULL, mesh.coords = NULL, ...) {
  if (is.null(effort)) { 
    stop("You have to provide an effort data frame.") 
  }
  
  if (is.null(geometry) | !any((geometry %in% c("euc", "geo")))) { 
    stop(paste0("You have to provide a geometry that is not supported: ",geometry)) 
  }
  
  if (is.null(mesh)) { 
    stop("You have to provide a mesh.")
  }
  
  if (is.null(mesh.coords)) { 
    stop("You have to provide mesh.coords.")
  }
  
  #
  # Seems like we have all we need. Let's assemble the dsdata object.
  #
  
  dset = list(effort = effort,
            geometry = geometry,
            mesh = mesh,
            mesh.coords = mesh.coords, ...)  
  
  class(dset) = c("dsdata", "list")
  
  #
  # Sanity check
  #
  
  sanity(dset)
  
  #
  # Done.
  # 
  
  return(dset)
  
}


#' Plot generic distance sampling data
#' 
#' @aliases plot.dsdata 
#' @export
#' @param dsdata Distance sampling data set
#' @name blah
#' @examples \dontrun{ toy=toy1() ; plot(toy) }
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 

plot.dsdata = function(data,rgl=FALSE,det.col="red",add=FALSE,R=1,col=NULL,asp=1,...){
  if (rgl==FALSE) {
    c1 = data$mesh.coords[1]
    c2 = data$mesh.coords[2]
    xlim=range(data$mesh$loc[,1])
    ylim=range(data$mesh$loc[,2])
    
    # Plot mesh
    plot(data$mesh, col = col, main="", asp = asp)
    
    # Plot transect lines
    spoint = startpoint(as.transect(data),data)
    epoint = endpoint(as.transect(data),data)
    lines(x=as.vector(t(cbind(spoint[,c1],epoint[,c1],rep(NA,dim(spoint)[1])))),
          y=as.vector(t(cbind(spoint[,c2],epoint[,c2],rep(NA,dim(spoint)[1])))),
          xlim=xlim,ylim=ylim,
          lwd=2)
    
    # Plot detections
    par(new=TRUE)
    plot(detdata(data)[,data$mesh.coords],xlim=xlim,ylim=ylim,col=det.col, asp = asp)
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



#' Check sanity of a data set
#' 
#' @aliases sanity.dsdata 
#' @export
#' @param dsdata Distance sampling data set
#' 
#' @examples \dontrun{ data(toy1) ; sanity(toy1) }
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

sanity.dsdata = function(dset){

  # Error handlers
  sanity.stop = function(ret) { if (!ret$sane) { stop(ret$message) } else { cat(" [OK]\n") } }
  sanity.warning = function(ret) { if (!ret$sane) { warning(ret$message) } else { cat(" [OK]\n") } }
  
  # Checks
  
  cat("Checking if data set has all the necessary fields ...")
  sanity.stop(sanity.fields(dset, message = TRUE))
  
  cat("Checking if $effort has all the necessary columns ...")
  sanity.stop(sanity.columns(dset, message = TRUE))
  
  #cat("Checking if effort$seg is sorted ...")
  #sanity.stop(sanity.unsorted(dset, message = TRUE))
  
  cat("Checking if all segments are inside the mesh boundary ...")
  sanity.stop(sanity.effort.inside(dset, message = TRUE))
  
  cat("Checking if all detections are inside the mesh boundary ...")
  sanity.stop(sanity.detection.inside(dset, message = TRUE))
  
  cat("Checking if any distances are NA or NaN...")
  sanity.stop(sanity.finite(dset, message = TRUE))
  
  #
  # The following checks will only throw warnings
  #
  
  cat("Checking if $effort has distance column ...")
  sanity.warning(sanity.distance(dset, message = TRUE))
  
}


# Finite data?
sanity.finite = function(dset, message = FALSE) {
  # cols = c(paste0("start.",dset$mesh.coords), paste0("end.",dset$mesh.coords), dset$mesh.coords, "distance")
  cols = "distance"
  sane = !any(!is.finite(as.matrix(detdata(dset)[,cols]))) 
  if ( message & !sane) { message = paste0("Non-finite values in one of the following columns: ", do.call(paste, c(as.list(cols),sep = ", "))) }
  return(list(sane = sane, message = message))
}


# Has all the fields? effort, mesh, mesh.coords, geometry
sanity.fields = function(dset, message = FALSE) {
  cols = c("mesh", "effort", "mesh.coords", "geometry")
  sane = !any(!(cols %in% names(dset)))
  if ( message & !sane) { message = paste0("The data set is missing one of the fields: ", do.call(paste, c(as.list(cols),sep = ", "))) }
  return(list(sane = sane, message = message))
}

# Has Distance column
sanity.distance = function(dset, message = FALSE) {
  cols = c("mesh", "effort", "mesh.coords", "geometry")
  sane = "distance" %in% names(dset$effort)
  if ( message & !sane) { message = "The effort has not distance column"  }
  return(list(sane = sane, message = message))
}


# Has all the columns: strat, trans, seg, det
sanity.columns = function(dset, message = FALSE) {
  cols = c("strat", "trans", "seg", "det")
  sane = !any(!(cols %in% names(dset$effort)))
  if ( message & !sane) { message = paste0("The data set is missing one of the following colums: ", do.call(paste, c(as.list(cols),sep = ", "))) }
  return(list(sane = sane, message = message))
}

# Segments sorted?
sanity.unsorted = function(dset, message = FALSE) {
  sane = is.sorted(segment.id(dset$effort))
  if ( message & !sane ) { message = "Segments are not sorted" }
  return(list(sane = sane, message = message))
}


# All detections inside mesh ?
sanity.detection.inside = function(dset, message = FALSE){
  sane = !any(!(is.inside(dset$mesh, detdata(dset), mesh.coords = dset$mesh.coords)))
  if ( message & !sane) { message = "There are detections outside the mesh." }
  return(list(sane = sane, message = message))
}

# All effort inside mesh?
sanity.effort.inside = function(dset, message = FALSE){
  sane = !any(!(is.inside(dset$mesh, dset$effort[is.na(dset$effort$det),], mesh.coords = paste0("start.",dset$mesh.coords))))
  if ( message & !sane) { message = "There is effort outside the mesh." }
  return(list(sane = sane, message = message))
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
    msk =  is.inside.polygon(whales$mesh, whales$inner.boundary,dsdata$mesh$loc, dsdata$mesh.coords)
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

sanity = function(...){UseMethod("sanity")}
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
  if (!is.sorted(transect.id(effort))) { stop("The transects of your data set are not sorted in increasing order. Fix that.") }
  # end.idx = findInterval(1:length(levels(effort$tr)),as.numeric(effort$tr))
  # end.idx = findInterval(1:max(transect.id(effort)), transect.id(effort))
  end.idx = which(diff(transect.id(effort))>0)
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
  if (!is.sorted(transect.id(effort))) { stop("The transects of your data set are not sorted in increasing order. Fix that.") }
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
  det.srt = as.numeric(unlist(lapply(strsplit(as.character(final$det),"[.]"), function(x) {return(x[length(x)])})))
  seg.srt = as.numeric(unlist(lapply(strsplit(as.character(final$seg),"[.]"), function(x) {return(x[length(x)])})))
  final = final[order(det.srt),]
  final = final[order(seg.srt),]
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
  cols = names(data$effort) %in% paste0("start.",data$mesh.coords)
  if (keep){ return(data$effort[tr$start,]) }
  else {
    pts = strip.coords(data$effort[tr$start,cols])
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
  cols = names(data$effort) %in% paste0("end.",data$mesh.coords)
  if (keep){ return(data$effort[tr$end,]) } 
  else {
    pts = strip.coords(data$effort[tr$end,cols])
    class(pts) = c(data$geometry,"data.frame")
    return(pts)
  } 
}

#' Number of transects
#' 
#' @aliases numtr.transect
#' @export
#' @param transects
#' @param data Data set with $effort
#' @return n
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

numtr.transect = function(tr,data=data) { return(dim(tr)[1]) }


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


#' Numeric statum id
#' 
#' @aliases stratum.id
#' @export
#' @param effort
#' @return id
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

stratum.id = function(effort) { as.numeric(effort$strat) }


#' Numeri transect id
#' 
#' @aliases transect.id
#' @export
#' @param effort
#' @return id
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

transect.id = function(effort) { as.numeric(gsub("^.*\\.","", as.character(effort$trans))) }

#' Numeric segment id
#'  
#' @aliases segment.id
#' @export
#' @param effort
#' @return id
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

segment.id = function(effort) { as.numeric(gsub("^.*\\.","", as.character(effort$seg))) }


#' Numeric detection id
#'  
#' @aliases detection.id
#' @export
#' @param effort
#' @return id
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

detection.id = function(effort) { as.numeric(gsub("^.*\\.","", as.character(effort$det))) }
