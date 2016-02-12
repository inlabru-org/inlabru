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
#' @param effort A data.frame object describing transects and detections
#' @param geometry Character describing the geometry of the data set. Either "euc" for euclidean or "geo" for geographic coordinates.
#' @param mesh An inla.mesh object modeling the domain.
#' @param mesh.coords Names of the effort columns that are interpreted as coordinates
#' @param mesh.args If no mesh if provided these arguments passed on to inla.mesh.create to construct the mesh. 
#' @examples \dontrun{ data(toy1) ; dset = make.dsdata(effort = toy1$effort) ; plot(dset) }
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 

make.dsdata = function(effort = NULL, geometry = NULL, mesh = NULL, mesh.coords = NULL, mesh.args = list()) {
  
  dset = list()
  class(dset) = c("dsdata", "list")
  
  # EFFORT
  if (is.null(effort)) { 
    stop("make.dsdata(): You have to provide an effort data frame.") 
  } else {
    class(effort) = c("effort", "data.frame")
    dset$effort = effort
    cat("make.dsdata(): Checking if the effort table has all the necessary columns ...")
    sanity.stop = function(ret) { if (!ret$sane) { stop(ret$message) } else { cat(" [OK]\n") } }
    sanity.stop(sanity.columns(dset, message = TRUE))
  }
  
  # GEOMETRY
  if ( is.null(geometry) ) { 
    cat("make.dsdata(): You did not provide a geometry. ")
    if ( "lon" %in% colnames(effort) & "lat" %in% colnames(effort)) {
      cat("However, your effort table has longitude and latitude data. Assuming geographic coordinates and setting geometry to 'geo'.\n")
      geometry = "geo"
    } else if ( "x" %in% colnames(effort) & "y" %in% colnames(effort)){
      cat("However, your effort table has 'x' and 'y' columns. Assuming euclidean coordinates and setting geometry to 'euc'.\n")
      geometry = "euc"
    }
  } else if ( !any((geometry %in% c("euc", "geo"))) ){
    stop(paste0("You have to provide a geometry that is not supported: ", geometry)) 
  }
  dset$geometry = geometry
  
  # MESH COORDS
  if ( is.null(mesh.coords) ) { 
    cat("make.dsdata(): You did not provide mesh.coords")
    if ( geometry=="geo" ) {
      cat(" but geometry is set to 'geo'. Will assume defaults 'lon' and 'lat' for mesh.coords.\n")
      mesh.coords = c("lon","lat")
      if (!( "lon" %in% colnames(effort) & "lat" %in% colnames(effort))) { 
        stop("make.dsdata(): Your effort data has not 'lon' and 'lat' columns.")
      }
    } else if ( geometry=="euc" ) {
      cat(" but geometry is set to 'euc'. Will assume defaults 'x' and 'y' for mesh.coords.\n")
      mesh.coords = c("x","y")
      if (!( "x" %in% colnames(effort) & "y" %in% colnames(effort))) { 
        stop("make.dsdata(): Your effort data has not 'x' and 'y' columns.")
      }
    } else {
      stop(".")
    }
  }
  dset$mesh.coords = mesh.coords
  
  if (is.null(mesh)) { 
    cat("make.dsdata(): You did not provide a mesh. A default mesh will be constructed.\n")
    mesh = do.call(make.mesh, c(list(dset), mesh.args))
  }
  dset$mesh = mesh
    
  #
  # Sanity check
  #
  cat("make.dsdata(): Done with constructing the data set. Running sanity checks...\n")
  sanity(dset)
  
  #
  # Done.
  # 
  
  return(dset)
  
}


#' Make a default mesh
#' 
#' @aliases make.mesh
#' @export
#' @examples \dontrun{  }
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 

make.mesh = function(dset, ...) {
  
  # Determine range of coordinates
  det.x1.range = range(dset$effort[,dset$mesh.coords[1]], na.rm = TRUE)
  det.x2.range = range(dset$effort[,dset$mesh.coords[2]], na.rm = TRUE)
  eff.x1.range = range(dset$effort[,paste0("start.",dset$mesh.coords[1])], na.rm = TRUE)
  eff.x2.range = range(dset$effort[,paste0("start.",dset$mesh.coords[2])], na.rm = TRUE)
  x1.range = c(min(det.x1.range[1],eff.x1.range[1]), max(det.x1.range[2],eff.x1.range[2]))
  x2.range = c(min(det.x2.range[1],eff.x2.range[1]), max(det.x2.range[2],eff.x2.range[2]))
  # Minimum range
  min.range = min(x1.range[2]-x1.range[1],x2.range[2]-x2.range[1])
  max.range = max(x1.range[2]-x1.range[1],x2.range[2]-x2.range[1])
  
  mesh.args = list(...)
  
  # Default mesh for geographic coordinates 
  if ( dset$geometry == "geo" ) {
    if ( length(mesh.args) == 0) {
      bnd = data.frame(rbind(c(x1.range[1],x2.range[1]),
                      c(x1.range[1],x2.range[2]),
                      c(x1.range[2],x2.range[1]),
                      c(x1.range[2],x2.range[2])))
      colnames(bnd) = dset$mesh.coords
      eloc = geo.to.euc(bnd, R=1)
      mesh = inla.mesh.create(loc = eloc,  extend=FALSE, refine = list(max.edge = min.range/360))
    } else {
      mesh = do.call(inla.mesh.create, mesh.args)
    }
  # Default mesh for coordinates in the plane
  } else if (dset$geometry == "euc") {
    if ( length(mesh.args) == 0 ) {
      lattice = inla.mesh.lattice(x = seq(x1.range[1], x1.range[2], length.out = 2),
                                  y = seq(x2.range[1], x2.range[2], length.out = 2))
      mesh = do.call(inla.mesh.create, list(lattice = lattice, refine = list(max.edge = min.range/5)))
    } else {
      mesh = do.call(inla.mesh.create, mesh.args)
    }
  }
  return(mesh)
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
    
    # 2D plot of spherical mesh
    if ( data$geometry == "geo" & data$mesh$manifold == "S2") {
      # Transform coordinates of mesh before plotting
      loc = data.frame(data$mesh$loc)
      colnames(loc)  = c("x","y","z")
      mesh = data$mesh
      mesh$manifold = "R2"
      mesh$loc = as.matrix(euc.to.geo(loc, R = 1)[,data$mesh.coords])
    }  else {
      mesh = data$mesh
    }
    
    c1 = data$mesh.coords[1]
    c2 = data$mesh.coords[2]
    xlim=range(mesh$loc[,1])
    ylim=range(mesh$loc[,2])
    
    # Plot mesh
    plot(mesh, col = col, main="", asp = asp)
    
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
      if( is.null(col) ) { col = "white" }
      # plot mesh
      plot(data$mesh, rgl = TRUE, col = col, edge.color = rgb(0.4,0.4,0.4))
      
      # plot detections
      detections = cbind(detdata(data)[,data$mesh.coords],z=rep(0.01,nrow(detdata(data))))
      rgl.points(x=detections, z=0.001,size=5,col="red")
      
      # plot transect lines
      sp = startpoint(as.transect(data),data)[,data$mesh.coords]
      ep = endpoint(as.transect(data),data)[,data$mesh.coords]
      sp$z = 0.001
      ep$z = 0.001
      lns = as.matrix(cbind(ep,sp,matrix(NaN, nrow = dim(sp)[1], ncol = 3)))
      qq = matrix(t(lns), ncol = 3, byrow = TRUE)
      rgl.linestrips(qq, lwd = 3, col = "black", alpha = 0.4)
    }
    else if (data$geometry == "geo"){
      # Elevation for plotting on top of earth surface
      R.delta = 0.003
      
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
      
      # do.call(par3d,data$par3d.args)
      
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
  if ( message & !sane) { message = paste0("The effort table is missing one of the following colums: ", do.call(paste, c(as.list(cols),sep = ", "))) }
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
  pts = detdata(dset)[, dset$mesh.coords]
  if ( dset$geometry=="geo" ) { pts = geo.to.euc(pts, R = 1) }
  sane = all(is.inside(dset$mesh, pts, mesh.coords = colnames(pts)))
  if ( message & !sane ) { message = "There are detections outside the mesh." }
  return(list(sane = sane, message = message))
}

# All effort inside mesh?
sanity.effort.inside = function(dset, message = FALSE){
  pts1 = dset$effort[is.na(dset$effort$det), paste0("start.",dset$mesh.coords)]; colnames(pts1) = dset$mesh.coords
  pts2 = dset$effort[is.na(dset$effort$det), paste0("end.",dset$mesh.coords)]; colnames(pts2) = dset$mesh.coords
  pts = rbind(pts1, pts2)
  if ( dset$geometry=="geo" ) { pts = geo.to.euc(pts, R = 1) }
  sane = all(is.inside(dset$mesh, pts, mesh.coords = colnames(pts)))
  if ( message & !sane ) { message = "There is effort outside the mesh." }
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
effort = function(...){UseMethod("effort")}
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

effort.dsdata = function(dsdata){ 
  return(dsdata$effort[is.na(dsdata$effort$det),])
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
  tr.id = transect.id(effort)
  if (tr.id[length(tr.id)]==1) {
    end.idx = length(tr.id)
    start.idx = 1
  } else {
    end.idx = which(diff(tr.id)>0)
    start.idx = c(1,end.idx[1:(length(end.idx)-1)]+1)
  }
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
