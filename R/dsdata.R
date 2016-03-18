#' iDistance data format for distance sampling surveys
#' 
#' For details on ho to construct such a data set see \link{make.dsdata}
#' 
#' Data sets included in iDistance:
#' \itemize{
#'  \item{\link{toy1}: }{A toy data set with a single transect line}
#'  \item{\link{weeds}: }{Johnson et al. (2010) weeds data (line transect)}
#'  \item{\link{whales}: }{ETP blue whales (line transect)}
#'  \item{\link{strdolphin}: }{ETP striped dolphins (line transect)}
#'  \item{\link{sbdolphin}: }{ETP short beaked dolphins (line transect)}
#'  \item{\link{mrsea}: }{Marine renewables strategic environmental assessment (line transect, see package MRSea)}
#' }
#' 
#' Data sets that can be imported from other packages
#' \itemize{
#'  \item{mexdolphin (dsm package):}{ Use \link{import.dsmdata} }
#' }
#' 
#' @name dsdata
#' 
NULL


#' Check if an object is of class dsdata
#' 
#' @aliases is.dsdata
#' @name is.dsdata
#' @export
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 

is.dsdata = function(data) { return( class(data)[1] == "dsdata") }


#' Make dsdata 
#' 
#' Constructs a distance sampling data set from transects, detections and model of the survey area
#'
#' The core item that this function needs is a \code{data.frame}, hereafter called \code{effort}
#' containing transect definitions, segments thereof, detections and possibly some covariates.
#' The effort table needs to have the following columns:
#' 
#'  \itemize{
#'  \item{strat: }{ An identifier of the stratum }
#'  \item{trans: }{ An identifier of the transect }
#'  \item{seg: }{An identifier of a segment}
#'  \item{det: }{An identifier of a detection}
#'  \item{start.x}{x-coordinate of a transect's/semgent's starting point}
#'  \item{start.y}{y-coordinate of a transect's/semgent's starting point}
#'  \item{end.x}{x-coordinate of a transect's/semgent's end point}
#'  \item{end.y}{y-coordinate of a transect's/semgent's end point}
#'  \item{x}{x-coordinate of a detection (NA for transects and segments)}
#'  \item{y}{y-coordinate of a detection (NA for transects and segments)}
#'  \item{distance}{distance of the detection from the observer}
#'  }
#'  
#'  Note that in pinciple the coordinates can have arbitrary names which can be set by the mesh.coords
#'  parameter explained below. For instance, if \code{mesh.coords = c("lon",("lat")}, the columns above
#'  must be called 'lon' and 'lat' for detections and start.lon etc. for trensects and segments.
#'  
#'  If no \link{inla.mesh} is provided via the \code{mesh} parameter, \code{make.dsdata} will try to
#'  construct a mesh from the effort table. By default, an inner boundary segment is constructed from the
#'  non-convex hull of all detections and transect/segment points. Thereafter, a convex hull is found
#'  to create the mesh. As this might not reflect the reality of your data please consider to construct
#'  a mesh manually.
#'  
#'  
#' @aliases make.dsdata
#' @name make.dsdata
#' @export
#' @param effort A data.frame object describing transects and detections
#' @param geometry Character describing the geometry of the data set. Either "euc" for euclidean or "geo" for geographic coordinates.
#' @param mesh An inla.mesh object modeling the domain.
#' @param mesh.coords Names of the effort columns that are interpreted as coordinates
#' @param mesh.args If no mesh if provided these arguments passed on to inla.mesh.create to construct the mesh. 
#' @return  a \link{dsdata} object
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


#' Plot distance sampling data (\link{dsdata})
#' 
#' @aliases plot.dsdata 
#' @name plot.dsdata
#' @export
#' @param data A \link{dsdata} object
#' @param rgl If TRUE, use RGL to plot the data
#' @param ggp if TRUE, use ggplot2 for plotting
#' @param add If TRUE, add the plot instead of starting a new one
#' @param segment If TRUE, plot the segment lines
#' @param segment.args Plotting argument for segments.
#' @param segment.colorize Colorize segments according to their transect. Only supported if rgl=FALSE
#' @param detection If TRUE, plot the detected animals
#' @param detection.args Plotting argument for detections.
#' @param add.mesh If TRUE, add the mesh to the plot (only ggplot)
#' @param col Color specification for the mesh. Requires rgl=TRUE
#' @param asp Apect ration of the plot
#' @param transect DEPRECATED. Use segment instead.
#' @param transect.args DEPRECATED. Use segment instead.
#' @examples \dontrun{ toy=toy1() ; plot(toy) }
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
plot.dsdata = function(data, 
                       rgl = FALSE,
                       ggp = TRUE,
                       add=FALSE,
                       transect = FALSE,
                       transect.args = list(lwd = 2, col = rgb(0.6, 0.6, 0.6)),
                       segment = TRUE,
                       segment.args = list(lwd = 1, color = "turquoise4"),
                       segment.colorize = FALSE,
                       detection = TRUE,
                       detection.args = list(color = "red2", size = 1),
                       add.mesh = TRUE,
                       col = NULL, 
                       asp = 1, ...){

  # Defaults
  
  if ( !("lwd" %in% names(segment.args)) ) { segment.args$lwd = 0.8 }
  if ( !("color" %in% names(segment.args)) ) { segment.args$color = "turquoise4" }
  if ( !("color" %in% names(detection.args)) ) { detection.args$color = "red2" }
  
  if (rgl==FALSE) {
    
    # 2D plot of spherical mesh
    if ( data$geometry == "geo" & data$mesh$manifold == "S2") {
      # Transform coordinates of mesh before plotting
      loc = data.frame(data$mesh$loc)
      R = sqrt(sum(loc[1,]^2))
      colnames(loc)  = c("x","y","z")
      mesh = data$mesh
      mesh$manifold = "R2"
      mesh$loc = as.matrix(euc.to.geo(loc, R = R)[,data$mesh.coords])
    }  else {
      mesh = data$mesh
    }
    
    if ( ggp ) {
      
      # Check for package
      if ( !requireNamespace("ggplot2", quietly = TRUE) ) { stop("This function requires the ggplot2 package.")}
      
      gg = ggplot()
      
      # Plot the mesh
      if ( add.mesh ) { gg = gg + gg.mesh(data) + gg.bnd(data) + gg.int(data) }
      
      # Plot segments
      if ( segment ) { gg = gg + do.call(gg.seg, c(list(data), segment.args)) }

      # Plot detections
      if ( detection ) { gg = gg + do.call(gg.det, c(list(data), detection.args)) }
      
      # Set labels
      gg = gg + xlab(data$mesh.coords[1]) + ylab(data$mesh.coords[2]) + coord_fixed()

      return(gg)
      
    } else {
      
      # Plot mesh
      plot(mesh, main="", asp = asp, add = add)
      
      # Axis
      if ( !add ) { axis(1); axis(2) ; box()}
      
      # Plot segments
      if ( segment | segment.colorize ) {
        if ( segment.colorize ){ segment.args$col = data$effort$trans }
        sp = startpoint(as.segment(data),data)
        ep = endpoint(as.segment(data),data)
        do.call(segments, c(list(sp[,1],sp[,2],ep[,1],ep[,2]), segment.args))
        
  
      }
      
      # Plot transect lines
      if ( transect ) {
        spoint = startpoint(as.transect(data),data)
        epoint = endpoint(as.transect(data),data)
        do.call(segments, c(list(spoint[, data$mesh.coords[1]], spoint[, data$mesh.coords[2]], 
                                 epoint[, data$mesh.coords[1]], epoint[, data$mesh.coords[2]]), 
                                 transect.args))
      }
    
      # Plot detections
      if ( detection ) {
        do.call(points, c(list(detdata(data)[,data$mesh.coords]), detection.args) )
      }
    }
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
      if ( segment ) {
        sp = startpoint(as.segment(data),data)[,data$mesh.coords]
        ep = endpoint(as.segment(data),data)[,data$mesh.coords]
        sp$z = 0.001
        ep$z = 0.001
        lns = as.matrix(cbind(ep,sp,matrix(NaN, nrow = dim(sp)[1], ncol = 3)))
        qq = matrix(t(lns), ncol = 3, byrow = TRUE)
        rgl.linestrips(qq, lwd = 3, col = "black", alpha = 0.4)
      }
    }
    else if (data$geometry == "geo"){
      # Elevation for plotting on top of earth surface
      if ( data$geometry == "geo" & data$mesh$manifold == "S2") {
        R = sqrt(sum(data$mesh$loc[1,]^2))
      } else { R = 1 }
      R.delta = 1.003
      
      # Draw earth
      rgl.earth(R = R)
      rgl.viewpoint(0,-90)

      # Draw detections
      if ( detection ) {
        detection.args = c(detection.args, list(size = 4*detection.args$cex))
        detections = detdata(data)[,data$mesh.coords]
        do.call(rgl.sphpoints, c(list(long = detections[,1] + 360, 
                                      lat = detections[,2],
                                      radius = R * R.delta)))
      }
      
      # Draw the mesh
      if (is.null(col)){
        rgl.sphmesh(data$mesh, add = TRUE, radius = R*R.delta, lwd = 1, edge.color = rgb(0,0,0), alpha = 0.3, draw.segments = TRUE)
      } else {
        rgl.sphmesh(data$mesh, add = TRUE, radius = R*R.delta, lwd = 1, edge.color = rgb(0,0,0), col = col)
      }
      
      # Draw segments
      if ( segment ) {
        if ( segment.colorize ){ segment.args$col = data$effort$trans }
        pseudo.transect = data.frame(start=1:nrow(data$effort), end=1:nrow(data$effort))
        class(pseudo.transect) = c("transect", "data.frame")
        do.call(plot, c(list(pseudo.transect, data = data, rgl = TRUE, add = TRUE, radius = R*R.delta), segment.args) )
      }
      
      # Draw transects
      if ( transect ) {
        do.call(plot, c(list(as.transect(data$effort), data = data, rgl = TRUE, add = TRUE, radius = R*(R.delta-0.001)), transect.args) )
      }
                
    }
  }
}

#' Change the coordinate system of \link{dsdata}
#' 
#' @aliases remap.dsdata 
#' @export
#' @param dsdata Distance sampling data set
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

remap.dsdata = function(data, p4s = "+proj=longlat", mesh.p4s = "+proj=longlat"){
  
  # We need the proj4 strings for the effort and the mesh
  if ( is.null(data$p4s) | is.null(data$mesh.p4s) ) { 
    stop("Your data set has no p4s or mesh.p4s entry. These are required for remapping") 
  }
  # rgdal is needed for the projection
  if ( !require("rgdal") ) { stop("This function requires the rgdal package")}
  
  # Our projection function
  projFun = function(df, from.p4s, to.p4s) {
    class(df) = "data.frame"
    msk = is.finite(df[,1]) & is.finite(df[,2])
    spts <-SpatialPoints(df[msk,], proj4string = CRS(from.p4s))
    geo.pts <-spTransform(spts, CRS(to.p4s))
    df[msk,] = as.data.frame(geo.pts)
    return(df)
  }
  
  # Project detections
  data$effort[,data$mesh.coords] = projFun(data$effort[,data$mesh.coords], data$p4s, p4s)
  
  # Project segment start end end points
  data$effort[,paste0("start.", data$mesh.coords)] = projFun(data$effort[,paste0("start.", data$mesh.coords)], data$p4s, p4s)
  data$effort[,paste0("end.", data$mesh.coords)] = projFun(data$effort[,paste0("end.", data$mesh.coords)], data$p4s, p4s)
  
  # If we are projecting to lon/lat change column names and mesh.coords entry
  if (p4s == "+proj=longlat") {
    names(data$effort)[which(names(data$effort) %in% data$mesh.coords)] = c("lon","lat")
    names(data$effort)[which(names(data$effort) %in% paste0("start.", data$mesh.coords))] = c("start.lon","start.lat")
    names(data$effort)[which(names(data$effort) %in% paste0("end.", data$mesh.coords))] = c("end.lon","end.lat")
    data$mesh.coords = c("lon","lat")
  }
  # Project mesh
  data$mesh$loc[,c(1,2)] = as.matrix(projFun(as.data.frame(data$mesh$loc[,c(1,2)]), data$mesh.p4s, mesh.p4s))
  
  # Set new proj4 string
  data$p4s = p4s
  data$mesh.p4s = p4s
  
  # Return
  return(data)
}

#' Sanity of a \link{dsdata} object
#' 
#' Analyzes if a distance sampling data object is well defined. 
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
  
  cat("Checking if transects are sorted ...")
  sanity.stop(sanity.unsorted.transect(dset, message = TRUE))
  
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
sanity.unsorted.transect = function(dset, message = FALSE) {
  sane = !is.unsorted(transect.id(dset$effort), na.rm=TRUE)
  if ( message & !sane ) { message = "Transects are not sorted" }
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



#' Summarize distance sampling data 
#' 
#' Gives a brief summary of the data and calls the \link{statistics.dsdata} function to print some statistics.
#' 
#' @aliases summary.dsdata 
#' @export
#' @param data A \link{dsdata} object
#' @param ... arguments passed on to \link{statistics.dsdata}
#' @examples \dontrun{ data(toy1) ; summary(toy1) }
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

summary.dsdata = function(data, ...){
  
  cat("#-------- data fields --------#\n")
  cat(paste0("Geometry: "),data$geometry,"\n")
  cat(paste0("Coordinate names: "),data$mesh.coords,"\n")
  cat(paste0("Effort columns names: "),"\n")
  cat(paste0(" ",names(data$effort)), "\n")
  
  cat("#-------- statistics --------#\n")

  statistics.dsdata(data, ...)
  
}

#' Calculate baseline statistics
#' 
#' A brief overview of important statistics like the area covered by the mesh and the number of observations
#' 
#' @aliases statistics.dsdata 
#' @name statistics.dsdata 
#' @export
#' @param dsdata A \link{dsdata} object
#' @param distance.truncation Set this to a numeric value to obtain more detailed statistics.
#' @examples \dontrun{ data(toy1) ; statistics(toy1) }
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

statistics.dsdata = function(data, distance.truncation = NULL){
  
  #
  # Calculate the stats
  #
  
  n.det = nrow(detdata(data))
  n.seg = nrow(as.segment(data))
  len.seg = sum(distance(startpoint(as.segment(data), data), endpoint(as.segment(data), data), geometry = data$geo))
  n.tr = nrow(as.transect(data))
  len.tr = sum(len.transect(as.transect(data),data))
  
  if (data$geometry == "geo") {
    if (data$mesh$manifold == "S2"){
      A.mesh = (6371^2)*sum(Matrix::diag(inla.mesh.fem(tmp.mesh, order=1)$c0))
      unit.mesh = "km"
    } else if (data$mesh$manifold == "R2") { # Special case: the mesh is R2 but coordinates are lon/lat
      tmp.mesh = data$mesh
      loc = data.frame(tmp.mesh$loc[,c(1,2)]); colnames(loc) = data$mesh.coords
      eloc = as.matrix(geo.to.euc(loc, R = 1))
      tmp.mesh$loc = as.matrix(geo.to.euc(loc, R = 1))
      tmp.mesh.manifold = "S2"
      A.mesh = (6371^2)*sum(Matrix::diag(inla.mesh.fem(tmp.mesh, order=1)$c0)) # in km^2
      unit.mesh = "km"
    }
  } else if ( data$geometry == "euc" ){
    A.mesh = sum(Matrix::diag(inla.mesh.fem(data$mesh, order=1)$c0))
    unit.mesh = "mesh units"  
  }
  
  #
  # Output
  #
  
  cat(paste0("Number of transects: N_tr ="), n.tr, "\n")
  cat(paste0("Total length of transects: L ="), len.tr, "units","\n")
  cat(paste0("Number of segments: N_s ="), n.seg, "\n")
  cat(paste0("Total length of segments: L_s ="), len.seg, "units","\n")
  cat(paste0("Number of detections N = "), n.det, "\n")
  cat(paste0("Total area of mesh A_m ="), A.mesh, " square ",unit.mesh,"\n")
  
  if ( !is.null(distance.truncation) ) {
    A.eff = 2 * distance.truncation * len.tr
    rate.eff = n.det/A.eff
    cat(paste0("Distance truncation: Z = "), distance.truncation, "\n")
    cat(paste0("Effort area (transect based): A_eff := 2 * Z * L = "), A.eff, "\n")
    cat(paste0("On-effort rate (transect based): R_e := N/A_eff ="), rate.eff, " =>  log(R_e) =", log(rate.eff),"\n")
    cat(paste0("Expected detections on mesh (approximate): N_m := A_m * R_e  = "), A.mesh * rate.eff, "per square", unit.mesh ,"\n")
    cat("--> Only rely on this if mesh units = transect units!")
  }
  
  #   # Survey area
  #   if ("inner.boundary" %in% names(data)) {
  #     msk =  is.inside.polygon(whales$mesh, whales$inner.boundary,data$mesh$loc, data$mesh.coords)
  #     A = sum(Matrix::diag(inla.mesh.fem(data$mesh, order=1)$c0) * msk)
  #     cat(paste0("Total survey area: A ="), A, " square mesh units","\n")
  #   }
  
}

remap = function(...){UseMethod("remap")}
sanity = function(...){UseMethod("sanity")}
statistics = function(...){UseMethod("statistics")}
as.transect = function(...){UseMethod("as.transect")}
as.segment = function(...){UseMethod("as.segment")}
as.detection = function(...){UseMethod("as.detection")}
startpoint = function(...){UseMethod("startpoint")}
endpoint = function(...){UseMethod("endpoint")}
detdata = function(...){UseMethod("detdata")}
segdata = function(...){UseMethod("segdata")}
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

segdata.dsdata = function(data,...){ 
 return(data$effort[is.na(data$effort$det),])
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
  if (is.unsorted(transect.id(effort),na.rm=TRUE)) { stop("The transects of your data set are not sorted in increasing order. Fix that.") }
  tr.id = transect.id(effort)
  if (tr.id[length(tr.id)]==1) {
    end.idx = length(tr.id)
    start.idx = 1
  } else {
    end.idx = c(which(diff(tr.id)>0), length(tr.id))
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
  if (is.unsorted(transect.id(effort), na.rm = TRUE)) { stop("The transects of your data set are not sorted in increasing order. Fix that.") }
  wh = which(is.na(effort$det))
  seg = data.frame(start = wh, end = wh)
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
