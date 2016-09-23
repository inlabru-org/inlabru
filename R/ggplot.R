# This file contains several helpers to plot distance sampling data using ggplot2
# gg.map() : Plot a ggmap covering the area some data sits in
# gg.point(): Point geom for Spatial* objects
# gg.segment(): Segment geom for Spatial* objects
# gg.mesh() : A geom for plotting the mesh
# gg.bnd() : A geom for plotting the boundary segment of a mesh
# gg.int() : A geom for plotting the interior segments of a mesh
# gg.seg() : A geom for plotting segment lines
# gg.det() : A geom for plotting detections
# gg.swath() : A geom for plotting the segment area


#' Plot a gg-map covering the area a spatial object resides in
#' 
#' @aliases gg.map
#' @name gg.map
#' @export
#' @param data A Spatial* object
#' @return A ggmap
#' 
gg.map = function(data, ...) {
  data = spTransform(data, CRS("+proj=longlat"))
  df = cbind(coordinates(data), data@data)
  
  # Figure out a sensible bounding box (range of data plus 30%)
  extend = 1.3
  cloc = apply(coordinates(data), MARGIN = 2, mean)
  crange = apply(coordinates(data), MARGIN = 2, range)
  lonlim = (extend*(crange[,1] - cloc[1])) + cloc[1]
  latlim = (extend*(crange[,2] - cloc[2])) + cloc[2]
  
  # Create map
  myMap = get_map(c(lonlim[1], latlim[1], lonlim[2], latlim[2]), ...)
  
  # Return map
  ggmap(myMap)
}

#' Point geom for Spatial* objects
#' 
#' @aliases gg.point
#' @name gg.point
#' @export
#' @param data A Spatial* object
#' @return geom_point
#' 
gg.point = function(data, ...) {
  # data = spTransform(data, CRS("+proj=longlat"))
  df = data.frame(coordinates(data))
  geom_point(data = df, aes_string(x = coordnames(data)[1], y = coordnames(data)[2]), color = "red", ...)
}

#' Segment geom for Spatial* objects
#' 
#' @aliases gg.segment
#' @name gg.segment
#' @export
#' @param data A Spatial* objectt
#' @return geom_segment
#' 
gg.segment = function(data, ...) {
  # data = spTransform(data, CRS("+proj=longlat"))
  qq = coordinates(data)
  cnames = coordnames(data)
  sp = do.call(rbind, lapply(qq, function(k) do.call(rbind, lapply(k, function(x) x[1:(nrow(x)-1),]))))
  ep = do.call(rbind, lapply(qq, function(k) do.call(rbind, lapply(k, function(x) x[2:(nrow(x)),]))))
  colnames(sp) = paste0("start.", cnames)
  colnames(ep) = paste0("end.", cnames)
  df = data.frame(cbind(sp, ep))
  geom_segment(data = df, aes_string(x = paste0("start.", cnames)[1], 
                                     y = paste0("start.", cnames)[2], 
                                     xend = paste0("end.", cnames)[1], 
                                     yend = paste0("end.", cnames)[2]),
               color = "green", ...)  
} 


#' Plot inla.mesh using ggplot2
#' 
#' @aliases gg.mesh
#' @name gg.mesh
#' @export
#' @param data A \link{dsdata} object or a \link{inla.mesh}
#' @param color Color of the mesh
#' @return A ggplot2 object
#' 

gg.mesh = function(data, color = rgb(0,0,0,0.1), ...) {
  if ( class(data)[1] == "inla.mesh" ) { mesh = data } else { mesh = data$mesh }
  if ( mesh$manifold == "S2" ) { stop("Geom not implemented for spherical meshes (manifold = S2)" ) }
  df = rbind(data.frame(a=mesh$loc[mesh$graph$tv[,1],c(1,2)],b=mesh$loc[mesh$graph$tv[,2],c(1,2)]),
               data.frame(a=mesh$loc[mesh$graph$tv[,2],c(1,2)],b=mesh$loc[mesh$graph$tv[,3],c(1,2)]),
               data.frame(a=mesh$loc[mesh$graph$tv[,1],c(1,2)],b=mesh$loc[mesh$graph$tv[,3],c(1,2)]))
  colnames(df) = c("x","y","xend","yend")
  gg = geom_segment(data = df, aes(x = x,y = y, xend = xend, yend = yend), color = color, ...)
  return(gg)
}


#' Plot inla.mesh using ggplot2
#' 
#' @aliases gg.col
#' @name gg.col
#' @export
#' @param data A \link{dsdata} object or a \link{inla.mesh}
#' @param color on mesh vertices
#' @param alpha on mesh vertices
#' @return A ggplot2 object
#' 

gg.col = function(data, color, alpha = 1, ...) {

  if ( class(data)[1] == "inla.mesh" ) { mesh = data } else { mesh = data$mesh }
  if ( length(alpha) == 1 ) { 
    alpha = rep(alpha,length(color))
    alph.fix = TRUE}
  else { alph.fix = FALSE }
  
  refine = TRUE
  if ( refine ) {
    omesh = mesh
    mesh = mesh.refine(mesh, refine = list(max.edge = diff(range(omesh$loc[,1]))/100))
    color = as.vector(inla.spde.make.A(omesh, loc = mesh$loc) %*% color)
    if ( !is.null(alpha) ) { alpha = as.vector(inla.spde.make.A(omesh, loc = mesh$loc) %*% alpha) }
  } 
  
  df = matrix(NA, nrow(mesh$graph$tv)*3, 2)
  df[seq(1, nrow(mesh$graph$tv)*3, by=3), ] = mesh$loc[mesh$graph$tv[,1],1:2]
  df[seq(2, nrow(mesh$graph$tv)*3, by=3), ] = mesh$loc[mesh$graph$tv[,2],1:2]
  df[seq(3, nrow(mesh$graph$tv)*3, by=3), ] = mesh$loc[mesh$graph$tv[,3],1:2]
  df = data.frame(df, tr = as.vector(matrix(1:nrow(mesh$graph$tv), 3, nrow(mesh$graph$tv),byrow=TRUE)))
  colnames(df) = c("x","y","tr")
  icol  = (color[mesh$graph$tv[,1]] + color[mesh$graph$tv[,2]] + color[mesh$graph$tv[,3]])/3
  ialpha  = (alpha[mesh$graph$tv[,1]] + alpha[mesh$graph$tv[,2]] + alpha[mesh$graph$tv[,3]])/3
  df$col = as.vector(matrix(icol, 3, nrow(mesh$graph$tv),byrow=TRUE))
  df$alph = as.vector(matrix(ialpha, 3, nrow(mesh$graph$tv),byrow=TRUE))
  
  aest = aes_string(x = "x", 
             y = "y",
             group = "tr", 
             fill = "col",
             alpha = df$alph)
  if (alph.fix) { 
    gg = geom_polygon(data = df, aest, linetype = "blank", show.legend = FALSE, alpha = alpha[1])
  } else {
    gg = geom_polygon(data = df, aest, linetype = "blank", show.legend = FALSE)
  }
  
  return(gg)
}


#' Plot exterior boundary of a mesh
#'
#' @aliases gg.bnd
#' @name gg.bnd
#' @export
#' @param data A \link{dsdata} object or a \link{inla.mesh}
#' @param mapping A set of aesthetics mappings created by \link{aes} or \link{aes_}.
#' @param ... Arguments passed on to \link{geom_segment}
#' @return A ggplot2 object
#' 

gg.bnd = function(data, mapping = NULL, ...) {
  if ( class(data)[1] == "inla.mesh" ) { mesh = data } else { mesh = data$mesh }
  # Make data frame
  df = data.frame(mesh$loc[mesh$segm$bnd$idx[,1],1:2], mesh$loc[mesh$segm$bnd$idx[,2],1:2])
  colnames(df) = c("x","y","xend","yend")
  # Make aesthetics
  mp = aes(x = x,y = y, xend = xend, yend = yend)
  # Add user aesthetics
  if ( !is.null(mapping) ) { mp = modifyList(mp, mapping) }
  # If user did not provide color, set color to red
  more.args = list(...)
  nms = c(names(mapping), names(more.args))
  if ( !("color" %in% nms | "colour" %in% nms | "col" %in% nms) )  { more.args$color = rgb(1,0,0,0.5) }
  # Make geom
  gg = do.call(geom_segment, c(list(data = df, mapping = mp), more.args))
  return(gg)
}


#' Plot interior boundaries of a mesh
#'
#' 
#' @aliases gg.int
#' @name gg.int
#' @export
#' @param data A \link{dsdata} object or a \link{inla.mesh}
#' @param mapping A set of aesthetics mappings created by \link{aes} or \link{aes_}.
#' @param ... Arguments passed on to \link{geom_segment}
#' @return A ggplot2 object
#' 

gg.int = function(data, mapping = NULL, ...) {
  if ( class(data)[1] == "inla.mesh" ) { mesh = data } else { mesh = data$mesh }
  # Make data frame
  df = data.frame(mesh$loc[mesh$segm$int$idx[,1],1:2], mesh$loc[mesh$segm$int$idx[,2],1:2])
  if ( nrow(df) == 0 ) { return(NULL) }
  colnames(df) = c("x","y","xend","yend")
  # Make aesthetics
  mp = aes(x = x,y = y, xend = xend, yend = yend)
  # Add user aesthetics
  if ( !is.null(mapping) ) { mp = modifyList(mp, mapping) }
  # If user did not provide color, set color to blue
  more.args = list(...)
  nms = c(names(mapping), names(more.args))
  if ( !("color" %in% nms | "colour" %in% nms | "col" %in% nms) )  { more.args$color = rgb(0,0,1,0.5) }
  # Make geom
  gg = do.call(geom_segment, c(list(data = df, mapping = mp), more.args))
  return(gg)
}


#' Plot distance sampling segments using ggplot2
#'
#' @aliases plot gg.seg
#' @export
#' @param data a \link{dsdata} object
#' @param mapping A set of aesthetics mappings created by \link{aes} or \link{aes_}.
#' @param ... Arguments passed on to \link{geom_segment}
#' @return A ggplot2 object
#' 

gg.seg = function(data, mapping = NULL, ...) {
  # Make the data frame
  df = data$effort[as.segment(data)[,1],]
  # Default aesthetics
  mp = aes_string(x = paste0("start.", data$mesh.coords[1]), 
                  y = paste0("start.", data$mesh.coords[2]), 
                  xend = paste0("end.", data$mesh.coords[1]), 
                  yend = paste0("end.", data$mesh.coords[2]))
  # Add user aesthetics
  if ( !is.null(mapping) ) { mp = modifyList(mp, mapping) }
  # Make geom
  gg = geom_segment(data = df, mapping = mp, ...)
  return(gg)
}


#' Plot detections using ggplot2
#'
#' 
#' @aliases gg.det
#' @export
#' @param data a \link{dsdata} object
#' @param mapping A set of aesthetics mappings created by \link{aes} or \link{aes_}.
#' @param ... Arguments passed on to \link{geom_point}
#' @return A ggplot2 object
#' 

gg.det = function(data, mapping = NULL, ...) {
  # Make data frame
  df = detdata(data)
  # Default aesthetics
  mp = aes_string(x = data$mesh.coords[1], y = data$mesh.coords[2])
  # Add user aesthetics
  if ( !is.null(mapping) ) { mp = modifyList(mp, mapping) }
  # If user did not provide color, set color to red
  more.args = list(...)
  nms = c(names(mapping), names(more.args))
  if ( !("color" %in% nms | "colour" %in% nms | "col" %in% nms) )  { more.args$color = rgb(1,0,0,1) }
  # Make geom
  gg = do.call(geom_point, c(list(mapping = mp, data = df), more.args))
  return(gg)
}

#' Plot segment swath
#'
#' 
#' @aliases gg.swath
#' @export
#' @param data a \link{dsdata} object
#' @param mapping A set of aesthetics mappings created by \link{aes} or \link{aes_}
#' @param width The width of the swath
#' @param ... Arguments passed on to \link{geom_polygon}
#' @return A ggplot2 object
#' 

gg.swath = function(data, mapping = NULL, width = NULL, ...) {
  
  if ( is.null(width) ) { 
    width = segdata(data)$width
    }
  
  # Make data frame
  sp = startpoint(as.segment(data),data)[,data$mesh.coords]
  ep = endpoint(as.segment(data),data)[,data$mesh.coords]
  v = ep-sp
  vo = data.frame(x = v[,2],y = -v[,1])
  v1 = data.frame(sp + width/2 * normalize.euc(vo), seg = 1:dim(vo)[1])
  v2 = data.frame(sp - width/2 * normalize.euc(vo), seg = 1:dim(vo)[1])
  v3 = data.frame(ep - width/2 * normalize.euc(vo), seg = 1:dim(vo)[1])
  v4 = data.frame(ep + width/2 * normalize.euc(vo), seg = 1:dim(vo)[1])
  df = rbind(v1, v2, v3, v4)
  # Default aesthetics
  mp = aes_string(x = data$mesh.coords[1], y = data$mesh.coords[2], group = "seg")
  # Add user aesthetics
  if ( !is.null(mapping) ) { mp = modifyList(mp, mapping) }
  more.args = list(...)
  
  if ( !("fill" %in% names(mapping)) )  { more.args$fill = rgb(0,0,0,0.1) }
  gg = do.call(geom_polygon, c(list(data = df, mapping = mp),list(...)))

  return(gg)
}



