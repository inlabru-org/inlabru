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
  extend = 2.5
  crange = apply(coordinates(data), MARGIN = 2, range)
  lonlim = (extend*(crange[,1] - mean(crange[,1]))) + mean(crange[,1])
  latlim = (extend*(crange[,2] - mean(crange[,2]))) + mean(crange[,2])
  
  # Create map
  myMap = get_map(c(lonlim[1], latlim[1], lonlim[2], latlim[2]), ...)
  
  # Return map
  ggmap(myMap)
}



#' ggplot2 geom for spatial objects
#' 
#' 
#' @aliases gg
#' @name gg
#' @export
#' @param data A Spatial* object or a mesh
#' 
gg = function(data, ...) {
  if (inherits(data, "SpatialPoints") | inherits(data, "SpatialPointsDataFrame")) gg.point(data, ...)
  else if (inherits(data, "SpatialLines") | inherits(data, "SpatialLinesDataFrame")) gg.segment(data, ...)
  else if (inherits(data, "SpatialPolygons") | inherits(data, "SpatialPolygonsDataFrame")) gg.polygon(data, ...)
  else if (inherits(data, "inla.mesh")) gg.mesh(data, ...)
}

#' ggmap geom for spatial objects
#' 
#' Coordinates are transformed to lon/lat for convenient plotting using ggmap
#' 
#' @aliases gm
#' @name gm
#' @export
#' @param data A Spatial* object or a mesh
#' 
gm = function(data, ...) { gg(data, crs = CRS("+proj=longlat"), ...) }



#' Point geom for Spatial* objects
#' 
#' @aliases gg.point
#' @name gg.point
#' @export
#' @param data A Spatial* object
#' @return geom_point
#' 
gg.point = function(data, crs = NULL, 
                    mapping = NULL, 
                    color = "#08519c", 
                    alpha = 0.5, ...) { 
  if ( !is.null(crs) ) { data = spTransform(data, crs) }
  df = data.frame(data)
  if (is.null(mapping)) { 
    mapping = aes_string(x = coordnames(data)[1], y = coordnames(data)[2])
    gp = geom_point(data = df, mapping = mapping, color = color, alpha = alpha,...)
  } else {
    gp = geom_point(data = df, mapping = mapping, ...)
  }

  return(gp)
}

#' Segment geom for Spatial* objects
#' 
#' @aliases gg.segment
#' @name gg.segment
#' @export
#' @param data A Spatial* objectt
#' @return geom_segment
#' 
gg.segment = function(data, crs = NULL,  color = "black", ...) {
  if ( !is.null(crs) ) { data = spTransform(data, crs) }
  
  qq = coordinates(data)
  cnames = coordnames(data)
  if (is.null(cnames)) { cnames = c("x","y") }
  sp = do.call(rbind, lapply(qq, function(k) do.call(rbind, lapply(k, function(x) x[1:(nrow(x)-1),]))))
  ep = do.call(rbind, lapply(qq, function(k) do.call(rbind, lapply(k, function(x) x[2:(nrow(x)),]))))
  colnames(sp) = paste0("start.", cnames)
  colnames(ep) = paste0("end.", cnames)
  df = data.frame(cbind(sp, ep))
  geom_segment(data = df, aes_string(x = paste0("start.", cnames)[1], 
                                     y = paste0("start.", cnames)[2], 
                                     xend = paste0("end.", cnames)[1], 
                                     yend = paste0("end.", cnames)[2]),
               color = color, ...)  
} 

#' Polygon geom for Spatial* objects
#' 
#' @aliases gg.polygon
#' @name gg.polygon
#' @export
#' @param data A SpatialPolygon* object
#' @return geom_polygon
#' 
gg.polygon = function(data, crs = NULL, colour = "black", alpha = 0.1, ...) {
  if ( !is.null(crs) ) { data = spTransform(data, crs) }
  geom_polygon(data= fortify(data), aes(x=long,y=lat,group=group), colour = colour, alpha = alpha, ...)
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

gg.mesh = function(data, crs = NULL, color = rgb(0,0,0,0.1), ...) {
  if ( class(data)[1] == "inla.mesh" ) { mesh = data } else { mesh = data$mesh }
  if ( mesh$manifold == "S2" ) { stop("Geom not implemented for spherical meshes (manifold = S2)" ) }
  if ( !is.null(crs) ) { mesh = inla.spTransform(mesh, CRSobj = crs)}
  
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


plot.prediction = function(..., property = "median") {
  args = list(...)
  pnames = sapply(substitute(list(...))[-1], deparse)
  
  ggopts = attr(args[[1]], "misc")
  data = args[[1]]
  
  if ( attr(data,"type") == "full" ) {
    
    df = do.call(rbind, lapply(1:length(args), function(k) {
      apr = approx(args[[k]]$x, args[[k]]$y, xout = seq(range(args[[k]]$x)[1],range(args[[k]]$x)[2],length.out = 1000))
      data.frame(effect = pnames[[k]], x = apr$x, y = 0.1*apr$y/max(apr$y))}))

    qtl = do.call(rbind, lapply(1:length(args), function(k) 
      data.frame(y = args[[k]]$quantiles,
                 yend = args[[k]]$quantiles,
                 x = k, xend = k + 0.5,
                 effect = pnames[[k]])))
    
    expec = do.call(rbind, lapply(1:length(args), function(k) 
      data.frame(y = args[[k]]$mean,
                 yend = args[[k]]$mean,
                 x = k, xend = k - 0.28,
                 sdy = args[[k]]$sd,
                 cvy = args[[k]]$cv,
                 effect = pnames[[k]])))
    
    sdev = do.call(rbind, lapply(1:length(args), function(k) 
      data.frame(y = args[[k]]$mean + args[[k]]$sd,
                 yend = args[[k]]$mean - args[[k]]$sd,
                 x = k - 0.28, xend = k - 0.28,
                 effect = pnames[[k]])))

    jit = do.call(rbind, lapply(1:length(args), function(k)
      data.frame(y = inla.rmarginal(5000,args[[k]]),
                 n = 500,
                 effect = pnames[[k]])))
    
    plt = ggplot() +  geom_violin(data = df, aes(x=as.numeric(effect),y = x, weight = y, fill=effect), width = 0.5, alpha = 0.7, adjust = 0.2) +
      geom_text(data = qtl, aes(x=xend, y=y, label = signif(y)), size = 3.5, family = "", vjust = -0.5, hjust = 1.1) + 
      geom_text(data = expec, aes(x=xend, y=y, label = paste0(signif(y)," ± ", signif(sdy), " [",round(100*cvy),"%]")), size = 3.5, family = "", vjust = -0.5, angle = 90) + 
      geom_segment(data = qtl, aes(x=x,xend=xend,y=y,yend=yend), linetype = 1, alpha = 0.2) +
      geom_segment(data = expec, aes(x=x,xend=xend,y=y,yend=yend), alpha = 0.5, linetype = 3) +
      geom_segment(data = sdev, aes(x=x,xend=xend,y=y,yend=yend), alpha = 0.5, linetype = 1) +
      geom_point(data = jit, aes(x=as.numeric(effect), y = y), shape = 95, size = 3, alpha = 0.05) +
      ylab(paste0("integral_{", paste(ggopts$idims, collapse = ","), "} (", ggopts$predictor, ")")) + 
      guides(fill=FALSE) + 
      scale_x_continuous(name = "", breaks = 1:length(pnames), labels = pnames) +
      scale_fill_brewer(palette = "YlOrRd", direction = -1) # YlOrRd, PuBu
    
    # This still causes a warning since the violin plot does not like the width argument
    suppressWarnings(print(plt))  
    
  } 
    else if ( attr(data,"type") == "1d" ) {
      
      data = do.call(rbind, lapply(1:length(args), function(k) data.frame(args[[k]], Prediction = pnames[[k]])))
    
      ggp = ggplot(data = data) + 
              geom_ribbon(aes_string(x = ggopts$dims, ymin = "lq", ymax = "uq",  fill = "Prediction"), alpha = 0.3) + 
              geom_line(aes_string(x = ggopts$dims, y = "mean", color = "Prediction"), size = 1.1) +
              xlab(ggopts$predictor) +
              ylab(ggopts$predictor[2])
      
      # Show guide only for multiple graphs
      if (length(pnames) == 1) { ggp = ggp + guides(color=FALSE, fill=FALSE) } 
      
      # Plot graph
      ggp
    }
      
     
    else if ( attr(data,"type") == "spatial" ) {
      # ggplot() + gg.col(ggopts$mesh, color = data$mean) + scale_fill_gradientn(colours = topo.colors(100))
      plot.spatial(mesh = ggopts$mesh, col = data[,property], add.mesh = FALSE)
  }
}


