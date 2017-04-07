#' Plot a map using extend of a spatial object
#' 
#' @aliases gmap
#' @name gmap
#' @export
#' @param data A Spatial* object
#' @param ... Arguments passed on to get_map
#' @return A ggmap
#' 
gmap = function(data, ...) {
  data = spTransform(data, CRS("+proj=longlat"))
  df = cbind(coordinates(data), data@data)
  
  # Figure out a sensible bounding box (range of data plus 30%)
  extend = 2.5
  crange = apply(coordinates(data), MARGIN = 2, range)
  lonlim = (extend*(crange[,1] - mean(crange[,1]))) + mean(crange[,1])
  latlim = (extend*(crange[,2] - mean(crange[,2]))) + mean(crange[,2])
  
  # Create map
  requireNamespace("ggmap")
  myMap = ggmap::get_map(c(lonlim[1], latlim[1], lonlim[2], latlim[2]), ...)
  
  # Return map
  ggmap(myMap)
}

#' ggplot2 geomes for spatial data
#' 
#' 
#' @aliases gg
#' @name gg
#' @export
#' 

gg = function(...){UseMethod("gg")}


#' ggmap geom for spatial objects
#' 
#' Coordinates are transformed to lon/lat for convenient plotting using ggmap
#' 
#' @aliases gm
#' @name gm
#' @export
#' @param data A Spatial* object or a mesh
#' @param ... Arguments passed on to \link{gg}
#' 
gm = function(data, ...) { gg(data, crs = CRS("+proj=longlat"), ...) }



#' Point geom for SpatialPoints objects
#' 
#' @aliases gg.SpatialPoints
#' @name gg.SpatialPoints
#' @export
#' @import ggplot2
#' @param data A SpatialPoints object
#' @param mapping Set of aesthetic mappings created by \link{aes} or \link{aes_}
#' @param crs A \link{CRS} object defining the coordinate system to project the data to before plotting
#' @param ... Arguments passed on to \link{geom_point}
#' @return geom_point
#' 
gg.SpatialPoints = function(data, mapping = NULL, crs = NULL, ...) {
  
  if ( !is.null(crs) ) { data = spTransform(data, crs) }
  
  df = data.frame(data)
  
  dmap = aes_string(x = coordnames(data)[1], 
                    y = coordnames(data)[2])
  
  # dmap = modifyList(dmap, aes(color = "#08519c")) # , alpha = 0.6))
  
  if ( !is.null(mapping) ) {dmap = modifyList(dmap, mapping)}
  
  geom_point(data = df, mapping = dmap, ...)
}

#' Segment geom for SpatialLines objects
#' 
#' @aliases gg.SpatialLines
#' @name gg.SpatialLines
#' @export
#' @import ggplot2
#' @param data A SpatialLines object
#' @param mapping Set of aesthetic mappings created by \link{aes} or \link{aes_}
#' @param crs A \link{CRS} object defining the coordinate system to project the data to before plotting
#' @param ... Arguments passed on to \link{geom_segment}
#' @return geom_segment
#' 
gg.SpatialLines = function(data, mapping = NULL, crs = NULL, ...) {
  
  if ( !is.null(crs) ) { data = spTransform(data, crs) }
  
  qq = coordinates(data)
  cnames = coordnames(data)
  if (is.null(cnames)) { cnames = c("x","y") }
  sp = do.call(rbind, lapply(qq, function(k) do.call(rbind, lapply(k, function(x) x[1:(nrow(x)-1),]))))
  ep = do.call(rbind, lapply(qq, function(k) do.call(rbind, lapply(k, function(x) x[2:(nrow(x)),]))))
  colnames(sp) = paste0("start.", cnames)
  colnames(ep) = paste0("end.", cnames)
  df = data.frame(cbind(sp, ep), data@data)
  
  dmap = aes_string(x = paste0("start.", cnames)[1], 
                    y = paste0("start.", cnames)[2], 
                    xend = paste0("end.", cnames)[1], 
                    yend = paste0("end.", cnames)[2])
  
  if ( !is.null(mapping) ) {dmap = modifyList(dmap, mapping)}
  
  geom_segment(data = df, mapping = dmap, ...)  
} 

#' Polygon geom for SpatialPolygons objects
#' 
#' @aliases gg.SpatialPolygons
#' @name gg.SpatialPolygons
#' @export
#' @import ggplot2
#' @param data A SpatialPolygons object
#' @param mapping Set of aesthetic mappings created by \link{aes} or \link{aes_}
#' @param crs A \link{CRS} object defining the coordinate system to project the data to before plotting
#' @param ... Arguments passed on to \link{geom_polygon}
#' @return geom_polygon
#' 
gg.SpatialPolygons = function(data, mapping = NULL, crs = NULL, alpha = NULL, color = color, ...) {
  
  if ( !is.null(crs) ) { data = spTransform(data, crs) }
  df = fortify(data)
  dmap = aes_string(x = "long", y = "lat", group = "group")
  if ( !is.null(mapping) ) { dmap = modifyList(dmap, mapping) }
  if (!("alpha" %in% names(dmap)) & is.null(alpha)) { alpha = 0.1 }
  if (!("color" %in% names(dmap)) & is.null(color)) { color = "black" }

  geom_polygon(data = df, mapping = dmap, alpha = alpha, color = color, ...)
}

#' Plot SpatialGrid using ggplot2
#' 
#' @aliases gg.SpatialGrid
#' @name gg.SpatialGrid
#' @export
#' @import ggplot2
#' @param sgdf A SpatialGrid object
#' @param ... Arguments passed on to \link{geom_tile}
#' @return A ggplot2 object
#' 

gg.SpatialGrid = function(sgdf, fill = names(sgdf)[[1]], ...) {
  df <- as.data.frame(sgdf)
  geom_tile(data = df, mapping = aes_string(x = coordnames(sgdf)[1], y = coordnames(sgdf)[2], fill = fill), ...)
}

#' Plot SpatialPixels using ggplot2
#' 
#' @aliases gg.SpatialPixels
#' @name gg.SpatialPixels
#' @export
#' @import ggplot2
#' @param sgdf A SpatialPixels object
#' @param ... Arguments passed on to \link{geom_tile}
#' @return A ggplot2 object
#' 

gg.SpatialPixels = function(sgdf, fill = names(sgdf)[[1]], ...) {
  df <- as.data.frame(sgdf)
  geom_tile(data = df, mapping = aes_string(x = coordnames(sgdf)[1], y = coordnames(sgdf)[2], fill = fill))
}



#' Plot inla.mesh using ggplot2
#' 
#' @aliases gg.inla.mesh
#' @name gg.inla.mesh
#' @export
#' @import ggplot2
#' @param mesh An \link{inla.mesh} object
#' @param color A vector of scalar values to fill the mesh with colors. The length of the vector mus correspond to the number of mesh vertices.
#' @param alpha A vector of scalar values setting the alpha value of the colors provided
#' @param edge.color Color of the mesh edges
#' @param interior If TRUE, plot the interior boundaries of the mesh
#' @param int.color Color used to plot the interior boundaries
#' @param exterior If TRUE, plot the exterior boundaries of the mesh
#' @param ext.color Color used to plot the interior boundaries
#' @param crs A \link{CRS} object defining the coordinate system to project the mesh to before plotting
#' @return geom_polygon or geom_segment
#' 

gg.inla.mesh = function(mesh, 
                        color = NULL, alpha = 1,
                        edge.color = "grey",
                        interior = TRUE, int.color = "blue", 
                        exterior = TRUE, ext.color = "black",
                        crs = NULL) {

  
if ( !is.null(color) ) {
  if ( length(alpha) == 1 ) { 
    alpha = rep(alpha,length(color))
    alph.fix = TRUE}
  else { alph.fix = FALSE }
  
  refine = TRUE
  if ( refine ) {
    omesh = mesh
    mesh = mesh.refine(mesh, refine = list(max.edge = diff(range(omesh$loc[,1]))/100))
    color = as.vector(INLA::inla.spde.make.A(omesh, loc = mesh$loc) %*% color)
    if ( !is.null(alpha) ) { alpha = as.vector(INLA::inla.spde.make.A(omesh, loc = mesh$loc) %*% alpha) }
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
  
} else {
  if ( mesh$manifold == "S2" ) { stop("Geom not implemented for spherical meshes (manifold = S2)" ) }
  if ( !is.null(crs) ) { mesh = inla.spTransform(mesh, CRSobj = crs)}
  
  df = rbind(data.frame(a=mesh$loc[mesh$graph$tv[,1],c(1,2)],b=mesh$loc[mesh$graph$tv[,2],c(1,2)]),
             data.frame(a=mesh$loc[mesh$graph$tv[,2],c(1,2)],b=mesh$loc[mesh$graph$tv[,3],c(1,2)]),
             data.frame(a=mesh$loc[mesh$graph$tv[,1],c(1,2)],b=mesh$loc[mesh$graph$tv[,3],c(1,2)]))
  
  colnames(df) = c("x","y","xend","yend")
  mp = aes(x = x,y = y, xend = xend, yend = yend)
  
  msh = geom_segment(data = df, mapping = mp, color = edge.color)
  
  # Outer boundary
  if ( exterior ) {
    df = data.frame(mesh$loc[mesh$segm$bnd$idx[,1],1:2], mesh$loc[mesh$segm$bnd$idx[,2],1:2])
    colnames(df) = c("x","y","xend","yend")
    bnd = geom_segment(data = df, mapping = mp, color = ext.color)
  } else { bnd = NULL }
  
  if ( interior ) {
    # Interior boundary
    df = data.frame(mesh$loc[mesh$segm$int$idx[,1],1:2], mesh$loc[mesh$segm$int$idx[,2],1:2])
    colnames(df) = c("x","y","xend","yend")
    if ( nrow(df) == 0 ) { int = NULL } else {
      int = geom_segment(data = df, mapping = mp, color = int.color)
    }
  } else { int = NULL }
  
  # Return combined geomes
  c(msh, bnd, int)
} 
  
  
}


#' Plot inla.mesh.1d using ggplot2
#' 
#' @aliases gg.inla.mesh.1d
#' @name gg.inla.mesh.1d
#' @export
#' @import ggplot2
#' @param mesh An inla.mesh.1d object
#' @param y Single numeric defining the y-coordinates of the mesh knots to plot
#' @param shape Shape of the knot markers
#' @param ... parameters passed on to \link{geom_point}
#' @return geom_point
#' 

gg.inla.mesh.1d = function(mesh, y = 0, shape = 4, ...) {
  
    df = data.frame(x = mesh$loc, y = y)
    geom_point(data = df, mapping = aes(x,y), shape = shape, ...)

}



#' Plot RasterLayer using ggplot2
#' 
#' @aliases gg.RasterLayer
#' @name gg.RasterLayer
#' @export
#' @import ggplot2
#' @param r A RasterLayer object
#' @param ... Arguments passed on to \link{geom_tile}
#' @return geom_tile
#' 

gg.RasterLayer = function(r) {
  
  library(raster)
  spdf <- as(r, "SpatialPixelsDataFrame")
  df <- as.data.frame(spdf)
  # head(r.df)
  # g <- ggplot(r.df, aes(x=x, y=y)) + geom_tile(aes(fill = layer)) + coord_equal()
  geom_tile(data = df, mapping = aes(x=x, y=y, fill = layer))
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
    
    txt_size = (5/14) * theme_get()$text$size
    
    # Function for formating numbers
    fmt = function(x) {
      th = 1E3
      if (abs(x)<th & abs(x)>0.01) {
        # sprintf("%.2f",x)
        as.character(signif(x,4))
      } else {
        sprintf("%.2E",x)
      }
    }
    
    plt = ggplot() +  geom_violin(data = df, aes(x=as.numeric(effect),y = x, weight = y, fill=effect), width = 0.5, alpha = 0.7, adjust = 0.2) +
      geom_text(data = qtl, aes(x=xend, y=y, label = fmt(y)), size = txt_size, family = "", vjust = -0.5, hjust = 1.1) + 
      geom_text(data = expec, aes(x=xend, y=y, label = paste0(fmt(y)," +/- ", fmt(sdy), " [",round(100*cvy),"%]")), size = txt_size, family = "", vjust = -0.5, angle = 90) + 
      geom_segment(data = qtl, aes(x=x,xend=xend,y=y,yend=yend), linetype = 1, alpha = 0.2) +
      geom_segment(data = expec, aes(x=x,xend=xend,y=y,yend=yend), alpha = 0.5, linetype = 3) +
      geom_segment(data = sdev, aes(x=x,xend=xend,y=y,yend=yend), alpha = 0.5, linetype = 1) +
      geom_point(data = jit, aes(x=as.numeric(effect), y = y), shape = 95, size = 3, alpha = 0.05) +
      ylab(paste0("integral_{", paste(ggopts$idims, collapse = ","), "} (", ggopts$predictor, ")")) + 
      guides(fill=FALSE) + 
      scale_x_continuous(name = "", breaks = 1:length(pnames), labels = pnames) +
      scale_fill_brewer(palette = "YlOrRd", direction = -1) # YlOrRd, PuBu
    
    # This still causes a warning since the violin plot does not like the width argument
    # suppressWarnings(print(plt))  
    return(plt)
    
  } 
  
    else if ( attr(data,"type") == "0d" ) {
      dnames = colnames(data)
      ggp = ggplot()
      
      if ( "summary" %in% names(attributes(data)) ){
        smy = attr(data, "summary")
        ggp = ggp + geom_ribbon(data = smy$inner.marg, ymin = 0, aes_string(x = dnames[1], ymax = dnames[2]), alpha = 0.1, colour = "grey")
      } 
      
      ggp = ggp + geom_path(data = data, aes_string(x = dnames[1], y = dnames[2]))
      
      if (length(pnames) == 1) { ggp = ggp + guides(color=FALSE, fill=FALSE) } 

      ggp
    }

    else if ( attr(data,"type") == "1d" ) {
      
      data = do.call(rbind, lapply(1:length(args), function(k) data.frame(args[[k]], Prediction = pnames[[k]])))
    
      ggp = ggplot(data = data) + 
              geom_ribbon(aes_string(x = ggopts$dims, ymin = "lq", ymax = "uq",  fill = "Prediction"), alpha = 0.3) + 
              geom_line(aes_string(x = ggopts$dims, y = "median", color = "Prediction"), size = 1.1) +
              xlab(ggopts$predictor) +
              ylab(ggopts$predictor[2])
      
      # Show guide only for multiple graphs
      if (length(pnames) == 1) { ggp = ggp + guides(color=FALSE, fill=FALSE) } 
      
      # Plot graph
      ggp
    }
      
     
    else if ( attr(data,"type") == "spatial" ) {
      # ggplot() + gg.col(ggopts$mesh, color = data$mean) + scale_fill_gradientn(colours = topo.colors(100))
      plot.spatial(mesh = ggopts$mesh, col = data[,property], add.mesh = FALSE, property = property, ...)
  }
}


#' @title Multiple ggplots on a page.
#'
#' @description
#' Renders multiple ggplots on a single page. 
#'
#' @param ... Comma-separated \code{ggplot} objects.
#' @param plotlist A list of \code{ggplot} objects - an alternative to the comma-separated argument above.
#' @param cols Number of columns of plots on the page.
#' @param layout A matrix specifying the layout. If present, 'cols' is ignored. If the layout is 
#' something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE), then plot 1 will go in the upper left, 
#' 2 will go in the upper right, and 3 will go all the way across the bottom.
#'
#' @source 
#' \url{http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/}
#'  
#' @examples 
#' df = data.frame(x=1:10,y=1:10,z=11:20)
#' pl1 = ggplot(data = df) + geom_line(mapping = aes(x,y), color = "red")
#' pl2 = ggplot(data = df) + geom_line(mapping = aes(x,z), color = "blue")
#' multiplot(pl1,pl2, cols = 2)
#' 
#' @export
#
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
