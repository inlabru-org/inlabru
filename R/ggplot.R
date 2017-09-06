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
  ggmap::ggmap(myMap)
}

#' ggplot2 geomes for spatial data
#' 
#' Depending on the class of the first argument one of the following functions is called:
#' \describe{
#' \item{}{\link{gg.data.frame}}
#' \item{}{\link{gg.prediction}}
#' \item{}{\link{gg.SpatialPoints}}
#' \item{}{\link{gg.SpatialLines}}
#' \item{}{\link{gg.SpatialPolygons}}
#' \item{}{\link{gg.SpatialPixels}}
#' \item{}{\link{gg.SpatialGrid}}
#' \item{}{\link{gg.RasterLayer}}
#' \item{}{\link{gg.inla.mesh}}
#' \item{}{\link{gg.inla.mesh.1d}}
#' }
#' @aliases gg
#' @name gg
#' @export
#' @param ... Arguments passed on the the geomes
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


#' Geom for data.frame
#' 
#' @aliases gg.data.frame
#' @name gg.data.frame
#' @export
#' @import ggplot2
#' @param data A data.frame
#' @param ... Arguments passed on to \link{gg.prediction}
#' @return c(geom_ribbon, geom_line)
#' 
gg.data.frame = function(...){ gg.prediction(...) }


#' Geom for predictions
#' 
#' @aliases gg.prediction
#' @name gg.prediction
#' @export
#' @import ggplot2
#' @importFrom utils modifyList
#' @param data A prediction object
#' @param mapping Set of aesthetic mappings created by \link{aes} or \link{aes_}
#' @param ribbon If TRUE, plot a ribbon around the line based on the upper and lower 2.5 percent quantiles
#' @param alpha Alpha level of the ribbon
#' @param bar If TRUE plot boxplot-style summary for each variable
#' @param ... Arguments passed on to \link{geom_line}
#' @return \code{c(geom_ribbon, geom_line)}
#' 
gg.prediction = function(data, mapping = NULL, ribbon = TRUE, alpha = 0.3, bar = FALSE, ...){
  
  if ( all(c("0.025quant", "0.975quant") %in% names(data)) ) {
    names(data)[names(data) == "0.975quant"] = "q0.975"
    names(data)[names(data) == "0.025quant"] = "q0.025"
    names(data)[names(data) == "0.5quant"] = "median"
  }
  lqname = "q0.025"
  uqname = "q0.975"

  
  if ( bar | ( nrow(data) == 1) ) {
    
    sz = 10 # bar width
    med_sz = 4 # median marker size
    
    data = cbind(data.frame(data), 
                 data.frame(variable = rownames(data), 
                            summary = data$mean[1],
                            sdmax = data$mean+data$sd,
                            sdmin = data$mean-data$sd))
    
    geom = c(geom_point(data = data, mapping = aes_string(x = "variable", y = "summary", color = "variable"), shape = 95, size = 0), # Fake ylab
             geom_segment(data = data, mapping = aes_string(y = lqname, yend = uqname, x = "variable", xend = "variable", color = "variable"), size = sz))
    
    # Min and max sample
    if ( all(c("smax", "smax") %in% names(data)) ) {
      geom = c(geom, 
               geom_segment(data = data, mapping = aes_string(y = "smin", yend = "smax", x = "variable", xend = "variable", color = "variable"), size = 1),
               geom_point(data = data, mapping = aes_string(x = "variable", y = "smax", color = "variable"), shape = 95, size = 5),
               geom_point(data = data, mapping = aes_string(x = "variable", y = "smin", color = "variable"), shape = 95, size = 5))
    }

    # Mean and median       
    geom = c(geom,       
             geom_point(data = data, mapping = aes_string(x = "variable", y = "mean"), color = "black", shape = 95, size = sz),
             geom_point(data = data, mapping = aes_string(x = "variable", y = "median"), color = "black", shape = 20, size = med_sz),
             coord_flip())
             
  } else {
    
    if ( "pdf" %in% names(data) ) { 
      y.str = "pdf"
      ribbon = FALSE
    } else if ( "mean" %in% names(data) ){ 
      y.str = "mean" 
    } else if ( "median" %in% names(data) ){
      y.str = "median"
    } else {
      stop("Prediction has neither mean nor median or pdf as column. Don't know what to plot.")
    }
    
    line.map = aes_string(x = names(data)[1], 
                          y = y.str)
    
    ribbon.map = aes_string(x = names(data)[1], 
                            ymin = lqname, 
                            ymax = uqname)
    
    if ( !is.null(mapping) ) { line.map = utils::modifyList(line.map, mapping) }
    
    # Use line color for ribbon filling
    if ( "colour" %in% names(line.map) ) { 
      ribbon.map = modifyList(ribbon.map, aes_string(fill = line.map[["colour"]])) 
      }
    
    geom = geom_line(data = data, line.map, ...)
    if ( ribbon ) {
      geom = c(geom, geom_ribbon(data = data, ribbon.map, alpha = alpha))
    }
  }
  geom
}

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
#' @param color Filling color for the polygons
#' @param alpha Alpha level for polygon filling
#' @param ... Arguments passed on to \link{geom_polygon}
#' @return geom_polygon
#' 
gg.SpatialPolygons = function(data, mapping = NULL, crs = NULL, color = "black", alpha = NULL, ...) {
  
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
#' @param fill Character identifying the data column used for plotting
#' @param ... Arguments passed on to \link{geom_tile}
#' @return A ggplot2 object
#' 

gg.SpatialGrid = function(sgdf, fill = names(sgdf)[[1]], ...) {
  df <- as.data.frame(sgdf)
  gm = geom_tile(data = df, mapping = aes_string(x = coordnames(sgdf)[1], y = coordnames(sgdf)[2], fill = fill), ...)
  
  # If data is not discrete (factor), add default color palette
  if (!inherits(sgdf[[fill]], "factor")) {
    gm = c(gm, scale_fill_gradientn(colours = bru.pal()))
  }
  gm
}

#' Plot SpatialPixelsDataFram using ggplot2
#' 
#' @aliases gg.SpatialPixelsDataFrame
#' @name gg.SpatialPixelsDataFrame
#' @export
#' @import ggplot2
#' @param sgdf A SpatialPixelsDataFrame object
#' @param fill Character array identifying the data column used for plotting
#' @param alpha Character array identifying the data column used for transparency
#' @param mask A SpatialPolygon defining the region that is plotted
#' @param ... Arguments passed on to \link{geom_tile}
#' @return A ggplot2 object
#' 

gg.SpatialPixelsDataFrame = function(sgdf, fill = names(sgdf)[[1]], alpha = NULL, mask = NULL, ...) {
  
  if ( !is.null(mask) ) {
    sgdf = sgdf[as.vector(!is.na(over(sgdf, mask))), ]
  }
  
  df <- as.data.frame(sgdf)
  dmap = aes_string(x = coordnames(sgdf)[1], y = coordnames(sgdf)[2], fill = fill)
  if ( !is.null(alpha) ) dmap = modifyList(dmap, aes_string(alpha = alpha))
  gm = geom_tile(data = df, mapping = dmap, ...)
  
  # If data is not discrete (factor), add default color palette
  if (!inherits(sgdf[[fill]], "factor")) {
    gm = c(gm, scale_fill_gradientn(colours = bru.pal()))
  }
  gm
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

gg.SpatialPixels = function(sgdf, ...) {
  gg(SpatialPoints(sgdf), ...) 
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
#' @param nx Number of pixels in x direction (when plotting using the color parameter)
#' @param ny Number of pixels in y direction (when plotting using the color parameter)
#' @param mask A SpatialPolygon defining the region that is plotted
#' @param ... ignored arguments (S3 generic compatibility)
#' @return geom_polygon or geom_segment
#' 

gg.inla.mesh = function(mesh, 
                        color = NULL, alpha = NULL,
                        edge.color = "grey",
                        interior = TRUE, int.color = "blue", 
                        exterior = TRUE, ext.color = "black",
                        crs = NULL,
                        mask = NULL,
                        nx = 500, ny = 500,
                        ...) {

  
if ( !is.null(color) ) {
  
  px = pixels(mesh, nx = nx, ny = ny)
  A = INLA::inla.spde.make.A(mesh, px)
  px$color = as.vector(A %*% color)
  if ( !is.null(alpha) ) { 
    px$alpha = as.vector(A %*% alpha)
    gg = gg(px, mask = mask, alpha = "alpha")
  } else {
    gg = gg(px, mask = mask)
  }
  
  return(gg)
  
} else {
  if ( mesh$manifold == "S2" ) { stop("Geom not implemented for spherical meshes (manifold = S2)" ) }
  if ( !is.null(crs) ) { mesh = INLA::inla.spTransform(mesh, CRSobj = crs)}
  
  df = rbind(data.frame(a=mesh$loc[mesh$graph$tv[,1],c(1,2)],b=mesh$loc[mesh$graph$tv[,2],c(1,2)]),
             data.frame(a=mesh$loc[mesh$graph$tv[,2],c(1,2)],b=mesh$loc[mesh$graph$tv[,3],c(1,2)]),
             data.frame(a=mesh$loc[mesh$graph$tv[,1],c(1,2)],b=mesh$loc[mesh$graph$tv[,3],c(1,2)]))
  
  colnames(df) = c("x","y","xend","yend")
  mp = aes_string(x = "x",y = "y", xend = "xend", yend = "yend")
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
    geom_point(data = df, mapping = aes_string("x","y"), shape = shape, ...)

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

gg.RasterLayer = function(r, ...) {
  
  requireNamespace("raster")
  spdf <- as(r, "SpatialPixelsDataFrame")
  df <- as.data.frame(spdf)
  # head(r.df)
  # g <- ggplot(r.df, aes(x=x, y=y)) + geom_tile(aes(fill = layer)) + coord_equal()
  geom_tile(data = df, mapping = aes_string(x="x", y="y", fill = "layer"),...)
}

#' Plot bru effects
#' 
#' @method plot bru
#' @export
#' @import ggplot2
#' @param x a \link{bru} object
#' @param ... A character naming the effect to plot
#' @return a gg object

plot.bru = function(x, ...) {
  plotmarginal.inla(x, ...)
}

#' Plot predction using ggplot2
#' 
#' @export
#' @method plot prediction
#' @import ggplot2
#' @param x a prediction object
#' @param y a mapping created by \link{aes} or \link{aes_string}
#' @param ... Arguments passed on to \link{gg.prediction}
#' @return a gg object

plot.prediction = function(x, y = NULL, ...) {
  ggplot() + gg(x, mapping = y, ...)
}

plot.prediction_old = function(..., property = "median") {
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
    
    # Workaround for non-visible bindings
    x = NULL
    y = NULL
    xend = NULL
    yend = NULL
    sdy = NULL
    cvy = NULL
    
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
      pixelplot.mesh(mesh = ggopts$mesh, col = data[,property], property = property, ...)
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
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



# Plot spatial field
#
# 
# @aliases plot plot.spatial
# @name plot.mesh
# @export
# @param mesh The mesh used for spatial plotting. If given, overrides meshes from data and model.
# @param property Which property of the spatial field to plot, e.g. 'mode' (default)
# @param mask A spatial mask restricting the region where we plot. A mask can be one of the following objects:
# (1) A data frame with x and y coordinated defining a polygon. (2) An inla.mesh.segment object. (3) A SpatialPolygonDataFrame (see sp package)
# @param col A numeric vector of values to plot. It's length must be equal to the number of mesh nodes. Overrides all data coming from the inla result.
# @return gg A ggplot2 object 
# 

pixelplot.mesh = function(mesh = NULL,
                          property = "mean",
                          mask = NULL,
                          col = NULL,
                          nx = 300, ...){
  
  data = NULL
  effect = NULL
  group = NULL
  model = NULL
  result = NULL 
  add.mesh = FALSE
  add.detection = FALSE
  add.segment = FALSE
  logscale = NULL
  rgl = FALSE
  add = FALSE  
  
  #
  # (0) Some defaults
  #
  if ( inherits(model, "inla") ) { result = model ; model = NULL}
  if ( is.null(model) && ("model" %in% names(result)) ) { model = result$model }
  if ( is.null(logscale) && !is.null(result$sppa$method) && (result$sppa$method=="lgcp" || result$sppa$method=="poiss")) { logscale = FALSE } else { logscale = TRUE } 
  
  #  
  # (1) The mesh and vertex location
  #
  
  if ( is.null(mesh) ) {
    # If no mesh is provided we first check if data is provided. If so, use the data mesh
    if ( !is.null(data) && ("mesh" %in% names(data)) ) { 
      mesh = data$mesh
      mloc = mesh$loc[,1:length(data$mesh.coords)]
      colnames(mloc) = data$mesh.coords
    }
    else {
      # No mesh and no data set provided. Check if the model has a mesh. If so, use the first mesh we find.
      if ( is.null(model) ) { stop("Your provided no mesh, no data set with a mesh and not model to take a mesh from") }
      if ( "mesh" %in% names(model) & length(model$mesh)>0) {
        mesh = model$mesh[[1]]
        if ("coordnames" %in% names(result$sppa)) {
          mloc = mesh$loc[,1:length(result$sppa$coordnames)]
          colnames(mloc) = result$sppa$coordnames
        } else {
          mloc = mesh$loc[,1:length(model$mesh.coords[[1]])]
          colnames(mloc) = model$mesh.coords[[1]]
        }
      } else {
        stop("Your provided no mesh, no data set with a mesh and your model does not have a mesh.")
      }
    }
  } else {
    # User provided mesh. Where to take the coordinate names from?
    if (all(mesh$loc[,2]==0)) { mdim = 1 } 
    else if(all(mesh$loc[,3]==0)) { mdim = 2 }
    else { mdim = 3 }
    mloc = mesh$loc[,1:mdim]
    colnames(mloc) = c("x","y","z")[1:mdim]
  }
  coords = colnames(mloc)
  
  #
  # (2) Projection locations
  #
  
  if ( rgl ) { ploc = mloc }
  else {
    xlim = range(mloc[,1])
    ylim = range(mloc[,2])
    aspect = diff(ylim)/diff(xlim)
    ny = round(nx * aspect)
    proj <- inla.mesh.projector(mesh, dims = c(nx,ny))
    x = seq(xlim[1], xlim[2], length.out = nx)
    y = seq(ylim[1], ylim[2], length.out = ny)
    grid = expand.grid(x=x,y=y)
    ploc = data.frame(cbind(grid$x,grid$y))
    colnames(ploc) = colnames(mloc)
  }
  if ( !is.null(group) ) { loc = merge(ploc, group) }
  
  #
  # (3) Values at the mesh location 
  #
  if ( logscale ) { link = function(x){x} } else { link = exp }
  
  if ( !is.null(col) ) {
    if ( is.data.frame(col) ) {
      col = do.call(c,lapply(property, 
                             function(prp) { 
                               col = inla.mesh.project(proj, field = col[, prp])
                               return(as.vector(col))}))
      col = data.frame(col = link(col), property = merge(rep(1,nrow(ploc)), property)[,2])
      ploc =  ploc[rep(seq_len(nrow(ploc)), length(property)), ]
      df = data.frame(col, ploc)
    } else {
      col = as.vector((inla.mesh.project(proj, field = col)))
      df = data.frame(col = link(col), ploc)
    }
  }
  else if ( !is.null(effect) ) {
    # We extract values from a tagged INLA result OR directly from a random effect
    if ( is.null(result) ) { stop("You are trying to plot without providing an INLA result.") } 
    if ( effect %in% names(result$summary.random)) {
      # Values come from a random effect
      ires = result$summary.random[[effect]] 
    }
    else {
      # Values come from a tagged effect, e.g. a prediction
      ires = result$summary.fitted.values[inla.stack.index(result$stack, effect)$data, ]
    }
    col = do.call(c,lapply(property, 
                           function(prp) { 
                             col = inla.mesh.project(proj, field = ires[, prp])
                             return(as.vector(col))}))
    col = data.frame(col = link(col), property = merge(rep(1,nrow(ploc)), property)[,2])
    ploc =  ploc[rep(seq_len(nrow(ploc)), length(property)), ]
    df = data.frame(col, ploc)
  }
  else if ( !is.null(model) ) {
    stop("Not supported anymore.")
  }
  else { stop("Not sure what to plot.") }
  
  #
  # (4) Select a region to plot and create alpha channel
  #
  df$col[!is.inside(mesh, ploc, colnames(ploc))] = NA
  df = cbind(df, alpha = is.inside(mesh, ploc, colnames(ploc)))
  
  if ( class(mask) == "inla.mesh.segment" ) {
    mask.mesh = inla.mesh.create(loc = mesh$loc, boundary = list(mask))
    df$col[!is.inside(mask.mesh, ploc, colnames(ploc))] = NA
    df$alpha = df$alpha & is.inside(mask.mesh, ploc, colnames(ploc))
  } else if ( class(mask) == "SpatialPolygonsDataFrame" ) {
    msk = is.na(sp::over(SpatialPoints(ploc, proj4string = CRS(proj4string(mask))), mask , fn = NULL)[,1])
    df$col[msk] = NA
    # df$alpha[msk] = 0 ; df$alpha[!msk] = 1
  } else if ( is.data.frame(mask) ){
    msk = point.in.polygon(ploc[,1], ploc[,2], mask[,1], mask[,2]) == 1
    df$col[!msk] = NA
    df$alpha = df$alpha & msk
  }
  
  #
  # (5) Create the plot using either rgl or ggplot2
  #
  

    gg = ggplot(df, aes_string(x = coords[1], y = coords[2]) )
    gg = gg + geom_raster(aes_string(fill = "col", alpha = "alpha"), interpolate = TRUE)
    gg = gg + scale_alpha_discrete(guide = 'none')
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      gg = gg + scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9,"YlOrRd"))
    } 
    gg = gg + theme(legend.title = element_blank()) + coord_fixed()
    
    # Facet properties
    if (length(property)>1) { gg = gg + facet_grid(. ~ property) }
    
    # If the data is grouped, use ggplot facets
    if ( !is.null(group) ) { gg = gg + facet_grid(as.formula(paste0(". ~ ", colnames(group)[1]))) }
    
    # Plot the mesh
    if ( add.mesh ) { gg = gg + gg.inla.mesh(mesh) }
    
    return(gg)
  
}

# Default color palette for plots using \link{gg}
# 
# @aliases bru.pal
# @name bru.pal
# @export
# @param 
# @return Color values as returned by \code{RColorBrewer::brewer.pal(9, "YlOrRd")}

bru.pal = function() {
  # library(RColorBrewer)
  # brewer.pal(9, "YlOrRd")
  pal = c("#FFFFCC","#FFEDA0","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#BD0026","#800026")
}

