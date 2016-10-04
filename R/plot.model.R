#' Plot a marginal posterior given an INLA result
#' 
#'
#' @aliases plot.marginal
#' @param effects Character array stating the name of the effect to be plotted
#' @export

plot.marginal = function(result, effects, ...){ 
  plot.marginal.inla(result, varname = effects, ...) 
}


#' Plot a half-normal detection function (DEPRECATED)
#' 
#' DEPRECATED. Use \link{plot.detfun}.
#'
#' @aliases plot plot.halfnormal

plot.halfnormal = function(...) { plot.detfun(...) }


#' Plot a detection function
#'
#' Plots a detection function from a given \link{model} and \link{INLA} estimate. 
#' 
#' This function is fairly general but in most cases it is applied to a model of the \link{model.detfun} family. The model is
#' evaluated for distances from 0 to a value picked by the \code{distance.truncation} parameter. In order
#' to give an idea about the uncertainty of the estimates, the upper an lower 2.5 percent quantiles are 
#' plotted as well. Note, however that their estimated are only plausible for some models and a better way
#' of caclulating the has yet to be implemented. Any of the plotted graphs can be replace by manually 
#' setting the \code{upper}, \code{lower} and \code{mode} as vectors of values. Note that these must have
#' length of 500 or otherwise the \code{loc} parameter has to be set to a data frame with a column
#' called 'distance' corresponding to the distances the graphs are representing.
#' 
#' For convenience multiple plotting options are available. If a \link{dsdata} object is provided, the
#' detection function is plotted against a histogram of the detection distances. In case the detection
#' function is only mean to model a subset of the detection, the \code{filter} parameter can be used
#' to select them. 
#' 
#' Plotting detection functions that are conditional on covariated is slightly more involved.
#' For this purpose please provide the \code{prior} parameter (function) and set \code{integration} to TRUE. 
#' For each of the location in \code{loc} the prior as well as the detection function are evaluated.
#' Thereafter all locations are integrated over except for the distance dimension.
#'
#' @aliases plot plot.detfun
#' @export
#' @param model A \link{model} object
#' @param result An \link{inla} object, i.e. the result of running INLA
#' @param distance.truncation Distance at which to truncate
#' @param data iDistance data structure. Used to plot a histogram of the detections.
#' @param loc Distance and possibly covariates at which to plot at (optional).
#' @param prior A function representiong a prior probability for each location in the \code{loc} parameter
#' @param do.ecdf If TRUE, plot an ECDF instead of the detection function
#' @param upper Vector of upper quantiles to plot
#' @param mode Vector of modes to plot (basically the detection function)
#' @param lower  Vector of lower quantiles to plot
#' @param integrate If TRUE, integrate detection function over all od loc's dimensions except for distance
#' @param covariate Function transforming distances into effect covariates
#' @param add.uncertainty If TRUE, plot uncertainty boundaries
#' @param add.histogram If TRUE, use \code{data} to plot a histogram of the detections.
#' @param filter A function used to filter the detections before calculating the histogram
#' @param col Color to plot mode of detection function
#' @param ucol Color to plot uncertainty polygon
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

plot.detfun = function(model = NULL,
                           result = NULL,
                           distance.truncation = NULL,
                           data = NULL,  
                           loc = NULL,
                           prior = NULL,
                           do.ecdf = FALSE,
                           upper = NULL,
                           mode = NULL,
                           lower = NULL,
                           integrate = FALSE,
                           covariate = model$covariates[[name]],
                           add = FALSE,
                           add.uncertainty = TRUE,
                           add.histogram = TRUE,
                           breaks = "sturges",
                           filter = NULL,
                           col = rgb(250/256, 83/256, 62/256, 1),
                           ucol = rgb(250/256, 83/256, 62/256, 0.3), 
                           ggp = TRUE,...) {
  
  if ( is.null(distance.truncation) ) {
    if ( is.null(data) & is.null (loc) ) {
      stop("Your have to provide at least one of the following parameters: distance.truncation, data or loc")
    } else {
      if (! is.null(loc) ) (distance.truncation = max(loc$distance))
      if (! is.null(data) ) (distance.truncation = max(detdata(data)$distance))
    }
  }

  if (is.null(loc)) { 
    x = data.frame(distance = seq(0, distance.truncation, length.out = 2000)) 

  }
  else { 
    x = loc 

  }
  
  if ( is.null(upper) ) { upper = evaluate.model(model = model, inla.result = result, loc = x, property = "0.975quant", link = exp) }
  if ( is.null(mode) ) { 
    dmean = evaluate.model(model = model, inla.result = result, loc = x, property = "mode", link = exp)
  } else {dmean = mode}
  if ( is.null(lower) ) { lower = evaluate.model(model = model, inla.result = result, loc = x, property = "0.025quant", link = exp) }
  
  if ( !is.null(prior) ) {
    upper = upper * prior(loc)
    dmean = dmean * prior(loc)
    lower = lower * prior(loc)
  }
  
  if ( integrate ) {
    upper = as.vector(by(upper * loc[,"weight"], loc$distance, sum))
    dmean = as.vector(by(dmean * loc[,"weight"], loc$distance, sum))
    lower = as.vector(by(lower * loc[,"weight"], loc$distance, sum))
    x = data.frame(distance = seq(0, distance.truncation, length.out = length(dmean))) 
  }
  
  # If data was provided, plot histogram
  if ( !is.null(data) & add.histogram) {
    
    # Filter
    dets = detdata(data)
    if ( !is.null(filter) ) { dets = filter(dets) } 
    
    if ( do.ecdf ) {
      edf = ecdf(dets$distance)
      x.plot = seq(0,distance.truncation,length.out=100)
      if ( ! ggp ) {
        if ( add ) { lines(x,1-edf(x$distance)) }
        else {
          plot(x.plot, 1-edf(x.plot), type = "l",
               xaxt = 'n', yaxt = 'n',
               ylab = "", xlab = "",
               main = "",)
        }
      }
      uy = 1
    } else {
      # Compute data histogram, replace values to plot by area normalized to 1
      if (is.numeric(breaks) & length(breaks)==1) { breaks = seq(0, distance.truncation,length.out = breaks) }
      hst = hist(dets$distance,plot=FALSE, breaks = breaks)
      hst$density = hst$density/mean(hst$density) # normalized 
      uy = max(hst$density[1],max(dmean/mean(dmean)))
      if (! ggp ) {
        if ( add ) {
          
        } else {
          plot(hst, 
               freq = FALSE,
               xaxt = 'n', yaxt = 'n',
               ylab = "", xlab = "",
               main = "",
               ylim = c(0,uy),
               xlim = c(0,distance.truncation))
        }
      }
    }
    scale = 1/mean(dmean)
    yaxt = 'n'
    par(new = TRUE)
  } else {
    uy = 1
    scale = 1/max(dmean)
    yaxt = NULL
  }
  
  # Are we plotting the ECDF?
  if ( do.ecdf ) { 
    y.plot = 1-cumsum(dmean)/sum(dmean)
    y.boundary = c(1-cumsum(lower)/sum(lower), rev(1-cumsum(upper)/sum(upper)) ) #1-cumsum(upper)/sum(upper)
  }
  else { 
    y.plot = scale * dmean
    y.boundary = scale*c(lower, rev(upper))
  }
  
  # Plot mode
  if ( ggp ) {
    if (!require(ggplot2)) {stop("ggplot2 not installed. Set ggp = FALSE to use default plotting method")}
    if ( do.ecdf ) {
      
      df = data.frame(distance=x$distance, 
                      mode = 1-cumsum(dmean)/sum(dmean), 
                      upper = 1-cumsum(upper)/sum(upper), 
                      lower = 1-cumsum(lower)/sum(lower),
                      edf = 1-edf(x$distance))
      
      gg = ggplot(data=df, aes(x=distance, y=mode))
      if ( !is.null(data) & add.histogram ) {
        gg = gg + geom_path(aes(y = edf),colour = "black" , alpha = 0.5)
      }
      gg = gg + geom_ribbon(aes(x=distance, ymax=upper, ymin=lower), fill = ucol, alpha = 0.2) + 
        geom_path(colour = col, size = 1) +
        ylab("detection probability ecdf")
      gg
      return(gg)
      
    } else {
      df = data.frame(distance=x$distance, mode = dmean, upper = upper, lower = lower)
      gg = ggplot(data=df, aes(x=distance, y=mode))
      
      if ( add.histogram & !do.ecdf ) {
        hstdf = data.frame(x=hst$mids, y = hst$density/scale)
        gg = gg + geom_bar(data = hstdf, stat="identity", aes(x=x, y=y), fill = rgb(0,0,0,0.2)) +
          coord_cartesian(ylim = c(0,max(hstdf$y,1)))
      }
      gg = gg + geom_ribbon(aes(x=distance, ymax=upper, ymin=lower), fill = ucol, alpha = 0.2) + 
        geom_path(colour = col, size = 1) +
        ylab("detection probability")
    
      return(gg)
    }
    

  } else {
    if (add) {
      lines(x$distance,y.plot, lwd = 3, col = col)
    } else {
      plot(x$distance,y.plot, lwd = 3, col = col,
           xlim = c(0,distance.truncation), 
           ylim = c(0,uy),
           type = "l",
           yaxt = yaxt,
           main = "",
           ylab = "1 - ECDF", 
           xlab = "Distance")
    }  
    # Plot uncertainty bounds
    if ( add.uncertainty ) {
      polygon(c(x$distance, rev(x$distance)), y.boundary,  col = ucol, border = NA, yaxt = "n")
    }
  }
  
}

#' Plot spatial field
#'
#' Plots spatial \link{model}s given an \link{inla} result object.
#' 
#' The main purpose of this function is to visualize spatial model posteriors using the parameters
#' \code{model}, \code{result} and \code{data}. The model is then evaluated on the data set's mesh
#' and using the \code{property} setting. The latter determins which property of the latent effects
#' are used of each effect to be combined. 
#' Options are: 'mode', 'mean', 'sd', '0.975quant', '0.025quant', 'kld'. Default is 'mean'.
#' 
#' The other paramters of this function mainly serve a purpose if one of the three main paramters
#' is not set. The \code{result} object is not needed if the \code{col} parameter is provided as
#' the latter already defines the spatial field to be plotted. Moreover, if no \code{model} is
#' provided, the \code{name} parameter is used to select an effect from the INLA \code{result}. 
#' 
#' @aliases plot plot.spatial
#' @name plot.spatial
#' @export
#' @param model A \link{model} object
#' @param result An \link{inla} object, i.e. the result of running INLA
#' @param data A \link{dsdata} object
#' @param effect Needed if no model or col is provided. This selects the random field to plot, e.g. "spde". For predictions this corresponds to the tag parameter.
#' @param property Which property of the spatial field to plot, e.g. 'mode' (default)
#' @param group A data frame that is merged with the spatial locations at which we plot the data. For instance, in temporal SPDE models, we can set group = data.frame(group=1:10).
#' @param mesh The mesh used for spatial plotting. If given, overrides meshes from data and model.
#' @param mask A spatial mask restricting the region where we plot. A mask can be one of the following objects:
#' (1) A data frame with x and y coordinated defining a polygon. (2) An inla.mesh.segment object. (3) A SpatialPolygonDataFrame (see sp package)
#' @param col A numeric vector of values to plot. It's length must be equal to the number of mesh nodes. Overrides all data coming from the inla result.
#' @param add.mesh If TRUE, plot the mesh
#' @param add.detection If TRUE, plot the detections (requires: data)
#' @param add.segment If TRUE, plot the segments (requires: data)
#' @param logscale If set to FALSE, exp() of the values is plotted
#' @param rgl Use rgl to plot
#' @param add Add to previous rgl plot 
#' @return gg A ggplot2 object (if rgl=FALSE)
#' 
plot.spatial = function(model = NULL,
                        result = NULL, 
                        data = NULL,
                        effect = NULL,
                        property = "mean",
                        group = NULL,
                        mesh = NULL,
                        mask = NULL,
                        col = NULL,
                        add.mesh = TRUE,
                        add.detection = FALSE,
                        add.segment = FALSE,
                        logscale = NULL,
                        rgl = FALSE,
                        add = FALSE,
                        nx = 300){
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
    # We extract values from a model and INLA result
    if ( is.null(result) ) { stop("You are trying to plot a model without providing an INLA result.") }
    if ( !is.null(group) ){ ploc = merge(ploc, group) }
    col = do.call(c,lapply(property, 
                         function(prp) { 
                           col = evaluate.model(model, inla.result = result, loc = ploc, do.sum = TRUE, property = prp) 
                           return(as.vector(col))}))
    col = data.frame(col = link(col), property = merge(rep(1,nrow(ploc)), property)[,2])
    ploc =  ploc[rep(seq_len(nrow(ploc)), length(property)), ]
    df = data.frame(col, ploc)
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
      msk = is.na(sp::over(SpatialPoints(ploc), mask , fn = NULL)[,1])
      df$col[msk] = NA
      df$alpha = df$alpha & msk
    } else if ( is.data.frame(mask) ){
      if ( !require(sp) ) { "You provided a data.frame as mask. The package sp is required to interpret it as a polygon."}
      msk = point.in.polygon(ploc[,1], ploc[,2], mask[,1], mask[,2]) == 1
      df$col[!msk] = NA
      df$alpha = df$alpha & msk
    }
    
#
# (5) Create the plot using either rgl or ggplot2
#
    
if ( !rgl ){ 
  # Use ggplot2
  if ( !requireNamespace("ggplot2", quietly = TRUE) ) { stop("This function requires the ggplot2 package")}
  
  gg = ggplot(df, aes_string(x = coords[1], y = coords[2]) )
  gg = gg + geom_raster(aes_string(fill = "col", alpha = "alpha"), interpolate = TRUE)
  gg = gg + scale_alpha_discrete(guide = 'none')
  gg = gg + scale_fill_gradientn(colours = topo.colors(100) ) 
  gg = gg + theme(legend.title = element_blank()) + coord_fixed()

  # Facet properties
  if (length(property)>1) { gg = gg + facet_grid(. ~ property) }
  
  # If the data is grouped, use ggplot facets
  if ( !is.null(group) ) { gg = gg + facet_grid(as.formula(paste0(". ~ ", colnames(group)[1]))) }
  
  # Plot the mesh
  if ( add.mesh ) { gg = gg + gg.mesh(mesh) }
  
  return(gg)
  
} else {
  # Use rgl to plot
  col = df$col
  if ( !requireNamespace("rgl", quietly = TRUE) ) { stop("This function requires the rgl package. Install rgl or set rgl=FALSE.")}
  
  if (data$geometry == "geo"){
    
    # Plot sphere
    if (!add){ rgl.open() }
    bg3d(color = "white")
    rgl.earth()
    
    # # Plot detections and integration points
    # if ( add.detections ) { rgl.sphpoints(long = det.points$lon+360, lat = det.points$lat, radius = 1.01*R, col="red", size = 5) }
    # if ( add.points ) { rgl.sphpoints(long = add.points$lon+360, lat = add.points$lat, radius=1.01*R, col = rgb(0,0,1), size = 3) }
    # 
    # Plot colorbar (This is a temporary workaround, rgl people are working on something like this)
    if ( TRUE ){
      cb.col = cm.colors(100)
      rgl.linestrips(x=-0.85,y=0.85,z=seq(-1,1,length.out=length(cb.col)),col=cb.col,lwd=50)
      rgl.texts(x=-0.79,y=0.79,z=c(-0.97,0.97),c(format(min(col),digits=3),format(max(col),digits=3)),col="black",cex=1)
    }
    rgl.pop()
    rgl.sphmesh(mesh, add=TRUE, radius=1.001, col = col)
    
  }
  else {
    plot(mesh, rgl = rgl, col = col,...)
    if ( add.detections ) { rgl.points(x = det.points$x, y = det.points$y, z =0.02, col = rgb(1,0,0), size=8) }
    if ( int.points ) { rgl.points(x = add.points$x, y = add.points$y, z = 0.01, col = rgb(0,0,0.6), size=4) }
  }
}
  
}

#' Play (animate) spatial field
#'
#' Animates a spatial field using RGL. 
#' 
#' @aliases play.spatial
#' @export
#' @param group Example: group = list(year = c(1,2)) animates the field for years 1 and 2
#' @param ... Parameters passed on to \link{plot.spatial}
#' 


play.spatial = function(group = list(), rgl, ...){
  if ( !rgl ) { rgl = TRUE}
  rgl.earth()
  sargs = list(...)
  myanim = function(time, ...) {
    par3d(skipRedraw = TRUE)
    grp = list()
    grp[[names(group)[[1]]]] = group[[1]][mod(floor(time),2)+1]
    do.call(plot.spatial, c(sargs, list(group = grp, add = TRUE, rgl = TRUE))) ;
    par3d(skipRedraw = FALSE)
    return("")
  }
  
  play3d(myanim, duration = 10, startTime = 0, fps = 1)
  
}