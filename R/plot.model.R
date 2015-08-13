#' Plot a half normal detection function
#'
#'
#' @aliases plot plot.halfnormal
#' @export
#' @param mdl A \link{model} object
#' @param result An \link{inla} object, i.e. the result of running INLA
#' @param data iDistance data structure. Used to plot a histogram of the detections.
#' @param name Name of the effect used to model the detection function. Default: "nhsd"
#' @param distance.truncation Distance at which to truncate
#' @param covariate Function transforming distances into effect covariates
#' @param add.uncertainty Plot uncertainty boundaries
#' @param add.histogram Use data to plot a histogram of the detections.
#' @param col Color to plot mode of detection function
#' @param ucol Color to plot uncertainty polygon
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

plot.halfnormal = function(mdl = NULL,
                           result = NULL, 
                           data = NULL,  
                           name = "nhsd", 
                           distance.truncation = environment(mdl$formula)$truncation,
                           covariate = mdl$covariates[[name]],
                           add.uncertainty = TRUE,
                           add.histogram = TRUE,
                           col = rgb(250/256, 83/256, 62/256, 1),
                           ucol = rgb(250/256, 83/256, 62/256, 0.3)) {

  hyper.name = paste("Beta for",name)
  
  if (hyper.name %in% rownames(result$summary.hyperpar)) {
    qtl.dist.lower = result$summary.hyperpar[hyper.name,"0.025quant"]
    qtl.dist.upper = result$summary.hyperpar[hyper.name,"0.975quant"]
    de = result$summary.hyperpar[hyper.name,"mean"]
  } 
  else {
    qtl.dist.lower = result$summary.fixed[name,"0.025quant"]
    qtl.dist.upper = result$summary.fixed[name,"0.975quant"]
    de = result$summary.fixed[name,"mode"]
  }
  
  x = data.frame(distance = seq(0,distance.truncation,0.1))
  
  upper = exp(qtl.dist.upper * covariate(x))
  dmean = exp(de*covariate(x))
  lower = exp(qtl.dist.lower * covariate(x))
  
  # If data was provided, plot histogram
  if ( !is.null(data) & add.histogram ) {
    # Compute data histogram, replace values to plot by area normalized to 1
    n.breaks = 10 # breaks=seq(0,mdl$int.args$truncation,length.out=n.breaks)
    hst = hist(detdata(data)$distance,plot=FALSE)
    hst$density = length(hst$density)*hst$density/sum(hst$density,feq=FALSE) # normalized
    plot(hst, 
         freq = FALSE,
         xaxt = 'n', yaxt = 'n',
         ylab = "", xlab = "",
         main = "",
         ylim = c(0,hst$density[1]))
    uy = hst$density[1]
    scale = 1/mean(dmean)
    yaxt = 'n'
    par(new = TRUE)
  } else {
    uy = 1
    scale = 1
    yaxt = NULL
  }
  
  # Plot mode
  plot(x$distance,scale * dmean,, lwd = 3, col = col,
       xlim = c(0,distance.truncation), 
       ylim = c(0,uy),
       type = "l",
       yaxt = yaxt,
       main = "",
       ylab = "Detection probability", 
       xlab = "Distance")
  
  # Plot uncertainty bounds
  if ( add.uncertainty ) {
    polygon(c(x$distance, rev(x$distance)), scale*c(lower, rev(upper)),  col = ucol, border = NA, yaxt = "n")
  }
  
}


#' Plot spatial field
#'
#' 
#' @aliases plot plot.spatial
#' @export
#' @param mdl A \link{model} object
#' @param result An \link{inla} object, i.e. the result of running INLA
#' @param data iDistance data structure. Used to plot a histogram of the detections.
#' @param name Name of the effect used to model the detection function. Default: "spde"
#' @param property Which property of the spatial field to plot, e.g. 'mode' (default)
#' @param stack Plot a combined effect via this stack \link{inla.stack}
#' @param mesh The mesh used for spatial plotting
#' @param mesh.coordinates Coordinates to plot over
#' @param rgl Use rgl to plot
#' @param add.detection Add detections to the plot (needs: data)
#' @param add.points Addtitional points to plot (data.frame). If Boolean and TRUE, plot integration points data$ips (needs: data$ips)
#' @param add Add to previous rgl plot
#' @param colorbar Plot a colorbar
#' @param logscale If set to FALSE, exp() of the predictor is plotted
#' @param col Plot this data instead of the model data
#' 

plot.spatial = function(mdl = NULL,
                     result = NULL, 
                     data = NULL,
                     name = "spde",
                     property = "mode",
                     stack = NULL,
                     mesh = mdl$mesh[[name]],
                     mesh.coords = mdl$mesh.coords[[name]],
                     rgl = FALSE,
                     add.detections = TRUE,
                     add.points = FALSE,
                     add = FALSE, 
                     colorbar = TRUE, 
                     logscale = TRUE, 
                     col = NULL, ...){
  
  geometry = data$geometry
  if ( !is.null(data)) { det.points = detdata(data) } else { add.detections = FALSE }
  if ( add.points == TRUE ) { add.points = data$ips }
  if ( is.data.frame(add.detections) ) { 
    det.points = add.detections
    add.detections = TRUE
  }

  if (rgl) {
    require(rgl)  
    
    # Get predicted values
    if (name %in% names(result$summary.ran)){
      col = result$summary.ran[[name]][[property]]
    } else {
      ind = inla.stack.index(stack, name)$data
      col = result$summary.fitted.values[ind,property]
    }
    
    if (geometry == "geo"){
      
      # Plot sphere
      if (!add){ rgl.open() }
      bg3d(color = "white")
      rgl.earth()
      rgl.sphmesh(mesh, add=TRUE, radius=1.001, col = col)
      
      # Plot detections and integration points
      if ( add.detections ) { rgl.sphpoints(long = det.points$lon+360, lat = det.points$lat, radius = 1.01, col="red", size = 5) }
      if ( add.points ) { rgl.sphpoints(long = add.points$lon+360, lat = add.points$lat, radius=1.01, col = rgb(0,0,1), size = 3) }
      
      # Plot colorbar (This is a temporary workaround, rgl people are working on something like this)
      if (colorbar){
        cb.col = cm.colors(100)
        rgl.linestrips(x=-0.85,y=0.85,z=seq(-1,1,length.out=length(cb.col)),col=cb.col,lwd=50)
        rgl.texts(x=-0.79,y=0.79,z=c(-0.97,0.97),c(format(min(col),digits=3),format(max(col),digits=3)),col="black",cex=1)
      }
      
    }
    else {
      plot(mesh, rgl = rgl, col = col,...)
      if ( add.detections ) { rgl.points(x = det.points$x, y = det.points$y, z =0.02, col = rgb(1,0,0), size=8) }
      if ( int.points ) { rgl.points(x = add.points$x, y = add.points$y, z = 0.01, col = rgb(0,0,0.6), size=4) }
    }
    if (any("par3d.args" %in% names(data))) { do.call(par3d,data$par3d.args) }
  }
  else {
    
    require(lattice)
    proj <- inla.mesh.projector(mesh, dims = c(300,300))
    xlim = range(mesh$loc[,1])
    ylim = range(mesh$loc[,2])
    x = seq(xlim[1], xlim[2], length.out = 300)
    y = seq(ylim[1], ylim[2], length.out = 300)
    grid = expand.grid(x=x,y=y)
    
    if ( !is.null(col) ) {
      A = inla.spde.make.A(mesh, loc = cbind(grid$x,grid$y))
      col = A%*%as.vector(col)
      msk = apply(abs(A), MARGIN = 1, sum) > 0
      col[!msk] = NA
    }
    else if (name %in% names(result$summary.ran)){
      col = inla.mesh.project(proj, field = result$summary.ran[[name]][[property]])
    } 
    else {
      ind = inla.stack.index(stack, name)$data
      col = inla.mesh.project(proj, field = result$summary.fitted.values[ind,property])
    }
    
    if (!logscale) { col = exp(col) }
    
    # print(paste0("min:", min(col,na.rm=TRUE),", max:",max(col,na.rm=TRUE), ", mean:",mean(col,na.rm=TRUE)))
    
    co1 = mesh.coords[1]
    co2 = mesh.coords[2]
    levelplot(col ~ grid$x + grid$y,
              col.regions = topo.colors(100),
              panel=function(...){
                panel.levelplot(...)
                if ( add.points ) { panel.points( x = add.points[,co1], 
                                                  y = add.points[,co2], 
                                                  pch = ".", 
                                                  col = "blue",
                                                  cex = 2) }
                if ( add.detections ) { panel.points( x = det.points[,co1],
                                                  y = det.points[,co2], 
                                                  col = rgb(0,0,0),
                                                  pch = ".",
                                                  cex = 2) }
              },
              xlab = co1, ylab = co2,...
    ) 
  }
}

