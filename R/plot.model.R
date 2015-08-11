#' Plot a half normal detection function
#'
#'
#' @aliases plot plot.halfnormal
#' @export
#' @param mdl A \link{model}
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
  
  x = data.frame(distance = seq(0,truncation,0.1))
  
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
       xlim = c(0,truncation), 
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