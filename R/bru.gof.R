#' 1D LGCP bin count simulation and comparison with data
#' 
#' A common procedure of analyzing the distribution of 1D points is to chose a binning
#' and plot the data's histogram with respect to this binning. This function compares the
#' counts that the histogram calculates to simulations from a 1D log Gaussian Cox process
#' conditioned on the number of data samples. For each bin this results in a median number 
#' of counts as well as a confidence interval. If the LGCP is a plausible model for the 
#' observed points then most of the histrogram counts (number of points within a bin) 
#' should be within the confidence intervals. Note that a proper comparison  is a multiple 
#' testing problem which the function does not solve for you.
#' 
#'
#' @aliases bincount
#' @export
#' @importFrom graphics hist
#' @param result A result object from a \link{bru} or \link{lgcp} call
#' @param predictor A formula describing the prediction of a 1D LGCP via \link{predict}.
#' @param observations A vector of observed values
#' @param breaks A vector of bin boundaries
#' @param nint Number of integration points per bin. Increase this if the bins are wide and
#' @param probs numeric vector of probabilities with values in [0,1]
#' @param ... arguments passed on to \link{predict}
#' @return An \link[INLA]{inla} object
#' 
#' @examples 
#' 
#' \dontrun{
#' 
#' # Load a point pattern
#' data(Poisson2_1D)
#' 
#' # Take a look at the point (and frequency) data
#' 
#' ggplot(pts2) + 
#'   geom_histogram(aes(x = x), binwidth = 55/20, boundary = 0, fill = NA, color = "black") +
#'   geom_point(aes(x), y = 0, pch = "|", cex = 4) + 
#'   coord_fixed(ratio = 1)
#' 
#' #' Fit an LGCP model 
#' x <- seq(0, 55, length = 50)
#' mesh1D <- inla.mesh.1d(x, boundary = "free")
#' mdl <- x ~ spde1D(map = x, model = inla.spde2.matern(mesh1D)) + Intercept # SOLUTION
#' fit.spde <- lgcp(mdl, pts2, domain = list(x = c(0,55)))
#' 
#' # Calculate bin statistics
#' bc <- bincount(result = fit.spde, 
#'                observations = pts2, 
#'                breaks = seq(0,max(pts2),length = 12), 
#'                predictor = x ~ exp(spde1D + Intercept))
#' 
#' 
#' # Plot them!
#' attributes(bc)$ggp
#' 
#' }


bincount = function(result, predictor, observations, breaks, nint = 20, probs = c(0.025, 0.5, 0.975), ...) {

  # Sort probabilities
  probs = sort(probs)
  
  # Filter out observations outside bins
  observations = observations[(observations>=min(breaks)) & (observations<=max(breaks))]
  
  # Number of ...
  nobs = length(observations) # observations
  nbins = length(breaks) - 1 # bins
  
  # Mid points of bins
  mid = breaks[1:(length(breaks)-1)] + diff(breaks)/2
  
  # Create integration points
  ip = int.quadrature(breaks[1:(length(breaks)-1)], breaks[2:(length(breaks))], scheme = "trapezoid", n.points = nint)
  points = data.frame(tmp = as.vector(ip$ips))
  colnames(points) = all.vars(update.formula(predictor, ~.0))
  points$bin = rep(1:nbins,each = nint)
  
  # Sampler
  smp = generate(result, points, predictor, ...)
  smp = do.call(cbind, smp)
  
  # Integrate per bin
  smp = ip$wl * smp
  qq = aggregate(smp, by = list(rep(1:nbins,each = nint)), FUN = sum, drop = TRUE)[,2:(ncol(smp)+1)]

  # Normalize bin probabilities
  for (s in 1:ncol(smp)) { qq[,s] = qq[,s]/sum(qq[,s]) }
  
  # For each bin calculate predictive interval
  pint = list()
  for (k in 1:nbins) {
    xx = 0:nobs
    cdff = function (p) pbinom(xx, size = nobs, prob = p)
    zz = apply(qq[k,,drop=FALSE], MARGIN = 2, cdff)
    zz = apply(zz, MARGIN = 1, mean)
    pint[[k]] = approx(zz, xx, xout = probs, rule = 2)$y
  }
  pint = data.frame(do.call(rbind, pint))
  colnames(pint) = paste0("q",probs)
  mxname = paste0("q", max(probs))
  miname = paste0("q", min(probs))
  mdname = paste0("q", probs[1+floor(length(probs)/2)])
  
  # Append more information
  pint = cbind(data.frame(
    mid = mid, 
    width = diff(breaks),
    counts = hist(observations, breaks = breaks, plot = FALSE)$counts),
    pint)
  pint$inside = ( pint$counts >= pint[[miname]] ) & ( pint$counts <= pint[[mxname]] ) 
  
  ggp = ggplot() + geom_crossbar(data = pint, mapping = aes_string(x="mid",y=mdname,ymin=miname,ymax=mxname,fill="inside",color="inside"), show.legend=FALSE) +
    geom_point(data = pint, mapping = aes_string(x="mid",y=mdname), shape = 95, size = 3, color = "blue") +
    geom_point(data = pint, mapping = aes_string(x="mid",y="counts"), shape = 20, size = 2) +
    xlab(all.vars(update.formula(predictor, ~.0))) +
    ylab("count")
  
  attr(pint, "ggp") = ggp
  
  # Return
  pint
}



#' Variance and correlations measures for prediction components
#' 
#' Calculates local and integrated variance and correlation measures as introduced by Yuan et al. (2017).
#'
#' @aliases devel.cvmeasure
#' @export
#' @param joint A joint \code{prediction} of two latent model components.
#' @param prediction1 A \code{prediction} of first component.
#' @param prediction2 A \code{prediction} of the first component.
#' @param samplers A SpatialPolygon object describing the area for which to compute the cummulative variance measure.
#' @param mesh The \code{inla.mesh} for which the prediction was performed (required for cummulative Vmeasure).
#' @return Variance and correlations measures.
#' 
#' @references 
#' Y. Yuan, F. E. Bachl, F. Lindgren, D. L. Brochers, J. B. Illian, S. T. Buckland, H. Rue, T. Gerrodette. 2017. 
#' Point process models for spatio-temporal distance sampling data from a large-scale survey of blue whales.
#' \url{https://arxiv.org/abs/1604.06013}
#' 
#' @example inst/examples/devel.cvmeasure.R

devel.cvmeasure = function(joint, prediction1, prediction2, samplers = NULL, mesh = NULL) {

  # Covariance
  joint$cov = (joint$var - prediction1$var - prediction2$var)/2
  
  # Correlation
  corr = function(joint, a, b) { ((joint - a - b) / (2 * sqrt(a * b)))}
  cor = corr(joint$var, prediction1$var, prediction2$var)
  if (any((cor>1) | (cor<(-1)))) {
    cor[(cor>1) | (cor<(-1))] = NA
  }
  joint$cor = cor
  
  if ( !is.null(samplers) ){
    
    #
    # PRESENTED TO YOU BY HACKY McHACKERSON
    #
 
    
    wips = ipoints(samplers)
    A = INLA::inla.spde.make.A(mesh, loc = wips)
    
    weights = wips$weight
    weights = weights/sum(weights) 
    vj = sum(as.vector(A %*% joint$var) * weights)
    v1 = sum(as.vector(A %*% prediction1$var) * weights)
    v2 = sum(as.vector(A %*% prediction2$var) * weights)
    cr = corr(vj, v1, v2)
    ret = data.frame(var.joint = vj, var1 = v1, var2 = v2, cor = cr)
    
  } else {
    tmp = attributes(joint)
    ret = joint[,c("cov","cor")]
    ret$var.joint = joint$var
    ret$var1 = prediction1$var
    ret$var2 = prediction2$var
    attr(ret, "misc") = tmp[["misc"]]
    attr(ret, "type") = tmp[["type"]]
  }
  
  return(ret)
}


