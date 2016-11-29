#' Bin checking
#'
#' @aliases bincount
#' @export
#' @param result A result object from a \link{bru} or \link{lgcp} call
#' @param predictor A formula describing what and where to predict. See \link{predict.lgcp} for details.
#' @param observations A vector of observed values
#' @param breaks A vector of bin boundaries
#' @param nint Number of integration points per bin
#' @param ... arguments passed on to \link{predict}
#' @return An \link{inla} object


bincount = function(result, predictor, observations, breaks, nint = 20, ...) {
  # Quantiles to comoute
  qtls = c(0.025, 0.5, 0.975)
  
  # Filter out observations outside bins
  observations = observations[(observations>=min(breaks)) & (observations<=max(breaks))]
  
  # Number of ...
  nobs = length(observations) # observations
  nbins = length(breaks) - 1 # bins
  
  # Mid points of bins
  mid = breaks[1:(length(breaks)-1)] + diff(breaks)/2
  
  # Create integration points
  ip = int.quadrature(breaks[1:(length(breaks)-1)], breaks[2:(length(breaks))], scheme = "trapezoid", n = nint)
  points = data.frame(tmp = as.vector(ip$ips))
  colnames(points) = all.vars(update.formula(predictor, ~.0))
  points$bin = rep(1:nbins,each = nint)
  
  # Predict
  prd = predict(r, predictor = predictor, points = points, ...)
  
  # Integrate per bin
  smp = attr(prd,"samples")
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
    pint[[k]] = approx(zz, xx, xout = qtls, rule = 2)$y
  }
  pint = data.frame(do.call(rbind, pint))
  colnames(pint) = c("lq","mq","uq")
  
  # Append more information
  pint = cbind(data.frame(
    mid = mid, 
    width = diff(breaks),
    counts = hist(observations, breaks = breaks, plot = FALSE)$counts),
    pint)
  pint$inside = ( pint$counts >= pint$lq ) & ( pint$counts <= pint$uq ) 
  
  ggp = ggplot(bc) + geom_crossbar(aes(x=mid,y=mq,ymin=lq,ymax=uq,fill=inside,color=inside), show.legend=FALSE) +
    geom_point(aes(x=mid,y=mq), shape = 95, size = 3, color = "blue") +
    geom_point(aes(x=mid,y=counts), shape = 20, size = 1) +
    xlab(all.vars(update.formula(predictor, ~.0))) +
    ylab("count")
  
  attr(pint, "ggp") = ggp
  
  # Return
  pint
}
