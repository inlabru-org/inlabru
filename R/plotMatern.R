#' @title Plots Matern correlation or covariance function.
#'
#' @description
#' Plots Matern correlation or covariance function when passed the output of a call to 
#' \code{inla.spde.result}. 
#'
#' @param spde.pars Output of a call to \code{inla.spde.result}.
#' @param corr If TRUE, plots correlation function, else plots covariance function.
#' @param ylim Limits of y-axis for plotting.
#' @param rangemult Number by which to multiply the mean posterior range parameter in 
#' order to calculate the maximum distance for which the Matern covariance or correlation 
#' function is to be plotted.
#' @param d The dimension of the domain in which the SPDE lives.
#'
#' @return A \code{SpatialPolygonsDataFrame}.
#'  
#' @export
plotMatern = function(spde.pars,corr=TRUE,ylim=NULL,rangemult=1,d=2) {
  
  materncov = function(dist, kappa,d,corr){
    inla.matern.cov(nu=1, kappa, x=dist,d=d, corr = corr)
  }
  
  kappaQ = inla.qmarginal(p=c(0.025,0.5,0.975), marginal=spde.pars$marginals.kappa[[1]])
  
  xmax = exp(spde.pars$summary.log.range.nominal$mean)*rangemult
  x = seq(0,xmax,length=200)
  # median
  materncov.kappaQ2 = apply(matrix(data=x, ncol=1), 1, 
                            function(X)materncov(dist=X, kappa =kappaQ[2],d=d,corr=corr ))
  # lower band
  materncov.kappaQ1 = apply(matrix(data=x, ncol=1), 1, 
                            function(X)materncov(dist=X, kappa =kappaQ[1],d=d,corr=corr))
  # upper band
  materncov.kappaQ3 = apply(matrix(data=x, ncol=1), 1,
                            function(X)materncov(dist=X, kappa =kappaQ[3],d=d,corr=corr))
  if(corr) ylab = "Matern Correlation"
  if(!corr) ylab = "Matern Covariance"
  if(is.null(ylim)) ylim = c(0,max(materncov.kappaQ1))
  plot(x,materncov.kappaQ2,type="n",ylab=ylab,ylim=ylim)
  polygon(c(rev(x),x),c(rev(materncov.kappaQ3), materncov.kappaQ1), col="gray", 
          border = NA)
  lines(x,materncov.kappaQ2,lwd=2)
  invisible(data.frame(x=x,median=materncov.kappaQ1,lower=materncov.kappaQ2,upper=materncov.kappaQ3))
}