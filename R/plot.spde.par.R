#' @title Plot SPDE range or variance parameter posterior.
#'
#' @description
#' Plots the posterior distribution of the range, log(range), variance, or log(variance) 
#' parameter of a model's SPDE component.
#'
#' @param result output from a call to \code{inla.spde.result}.
#' @param varname One of "range", "log.range", "variance" or "log.variance".
#' @param add If TRUE, adds to existing plot.
#' @param ggp If TRUE, plots with ggplot2.
#' @param lwd Line width.
#' @param ... Other parameters passed to \code{plot()}, \code{points()} or \code{lines()}.
#'
#' @return Nothing.
#' 
#' @examples 
#' library(inlabru)
#' data(Poisson3_1D)
#' cd = countdata3b
#' #cd = Poisson3_1D$countdata3b
#' x = seq(0, 55, length.out = 50)
#' mesh1D <- inla.mesh.1d(x, boundary="free")
#' ggplot() + gg(mesh1D) + xlim(0,55)
#' mdl = count ~ spde1D(map=x, model=inla.spde2.matern(mesh1D), mesh=mesh1D) + Intercept
#' fit3b.bru = bru(cd,model=mdl,family="poisson",E=cd$exposure)
#' spdespec<-inla.spde2.matern(mesh=mesh1D)
#' spde.pars <- inla.spde.result(fit3b.bru, "spde1D", spdespec)
#' plot.spde.par(spde.pars)
#' plot.spde.par(spde.pars,"variance")
#' plot.spde.par(spde.pars,"log.range")
#' plot.spde.par(spde.pars,"log.variance")
#'  
#' @export
plot.spde.par = function(result,varname="range", add = FALSE, ggp = TRUE, lwd=3,...) {
  marg = switch(varname,
                range = result$marginals.range.nominal[[1]],
                log.range = result$marginals.log.range.nominal[[1]],
                variance = result$marginals.variance.nominal[[1]],
                log.variance = result$marginals.log.variance.nominal[[1]]
  )
  if(is.null(marg)) stop("Invalid varname: ",varname,". must be one of 'range', 
                         'log.range',  'variance',  'log.variance'")
  uq = inla.qmarginal(0.975, marg)
  uqy = inla.dmarginal(uq, marg)
  lq = inla.qmarginal(0.025, marg)
  lqy = inla.dmarginal(lq, marg)
  inner.x = seq(lq, uq, length.out = 100)
  inner.marg = data.frame(x = inner.x, y = inla.dmarginal(inner.x, marg))
  if ( ggp ) {
    require( ggplot2 )
    df = data.frame(marg)
    ggplot(data = df, aes(x=x,y=y)) + geom_path() + geom_ribbon(ymin = 0,aes(ymax = y), alpha = 0.1) +
      geom_segment(x = lq, y = 0, xend = lq, yend = lqy) +
      geom_segment(x = uq, y = 0, xend = uq, yend = uqy) +
      geom_ribbon(data = inner.marg, ymin = 0, aes(ymax = y), alpha = 0.1) +
      xlab(varname) + ylab("pdf")
  } else {
    if (!add) {
      plot(marg,type='l',xlab=varname,ylab="Posterior density",...)
    } else {
      points(marg, type='l', xlab="",ylab="", ...)
    }
    lheight = max(marg[,"y"])
    lines(x=c(vars[varname,"mode"],vars[varname,"mode"]),y=c(0,lheight),col="blue",lwd=lwd,...)
    lines(x=c(vars[varname,"mean"],vars[varname,"mean"]),y=c(0,lheight),col="red",lwd=lwd,...)
    lines(x=c(vars[varname,"0.025quant"],vars[varname,"0.025quant"]),y=c(0,lheight),col=rgb(0,0.6,0),lwd=lwd,...)
    lines(x=c(vars[varname,"0.975quant"],vars[varname,"0.975quant"]),y=c(0,lheight),col=rgb(0,0.6,0),lwd=lwd,...)
  }
} 