#' Plot a generic detection function
#'
#'
#' @aliases plot plot.dfun
#' @export
#' @param mdl A \link{dfun} model
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

plot.dfun = function(dfun,...){
  # temporarily we divert to half normal case
  plot.dfun_halfnormal(dfun,...)
}

#' Plot a half normal detection function
#'
#'
#' @aliases plot plot.dfun_halfnormal
#' @export
#' @param mdl A \link{dfun} model
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>


plot.dfun_halfnormal <- function(mdl,truth=FALSE,...){

  # colors
  
  r=250/256
  g=83/256
  b=62/256

  if ("Beta for distance" %in% rownames(mdl$INLA$result$summary.hyperpar)) {
    qtl.dist.lower = mdl$INLA$result$summary.hyperpar["Beta for distance","0.025quant"]
    qtl.dist.upper = mdl$INLA$result$summary.hyperpar["Beta for distance","0.975quant"]
    de = mdl$INLA$result$summary.hyperpar["Beta for distance","mean"]
  } 
  else {
    qtl.dist.lower = mdl$INLA$result$summary.fixed["distance","0.025quant"]
    qtl.dist.upper = mdl$INLA$result$summary.fixed["distance","0.975quant"]
    de = mdl$INLA$result$summary.fixed["distance","mode"]
  }


  x = seq(0,mdl$int.args$truncation,0.1)
  
  upper = exp(qtl.dist.upper*mdl$dist.trafo(x))
  dmean = exp(de*mdl$dist.trafo(x))
  lower = exp(qtl.dist.lower*mdl$dist.trafo(x))

  # Compute data histogram, replace values to plot by area normalized to 1
  n.breaks = 10 # breaks=seq(0,mdl$int.args$truncation,length.out=n.breaks)
  hst = hist(detdata(mdl$data)$distance,plot=FALSE)
  hst$density = length(hst$density)*hst$density/sum(hst$density,feq=FALSE) # normalized
  plot(hst,freq=FALSE,xaxt='n', yaxt='n', ylab="", xlab="",main="",ylim=c(0,hst$density[1]))
  
  # Plot mean
  par(new=TRUE)
  scale = 1/mean(dmean)
  plot(x,scale * dmean,, lwd=3, col=rgb(r,g,b,1),
       xlim=c(0,mdl$int.args$truncation), 
       ylim=c(0,hst$density[1]),type="l",
       yaxt='n',main="",
       ylab="Detection probability",xlab="Distance")
  
  # Plot uncertainty bounds
  axis(2,at=c(0,scale),labels=c(0,1))
  polygon(c(x, rev(x)), scale*c(lower, rev(upper)),  col=rgb(r,g,b,0.3), border=NA)
  
  # title(main="Half-normal detection function mean with 2.5% quantiles")
}



#' Plot empirical CDF for half-normal detection function
#'
#'
#' @aliases plot plot.ecdf.dfun_logconcave
#' @export
#' @param mdl A \link{dfun_logconcave} model
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'
plot.ecdf.dfun_halfnormal = function(mdl){
  x = seq(0,mdl$int.args$truncation,0.01)
  
  if ("Beta for distance" %in% rownames(mdl$INLA$result$summary.hyperpar)) {
    de = mdl$INLA$result$summary.hyperpar["Beta for distance","mean"]
  } 
  else {
    de = mdl$INLA$result$summary.fixed["distance","mode"]
  }
  
  
  mean = exp(de*mdl$dist.trafo(x))
  
  # data ecdf
  xecdf = ecdf(detdata(mdl$data)$distance)
  plot(x,xecdf(x),col=rgb(0.3,0.3,0.3,1),ylim=c(1,0),xlim=c(0,mdl$int.args$truncation),type='l',lwd=3,xaxt='n', yaxt='n', ylab="", xlab="",main="")
  
  # dfun ecdf
  par(new=TRUE)
  plot(seq(0,mdl$int.args$truncation,length.out=length(mean)),cumsum(mean)/sum(mean),col=rgb(1,0,0,1),
       ylim=c(1,0),xlim=c(0,mdl$int.args$truncation),
       type='l',lwd=3,
       xlab = "Distance",
       ylab = "ecdf")
  
}


    
#' Plot a log concave detection function
#'
#'
#' @aliases plot plot.dfun_logconcave
#' @export
#' @param mdl A \link{dfun_logconcave} model
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#' 
plot.dfun_logconcave <- function(mdl,truth=FALSE,...){
  
  r=250/256
  g=83/256
  b=62/256
  
  x = seq(0,mdl$int.args$truncation,0.01)
  dmean = dfun_logconcave.value(mdl,data=data.frame(distance=x),field="mean")
  
  # compute quantiles
  smp = dfun_logconcave.sample.value(mdl,x,n=100)
  smp.matrix = do.call(cbind,smp)
  lower = apply(smp.matrix,1,quantile,probs=c(0.025))
  upper = apply(smp.matrix,1,quantile,probs=c(0.975))
  
  # Compute data histogram, replace values to plot by area normalized to 1
  n.breaks = 10 # breaks=seq(0,mdl$int.args$truncation,length.out=n.breaks)
  hst = hist(detdata(mdl$data)$distance,plot=FALSE)
  hst$density = length(hst$density)*hst$density/sum(hst$density,feq=FALSE) # normalized
  plot(hst,freq=FALSE,xaxt='n', yaxt='n', ylab="", xlab="",main="",ylim=c(0,hst$density[1]))
  
  # Plot mean
  par(new=TRUE)
  scale = 1/mean(dmean)
  plot(x,scale * dmean,, lwd=3, col=rgb(r,g,b,1),
       xlim=c(0,mdl$int.args$truncation), 
       ylim=c(0,hst$density[1]),type="l",
       yaxt='n',main="",
       ylab="Detection probability",xlab="Distance")
  
  # Plot uncertainty bounds
  axis(2,at=c(0,scale),labels=c(0,1))
  polygon(c(x, rev(x)), scale*c(lower, rev(upper)),  col=rgb(r,g,b,0.3), border=NA)
  
}


#' Plot empirical CDF for log concave detection function
#'
#'
#' @aliases plot plot.ecdf.dfun_logconcave
#' @export
#' @param mdl A \link{dfun_logconcave} model
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'
plot.ecdf.dfun_logconcave = function(mdl){
  x = seq(0,mdl$truncation,0.01)
  mean = dfun_logconcave.value(mdl,data=data.frame(distance=x),field="mean")
  
  # data ecdf
  xecdf = ecdf(detdata(mdl$data)$distance)
  plot(x,xecdf(x),col=rgb(0.3,0.3,0.3,1),ylim=c(1,0),xlim=c(0,mdl$truncation),type='l',lwd=3,xaxt='n', yaxt='n', ylab="", xlab="",main="")
  
  # dfun ecdf
  par(new=TRUE)
  plot(seq(0,mdl$truncation,length.out=length(mean)),cumsum(mean)/sum(mean),col=rgb(1,0,0,1),
       ylim=c(1,0),xlim=c(0,mdl$truncation),
       type='l',lwd=3,
       xlab = "Distance",
       ylab = "ecdf")
  
}


#' Plot a 1D SPDE detection function
#'
#'
#' @aliases plot plot.dfun_spde1d
#' @export
#' @param mdl A \link{dfun_spde1d} model
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'  
plot.dfun_spde1d <- function(mdl,truth=FALSE,...){
  warning("This function is is very unfinished...")
  lines(mdl$INLA$result$summary.fitted.values$mean[index])
}
