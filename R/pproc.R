#' Point processes for wildlife abundance estimation using INLA.
#'
#' This method provides an interface to the class constructors
#' \code{\link{pproc_lgcp}}
#'
#' @aliases pproc
#' @export
#' @param family A string determining the point process model.
#' Choices are: "\code{lgcp}",
#' @param data A data set containing sightings, transects and a mesh of the survey area,
#' for example the \code{\link{whales}} data set.
#' @return \code{pproc} point process model.
#' @examples \\dontrun{data(whales) ; pp = pproc("lgcp",whales); plot(pp)}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

pproc= function(family,data,...){
  if (family=="lgcp") {return(pproc_lgcp(data,...))}
}

lgcp = function(...){ return(pproc_lgcp(...)) }

# GENERICS
simulate = function(...){UseMethod("simulate")}
logintensity = function(...){UseMethod("logintensity")}

summary.pproc = function(pproc, include.random=FALSE){
  vars = variables.inla(pproc$INLA$result, include.random=include.random)
  tbl = format(vars,scientific=FALSE,digits=2)[,c("type","mean","sd","0.025quant","0.5quant","0.975quant","mode")]
  colnames(tbl) = c("type","mean","sd","0.025q","0.5q","0.975q","mode")
  print(tbl)
  cat("\n")
  cat("Exponential transformation of quantiles: \n")
  print(exp(vars[,c("0.025quant","0.5quant","0.975quant")]))
  cat("\n")
  cat("Other statistics: \n")
  cat(paste0("Truncation K = ", pproc$int.args$truncation),"\n")
  cat(paste0("Number of integration points = ", dim(pproc$int.points)[1]),"\n")
  cat(paste0("Total weight of integration points W = ", sum(pproc$int.points$weight)),"\n")
  cat(paste0("W/K = ", sum(pproc$int.points$weight)/pproc$int.args$truncation),"\n")
  #mrate = exp(vars["Intercept","mean"])
  #cat(paste0("Rate mean at perfect sight: ", mrate),"\n")
  #trlen = sum(len(as.transect(pproc$data$effort),pproc$data))
  #cat(paste0("Total length of transects: ", trlen),"\n")
  #cat(paste0("Rate mean per transect length unit: ", mrate*trlen)," -log-> ", log(mrate*trlen))
}

blah = function(loc,data){
  d = rep(Inf,nrow(loc))
  for (k in 1:nrow(data$sighting)){
    dst = sqrt((data$sighting$lat[k]-loc$lat)^2+(data$sighting$lon[k]-loc$lon)^2)
    d[dst<d] = dst[dst<d]
  }
  return(d)
}
