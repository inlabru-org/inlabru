#' Detection function fitting using INLA.
#' 
#' This method provides an interface to the class constructors 
#' \code{\link{dfun_halfnormal}},\code{\link{dfun_logconcave}} and
#' \code{\link{dfun_spde1d}}. 
#'
#' @aliases dfun
#' @export 
#' @param family A string determining the detection function family. 
#' Choices are: "\code{halfnormal}", "\code{logconcave}" and "\code{spde1d}".
#' @param data A data set containing a list of sightings and their distances,
#' for example the \code{\link{whales}} data set.
#' @return \code{dfun}, detection function object.
#' @examples \\dontrun{data(whales) ; fun = dfun("halfnormal",whales); plot(fun)}
#' @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>
#'

dfun = function(family,data,...){
  if (family=="halfnormal") {return(dfun_halfnormal(data,...))}
  else if (family=="logconcave") {return(dfun_logconcave(data,...))}
  else if (family=="spde1d") {return(dfun_spde1d(data,...))}
}


summary.dfun = function(dfun){
  cat("INLA variables:\n")
  vars = variables.inla(dfun$INLA$result)
  tbl = format(vars,scientific=FALSE,digits=2)[,c("type","mean","sd","0.025quant","0.5quant","0.975quant","mode")]
  colnames(tbl) = c("type","mean","sd","0.025q","0.5q","0.975q","mode")
  print(tbl)
  cat("\n")
  cat("Other statistics: \n")
  cat(paste0("Truncation: ", dfun$int.args$truncation),"\n")
  mrate = exp(vars["(Intercept)","mean"])
  cat(paste0("Rate mean at perfect sight: ", mrate),"\n")
  trlen = sum(len(as.transect(dfun$data$effort),dfun$data))
  cat(paste0("Total length of transects: ", trlen),"\n")
  cat(paste0("Rate mean per transect length unit: ", mrate/trlen)," -log-> ", log(mrate/trlen))
}

#' Evaluate detection function at given distances
#' 
#' @aliases value
#' @export
#' @return A vector with detection function evaluations.
#' 
value <- function(x) UseMethod("value")
value.dfun_halfnormal <- function(...) dfun_halfnormal.value(...)
value.dfun_logconcave <- function(...) dfun_logconcave.value(...)
value.dfun_spde1d <- function(...) dfun_spde1d.value(...)

#' Sample detection function at given distances
#'
sample.value = function(...){UseMethod("sample.value")}
sample.logvalue = function(...){UseMethod("sample.logvalue")}



#' Check if an object is a detection function
#' 
#' @aliases is.dfun
#' @export
#' @return A vector with detection function evaluations.
#' 

is.dfun = function(obj){
  if (class(obj)[[1]] == "dfun") { return(TRUE) }
  else { return(FALSE) }
}

