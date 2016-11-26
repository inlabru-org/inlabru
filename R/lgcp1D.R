#' @name lgcp1D
#' @title 1-Dimensional LGCP example utility functions.
#' @docType data
#' @description Utility functions for the illustrative 1D LGCP example. 
#' 
#' @usage data(lgcp1D)
#' 
#' @format The data contain two \code{R} functions:
#'  \describe{
#'    \item{\code{lambda_1D}:}{ A function defining the intensity function of a 
#'    nonhomogeneous Poisson process. Note that this function is only defined on
#'    the interval (1,55).}
#'    \item{\code{cov_1D}:}{ A function that gives what we will call a 
#'    'habitat suitability' covariate in 1D space.}
#'  }
#' 
#' @examples
#' data(lgcp1D)
#' x = seq(1,55,by=2)
#' plot(x,lambda_1D(x),type="l",ylab=expression(lambda(x)))
#' lines(x,cov_1D(x),lty=2)
#' legend(legend=c(expression(lambda(s)),"habitat suitability"),lty=1:2,x="bottom",bty="n")
#'  
NULL