#' @name Poisson2_1D
#' @title 1-Dimensional NonHomogeneous Poisson example.
#' @docType data
#' @description Point data and count data, together with intensity function and expected counts for 
#' a unimodal nonhomogeneous 1-dimensional Poisson process example. 
#' 
#' @aliases lambda2_1D cov2_1D E_nc2 pts2 countdata2
#' 
#' @usage data(Poisson2_1D)
#' 
#' @format The data contain the following \code{R} objects:
#'  \describe{
#'    \item{\code{lambda2_1D}:}{ A function defining the intensity function of a 
#'    nonhomogeneous Poisson process. Note that this function is only defined on
#'    the interval (0,55).}
#'    \item{\code{cov2_1D}:}{ A function that gives what we will call a 
#'    'habitat suitability' covariate in 1D space.}
#'    \item{\code{E_nc2}}{ The expected counts of the gridded data.}
#'    \item{\code{pts2}}{ The locations of the observed points (a data frame with one column, named \code{x}).}
#'    \item{\code{countdata2}}{ A data frame with three columns, containing the count data:}
#'    \describe{
#'      \item{\code{x}}{ The grid cell midpoint.}
#'      \item{\code{count}}{ The number of detections in the cell.}
#'      \item{\code{exposure}}{ The width of the cell.}
#'    }
#'  }
#' 
#' @examples
#' \donttest{
#' library(ggplot2)
#' data(Poisson2_1D)
#' p1 = ggplot(countdata2) + 
#' geom_point(data = countdata2, aes(x=x,y=count),col="blue") +ylim(0,max(countdata2$count,E_nc2)) +
#'   geom_point(data = countdata2, aes(x=x), y = 0, shape = "+",col="blue",cex=4) +
#'   geom_point(data=data.frame(x=countdata2$x,y=E_nc2), aes(x=x), y = E_nc2, shape = "_",cex=5) +
#'   xlab(expression(bold(s))) +ylab("count")
#' ss = seq(0,55,length=200)
#' lambda = lambda2_1D(ss)
#' p2 = ggplot() + 
#'   geom_line(data=data.frame(x=ss,y=lambda), aes(x=x,y=y),col="blue") +ylim(0,max(lambda)) +
#'   geom_point(data = pts2, aes(x=x), y = 0.2, shape = "|",cex=4) +
#'   xlab(expression(bold(s))) +ylab(expression(lambda(bold(s))))
#' multiplot(p1,p2,cols=1)
#' }
NULL