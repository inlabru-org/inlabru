#' @name Poisson1_1D
#' @title 1-Dimensional Homogeneous Poisson example.
#' @docType data
#' @description Point data and count data, together with intensity function and expected counts for 
#' a homogeneous 1-dimensional Poisson process example. 
#' 
#' @aliases lambda1_1D E_nc1 pts1 countdata1
#' 
#' @usage data(Poisson1_1D)
#' 
#' @format The data contain the following \code{R} objects:
#'  \describe{
#'    \item{\code{lambda1_1D}:}{ A function defining the intensity function of a 
#'    nonhomogeneous Poisson process. Note that this function is only defined on
#'    the interval (0,55).}
#'    \item{\code{E_nc1}}{ The expected counts of the gridded data.}
#'    \item{\code{pts1}}{ The locations of the observed points (a data frame with one column, named \code{x}).}
#'    \item{\code{countdata1}}{ A data frame with three columns, containing the count data:}
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
#' data(Poisson1_1D)
#' ggplot(countdata1) + 
#' geom_point(data = countdata1, aes(x=x,y=count),col="blue") +ylim(0,max(countdata1$count)) +
#'   geom_point(data = pts1, aes(x=x), y = 0.2, shape = "|",cex=4) +
#'   geom_point(data = countdata1, aes(x=x), y = 0, shape = "+",col="blue",cex=4) +
#'   xlab(expression(bold(s))) +ylab("count")
#' }
NULL