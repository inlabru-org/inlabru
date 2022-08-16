#' @name robins_subset
#' @title Robins
#' @docType data
#' @description This is the `robins_subset` dataset which is a subset of the
#'   full Robins data set from ...
#'
#' @usage data(robins_subset)
#'
#' @format The data are a data.frame with variables
#'  \describe{
#'    \item{`variablename`:}{ The description}
#'  }
#' @source
#' github ...
#'
#'
#' @references
#' Citation?
#'
#' @examples
#' if (require(ggplot2, quietly = TRUE)) {
#'   data(robins_subset, package = "inlabru") # get the data
#'
#'   # plot the counts for one year of data
#'   ggplot(robins_subset[robins_subset$std_yr == 0, ]) +
#'     geom_point(aes(lon, lat, colour = count + 1)) +
#'     scale_colour_gradient(low = "blue", high = "red", trans = "log")
#' }
NULL
