#' @name toypoints
#' @title Simulated 2D point process data
#' @docType data
#' @description
#'
#' This data set serves as an example for basic inlabru.
#'
#' @usage data(toypoints)
#'
#' @format The data are a list that contains these elements:
#'  \describe{
#'    \item{`points`}{An `sf` object of point locations and and `z` measurements}
#'    \item{`mesh`}{An `fm_mesh_2d` object}
#'    \item{`boundary`}{An `sf` polygon denting the region of interest}
#'    \item{`pred_locs`}{A `sf` object with prediction point locations}
#'  }
#'
#'
#' @examples
#' if (require("ggplot2")) {
#'   ggplot() +
#'     fmesher::geom_fm(data = toypoints$mesh, alpha = 0) +
#'     geom_sf(data = toypoints$boundary, fill = "blue", alpha = 0.1) +
#'     geom_sf(data = toypoints$points, aes(color = z)) +
#'     scale_color_viridis_c()
#' }
"toypoints"
