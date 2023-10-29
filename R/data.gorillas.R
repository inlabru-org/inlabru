#' @title Gorilla nesting sites
#' @docType data
#' @description This is the `gorillas` dataset from the package `spatstat.data`, reformatted
#' as point process data for use with `inlabru`.
#'
#' @usage
#' gorillas
#' # To avoid the name clash with spatstat.data::gorillas, use
#' data(gorillas, package = "inlabru")
#'
#' @format The data are a list that contains these elements:
#'  \describe{
#'    \item{`nests`:}{ A `SpatialPointsDataFrame` object containing the locations of
#'    the gorilla nests.}
#'    \item{`boundary`:}{ An `SpatialPolygonsDataFrame` object defining the boundary
#'    of the region that was searched for the nests.}
#'    \item{`mesh`:}{ An `inla.mesh` object containing a mesh that can be used
#'    with function `lgcp` to fit a LGCP to the nest data.}
#'    \item{`gcov`:}{ A list of SpatialGridDataFrame objects, one for each of these spatial covariates:
#'     \describe{
#'       \item{`aspect`}{ Compass direction of the terrain slope. Categorical, with levels
#'       N, NE, E, SE, S, SW, W and NW, which are coded as integers 1 to 8.}
#'       \item{`elevation`}{ Digital elevation of terrain, in metres.}
#'       \item{`heat`}{ Heat Load Index at each point on the surface (Beer's aspect),
#'       discretised. Categorical with values Warmest (Beer's aspect between 0 and 0.999),
#'       Moderate (Beer's aspect between 1 and 1.999), Coolest (Beer's aspect equals 2). These are
#'       coded as integers 1, 2 and 3, in that order.}
#'       \item{`slopangle`}{ Terrain slope, in degrees.}
#'       \item{`slopetype`}{ Type of slope. Categorical, with values Valley, Toe (toe slope),
#'       Flat, Midslope, Upper and Ridge. These are coded as integers 1 to 6.}
#'       \item{`vegetation`}{ Vegetation type: a categorical variable with 6 levels coded as
#'       integers 1 to 6 (in order of increasing expected habitat suitability)}
#'       \item{`waterdist`}{ Euclidean distance from nearest water body, in metres.}
#'     }
#'    }
#'    \item{`plotsample`}{Plot sample of gorilla nests, sampling 9x9 over the region, with 60\% coverage. Components:
#'    \describe{
#'      \item{`counts`}{ A SpatialPointsDataFrame frame with elements `x`, `y`, `count`,
#'      `exposure`, being the x- and y-coordinates of the centre of each plot, the count in each
#'      plot and the area of each plot.}
#'      \item{`plots`}{ A `SpatialPolygonsDataFrame` defining the individual plot boundaries.}
#'      \item{`nests`}{ A `SpatialPointsDataFrame` giving the locations of each detected nest.}
#'    }
#'    }
#'  }
#' @source
#' Library `spatstat.data`.
#'
#'
#' @references
#' Funwi-Gabga, N. (2008) A pastoralist survey and fire impact assessment in the Kagwene Gorilla
#' Sanctuary, Cameroon. M.Sc. thesis, Geology and Environmental Science, University of Buea,
#' Cameroon.
#'
#' Funwi-Gabga, N. and Mateu, J. (2012) Understanding the nesting spatial behaviour of gorillas
#' in the Kagwene Sanctuary, Cameroon. Stochastic Environmental Research and Risk Assessment
#' 26 (6), 793-811.
#'
#' @examples
#' if (bru_safe_inla() &&
#'   bru_safe_sp() &&
#'   require("sp") &&
#'   require(ggplot2, quietly = TRUE)) {
#'   data(gorillas, package = "inlabru") # get the data
#'
#'   # plot all the nests, mesh and boundary
#'   ggplot() +
#'     gg(gorillas$mesh) +
#'     gg(gorillas$boundary) +
#'     gg(gorillas$nests)
#'
#'   # Plot the elevation covariate
#'   plot(gorillas$gcov$elevation)
#'
#'   # Plot the plot sample
#'   ggplot() +
#'     gg(gorillas$plotsample$plots) +
#'     gg(gorillas$plotsample$nests)
#' }
"gorillas"

#' @title Gorilla nesting sites in sf format
#' @docType data
#' @description This is the `gorillas` dataset from the package `spatstat.data`, reformatted
#' as point process data for use with `inlabru`.
#'
#' @usage
#' gorillas_sf
#' data(gorillas_sf, package = "inlabru")
#'
#' @format The data are a list that contains these elements:
#'  \describe{
#'    \item{`nests`:}{ An `sf` object containing the locations of
#'    the gorilla nests.}
#'    \item{`boundary`:}{ An `sf` object defining the boundary
#'    of the region that was searched for the nests.}
#'    \item{`mesh`:}{ An `fm_mesh_2d` object containing a mesh that can be used
#'    with function `lgcp` to fit a LGCP to the nest data.}
#'    \item{`gcov_file`:}{ The in-package filename of a `terra::SpatRaster` object,
#'    with one layer for each of these spatial covariates:
#'     \describe{
#'       \item{`aspect`}{ Compass direction of the terrain slope. Categorical, with levels
#'       N, NE, E, SE, S, SW, W and NW, which are coded as integers 1 to 8.}
#'       \item{`elevation`}{ Digital elevation of terrain, in metres.}
#'       \item{`heat`}{ Heat Load Index at each point on the surface (Beer's aspect),
#'       discretised. Categorical with values Warmest (Beer's aspect between 0 and 0.999),
#'       Moderate (Beer's aspect between 1 and 1.999), Coolest (Beer's aspect equals 2). These are
#'       coded as integers 1, 2 and 3, in that order.}
#'       \item{`slopangle`}{ Terrain slope, in degrees.}
#'       \item{`slopetype`}{ Type of slope. Categorical, with values Valley, Toe (toe slope),
#'       Flat, Midslope, Upper and Ridge. These are coded as integers 1 to 6.}
#'       \item{`vegetation`}{ Vegetation type: a categorical variable with 6 levels coded as
#'       integers 1 to 6 (in order of increasing expected habitat suitability)}
#'       \item{`waterdist`}{ Euclidean distance from nearest water body, in metres.}
#'     }
#'     Loading of the covariates can be done with `gorillas_sf_gcov()` or
#'
#'     gorillas_sf$gcov <- terra::rast(
#'       system.file(gorillas_sf$gcov_file, package = "inlabru")
#'     )
#'
#'    }
#'    \item{`plotsample`}{Plot sample of gorilla nests, sampling 9x9 over the region, with 60\% coverage. Components:
#'    \describe{
#'      \item{`counts`}{ A SpatialPointsDataFrame frame with elements `x`, `y`, `count`,
#'      `exposure`, being the x- and y-coordinates of the centre of each plot, the count in each
#'      plot and the area of each plot.}
#'      \item{`plots`}{ A `SpatialPolygonsDataFrame` defining the individual plot boundaries.}
#'      \item{`nests`}{ A `SpatialPointsDataFrame` giving the locations of each detected nest.}
#'    }
#'    }
#'  }
#' @source
#' Library `spatstat.data`.
#'
#'
#' @references
#' Funwi-Gabga, N. (2008) A pastoralist survey and fire impact assessment in the Kagwene Gorilla
#' Sanctuary, Cameroon. M.Sc. thesis, Geology and Environmental Science, University of Buea,
#' Cameroon.
#'
#' Funwi-Gabga, N. and Mateu, J. (2012) Understanding the nesting spatial behaviour of gorillas
#' in the Kagwene Sanctuary, Cameroon. Stochastic Environmental Research and Risk Assessment
#' 26 (6), 793-811.
#'
#' @examples
#' if (interactive() &&
#'   bru_safe_inla() &&
#'   bru_safe_sp() &&
#'   require("sp") &&
#'   require(ggplot2, quietly = TRUE) &&
#'   requireNamespace("terra")) {
#'   # plot all the nests, mesh and boundary
#'   ggplot() +
#'     gg(gorillas_sf$mesh) +
#'     geom_sf(
#'       data = gorillas_sf$boundary,
#'       alpha = 0.1, fill = "blue"
#'     ) +
#'     geom_sf(data = gorillas_sf$nests)
#'
#'   # Plot the elevation covariate
#'   gorillas_sf$gcov <- terra::rast(
#'     system.file(gorillas_sf$gcov_file, package = "inlabru")
#'   )
#'   plot(gorillas_sf$gcov$elevation)
#'
#'   # Plot the plot sample
#'   ggplot() +
#'     geom_sf(data = gorillas_sf$plotsample$plots) +
#'     geom_sf(data = gorillas_sf$plotsample$nests)
#' }
"gorillas_sf"


#' @describeIn gorillas_sf Access the `gorillas_sf` covariates data as a
#' `terra::rast()` object.
#' @export
#' @examples
#' \dontrun{
#' gorillas_sf$gcov <- gorillas_sf_gcov()
#' }
gorillas_sf_gcov <- function() {
  terra::rast(system.file(inlabru::gorillas_sf$gcov_file, package = "inlabru"))
}
