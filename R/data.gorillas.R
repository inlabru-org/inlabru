#' @title Gorilla nesting sites in sf format
#' @docType data
#' @description This is the `gorillas` dataset from the package `spatstat.data`,
#'   reformatted as point process data for use with `inlabru`.
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
#'    \item{`gcov_file`:}{ The in-package filename of a `terra::SpatRaster`
#'    object, with one layer for each of these spatial covariates:
#'     \describe{
#'       \item{`aspect`}{ Compass direction of the terrain slope. Categorical,
#'       with levels
#'       N, NE, E, SE, S, SW, W and NW, which are coded as integers 1 to 8.}
#'       \item{`elevation`}{ Digital elevation of terrain, in metres.}
#'       \item{`heat`}{ Heat Load Index at each point on the surface
#'       (Beer's aspect),
#'       discretised. Categorical with values Warmest (Beer's aspect between
#'       0 and 0.999),
#'       Moderate (Beer's aspect between 1 and 1.999), Coolest (Beer's aspect
#'       equals 2). These are
#'       coded as integers 1, 2 and 3, in that order.}
#'       \item{`slopangle`}{ Terrain slope, in degrees.}
#'       \item{`slopetype`}{ Type of slope. Categorical, with values Valley,
#'       Toe (toe slope),
#'       Flat, Midslope, Upper and Ridge. These are coded as integers 1 to 6.}
#'       \item{`vegetation`}{ Vegetation type: a categorical variable with 6
#'       levels coded as
#'       integers 1 to 6 (in order of increasing expected habitat suitability)}
#'       \item{`waterdist`}{ Euclidean distance from nearest water body, in
#'       metres.}
#'     }
#'     Loading of the covariates can be done with `gorillas_sf_gcov()` or
#'
#'     gorillas_sf$gcov <- terra::rast(
#'       system.file(gorillas_sf$gcov_file, package = "inlabru")
#'     )
#'
#'    }
#'    \item{`plotsample`}{Plot sample of gorilla nests, sampling 9x9 over the
#'    region, with 60\% coverage. Components:
#'    \describe{
#'      \item{`counts`}{An `sf` object with elements
#'      `count`, `exposure`, and `geometry`, holding the point geometry for the
#'      centre of each plot, the count in each
#'      plot and the area of each plot.}
#'      \item{`plots`}{An `sf` object with `MULTIPOLYGON` objects defining the
#'      individual plot boundaries and an all-ones `weight` column.}
#'      \item{`nests`}{An `sf` giving the locations of
#'      each detected nests, `group` ("minor" or "major"),
#'      `season` ("dry" or "rainy"), and `date` (in `Date` format).}
#'    }
#'    }
#'  }
#' @source
#' Library `spatstat.data`.
#'
#'
#' @references
#' Funwi-Gabga, N. (2008) A pastoralist survey and fire impact assessment in the
#' Kagwene Gorilla Sanctuary, Cameroon. M.Sc. thesis, Geology and Environmental
#' Science, University of Buea, Cameroon.
#'
#' Funwi-Gabga, N. and Mateu, J. (2012) Understanding the nesting spatial
#' behaviour of gorillas in the Kagwene Sanctuary, Cameroon. Stochastic
#' Environmental Research and Risk Assessment 26 (6), 793-811.
#'
#' @examples
#' if (interactive() &&
#'   require(ggplot2, quietly = TRUE) &&
#'   requireNamespace("terra", quietly = TRUE) &&
#'   requireNamespace("tidyterra", quietly = TRUE)) {
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
#'   gorillas_sf$gcov <- gorillas_sf_gcov()
#'   ggplot() +
#'     tidyterra::geom_spatraster(data = gorillas_sf$gcov$elevation)
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
#' if (interactive() &&
#'   requireNamespace("terra", quietly = TRUE)) {
#'   gorillas_sf$gcov <- gorillas_sf_gcov()
#' }
gorillas_sf_gcov <- function() {
  requireNamespace("terra")
  terra::rast(system.file(inlabru::gorillas_sf$gcov_file, package = "inlabru"))
}

#' @describeIn gorillas_sf Access the `gorillas_sf` data in `sp` format.
#' The covariate data is added as `gcov`, a list of `sp::SpatialPixelsDataFrame`
#' objects. Requires the `sp`, `sf`, and `terra` packages to be installed.
#' @export
gorillas_sp <- function() {
  requireNamespace("sf")
  requireNamespace("terra")
  stopifnot(bru_safe_sp())
  dat <- inlabru::gorillas_sf
  gcov <- gorillas_sf_gcov()
  gcov_ <- as.data.frame(gcov, xy = TRUE)
  gcov_sp <- list()
  for (nm in names(gcov)) {
    gcov_sp[[nm]] <- sp::SpatialPixelsDataFrame(
      cbind(gcov_$x, gcov_$y),
      data = gcov_[, nm, drop = FALSE],
      proj4string = fm_CRS(gcov)
    )
  }

  out <- list(
    nests = sf::as_Spatial(dat$nests),
    mesh = dat$mesh,
    boundary = sf::as_Spatial(dat$boundary),
    plotsample = list(
      counts = sf::as_Spatial(dat$plotsample$counts),
      plots = sf::as_Spatial(dat$plotsample$plots),
      nests = sf::as_Spatial(dat$plotsample$nests)
    ),
    gcov = gcov_sp
  )
  out
}

#' @title Deprecated alias for sp version of the gorillas dataset
#' @name gorillas
#' @rdname gorillas
#' @description
#' Deprecated dataset name for the `sp` version
#' of [gorillas_sf]. Use [gorillas_sp()] instead.
#' @seealso [gorillas_sf]
#' @keywords internal
NULL
