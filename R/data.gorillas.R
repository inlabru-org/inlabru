#' @name gorillas
#' @title Gorilla nesting sites
#' @docType data
#' @description This is the `gorillas` dataset from the package `spatstat.data`, reformatted
#' as point process data for use with `inlabru`.
#'
#' @usage data(gorillas)
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
#'   require(ggplot2, quietly = TRUE) &&
#'   require(ggpolypath, quietly = TRUE)) {
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
NULL

#' Gorilla data import
#'
#' @aliases import.gorillas
#' @keywords internal
#' @importFrom utils data
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}, David L. Borchers \email{dlb@@st-andrews.ac.uk}
#' @return gorilla data

import.gorillas <- function() {
  if (!check_spatstat("spatstat.data")) {
    # Load Gorilla data from spatstat
    gorillas <- NULL
  }
  gorillas.extra <- NULL
  data(gorillas, package = "spatstat.data", envir = environment())

  # Create SpatialPoints representing nest locations
  requireNamespace("spatstat.geom")
  nests <- as.data.frame(gorillas)
  coordinates(nests) <- c("x", "y")
  crs <- sp::CRS("+proj=utm +zone=32 N +datum=WGS84") # from the Gorillas help file
  crs_km <- sp::CRS("+proj=utm +zone=32 N +datum=WGS84 +units=km")
  proj4string(nests) <- crs

  # Turn the observation window into spatial polygon
  boundary <- spoly(as.data.frame(gorillas$window$bdry[[1]]),
    crs = crs
  )

  # Build mesh
  bnd <- fm_as_inla_mesh_segment(boundary)
  mesh <- INLA::inla.mesh.2d(
    interior = bnd, max.edge = 222,
    crs = crs
  ) # With higher max.edge we run into various INLA errors/warnings

  # Turn covariates int SpatialGridDataFrame
  gcov <- list()
  for (nm in names(gorillas.extra)) {
    gcov[[nm]] <- as(gorillas.extra[[nm]], "SpatialGridDataFrame")
    proj4string(gcov[[nm]]) <- proj4string(nests)
    coordnames(gcov[[nm]]) <- c("x", "y")
    names(gcov[[nm]]) <- nm
    if (is.character(gcov[[nm]][[nm]])) {
      gcov[[nm]][[nm]] <- as.factor(gcov[[nm]][[nm]])
    }
  }

  # Hack: Change CRS units of the covariates to km
  for (nm in names(gcov)) {
    ga <- attributes(gcov[[nm]])$grid
    attributes(ga)$cellcentre.offset <- attributes(ga)$cellcentre.offset / 1000
    attributes(ga)$cellsize <- attributes(ga)$cellsize / 1000
    attributes(gcov[[nm]])$grid <- ga
    attributes(gcov[[nm]])$bbox <- attributes(gcov[[nm]])$bbox / 1000
    attributes(gcov[[nm]])$proj4string <- crs_km
  }


  ### Make final gorilla data set
  gorillas <- list(
    nests = nests,
    mesh = mesh,
    boundary = boundary,
    gcov = gcov
  )

  gorillas$nests <- fm_transform(gorillas$nests, crs_km)
  gorillas$mesh <- fm_transform(gorillas$mesh, crs_km)
  gorillas$boundary <- fm_transform(gorillas$boundary, crs_km)

  # Create a plot sampling data set
  set.seed(121)
  plotpts <- plotsample(gorillas$nests, gorillas$boundary,
    x.ppn = 0.6, y.ppn = 0.6, nx = 5.4, ny = 5.4
  )
  counts <- point2count(plotpts$plots, plotpts$dets)
  x <- coordinates(counts)[, 1]
  y <- coordinates(counts)[, 2]
  count <- counts@data$n

  # Make gam data frame
  gnestcount_9x9_60pc <-
    data.frame(x = x, y = y, count = count, exposure = counts$area)
  gnestplots_9x9_60pc <- plotpts$plots
  gnestdets_9x9_60pc <- plotpts$dets
  sample_9x9_60pc <- list(
    counts = gnestcount_9x9_60pc,
    plots = gnestplots_9x9_60pc,
    nests = gnestdets_9x9_60pc
  )

  # Attach plotsample to gorilla data
  gorillas$plotsample <- sample_9x9_60pc

  # Make the count data frame spatial
  coordinates(gorillas$plotsample$counts) <- c("x", "y")
  proj4string(gorillas$plotsample$counts) <- crs

  # Extrapolate covariate
  pxl <- fm_pixels(gorillas$mesh,
    mask = FALSE, nx = 220, ny = 180,
    format = "sp"
  )
  pxl <- fm_transform(pxl, fm_crs(gorillas$gcov[[1]]))
  for (k in names(gorillas$gcov)) {
    NA_value <- gorillas$gcov[[k]][[1]][1]
    is.na(NA_value) <- NA
    pxl[[k]] <- NA_value
    pxl[[k]] <- bru_fill_missing(gorillas$gcov[[k]], pxl, values = pxl[[k]])
  }
  gorillas$gcov <- pxl

  return(gorillas)
}


#' @describeIn import.gorillas Convert gorillas to `sf` and `terra` format
import.gorillas.sf <- function() {
  gorillas <- NULL
  data(gorillas, package = "inlabru", envir = environment())

  gorillas_sf <- list()
  gorillas_sf$nests <- sf::st_as_sf(gorillas$nests)
  gorillas_sf$mesh <- gorillas$mesh
  gorillas_sf$boundary <- sf::st_as_sf(gorillas$boundary)
  gorillas_sf$gcov <- terra::rast(gorillas$gcov[[1]])
  for (k in seq_len(length(gorillas$gcov) - 1L) + 1L) {
    terra::add(gorillas_sf$gcov) <- terra::rast(gorillas$gcov[[k]])
  }
  gorillas_sf$plotsample <- lapply(gorillas$plotsample, sf::st_as_sf)

  gorillas_sf
}


# gorillas <- import.gorillas()
# use_data(gorillas)
