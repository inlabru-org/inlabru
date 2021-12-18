#' @name mrsea
#' @title Marine renewables strategic environmental assessment
#' @docType data
#' @description Data imported from package MRSea, see <http://creem2.st-andrews.ac.uk/software/>
#'
#' @usage data(mrsea)
#'
#' @format A list of objects:
#'  \describe{
#'    \item{`points`:}{ A `SpatialPointsDataFrame` object containing the locations of
#'    XXXXX.}
#'    \item{`samplers`:}{ A `SpatialLinesDataFrame` object containing the transect lines
#'    that were surveyed.}
#'    \item{`mesh`:}{ An `inla.mesh` object containing a Delaunay triangulation
#'    mesh (a type of discretization of continuous space) covering the survey region.}
#'    \item{`boundary`:}{ An `SpatialPolygonsDataFrame` object defining the boundary of the
#'    survey region.}
#'    \item{`covar`:}{ An `SpatialPointsDataFrame` containing sea depth estimates.}
#'  }
#' @source
#' Library `MRSea`.
#'
#' @references
#'
#' NONE YET
#'
#' @examples
#' if (require(ggplot2, quietly = TRUE)) {
#'   data(mrsea)
#'   ggplot() +
#'     gg(mrsea$mesh) +
#'     gg(mrsea$samplers) +
#'     gg(mrsea$points) +
#'     gg(mrsea$boundary)
#' }
NULL

#' MRSea data import
#'
#' Load mrsea survey data from the MRSea package and convert to spatial formats defined by the [sp] package.
#' Requires the MRSea package from <https://github.com/lindesaysh/MRSea>, and
#' is normally only run by the package maintainer. For regular inlabru use of the data,
#' use `data("MRSea", package = "inlabru")`, which does not require the MRSea package.
#'
#' @aliases import.mrsea
#' @keywords internal
#' @return The [mrsea] data set
#' @examples
#' \dontrun{
#' mrsea <- import.mrsea()
#' }
#' @author Lindesay Scott-Hayward \email{lass@@st-andrews.ac.uk}
#'

import.mrsea <- function() {
  pkg <- "MRSea"
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("This package development function require the MRSea package from https://github.com/lindesaysh/MRSea")
  }

  # library(MRSea)
  predict.data.re <- NULL
  dis.data.re <- NULL
  data("dis.data.re", package = "MRSea", envir = environment())
  data("predict.data.re", package = "MRSea", envir = environment())
  impact <- dis.data.re$impact
  segment.id <- dis.data.re$segment.id

  preddata_gpseas <- dplyr::group_by(predict.data.re, impact, segment.id)

  # Some housekeeping to change the name labels to the right format and correct the effort unit to match
  # the coordinate information.
  # Effort should have same units as coordinates:
  # 2021-02-16 the original data was stored in km units in MRSea
  dis.data.re <- mrseanames2dsmnames(dis.data.re)

  segdata <- dis.data.re[, c(
    "Transect.Label", "Transect.label", "season", "impact", "depth", "Sample.Label",
    "segment.label", "length", "Effort", "x", "y"
  )]
  segdata <- dplyr::distinct(segdata, segdata$Sample.Label, .keep_all = TRUE)

  # effort, object and distance.
  # Not taken x and y as these are segement mid points not detection locations
  distdata <- makedistdata(dis.data.re)

  # obsdata
  obsdata <- makeobsdata(dis.data.re)

  # 2021-02-16 the original data was stored in metres units in MRSea; convert
  preddata <- makepreddata(predict.data.re)

  dsmdata <- list(obsdata = obsdata, distdata = distdata, segdata = segdata, preddata = preddata)
  dset <- import.dsmdata(dsmdata, covar.col = 5)

  # Depth data to the data set

  depth <- makecovardata(predict.data.re)

  ############ NEW FORMAT USING sp objects ##############

  crs <- CRS("+proj=utm +zone=32 +units=km")

  # Transects lines
  lns <- subset(dset$effort, is.na(det))
  class(lns) <- "data.frame"
  lns <- sline(lns, c("start.x", "start.y"), c("end.x", "end.y"), crs = crs)
  lns$weight <- lns$Effort

  # Detections
  pts <- subset(dset$effort, !is.na(det))
  class(pts) <- "data.frame"
  coordinates(pts) <- c("x", "y")
  proj4string(pts) <- crs

  # Mesh
  mesh <- dset$mesh
  mesh$crs <- crs


  # Boundary
  boundary <- spoly(dset$mesh$loc[dset$mesh$segm$int$idx[, 1], 1:2],
    cols = c(1, 2),
    crs = crs
  )

  # Covariates
  covar <- SpatialPointsDataFrame(depth[, 1:2],
    data = depth[, 3, drop = FALSE],
    proj4string = crs
  )

  # Remove `distance` column from transects
  lns$distance <- NULL

  mrsea <- list(
    points = pts,
    samplers = lns,
    boundary = boundary,
    mesh = mesh,
    covar = covar
  )
}

# mrsea <- import.mrsea()
# use_data(mrsea)

mrseanames2dsmnames <- function(data) {
  nam <- names(data)
  cols2change <- c("transect.id", "transect.label", "x.pos", "y.pos", "segment.id")
  id <- NULL
  for (i in seq_along(cols2change)) {
    id <- c(id, grep(cols2change[i], nam))
  }
  names(data)[id] <- c("Transect.Label", "Transect.label", "x", "y", "Sample.Label")
  # make sure effort column is same unit as coordinate system
  data$Effort <- data$length
  # At 2021-01-16, the distances were stored in metres, so need to convert
  data$distance <- data$distance / 1000
  return(data)
}
# ---------------------------------------------------------------------
makepreddata <- function(data) {
  id <- c(grep("x.pos", colnames(data)), grep("y.pos", colnames(data)))
  names(data)[id] <- c("x", "y")
  # 2021-02-16 the original data was stored in metres units in MRSea; convert:
  data$x <- data$x / 1000
  data$y <- data$y / 1000
  return(data)
}
makecovardata <- function(data) {
  depth <- data[
    (data$season == 1) &
      (data$impact == 0),
    c("x.pos", "y.pos", "depth")
  ]
  colnames(depth) <- c("x", "y", "depth")
  data$x <- data$x / 1000
  data$y <- data$y / 1000
  data
}
# ---------------------------------------------------------------------
makedistdata <- function(data) {
  # effort, object and distance.
  # Not taken x and y as these are segement mid points not detection locations
  distdata <- data[, c("object", "distance", "Effort")]
  if (is.null(data$size)) {
    distdata$size <- rep(1, nrow(distdata))
  } else {
    distdata$size <- data$size
  }
  distdata <- na.omit(distdata)
  return(distdata)
}
# ---------------------------------------------------------------------
makeobsdata <- function(data) {
  obsdata <- na.omit(data)
  if (is.null(data$size)) {
    obsdata$size <- rep(1, nrow(obsdata))
  } else {
    obsdata$size <- obsdata$size
  }
  obsdata <- obsdata[, c("object", "Sample.Label", "distance", "Effort", "size")]
  return(obsdata)
}
# ---------------------------------------------------------------------
