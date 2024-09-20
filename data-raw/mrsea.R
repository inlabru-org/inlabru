source(here::here("data-raw", "dsmdata.tools.R"))
source(here::here("data-raw", "dsmdata.R"))

#' MRSea data import
#'
#' Load mrsea survey data from the MRSea package and convert to spatial formats
#' defined by the `sf` or `sp` packages.
#' Requires the MRSea package from <https://github.com/lindesaysh/MRSea>, and is
#' normally only run by the package maintainer. For regular inlabru use of the
#' data, use `data("mrsea", package = "inlabru")`, which does not require the
#' MRSea package.
#'
#' @keywords internal
#' @return The [mrsea] data set
#' @examples
#' \dontrun{
#' mrsea <- import.mrsea()
#' }
#' @author Lindesay Scott-Hayward \email{lass@@st-andrews.ac.uk}
#'

import.mrsea <- function(format = c("sf", "sp")) {
  format <- match.arg(format)
  pkg <- "MRSea"
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0(
      "This package development function require the MRSea package ",
      "from https://github.com/lindesaysh/MRSea"
    ))
  }

  # library(MRSea)
  predict.data.re <- NULL
  dis.data.re <- NULL
  data("dis.data.re", package = "MRSea", envir = environment())
  data("predict.data.re", package = "MRSea", envir = environment())
  impact <- dis.data.re$impact
  segment.id <- dis.data.re$segment.id

  preddata_gpseas <- dplyr::group_by(predict.data.re, impact, segment.id)

  # Some housekeeping to change the name labels to the right format and correct
  # the effort unit to match the coordinate information.
  # Effort should have same units as coordinates:
  # 2021-02-16 the original data was stored in km units in MRSea
  dis.data.re <- mrseanames2dsmnames(dis.data.re)

  segdata <- dis.data.re[, c(
    "Transect.Label",
    "Transect.label",
    "season",
    "impact",
    "depth",
    "Sample.Label",
    "segment.label",
    "length",
    "Effort",
    "x",
    "y"
  )]
  segdata <- dplyr::distinct(segdata, segdata$Sample.Label, .keep_all = TRUE)

  # effort, object and distance.
  # Not taken x and y as these are segement mid points not detection locations
  distdata <- makedistdata(dis.data.re)

  # obsdata
  obsdata <- makeobsdata(dis.data.re)

  # 2021-02-16 the original data was stored in metres units in MRSea; convert
  preddata <- makepreddata(predict.data.re)

  dsmdata <- list(
    obsdata = obsdata,
    distdata = distdata,
    segdata = segdata,
    preddata = preddata
  )
  dset <- import.dsmdata(dsmdata, covar.col = 5)

  # Depth data to the data set

  depth <- makecovardata(predict.data.re)

  if (format == "sp") {
    ############ FORMAT USING sp objects ##############

    crs <- fm_CRS("+proj=utm +zone=32 +units=km")

    # Transects lines
    lns <- subset(dset$effort, is.na(det))
    class(lns) <- "data.frame"
    lns <- sline(lns, c("start.x", "start.y"), c("end.x", "end.y"), crs = crs)
    lns$weight <- lns$Effort

    # Detections
    pts <- subset(dset$effort, !is.na(det))
    class(pts) <- "data.frame"
    sp::coordinates(pts) <- c("x", "y")
    sp::proj4string(pts) <- crs

    # Mesh
    mesh <- dset$mesh
    mesh$crs <- fm_crs(crs)


    # Boundary
    boundary <- spoly(dset$mesh$loc[dset$mesh$segm$int$idx[, 1], 1:2],
      cols = c(1, 2),
      crs = crs
    )

    # Covariates
    covar <- sp::SpatialPointsDataFrame(depth[, c("x", "y")],
      data = depth[, "depth", drop = FALSE],
      proj4string = crs
    )

    # Remove `distance` column from transects
    lns$distance <- NULL
  } else {
    ############ FORMAT USING sf objects ##############

    crs <- fm_crs("+proj=utm +zone=32 +units=km")

    # Transects lines
    lns <- subset(dset$effort, is.na(det))
    class(lns) <- "data.frame"
    lns <- sline(lns, c("start.x", "start.y"), c("end.x", "end.y"),
      crs = crs,
      format = "sf"
    )
    lns$weight <- lns$Effort

    # Detections
    pts <- subset(dset$effort, !is.na(det))
    class(pts) <- "data.frame"
    pts <- sf::st_as_sf(pts, coords = c("x", "y"), crs = crs)

    # Mesh
    mesh <- dset$mesh
    mesh$crs <- crs

    # Boundary
    boundary <- spoly(dset$mesh$loc[dset$mesh$segm$int$idx[, 1], 1:2],
      cols = c(1, 2),
      crs = crs,
      format = "sf"
    )

    # Covariates
    covar <- sf::st_as_sf(depth, coords = c("x", "y"), crs = crs)

    # Remove `distance` column from transects
    lns$distance <- NULL
  }

  mrsea <- list(
    points = pts,
    samplers = lns,
    boundary = boundary,
    mesh = mesh,
    covar = covar
  )
}

mrseanames2dsmnames <- function(data) {
  nam <- names(data)
  cols2change <- c("transect.id",
                   "transect.label",
                   "x.pos",
                   "y.pos",
                   "segment.id")
  id <- NULL
  for (i in seq_along(cols2change)) {
    id <- c(id, grep(cols2change[i], nam))
  }
  names(data)[id] <- c(
    "Transect.Label",
    "Transect.label",
    "x",
    "y",
    "Sample.Label"
  )
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
  depth$x <- depth$x / 1000
  depth$y <- depth$y / 1000
  depth
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
  obsdata <- obsdata[, c("object",
                         "Sample.Label",
                         "distance",
                         "Effort",
                         "size")]
  return(obsdata)
}
# ---------------------------------------------------------------------

# mrsea <- import.mrsea(format = "sp")
# usethis::use_data(mrsea, overwrite = TRUE)
# mrsea <- import.mrsea(format = "sf")
# use_data(mrsea, compress = "xz", overwrite = TRUE)
