#' Gorilla data import
#'
#' @aliases import.gorillas
#' @keywords internal
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}, David L. Borchers \email{dlb@@st-andrews.ac.uk}
#' @return gorilla data

import.gorillas <- function() {
  stopifnot(check_spatstat("spatstat.data"))
  gorillas <- spatstat.data::gorillas
  gorillas.extra <- spatstat.data::gorillas.extra

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
  bnd <- fm_as_segm(boundary)
  mesh <- fm_mesh_2d_inla(
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
    mask = FALSE, dims = c(220, 180),
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
import.gorillas.sf <- function(overwrite = FALSE) {
  gorillas <- inlabru::gorillas

  gorillas_sf <- list()
  gorillas_sf$nests <- sf::st_as_sf(gorillas$nests)
  gorillas_sf$mesh <- fmesher::fm_as_fm(gorillas$mesh)
  gorillas_sf$boundary <- sf::st_as_sf(gorillas$boundary)
  gcov <- terra::rast(gorillas$gcov[[1]])
  for (k in seq_len(length(gorillas$gcov) - 1L) + 1L) {
    terra::add(gcov) <- terra::rast(gorillas$gcov[[k]])
  }
  gorillas_sf$plotsample <- lapply(gorillas$plotsample, sf::st_as_sf)

  gorillas_sf$gcov_file <- "extdata/gorillas_cov.tif"
  filename <- file.path("inst", gorillas_sf$gcov_file)
  if (!file.exists(filename) || overwrite) {
    if (file.exists(filename)) {
      message(paste0("Overwriting existing ", filename))
    }
    terra::writeRaster(
      gcov,
      filename = filename,
      overwrite = overwrite,
      gdal = c("COMPRESS=LZW")
    )
  }

  gorillas_sf
}


# gorillas <- import.gorillas()
# use_data(gorillas)
# gorillas_sf <- import.gorillas.sf()
# usethis::use_data(gorillas_sf, overwrite = TRUE)
