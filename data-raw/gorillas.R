#' Gorilla data import
#'
#' @keywords internal
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com},
#' David L. Borchers \email{dlb@@st-andrews.ac.uk}
#' @return gorilla data

import_gorillas_sp <- function() {
  stopifnot(inlabru:::check_spatstat("spatstat.data"))
  gorillas <- spatstat.data::gorillas
  gorillas.extra <- spatstat.data::gorillas.extra

  # Create SpatialPoints representing nest locations
  requireNamespace("spatstat.geom")
  nests <- as.data.frame(gorillas)
  sp::coordinates(nests) <- c("x", "y")
  # from the Gorillas help file:
  crs <- sp::CRS("+proj=utm +zone=32 N +datum=WGS84")
  crs_km <- sp::CRS("+proj=utm +zone=32 N +datum=WGS84 +units=km")
  sp::proj4string(nests) <- crs

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
    sp::proj4string(gcov[[nm]]) <- sp::proj4string(nests)
    sp::coordnames(gcov[[nm]]) <- c("x", "y")
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
  x <- sp::coordinates(counts)[, 1]
  y <- sp::coordinates(counts)[, 2]
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
  sp::coordinates(gorillas$plotsample$counts) <- c("x", "y")
  sp::proj4string(gorillas$plotsample$counts) <- crs_km

  # Extrapolate covariate
  pxl_ <- fm_pixels(gorillas$mesh,
                    mask = FALSE, dims = c(220, 180),
                    format = "sp"
  )
  pxl_ <- fm_transform(pxl_, fm_crs(gorillas$gcov[[1]]))
  for (k in names(gorillas$gcov)) {
    pxl <- pxl_
    NA_value <- gorillas$gcov[[k]][[1]][1]
    is.na(NA_value) <- TRUE
    pxl[[k]] <- NA_value
    pxl[[k]] <- bru_fill_missing(gorillas$gcov[[k]], pxl, values = pxl[[k]])
    gorillas$gcov[[k]] <- pxl
  }

  return(gorillas)
}


#' @describeIn import_gorillas_sp Import gorillas in `sf` and `terra` format
import_gorillas_sf <- function(overwrite = FALSE) {
  stopifnot(inlabru:::check_spatstat("spatstat.data"))
  gorillas <- spatstat.data::gorillas
  gorillas.extra <- spatstat.data::gorillas.extra

  requireNamespace("spatstat.geom")
  # from the Gorillas help file:
  crs <- fm_crs("+proj=utm +zone=32 N +datum=WGS84")
  crs_km <- fm_crs("+proj=utm +zone=32 N +datum=WGS84 +units=km")

  # Create sf representing nest locations
  nests <- sf::st_as_sf(as.data.frame(gorillas),
                        coords = c("x", "y"),
                        crs = crs)

  # Turn the observation window into spatial polygon
  boundary <- spoly(as.data.frame(gorillas$window$bdry[[1]]),
                    crs = crs,
                    format = "sf")

  # Build mesh
  bnd <- fm_as_segm(boundary)
  mesh <- fm_mesh_2d_inla(
    interior = bnd,
    max.edge = 222,
    crs = crs
  ) # With higher max.edge we run into various INLA errors/warnings

  # Turn covariates into SGDF
  gcov <- NULL
  for (nm in names(gorillas.extra)) {
    gcov_ <- terra::rast(gorillas.extra[[nm]])
    terra::crs(gcov_) <- fm_proj4string(crs)
    names(gcov_) <- nm
    if (is.character(terra::values(gcov_,
                                   mat = FALSE,
                                   dataframe = TRUE)[[1]])) {
      terra::values(gcov_) <- as.factor(terra::values(gcov_))
    }
    if (is.null(gcov)) {
      gcov <- gcov_
    } else {
      terra::add(gcov) <- gcov_
    }
  }

  ### Make final gorilla data set
  gorillas_sf <- list(
    nests = fm_transform(nests, crs_km),
    mesh = fm_transform(mesh, crs_km),
    boundary = fm_transform(boundary, crs_km)
  )

  # Create a plot sampling data set
  set.seed(121)
  plotpts <- plotsample(
    sf::as_Spatial(gorillas_sf$nests),
    sf::as_Spatial(gorillas_sf$boundary),
    x.ppn = 0.6, y.ppn = 0.6, nx = 5.4, ny = 5.4
  )
  counts <- point2count(plotpts$plots, plotpts$dets)
  counts <- sf::st_as_sf(counts, crs = crs_km)
  x <- sf::st_coordinates(counts)[, 1]
  y <- sf::st_coordinates(counts)[, 2]
  count <- counts$n

  # Make gam data frame
  gnestcount_9x9_60pc <-
    sf::st_as_sf(
      data.frame(
        x = x,
        y = y,
        count = count,
        exposure = counts$area
      ),
      coords = c("x", "y"),
      crs = crs_km
    )
  gnestplots_9x9_60pc <- sf::st_as_sf(plotpts$plots, crs = crs_km)
  gnestdets_9x9_60pc <- sf::st_as_sf(plotpts$dets, crs = crs_km)
  sample_9x9_60pc <- list(
    counts = gnestcount_9x9_60pc,
    plots = gnestplots_9x9_60pc,
    nests = gnestdets_9x9_60pc
  )

  # Attach plotsample to gorilla data
  gorillas_sf$plotsample <- sample_9x9_60pc

  # Extrapolate covariate
  pxl_ <- fm_pixels(
    gorillas_sf$mesh,
    mask = FALSE,
    dims = c(220, 180),
    format = "sf"
  )
  pxl_c <- sf::st_coordinates(pxl_)
  ex <- c(range(pxl_c[, 1]), range(pxl_c[, 2]))
  ex[1:2] <- ex[1:2] + c(-1, 1) * 0.5 * diff(ex[1:2]) / (220 - 1)
  ex[3:4] <- ex[3:4] + c(-1, 1) * 0.5 * diff(ex[3:4]) / (180 - 1)
  gcov_ext <- terra::rast(
    nrows = 180,
    ncols = 220,
    extent = ex,
    crs = fm_proj4string(crs_km),
    nlyrs = terra::nlyr(gcov),
    names = names(gcov)
  )
  pxl <- fm_transform(pxl_, fm_crs(gcov))
  for (k in names(gcov)) {
    NA_value <- terra::values(gcov[[k]], mat = FALSE, dataframe = TRUE)[[1]][1]
    if (is.numeric(NA_value)) {
      NA_value <- NA_real_
    } else {
      is.na(NA_value) <- TRUE
    }
    vals <- rep(NA_value, 180 * 220)
    vals <- bru_fill_missing(gcov[[k]], pxl, values = vals)
    vals <- matrix(vals, nrow = 220, ncol = 180)
    vals <- vals[, rev(seq_len(180))]
    terra::values(gcov_ext[[k]]) <- as.vector(vals)
  }

  gorillas_sf$gcov_file <- "extdata/gorillas_cov.tif"
  filename <- file.path("inst", gorillas_sf$gcov_file)
  if (!file.exists(filename) || overwrite) {
    if (file.exists(filename)) {
      message(paste0("Overwriting existing ", filename))
    }

    terra::writeRaster(
      gcov_ext,
      datatype = "FLT8S", # Avoid integer storage of doubles
      filename = filename,
      overwrite = overwrite,
      gdal = c("COMPRESS=LZW")
    )
  }

  gorillas_sf
}


#' @describeIn import_gorillas_sp Convert gorillas to `sf` and `terra` format
import_gorillas_sf_old <- function(gorillas = NULL, overwrite = FALSE) {
  if (is.null(gorillas)) {
    gorillas <- inlabru::gorillas
  }

  gorillas_sf <- list()
  gorillas_sf$nests <- sf::st_as_sf(gorillas$nests)
  gorillas_sf$mesh <- fmesher::fm_as_fm(gorillas$mesh)
  gorillas_sf$boundary <- sf::st_as_sf(gorillas$boundary)
  gorillas_sf$plotsample <- lapply(gorillas$plotsample, sf::st_as_sf)

  gorillas_sf$gcov_file <- "extdata/gorillas_cov.tif"
  filename <- file.path("inst", gorillas_sf$gcov_file)
  if (!file.exists(filename) || overwrite) {
    if (file.exists(filename)) {
      message(paste0("Overwriting existing ", filename))
    }
    gcov <- terra::rast(gorillas$gcov[[1]])
    for (k in seq_len(length(gorillas$gcov) - 1L) + 1L) {
      terra::add(gcov) <- terra::rast(gorillas$gcov[[k]])
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


# gorillas_sf <- import_gorillas_sf(overwrite = TRUE)
# usethis::use_data(gorillas_sf, overwrite = TRUE, compress = "xz")
