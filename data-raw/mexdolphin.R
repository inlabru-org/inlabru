source(here::here("data-raw", "dsmdata.tools.R"))
source(here::here("data-raw", "dsmdata.R"))


# Old methods needed to import mexdolphins data

as.spatial.dsdata <- function(dset, cnames, crs) {
  as.detection <- function(...) {
    UseMethod("as.detection")
  }
  detdata <- function(...) {
    UseMethod("detdata")
  }
  segdata <- function(...) {
    UseMethod("segdata")
  }

  as.detection.dsdata <- function(dsdata, ...) {
    as.detection.effort(dsdata$effort, ...)
  }

  detdata.dsdata <- function(data, detection = NULL, ...) {
    if (is.null(detection)) {
      data$effort[as.detection(data)[, "start"], ]
    } else {
      data$effort[detection[, "start"], ]
    }
  }

  segdata.dsdata <- function(data, ...) {
    data$effort[is.na(data$effort$det), ]
  }

  # Extract detection pointers from effort data
  #
  # @export
  # @param effort Effort data set
  # @return sighting sightings
  # @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

  as.detection.effort <- function(effort) {
    idx <- which(!(is.na(effort$det)))
    det <- data.frame(start = idx, end = idx)
    class(det) <- c("detection", "data.frame")
    det
  }

  # Convert dsdata into sf objects
  #
  # @export
  # @param dsdata dsdata object
  # @param cnames Column names of the coordinates, e.g. c("lon","lat")
  # @param crs Coordinate reference system, e.g. CRS("+proj=longlat")
  # @return list with spatial points, spatial lines and a mesh
  # @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

  # as.spatial.dsdata = function(dset, cnames, crs)

  # Extract detections and convert to sf points
  not_cnames <- setdiff(names(detdata(dset)), cnames)
  points <- sf::st_as_sf(
    tibble::as_tibble(detdata(dset)),
    coords = cnames,
    crs = fmesher::fm_crs(crs)
  )

  # Extract effort and convert to sf lines
  sp <- as.data.frame(segdata(dset)[, paste0("start.", cnames)])
  ep <- as.data.frame(segdata(dset)[, paste0("end.", cnames)])
  colnames(sp) <- cnames
  colnames(ep) <- cnames
  lilist <- lapply(seq_len(nrow(sp)), function(k) {
    sf::st_linestring(as.matrix(rbind(sp[k, ], ep[k, ])))
  })

  spl <- sf::st_as_sfc(lilist, crs = fmesher::fm_crs(crs))
  tmp_names <- setdiff(
    names(segdata(dset)),
    c(
      paste0("start.", cnames),
      paste0("end.", cnames)
    )
  )
  df <- as.data.frame(segdata(dset))[, tmp_names]
  rownames(df) <- seq_len(nrow(df))
  sf::st_geometry(df) <- spl
  samplers <- df

  # Return detections ("points") and effort ("samplers")
  list(points = points, samplers = samplers, mesh = dset$mesh)
}


#' Mexdolphin data import
#'
#'
#' Load `mexdolphins` survey data from `dsm` package and convert to spatial
#' formats defined by the `sf` package.
#'
#' @keywords internal
#' @return The [mexdolphin] data set
#' @examples
#' \dontrun{
#' mexdolphin <- import_mexdolphin_sf()
#' }
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com} and Finn Lindgren
#' \email{finn.lindgren@@gmail.com}
#'
import_mexdolphin_sf <- function() {
  envir <- new.env()
  data("mexdolphins", package = "dsm", envir = envir)
  mexdolphins <- as.list(envir)
  data.p4s <- paste0(
    "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 ",
    "+ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  )
  data_crs <- fmesher::fm_crs(data.p4s)

  dset <- import.dsmdata(mexdolphins, covar.col = 8)
  dset$mesh$crs <- fmesher::fm_crs(data_crs)
  mexdolphin <- as.spatial.dsdata(
    dset,
    cnames = c("x", "y"),
    crs = fmesher::fm_crs(data_crs)
  )

  # Target CRS
  target.p4s <- paste0(
    "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 ",
    "+ellps=GRS80 +datum=NAD83 +units=km +no_defs +towgs84=0,0,0"
  )
  target_crs <- fmesher::fm_crs(target.p4s)

  # Units to km
  mexdolphin$points <- fmesher::fm_transform(
    mexdolphin$points,
    crs = target_crs
  )
  mexdolphin$samplers <- fmesher::fm_transform(
    mexdolphin$samplers,
    crs = target_crs
  )
  mexdolphin$points$distance <- mexdolphin$points$distance / 1000
  mexdolphin$points$Effort <- mexdolphin$points$Effort / 1000
  mexdolphin$points$mid.x <- mexdolphin$points$mid.x / 1000
  mexdolphin$points$mid.y <- mexdolphin$points$mid.y / 1000
  mexdolphin$points$start.x <- mexdolphin$points$start.x / 1000
  mexdolphin$points$start.y <- mexdolphin$points$start.y / 1000
  mexdolphin$points$end.x <- mexdolphin$points$end.x / 1000
  mexdolphin$points$end.y <- mexdolphin$points$end.y / 1000
  mexdolphin$samplers$mid.x <- mexdolphin$samplers$mid.x / 1000
  mexdolphin$samplers$mid.y <- mexdolphin$samplers$mid.y / 1000
  mexdolphin$samplers$Effort <- mexdolphin$samplers$Effort / 1000
  mexdolphin$mesh <- fmesher::fm_transform(mexdolphin$mesh, crs = target_crs)
  mexdolphin$mesh$loc <- fmesher::fm_unify_coords(mexdolphin$mesh$loc)

  # Remove all-NA columns such as "distance" from samplers, since they may
  # interfere with normal usage.
  all_na <- vapply(
    names(mexdolphin$samplers),
    function(nm) {
      all(is.na(mexdolphin$samplers[[nm]]))
    },
    TRUE
  )
  mexdolphin$samplers <- mexdolphin$samplers[!all_na]

  ##### Prediction polygon
  poly_pred <- fmesher::fm_segm(mexdolphin$mesh, boundary = FALSE)
  fmesher::fm_is_bnd(poly_pred) <- TRUE
  fmesher::fm_crs(poly_pred) <- fmesher::fm_crs(target_crs)
  poly_pred <- sf::st_as_sf(
    tibble::tibble(weight = 1, geometry = fmesher::fm_as_sfc(poly_pred))
  )
  poly_pred <- fmesher::fm_zm(poly_pred, target = "XY")
  mexdolphin$ppoly <- poly_pred

  ##### Simulate a whole population #####
  distance <- seq(0, 8, length.out = 20)
  matern <- INLA::inla.spde2.pcmatern(
    mexdolphin$mesh,
    prior.range = c(50, 0.01),
    prior.sigma = c(2, 0.01)
  )
  cmps <-
    geometry + distance ~
      Intercept(1) +
      df.lsigma(1) +
      spat(geometry, model = matern)
  pred <- ~ log(1 - exp(-(distance / exp(df.lsigma))^-1)) +
    spat +
    Intercept +
    log(2)
  r <- lgcp(
    data = mexdolphin$points,
    samplers = cbind(mexdolphin$samplers, data.frame(weight = 1)),
    components = cmps,
    formula = pred,
    domain = list(
      geometry = mexdolphin$mesh,
      distance = distance
    ),
    options = list(
      bru_verbose = 3,
      control.inla = list(int.strategy = "eb")
    )
  )
  llambda <- r$summary.random$spat$mean + r$summary.fixed["Intercept", "mean"]

  smexdolphin <- mexdolphin
  smexdolphin$llambda <- llambda
  mexdolphin$lambda <- exp(llambda)

  set.seed(1234L)
  pts <- sample.lgcp(
    mexdolphin$mesh,
    loglambda = llambda,
    samplers = mexdolphin$ppoly
  )

  mexdolphin$simulated <- sf::st_as_sf(pts)
  mexdolphin$simulated$size <- 1
  #  sf::st_crs(mexdolphin$simulated) <- fmesher::fm_CRS(target_crs)

  #### Depth covariate #####
  depth <- sf::st_as_sf(
    mexdolphins$preddata[, c("x", "y", "depth"), drop = FALSE],
    coords = c("x", "y"),
    crs = fmesher::fm_crs(data_crs)
  )
  mexdolphin$depth <- fmesher::fm_transform(depth, target_crs)

  # return
  mexdolphin
}

# mexdolphin_sf <- import_mexdolphin_sf()
# use_data(mexdolphin_sf, compress = "xz", overwrite = TRUE)
