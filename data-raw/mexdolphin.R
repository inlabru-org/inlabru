source(here::here("data-raw", "dsmdata.tools.R"))
source(here::here("data-raw", "dsmdata.R"))

#' Mexdolphin data import
#'
#'
#' Load `mexdolphins` survey data from `dsm` package and convert to spatial formats defined by the [sp] package.
#'
#' @aliases import.mexdolphin
#' @keywords internal
#' @return The [mexdolphin] data set
#' @examples
#' \dontrun{
#' mexdolphin <- import.mexdolphin()
#' }
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'


import.mexdolphin <- function() {
  envir <- new.env()
  data("mexdolphins", package = "dsm", envir = envir)
  mexdolphins <- as.list(envir)
  data.p4s <- "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  data_crs <- sp::CRS(data.p4s)

  dset <- import.dsmdata(mexdolphins, covar.col = 8)
  mexdolphin <- as.spatial.dsdata(dset, cnames = c("x", "y"), crs = data_crs)

  # Target CRS
  target.p4s <- "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs +towgs84=0,0,0"
  target_crs <- sp::CRS(target.p4s)

  # Units to km
  mexdolphin$points <- fm_transform(mexdolphin$points, crs = target_crs)
  mexdolphin$samplers <- fm_transform(mexdolphin$samplers, crs = target_crs)
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
  mexdolphin$mesh$loc <- mexdolphin$mesh$loc / 1000

  # Set mesh crs
  mexdolphin$mesh$crs <- fm_crs(target_crs)

  coordnames(mexdolphin$samplers) <- coordnames(mexdolphin$points)

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
  polyloc <- as.data.frame(
    mexdolphin$mesh$loc[mexdolphin$mesh$segm$int$idx[, 1], c(1, 2)]
  )
  colnames(polyloc) <- c("x", "y")
  po <- sp::Polygon(polyloc, hole = FALSE)
  pos <- sp::Polygons(list(po), ID = "c")
  predpoly <- sp::SpatialPolygons(list(pos), proj4string = target_crs)
  df <- data.frame(weight = 1)
  rownames(df) <- "c"
  predpolyd <- sp::SpatialPolygonsDataFrame(predpoly, data = df)
  # plot(predpolyd)
  mexdolphin$ppoly <- predpolyd

  ##### Simulate a whole population #####
  distance <- seq(0, 8, length.out = 20)
  matern <- INLA::inla.spde2.pcmatern(
    mexdolphin$mesh,
    prior.range = c(50, 0.01),
    prior.sigma = c(2, 0.01)
  )
  cmps <-
    coordinates + distance ~
    Intercept(1) +
    df.lsigma(1) +
    spat(coordinates, model = matern)
  pred <- ~ log(1 - exp(-(distance / exp(df.lsigma))^-1)) + spat + Intercept + log(2)
  r <- lgcp(
    data = mexdolphin$points,
    samplers = cbind(mexdolphin$samplers, data.frame(weight = 1)),
    components = cmps,
    formula = pred,
    domain = list(
      coordinates = mexdolphin$mesh,
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
  coordnames(pts) <- c("x", "y", "z")

  mexdolphin$simulated <- sp::SpatialPointsDataFrame(
    coords = pts,
    data = data.frame(size = rep(1, length(pts))),
    proj4string = target_crs
  )

  #### Depth covariate #####
  depth <- sp::SpatialPointsDataFrame(
    coords = mexdolphins$preddata[, c("x", "y")],
    data = mexdolphins$preddata[, "depth", drop = FALSE],
    proj4string = data_crs
  )
  mexdolphin$depth <- fm_transform(depth, target_crs)

  # return
  mexdolphin
}

#' @describeIn import.mexdolphin Import mexdolphin data as `sf`
import.mexdolphin.sf <- function() {
  envir <- new.env()
  data("mexdolphins", package = "dsm", envir = envir)
  mexdolphins <- as.list(envir)
  data.p4s <- "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  data_crs <- fm_crs(data.p4s)

  dset <- import.dsmdata(mexdolphins, covar.col = 8)
  dset$mesh$crs <- fm_CRS(data_crs)
  mexdolphin <- as.spatial.dsdata(dset, cnames = c("x", "y"), crs = fm_CRS(data_crs))
  mexdolphin$points <- sf::st_as_sf(mexdolphin$points)
  mexdolphin$samplers <- sf::st_as_sf(mexdolphin$samplers)

  # Target CRS
  target.p4s <- "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs +towgs84=0,0,0"
  target_crs <- fm_crs(target.p4s)

  # Units to km
  mexdolphin$points <- fm_transform(mexdolphin$points, crs = target_crs)
  mexdolphin$samplers <- fm_transform(mexdolphin$samplers, crs = target_crs)
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
  mexdolphin$mesh <- fm_transform(mexdolphin$mesh, crs = target_crs)
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
  polyloc <- as.data.frame(
    mexdolphin$mesh$loc[mexdolphin$mesh$segm$int$idx[, 1], c(1, 2)]
  )
  colnames(polyloc) <- c("x", "y")
  po <- sp::Polygon(polyloc, hole = FALSE)
  pos <- sp::Polygons(list(po), ID = "c")
  predpoly <- sp::SpatialPolygons(list(pos), proj4string = fm_CRS(target_crs))
  df <- data.frame(weight = 1)
  rownames(df) <- "c"
  predpolyd <- sp::SpatialPolygonsDataFrame(predpoly, data = df)
  # plot(predpolyd)
  mexdolphin$ppoly <- sf::st_as_sf(predpolyd)

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
  pred <- ~ log(1 - exp(-(distance / exp(df.lsigma))^-1)) + spat + Intercept + log(2)
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
  #  sf::st_crs(mexdolphin$simulated) <- fm_CRS(target_crs)

  #### Depth covariate #####
  depth <- sf::st_as_sf(
    mexdolphins$preddata[, c("x", "y", "depth"), drop = FALSE],
    coords = c("x", "y"),
    crs = fm_crs(data_crs)
  )
  mexdolphin$depth <- fm_transform(depth, target_crs)

  # return
  mexdolphin
}

# use_data(mexdolphin, compress = "xz", overwrite = TRUE)
# use_data(mexdolphin_sf, compress = "xz", overwrite = TRUE)
