#' @name mexdolphin
#' @title Pan-tropical spotted dolphins in the Gulf of Mexico
#' @docType data
#' @description This a version of the `mexdolphins` dataset from the package `dsm`, reformatted
#' as point process data for use with `inlabru`. The data are from a combination of several NOAA
#' shipboard surveys conducted on pan-tropical spotted dolphins in the Gulf of Mexico. 47 observations
#' of groups of dolphins wre detected. The group size was recorded, as well as the Beaufort sea state at
#' the time of the observation. Transect width is 16 km, i.e. maximal detection
#' distance 8 km (transect half-width 8 km).
#'
#' @usage data(mexdolphin)
#'
#' @format A list of objects:
#'  \describe{
#'    \item{`points`:}{ A `SpatialPointsDataFrame` object containing the locations of
#'    detected dolphin groups, with their size as an attribute.}
#'    \item{`samplers`:}{ A `SpatialLinesDataFrame` object containing the transect lines
#'    that were surveyed.}
#'    \item{`mesh`:}{ An `inla.mesh` object containing a Delaunay triangulation
#'    mesh (a type of discretization of continuous space) covering the survey region.}
#'    \item{`ppoly`:}{ An `SpatialPolygonsDataFrame` object defining the boundary of the
#'    survey region.}
#'    \item{`simulated`:}{ A `SpatialPointsDataFrame` object containing the locations of a
#'    *simulated* population of dolphin groups. The population was simulated from a 'code{inlabru}
#'    model fitted to the actual survey data. Note that the simulated data do not have any associated
#'    size information.}
#'  }
#' @source
#' Library `dsm`.
#'
#' @references
#' Halpin, P.N., A.J. Read, E. Fujioka, B.D. Best, B. Donnelly, L.J. Hazen, C. Kot, K. Urian,
#' E. LaBrecque, A. Dimatteo, J. Cleary, C. Good, L.B. Crowder, and K.D. Hyrenbach. 2009.
#' OBIS-SEAMAP: The world data center for marine mammal, sea bird, and sea turtle distributions.
#' Oceanography 22(2):104-115
#'
#' NOAA Southeast Fisheries Science Center. 1996. Report of a Cetacean Survey of Oceanic and
#' Selected Continental Shelf Waters of the Northern Gulf of Mexico aboard NOAA Ship Oregon II
#' (Cruise 220)
#'
#' @examples
#' \donttest{
#' if (require("ggplot2", quietly = TRUE)) {
#'   data(mexdolphin, package = "inlabru")
#'   ggplot() +
#'     gg(mexdolphin$mesh) +
#'     gg(mexdolphin$ppoly, color = "blue") +
#'     gg(mexdolphin$samplers) +
#'     gg(mexdolphin$points, aes(size = size), color = "red") +
#'     coord_equal()
#'
#'   ggplot() +
#'     gg(mexdolphin$mesh, col = mexdolphin$lambda, mask = mexdolphin$ppoly) +
#'     coord_equal()
#' }
#' }
#' \dontrun{
#' if (requireNamespace("ggmap", quietly = TRUE) &&
#'   require("ggplot2", quietly = TRUE)) {
#'   gmap(mexdolphin$depth) +
#'     gm(mexdolphin$ppoly, color = "blue") +
#'     gm(mexdolphin$samplers) +
#'     gm(mexdolphin$points, aes(size = size), color = "red")
#'
#'   gmap(mexdolphin$depth) +
#'     gm(mexdolphin$depth, aes(col = depth)) +
#'     gm(mexdolphin$ppoly)
#' }
#' }
NULL

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
  mexdolphin$points <- sp::spTransform(mexdolphin$points, CRSobj = target_crs)
  mexdolphin$samplers <- sp::spTransform(mexdolphin$samplers, CRSobj = target_crs)
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
  mexdolphin$mesh$crs <- target_crs

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
  mexdolphin$depth <- sp::spTransform(depth, target_crs)

  # return
  mexdolphin
}

# use_data(mexdolphin, compress = "xz")
