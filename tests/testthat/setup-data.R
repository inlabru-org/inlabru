disable_PROJ6_warnings <- function() {
  library(rgdal)
  if (inla.has_PROJ6()) {
    rgdal::set_rgdal_show_exportToProj4_warnings(FALSE)
    rgdal::set_thin_PROJ6_warnings(TRUE)
  }
}


basic_intercept_testdata <- function() {
  set.seed(123)
  data.frame(
    Intercept = 1,
    y = rnorm(100)
  )
}

basic_fixed_effect_testdata <- function() {
  cbind(
    basic_intercept_testdata(),
    data.frame(x1 = rnorm(100))
  )
}


gorillas_update_CRS <- function(gorillas) {
  if (inla.has_PROJ6()) {
    gorillas$nests <- rebuild_CRS(gorillas$nests)
    gorillas$boundary <- rebuild_CRS(gorillas$boundary)
    gorillas$mesh$crs <- rebuild_CRS(gorillas$mesh$crs)
    for (name in names(gorillas$gcov)) {
      gorillas$gcov[[name]] <- rebuild_CRS(gorillas$gcov[[name]])
    }
    for (name in names(gorillas$plotsample)) {
      gorillas$plotsample[[name]] <- rebuild_CRS(gorillas$plotsample[[name]])
    }
    crs <- sp::CRS(SRS_string = wkt(gorillas$mesh$crs))
    gorillas$nests <- spTransform(gorillas$nests, crs)
    gorillas$boundary <- spTransform(gorillas$boundary, crs)
    gorillas$mesh <- inla.spTransform(gorillas$mesh, crs)
    for (name in names(gorillas$gcov)) {
      gorillas$gcov[[name]] <-
        SpatialPixelsDataFrame(coordinates(gorillas$gcov[[name]]),
                               data = as.data.frame(gorillas$gcov[[name]]),
                               proj4string = crs)
    }
    for (name in names(gorillas$plotsample)) {
      gorillas$plotsample[[name]] <- spTransform(gorillas$plotsample[[name]], crs)
    }
  }
  gorillas
}


gorillas_lgcp_2d_testdata <- function() {
  data(gorillas, package = "inlabru")
  gorillas <- gorillas_update_CRS(gorillas)
  
  matern <- inla.spde2.pcmatern(gorillas$mesh, prior.sigma = c(0.1, 0.01), prior.range = c(5, 0.01))
  cmp <- coordinates ~ mySmooth(map = coordinates, model = matern) + Intercept

  fit <- lgcp(cmp, gorillas$nests,
    samplers = gorillas$boundary,
    domain = list(coordinates = gorillas$mesh),
    options = list(control.inla = list(int.strategy = "eb",
                                       h = 0.005),
                   num.threads = "1:1")
  )

  list(
    gorillas = gorillas,
    matern = matern,
    cmp = cmp,
    fit = fit
  )
}


mrsea_rebuild_CRS <-function(x, use_km = FALSE) {
  if (inla.has_PROJ6()) {
    x$points <- rebuild_CRS(x$points)
    x$samplers <- rebuild_CRS(x$samplers)
    x$mesh$crs <- rebuild_CRS(x$mesh$crs)
    x$boundary <- rebuild_CRS(x$boundary)
    x$covar <- rebuild_CRS(x$covar)
  }
  if (use_km) {
    # The estimation is numerically unreliable when the spatial
    # domain is represented in metres, and has been seen to produce
    # different results on different systems (e.g. Travis CI).
    # Transform m to km:
    crs_km <- inla.crs_set_lengthunit(x$mesh$crs, "km")
    x$mesh <- inla.spTransform(x$mesh, crs_km)
    x$samplers <- sp::spTransform(x$samplers, crs_km)
    x$points <- sp::spTransform(x$points, crs_km)
    x$boundary <- sp::spTransform(x$boundary, crs_km)
    x$covar <- sp::spTransform(x$covar, crs_km)
  }  
  x
}

