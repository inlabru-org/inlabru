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
  }
  gorillas
}


gorillas_lgcp_2d_testdata <- function() {
  data(gorillas, package = "inlabru")
  gorillas <- gorillas_update_CRS(gorillas)
  
  matern <- inla.spde2.pcmatern(gorillas$mesh, prior.sigma = c(0.1, 0.01), prior.range = c(5, 0.01))
  cmp <- coordinates ~ mySmooth(main = coordinates, model = matern) + Intercept(1)

  fit <- lgcp(cmp, gorillas$nests,
    samplers = gorillas$boundary,
    domain = list(coordinates = gorillas$mesh),
    options = list(control.inla = list(int.strategy = "eb",
                                       h = 0.005),
                   num.threads = 1)
  )

  list(
    gorillas = gorillas,
    matern = matern,
    cmp = cmp,
    fit = fit
  )
}


mrsea_rebuild_CRS <-function(x) {
  disable_PROJ6_warnings()
  if (inla.has_PROJ6()) {
    x$points <- rebuild_CRS(x$points)
    x$samplers <- rebuild_CRS(x$samplers)
    x$mesh$crs <- rebuild_CRS(x$mesh$crs)
    x$boundary <- rebuild_CRS(x$boundary)
    x$covar <- rebuild_CRS(x$covar)
  }
  x
}

