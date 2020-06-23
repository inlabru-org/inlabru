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
  cmp <- coordinates ~ mySmooth(map = coordinates, model = matern) + Intercept

  fit <- lgcp(cmp, gorillas$nests,
    samplers = gorillas$boundary,
    options = list(control.inla = list(int.strategy = "eb"))
  )

  list(
    gorillas = gorillas,
    matern = matern,
    cmp = cmp,
    fit = fit
  )
}

mrsea_rebuild_CRS <-function(x) {
  if (inla.has_PROJ6()) {
    x$points <- rebuild_CRS(x$points)
    x$samplers <- rebuild_CRS(x$samplers)
    x$mesh$crs <- rebuild_CRS(x$mesh$crs)
    x$boundary <- rebuild_CRS(x$boundary)
    x$covar <- rebuild_CRS(x$covar)
  }
  x
}

