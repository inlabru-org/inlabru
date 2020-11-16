#' Disable PROJ4/6 warnings
#' 
#' To be used within package tests. Restores state on exit.
#' 
#' @param envir environment for exit handlers
#' @export
#' @keywords internal
disable_PROJ6_warnings <- function(envir = parent.frame()) {
  withr::local_options(
    list("rgdal_show_exportToProj4_warnings" = "none"),
    envir
  )
  requireNamespace("rgdal", quietly = TRUE)
  if (fm_has_PROJ6()) {
    old1 <- rgdal::get_rgdal_show_exportToProj4_warnings()
    rgdal::set_rgdal_show_exportToProj4_warnings(FALSE)
    withr::defer(
      rgdal::set_rgdal_show_exportToProj4_warnings(old1),
      envir = envir
    )

    old2 <- rgdal::get_thin_PROJ6_warnings()
    rgdal::set_thin_PROJ6_warnings(TRUE)
    withr::defer(
      rgdal::set_thin_PROJ6_warnings(old2),
      envir = envir
    )
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
  set.seed(123)
  cbind(
    basic_intercept_testdata(),
    data.frame(x1 = rnorm(100))
  )
}






mrsea_rebuild_CRS <- function(x, use_km = FALSE) {
  if (fm_has_PROJ6()) {
    x$points <- sp::rebuild_CRS(x$points)
    x$samplers <- sp::rebuild_CRS(x$samplers)
    x$mesh$crs <- sp::rebuild_CRS(x$mesh$crs)
    x$boundary <- sp::rebuild_CRS(x$boundary)
    x$covar <- sp::rebuild_CRS(x$covar)
  }
  if (use_km) {
    # The estimation is numerically unreliable when the spatial
    # domain is represented in metres, and has been seen to produce
    # different results on different systems (e.g. Travis CI).
    # Transform m to km:
    crs_km <- fm_crs_set_lengthunit(x$mesh$crs, "km")
    x$mesh <- fm_spTransform(x$mesh, crs_km)
    x$samplers <- sp::spTransform(x$samplers, crs_km)
    x$points <- sp::spTransform(x$points, crs_km)
    x$boundary <- sp::spTransform(x$boundary, crs_km)
    x$covar <- sp::spTransform(x$covar, crs_km)
  }
  x
}
