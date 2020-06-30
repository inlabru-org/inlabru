sp_get_crs <- function(...) {
  INLA::inla.sp_get_crs(...)
}

crs_is_null <- function(crs) {
  if (is.null(crs)) {
    TRUE
  } else if (INLA::inla.has_PROJ6()) {
    wkt <- INLA::inla.crs_get_wkt(crs)
    is.null(wkt)
  } else {
    is.na(crs)
  }
}

ensure_crs <- function(crs) {
  if (crs_is_null(crs)) {
    crs <- CRS(NA_character_)
  }
  crs
}