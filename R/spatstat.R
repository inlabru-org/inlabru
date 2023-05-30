# Internal checker for spatstat packages
# Adapted from Adrian Baddeley and Ege Rubak
check_spatstat <- function(pkg = "spatstat.geom") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    caller_name <- fm_caller_name(1)
    stop(paste0(
      "package '",
      pkg,
      "' required",
      if (identical(caller_name, "")) {
        ""
      } else {
        paste0(" by '", caller_name, "'")
      },
      "; please install it (or the full spatstat package) first"
    ))
  } else {
    spst_ver <- try(utils::packageVersion("spatstat"), silent = TRUE)
    if (!inherits(spst_ver, "try-error") && spst_ver < "2.0-0") {
      warning(paste0(
        "You have an old version of 'spatstat' installed which is incompatible with '",
        pkg,
        "'. Please update 'spatstat' (or uninstall it)."
      ))
      return(FALSE)
    }
  }
  return(TRUE)
}

#' Convert SpatialPoints and boundary polygon to spatstat ppp object
#'
#' Spatstat point pattern objects consist of points and an observation windows. This
#' function uses a SpatialPoints object and a SpatialPolygon object to generate the points
#' and the window. Lastly, the ppp() function is called to create the `ppp` object.
#'
#' @aliases spatial.to.ppp
#' @export
#' @param points A `SpatialPoints[DataFrame]` object describing the point pattern.
#' @param samplers A `SpatialPolygons[DataFrame]` object describing the observation window.
#' @return A spatstat `spatstat` `ppp` object
#'
#' @examples
#' \donttest{
#' if (require("spatstat.geom") &&
#'   bru_safe_sp()) {
#'   # Load Gorilla data
#'
#'   data("gorillas", package = "inlabru")
#'
#'   # Use nest locations and survey boundary to create a spatstat ppp object
#'
#'   gp <- spatial.to.ppp(gorillas$nests, gorillas$boundary)
#'   class(gp)
#'
#'   # Plot it
#'
#'   plot(gp)
#' }
#' }
#'
spatial.to.ppp <- function(points, samplers) {
  check_spatstat("spatstat.geom")

  bnd <- samplers@polygons[[1]]@Polygons[[1]]@coords
  bnd <- bnd[1:(nrow(bnd) - 1), ]
  gp <- spatstat.geom::ppp(
    x = coordinates(points)[, 1],
    y = coordinates(points)[, 2],
    window = spatstat.geom::owin(poly = list(x = rev(bnd[, 1]), y = rev(bnd[, 2])))
  )
}
