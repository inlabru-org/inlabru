#' @importFrom sp coordinates proj4string `proj4string<-`

#' @describeIn inlabru-deprecated Detect whether PROJ6 is available
#'
#' @export

fm_has_PROJ6 <- function() {
  lifecycle::deprecate_warn("2.8.0",
    "fm_has_PROJ6()",
    details = c(
      i = "Since inlabru 2.7.1, fm_has_PROJ6() always returns TRUE",
      i = "rgdal/PROJ4 is no longer supported."
    )
  )
  TRUE
}

#' @describeIn inlabru-deprecated `fm_not_for_PROJ6` is called to warn about using old PROJ4
#' features even though PROJ6 is available

fm_not_for_PROJ6 <- function(fun = NULL) {
  lifecycle::deprecate_stop("2.8.0",
    "fm_not_for_PROJ6()",
    details = c(x = "rgdal/PROJ4 is no longer supported.")
  )
}

#' @describeIn inlabru-deprecated `fm_not_for_PROJ4` is called to give an error when
#' calling methods that are only available for PROJ6

fm_not_for_PROJ4 <- function(fun = NULL) {
  lifecycle::deprecate_stop("2.8.0",
    "fm_not_for_PROJ4()",
    details = c(x = "rgdal/PROJ4 is no longer supported.")
  )
}

#' @describeIn inlabru-deprecated Called to warn about falling back
#' to using old PROJ4 methods when a PROJ6 method hasn't been implemented

fm_fallback_PROJ6 <- function(fun = NULL) {
  lifecycle::deprecate_stop("2.8.0",
    "fm_not_for_PROJ4()",
    details = c(x = "rgdal/PROJ4 requested by PROJ4 is no longer supported.")
  )
}


#' @param fun The name of the function that requires PROJ6. Default: NULL,
#' which uses the name of the calling function.
#' @describeIn inlabru-deprecated Called to give an error when PROJ6
#' is required but not available

fm_requires_PROJ6 <- function(fun = NULL) {
  lifecycle::deprecate_stop("2.8.0",
    "fm_requires_PROJ6()",
    details = c(x = "rgdal/PROJ4 is no longer supported.")
  )
}


#' @rdname fm_as
#' @export
fm_as_sp_crs <- function(x, ...) {
  fm_CRS(x, ...)
}




#' @describeIn inlabru-deprecated Wrapper for CRS(projargs) (PROJ4) and CRS(wkt) for
#' `sp::Spatial` objects.
#' @param x A `sp::Spatial` object
#' @return A `CRS` object, or NULL if no valid CRS identified
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @details This function is a convenience method to workaround PROJ4/PROJ6
#' differences, and the lack of a crs extraction method for Spatial objects.
#' For newer code, use [fm_crs()] instead, that returns `crs` objects,
#' and use [fm_as_sp_crs()] to convert to old style `sp::CRS` objects.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   s <- sp::SpatialPoints(matrix(1:6, 3, 2), proj4string = fm_CRS("sphere"))
#'   fm_sp_get_crs(s)
#' }
#' }
#' @export

fm_sp_get_crs <- function(x) {
  lifecycle::deprecate_warn("2.8.0", "fm_sp_get_crs()", "fm_CRS()")
  fm_CRS(x)
}


#' @describeIn fm_crs Check if an object is or has `NULL` or `NA` CRS information
#' @export
fm_crs_is_null <- function(x) {
  if (is.null(x)) {
    return(TRUE)
  }
  is.na(fm_crs(x))
}





#' @title Handling CRS/WKT
#' @description Get and set CRS object or WKT string properties.
#' @export
#' @name fm_crs_wkt
#' @rdname fm_crs_wkt

fm_wkt_is_geocent <- function(wkt) {
  if (is.null(wkt) || identical(wkt, "") || is.na(wkt)) {
    return(FALSE)
  }
  # See https://proceedings.esri.com/library/userconf/proc17/tech-workshops/tw_2588-212.pdf
  geo_crs_items <- c(
    "GEODCRS", "GEOGCRS",
    "BASEGEODCRS", "BASEGEOGCRS"
  )
  wt <- fm_wkt_as_wkt_tree(wkt)
  if (identical(wt[["label"]], "BOUNDCRS")) {
    wt <- fm_wkt_tree_get_item(wt, "SOURCECRS")
    wt <- fm_wkt_tree_get_item(wt, c("PROJCRS", geo_crs_items))
  }
  if (identical(wt[["label"]], "PROJCRS")) {
    wt <- fm_wkt_tree_get_item(wt, geo_crs_items)
  }
  if (!(wt[["label"]] %in% geo_crs_items)) {
    return(FALSE)
  }
  cs <- fm_wkt_tree_get_item(wt, "CS")
  if (is.null(cs)) {
    return(FALSE)
  }
  cart <- ((cs[["params"]][[1]] == "Cartesian") &&
    (cs[["params"]][[2]] == "3"))
  if (!cart) {
    return(FALSE)
  }
  axis_names <- c('"(X)"', '"(Y)"', '"(Z)"')
  axis_types <- c("geocentricX", "geocentricY", "geocentricZ")
  for (k in seq_len(3)) {
    axis <- fm_wkt_tree_get_item(wt, "AXIS", k)
    if (!((axis[["params"]][[1]] == axis_names[k]) &&
      (axis[["params"]][[2]] == axis_types[k]))) {
      return(FALSE)
    }
  }
  TRUE
}

#' @export
#' @rdname fm_crs_wkt

fm_crs_is_geocent <- function(crs) {
  wkt <- fm_wkt(crs)
  result <- fm_wkt_is_geocent(wkt)
  result
}


# ellipsoid_radius ----

#' @rdname fm_crs_wkt
#' @export

fm_wkt_get_ellipsoid_radius <- function(wkt) {
  geo_crs_items <- c(
    "GEODCRS", "GEOGCRS",
    "BASEGEODCRS", "BASEGEOGCRS"
  )
  wt <- fm_wkt_as_wkt_tree(wkt)

  if (identical(wt[["label"]], "BOUNDCRS")) {
    wt <- fm_wkt_tree_get_item(wt, "SOURCECRS")
    wt <- fm_wkt_tree_get_item(wt, c("PROJCRS", geo_crs_items))
  }
  if (identical(wt[["label"]], "PROJCRS")) {
    wt <- fm_wkt_tree_get_item(wt, geo_crs_items)
  }
  if (is.null(wt) || !(wt[["label"]] %in% geo_crs_items)) {
    stop("Ellipsoid settings not found (no geo crs)")
  }

  datum <- fm_wkt_tree_get_item(wt, c("DATUM", "ENSEMBLE"))
  if (is.null(datum)) {
    stop("Ellipsoid settings not found (no datum/ensemble)")
  }
  ellipsoid <- fm_wkt_tree_get_item(datum, "ELLIPSOID")
  if (is.null(ellipsoid)) {
    stop("Ellipsoid settings not found (no ellipsoid)")
  }
  as.numeric(ellipsoid[["params"]][[2]])
}

#' @rdname fm_crs_wkt
#' @export

fm_crs_get_ellipsoid_radius <- function(crs) {
  fm_wkt_get_ellipsoid_radius(fm_wkt(crs))
}

#' @param x crs object to extract value from or assign values in
#' @rdname fm_crs_wkt
#' @export

fm_ellipsoid_radius <- function(x) {
  UseMethod("fm_ellipsoid_radius")
}

#' @rdname fm_crs_wkt
#' @export

fm_ellipsoid_radius.default <- function(x) {
  fm_wkt_get_ellipsoid_radius(fm_wkt(x))
}

#' @rdname fm_crs_wkt
#' @export

fm_ellipsoid_radius.character <- function(x) {
  fm_wkt_get_ellipsoid_radius(x)
}


#' @param radius numeric; The new radius value
#' @rdname fm_crs_wkt
#' @export

fm_wkt_set_ellipsoid_radius <- function(wkt, radius) {
  geo_crs_items <- c(
    "GEODCRS", "GEOGCRS",
    "BASEGEODCRS", "BASEGEOGCRS"
  )

  set_radius <- function(wt) {
    if (is.null(wt)) {
      stop("Ellipsoid settings not found")
    } else if (wt[["label"]] %in% geo_crs_items) {
      datum <- fm_wkt_tree_get_item(wt, c("DATUM", "ENSEMBLE"))
      null_datum <- is.null(datum)
      if (null_datum) {
        stop("Ellipsoid settings not found (no datum/ensemble)")
      }
      ellipsoid <- fm_wkt_tree_get_item(datum, "ELLIPSOID")
      if (is.null(ellipsoid)) {
        stop("Ellipsoid settings not found")
      }
      ellipsoid[["params"]][[2]] <- as.character(radius)
      datum <- fm_wkt_tree_set_item(datum, ellipsoid)
      wt <- fm_wkt_tree_set_item(wt, datum)
    } else if (wt[["label"]] %in% c("BOUNDCRS", "SOURCECRS", "PROJCRS")) {
      wt_sub <- fm_wkt_tree_get_item(
        wt,
        c(
          "BOUNDCRS", "SOURCECRS", "PROJCRS",
          geo_crs_items
        )
      )
      if (is.null(wt_sub)) {
        stop("Ellipsoid settings not found")
      }
      wt_sub_new <- set_radius(wt_sub)
      wt <- fm_wkt_tree_set_item(wt, wt_sub_new)
    } else {
      stop("Ellipsoid settings not found")
    }
    wt
  }

  wt <- fm_wkt_as_wkt_tree(wkt)
  wt <- set_radius(wt)
  fm_wkt_tree_as_wkt(wt)
}

#' @rdname fm_crs_wkt
#' @param value Value to assign
#' @export

`fm_ellipsoid_radius<-` <- function(x, value) {
  UseMethod("fm_ellipsoid_radius<-")
}

#' @rdname fm_crs_wkt
#' @export
`fm_ellipsoid_radius<-.character` <- function(x, value) {
  x <- fm_wkt_set_ellipsoid_radius(x, value)
  invisible(x)
}

#' @rdname fm_crs_wkt
#' @export
`fm_ellipsoid_radius<-.CRS` <- function(x, value) {
  crs <- fm_crs(x)
  fm_ellipsoid_radius(crs) <- value
  new_crs <- fm_CRS(crs)

  new_crs
}

#' @rdname fm_crs_wkt
#' @export
`fm_ellipsoid_radius<-.inla.CRS` <- function(x, value) {
  crs <- fm_crs(x)
  fm_ellipsoid_radius(crs) <- value
  new_crs <- fm_CRS(crs)

  new_crs
}

#' @rdname fm_crs_wkt
#' @export
`fm_ellipsoid_radius<-.crs` <- function(x, value) {
  wkt <- fm_wkt(x)
  fm_ellipsoid_radius(wkt) <- value
  new_crs <- fm_crs(wkt)
  new_crs
}

#' @rdname fm_crs_wkt
#' @export
`fm_ellipsoid_radius<-.fm_crs` <- function(x, value) {
  wkt <- fm_wkt(x)
  fm_ellipsoid_radius(wkt) <- value
  x$crs <- fm_crs(wkt)
  x
}


#' @rdname fm_crs_wkt
#' @export

fm_crs_set_ellipsoid_radius <- function(crs, radius) {
  fm_ellipsoid_radius(crs) <- radius
  crs
}

# Length unit ----



#' @param crs An `sf::crs`, `sp::CRS`, `fm_crs` or `inla.CRS` object
#' @param wkt A WKT2 character string
#' @param unit character, name of a unit. Supported names are
#' "metre", "kilometre", and the aliases "meter", "m", International metre",
#' "kilometer", and "km", as defined by `fm_wkt_unit_params` or the
#' `params` argument. (For legacy PROJ4 use, only "m" and "km" are
#' supported)
#' @param params Length unit definitions, in the list format produced by
#' `fm_wkt_unit_params()`, Default: NULL, which invokes
#' `fm_wkt_unit_params()`
#' @return For `fm_wkt_unit_params`, a
#' list of named unit definitions
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @examples
#' \dontrun{
#' c1 <- fm_crs("globe")
#' fm_crs_get_lengthunit(c1)
#' c2 <- fm_crs_set_lengthunit(c1, "m")
#' fm_crs_get_lengthunit(c2)
#' }
#' @export
#' @seealso [fm_crs()]
#' @aliases fm_crs_wkt
#' @rdname fm_crs_wkt

fm_wkt_unit_params <- function() {
  params <- list(
    "metre" =
      list(
        '"metre"',
        "1",
        list(
          label = "ID",
          params = list('"EPSG"', "9001")
        )
      ),
    "kilometre" =
      list(
        '"kilometre"',
        "1000",
        list(
          label = "ID",
          params = list('"EPSG"', "9036")
        )
      )
  )
  params[["meter"]] <- params[["metre"]]
  params[["m"]] <- params[["metre"]]
  params[["International metre"]] <- params[["metre"]]
  params[["kilometer"]] <- params[["kilometre"]]
  params[["km"]] <- params[["kilometre"]]
  params
}

#' @export
#' @rdname fm_crs_wkt
#' @return For `fm_wkt_get_lengthunit`, a
#' list of length units used in the wkt string, excluding the ellipsoid radius
#' unit.

fm_wkt_get_lengthunit <- function(wkt) {
  extract <- function(wt) {
    # 1. Recursively find LENGTHUNIT, except within ELLIPSOID
    # 2. Return unit

    if (wt[["label"]] == "LENGTHUNIT") {
      result <- list(wt[["params"]])
    } else if (wt[["label"]] != "ELLIPSOID") {
      result <- list()
      for (k in seq_along(wt$param)) {
        if (is.list(wt[["params"]][[k]])) {
          result <- c(result, extract(wt[["params"]][[k]]))
        }
      }
    } else {
      result <- list()
    }
    result
  }

  wt <- fm_wkt_as_wkt_tree(wkt)
  params <- unique(extract(wt))
  names(params) <-
    vapply(
      params,
      function(x) {
        gsub('"', "", x[[1]])
      },
      ""
    )
  params
}

#' @export
#' @rdname fm_crs_wkt
#' @return For `fm_wkt_set_lengthunit`, a
#' WKT2 string with altered length units.
#' Note that the length unit for the ellipsoid radius is unchanged.

fm_wkt_set_lengthunit <- function(wkt, unit, params = NULL) {
  convert <- function(wt, unit) {
    # 1. Recursively find LENGTHUNIT, except within ELLIPSOID
    # 2. Change unit

    if (wt[["label"]] == "LENGTHUNIT") {
      wt[["params"]] <- unit
    } else if (wt[["label"]] != "ELLIPSOID") {
      if ((wt[["label"]] == "PARAMETER") &&
        (wt[["params"]][[1]] %in% c('"False easting"', '"False northing"'))) {
        orig_unit <- (wt[["params"]][[3]][["params"]][[1]])
        new_unit <- (unit[[1]])
        if (orig_unit != new_unit) {
          if (orig_unit == '"metre"' && new_unit == '"kilometre"') {
            wt[["params"]][[2]] <-
              as.character(as.numeric(wt[["params"]][[2]]) / 1000)
          } else if (orig_unit == '"kilometre"' && new_unit == '"metre"') {
            wt[["params"]][[2]] <-
              as.character(as.numeric(wt[["params"]][[2]]) * 1000)
          } else {
            warning("False easting/northing could not be properly converted.")
          }
        }
      }
      for (k in seq_along(wt$param)) {
        if (is.list(wt[["params"]][[k]])) {
          wt[["params"]][[k]] <- convert(wt[["params"]][[k]], unit)
        }
      }
    }
    wt
  }

  if (is.null(params)) {
    params <- fm_wkt_unit_params()
  }
  if (!(unit %in% names(params))) {
    warning(paste0(
      "'fm_wkt_set_lengthunit' unit conversion to '",
      unit,
      "' not supported. Unit left unchanged."
    ))
    return(wkt)
  }

  wt <- fm_wkt_as_wkt_tree(wkt)
  wt <- convert(wt, params[[unit]])
  fm_wkt_tree_as_wkt(wt)
}


#' @return For `fm_crs_get_lengthunit`, a
#' list of length units used in the wkt string, excluding the ellipsoid radius
#' unit. (For legacy PROJ4 code, the raw units from the proj4string are
#' returned, if present.)
#' @export
#' @rdname fm_crs_wkt

fm_crs_get_lengthunit <- function(crs) {
  fm_wkt_get_lengthunit(fm_wkt(crs))
}


#' @rdname fm_crs_wkt
#' @export

fm_crs_set_lengthunit <- function(crs, unit) {
  fm_length_unit(crs) <- unit
  crs
}


#' @rdname fm_crs_wkt
#' @export

fm_length_unit <- function(x) {
  UseMethod("fm_length_unit")
}

#' @rdname fm_crs_wkt
#' @export

fm_length_unit.default <- function(x) {
  fm_wkt_get_lengthunit(fm_wkt(x))
}

#' @rdname fm_crs_wkt
#' @export

fm_length_unit.character <- function(x) {
  fm_wkt_get_lengthunit(x)
}




#' @return For `fm_length_unit<-`, a crs object with
#' altered length units.
#' Note that the length unit for the ellipsoid radius is unchanged.
#' @rdname fm_crs_wkt
#' @export

`fm_length_unit<-` <- function(x, value) {
  UseMethod("fm_length_unit<-")
}

#' @rdname fm_crs_wkt
#' @export
`fm_length_unit<-.character` <- function(x, value) {
  x <- fm_wkt_set_lengthunit(x, value)
  invisible(x)
}

#' @rdname fm_crs_wkt
#' @export
`fm_length_unit<-.CRS` <- function(x, value) {
  crs <- fm_crs(x)
  fm_length_unit(crs) <- value
  new_crs <- fm_CRS(crs)

  new_crs
}

#' @rdname fm_crs_wkt
#' @export
`fm_length_unit<-.inla.CRS` <- function(x, value) {
  crs <- fm_crs(x)
  fm_length_unit(crs) <- value
  new_crs <- fm_CRS(crs)

  new_crs
}

#' @rdname fm_crs_wkt
#' @export
`fm_length_unit<-.crs` <- function(x, value) {
  wkt <- fm_wkt(x)
  fm_length_unit(wkt) <- value
  new_crs <- fm_crs(wkt)
  new_crs
}

#' @rdname fm_crs_wkt
#' @export
`fm_length_unit<-.fm_crs` <- function(x, value) {
  wkt <- fm_wkt(x)
  fm_length_unit(wkt) <- value
  x$crs <- fm_crs(wkt)
  x
}



# fm_crs ----

#' @title Obtain coordinate reference system object
#'
#' @description Obtain an `sf::crs` or `fm_crs` object from a spatial object, or
#' convert crs information to construct a new `sf::crs` object.
#'
#' @param x Object to convert to `crs` or  to extract `crs` information from.
#' If `character`, a string suitable for `sf::st_crs(x)`, or the name of a
#' predefined `wkt` string from ``names(fm_wkt_predef())`.
#' @param \dots Additional arguments passed on the `sf::st_crs()`
#' @param crsonly logical; if TRUE, remove `oblique` information from `fm_crs`
#' objects and return a plain `crs` object instead.
#' @param value Vector of length at most 4 of rotation angles (in degrees)
#' for an oblique projection, all values defaulting to zero. The values
#' indicate (longitude, latitude, orientation, orbit), as explained in the
#' Details section below.
#'
#' @details The first two
#' elements of the `oblique` vector are the (longitude, latitude)
#' coordinates for the oblique centre point. The third value (orientation) is a
#' counterclockwise rotation angle for an observer looking at the centre point
#' from outside the sphere. The fourth value is the quasi-longitude (orbit
#' angle) for a rotation along the oblique observers equator.
#'
#' Simple oblique: `oblique=c(0, 45)`
#'
#' Polar: `oblique=c(0, 90)`
#'
#' Quasi-transversal: `oblique=c(0, 0, 90)`
#'
#' Satellite orbit viewpoint: `oblique=c(lon0-time*v1, 0, orbitangle,
#' orbit0+time*v2)`, where `lon0` is the longitude at which a satellite
#' orbit crosses the equator at `time=0`, when the satellite is at an
#' angle `orbit0` further along in its orbit.  The orbital angle relative
#' to the equatorial plane is `orbitangle`, and `v1` and `v2`
#' are the angular velocities of the planet and the satellite, respectively.
#' Note that "forward" from the satellite's point of view is "to the right" in
#' the projection.
#'
#' When `oblique[2]` or `oblique[3]` are non-zero, the resulting
#' projection is only correct for perfect spheres.
#' @param \dots Additional parameters. Not currently in use.
#' @return Either an `sf::crs` object or an `fm_crs` object,
#' depending on if the coordinate reference system described by the parameters
#' can be expressed with a pure `crs` object or not.
#'
#' @returns A `crs` object ([sf::st_crs()]) or a `fm_crs` object.
#' An S3 `fm_crs` object is a list with elements `crs` and `oblique`.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [sf::st_crs()], [`fm_crs_wkt`]
#' @examples
#' crs1 <- fm_crs("longlat_globe")
#' crs2 <- fm_crs("lambert_globe")
#' crs3 <- fm_crs("mollweide_norm")
#' crs4 <- fm_crs("hammer_globe")
#' crs5 <- fm_crs("sphere")
#' crs6 <- fm_crs("globe")
#' @export
#' @rdname fm_crs
fm_crs <- function(x, ..., crsonly = FALSE) {
  UseMethod("fm_crs")
}

#' @export
#' @rdname fm_crs
fm_crs_oblique <- function(x) {
  stopifnot(inherits(x, c("crs", "fm_crs", "inla.CRS")))
  if (inherits(x, c("fm_crs", "inla.CRS"))) {
    x$oblique
  } else {
    NULL
  }
}

#' @export
#' @rdname fm_crs
`fm_crs_oblique<-` <- function(x, value) {
  stopifnot(inherits(x, c("crs", "fm_crs")))
  if (is.null(value)) {
    if (inherits(x, "fm_crs")) {
      x <- x[["crs"]]
    }
  } else {
    if (!inherits(x, "fm_crs")) {
      stopifnot(is.vector(value))
      if (length(value) < 4) {
        value <- c(value, rep(0, 4 - length(value)))
      }
      x <- list(crs = x, oblique = value)
      class(x) <- "fm_crs"
    }
  }
  x
}

#' @export
#' @rdname fm_crs
print.fm_crs <- function(x, ...) {
  print(x[["crs"]])
  cat(paste0("Oblique: c(", paste0(x[["oblique"]], collapse = ", "), ")\n"))
}



#' @export
#' @rdname fm_crs
fm_crs.default <- function(x, ..., crsonly = FALSE) {
  if (missing(x)) {
    return(sf::NA_crs_)
  }
  sf::st_crs(x, ...)
}

#' @importFrom sf st_crs
#' @exportS3Method sf::st_crs fm_crs
#' @describeIn fm_crs `st_crs(x, ...)` is equivalent to `fm_crs(x, ..., crsonly = TRUE)`
#' when `x` is a `fm_crs` object.
st_crs.fm_crs <- function(x, ...) {
  fm_crs(x, ..., crsonly = TRUE)
}

#' @rawNamespace S3method("$", fm_crs)
#' @describeIn fm_crs For a `fm_crs` object `x`, `x$name` calls the accessor method for the
#' `crs` object inside it. If `name` is "crs", the internal crs object itself is returned.
#' If `name` is "oblique", the internal oblique angle parameter vector is returned.
#' @param name element name
`$.fm_crs` <- function(x, name) {
  if (name %in% c("crs", "oblique")) {
    x[[name]]
  } else {
    eval(parse(text = paste0('x[["crs"]]$', name)))
  }
}



#' @export
#' @param crsonly logical; if `TRUE`, remove any `oblique` information
#' for `fm_crs` class objects and return a pure `crs` class object. Default: `FALSE`.
#' @rdname fm_crs
fm_crs.fm_crs <- function(x, ..., crsonly = FALSE) {
  if (crsonly) {
    fm_crs_oblique(x) <- NULL
  }
  x
}

#' @export
#' @rdname fm_crs
fm_crs.inla.CRS <- function(x, ..., crsonly = FALSE) {
  y <- fm_crs(x$crs, ..., crsonly = TRUE)
  if (!crsonly) {
    fm_crs_oblique(y) <- x$oblique
  }
  y
}

#' @export
#' @rdname fm_crs
fm_crs.character <- function(x, ..., crsonly = FALSE) {
  predef <- fm_wkt_predef()
  if (x %in% names(predef)) {
    x <- predef[[x]]
  }
  if (identical(x, "")) {
    x <- NA_character_
  }
  y <- sf::st_crs(x, ...)
  # Would like nicer proj4string/input display.
  # Possible approach: sf::st_crs(as(sf::st_crs(x), "CRS"))
  # Borrowing from the sf::CRS_from_crs code:
  x <- y$proj4string
  if (is.na(x) || identical(x, y$input)) {
    return(y)
  }
  sf::st_crs(x, ...)
}


#' @rdname fm_crs
#' @export
fm_crs.Spatial <- function(x, ..., crsonly = FALSE) {
  if (is.null(x)) {
    return(sf::NA_crs_)
  }
  suppressWarnings(crs <- fm_crs(sp::wkt(x)))
  crs
}



#' @rdname fm_crs
#' @export
fm_crs.SpatVector <- function(x, ..., crsonly = FALSE) {
  tcrs <- terra::crs(x)
  if (is.null(tcrs) || is.na(tcrs) || identical(tcrs, "")) {
    fm_crs()
  } else {
    fm_crs(tcrs, ...)
  }
}

#' @rdname fm_crs
#' @export
fm_crs.SpatRaster <- function(x, ..., crsonly = FALSE) {
  tcrs <- terra::crs(x)
  if (is.null(tcrs) || is.na(tcrs) || identical(tcrs, "")) {
    fm_crs()
  } else {
    fm_crs(tcrs, ...)
  }
}

#' @rdname fm_crs
#' @export
fm_crs.sf <- function(x, ..., crsonly = FALSE) {
  sf::st_crs(x, ...)
}

#' @rdname fm_crs
#' @export
fm_crs.sfc <- function(x, ..., crsonly = FALSE) {
  sf::st_crs(x, ...)
}

#' @rdname fm_crs
#' @export
fm_crs.sfg <- function(x, ..., crsonly = FALSE) {
  sf::st_crs(x, ...)
}

#' @rdname fm_crs
#' @export
fm_crs.inla.mesh <- function(x, ..., crsonly = FALSE) {
  fm_crs(x[["crs"]], ..., crsonly = crsonly)
}

#' @rdname fm_crs
#' @export
fm_crs.inla.mesh.lattice <- function(x, ..., crsonly = FALSE) {
  fm_crs(x[["crs"]], ..., crsonly = crsonly)
}

#' @rdname fm_crs
#' @export
fm_crs.inla.mesh.segment <- function(x, ..., crsonly = FALSE) {
  fm_crs(x[["crs"]], ..., crsonly = crsonly)
}



# fm_CRS ----

#' Create a coordinate reference system object
#'
#' Creates either a CRS object or an inla.CRS object, describing a coordinate
#' reference system
#'
#' The first two
#' elements of the `oblique` vector are the (longitude, latitude)
#' coordinates for the oblique centre point. The third value (orientation) is a
#' counterclockwise rotation angle for an observer looking at the centre point
#' from outside the sphere. The fourth value is the quasi-longitude (orbit
#' angle) for a rotation along the oblique observers equator.
#'
#' Simple oblique: `oblique=c(0, 45)`
#'
#' Polar: `oblique=c(0, 90)`
#'
#' Quasi-transversal: `oblique=c(0, 0, 90)`
#'
#' Satellite orbit viewpoint: `oblique=c(lon0-time*v1, 0, orbitangle,
#' orbit0+time*v2)`, where `lon0` is the longitude at which a satellite
#' orbit crosses the equator at `time=0`, when the satellite is at an
#' angle `orbit0` further along in its orbit.  The orbital angle relative
#' to the equatorial plane is `orbitangle`, and `v1` and `v2`
#' are the angular velocities of the planet and the satellite, respectively.
#' Note that "forward" from the satellite's point of view is "to the right" in
#' the projection.
#'
#' When `oblique[2]` or `oblique[3]` are non-zero, the resulting
#' projection is only correct for perfect spheres.
#'
#' @param projargs Either 1) a projection argument string suitable as input to
#' `sp::CRS`, or 2) an existing `CRS` object, or 3) a shortcut
#' reference string to a predefined projection; run
#' `names(fm_wkt_predef())` for valid predefined projections.
#' @param doCheckCRSArgs ignored.
#' @param args An optional list of name/value pairs to add to and/or override
#' the PROJ4 arguments in `projargs`.  `name=value` is converted to
#' `"+name=value"`, and `name=NA` is converted to `"+name"`.
#' @param oblique Vector of length at most 4 of rotation angles (in degrees)
#' for an oblique projection, all values defaulting to zero. The values
#' indicate (longitude, latitude, orientation, orbit), as explained in the
#' Details section below.
#' @param SRS_string a WKT2 string defining the coordinate system;
#' see `sp::CRS`. This takes precedence over `projargs`.
#' @param \dots Additional parameters. Not currently in use.
#' @return Either an `sp::CRS` object or an `inla.CRS` object,
#' depending on if the coordinate reference system described by the parameters
#' can be expressed with a pure `sp::CRS` object or not.
#'
#' An S3 `inla.CRS` object is a list, usually (but not necessarily)
#' containing at least one element: \item{crs }{The basic `sp::CRS`
#' object}
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [sp::CRS()], [`fm_crs_wkt`],
#' [fm_sp_get_crs()], [fm_identical_CRS()]
#' @examples
#' crs1 <- fm_CRS("longlat_globe")
#' crs2 <- fm_CRS("lambert_globe")
#' crs3 <- fm_CRS("mollweide_norm")
#' crs4 <- fm_CRS("hammer_globe")
#' crs5 <- fm_CRS("sphere")
#' crs6 <- fm_CRS("globe")
#' @export
#' @rdname fm_CRS_sp
fm_CRS <- function(...) {
  UseMethod("fm_CRS")
}

#' @export
#' @param x Object to convert to CRS or to extract CRS information from.
#' @rdname fm_CRS_sp
fm_CRS.crs <- function(x, ...) {
  if (is.na(x)) {
    return(fm_CRS(NA_character_))
  }
  fm_CRS(SRS_string = x$wkt)
}

#' @export
#' @rdname fm_CRS_sp
fm_CRS.fm_crs <- function(x, ...) {
  if (fm_crs_is_null(x)) {
    return(fm_CRS(NA_character_))
  }
  fm_CRS(SRS_string = x$crs$wkt, oblique = fm_crs_oblique(x))
}

#' @rdname fm_CRS_sp
#' @export
fm_CRS.Spatial <- function(x, ...) {
  suppressWarnings(crs <- sp::CRS(SRS_string = sp::wkt(x)))
  crs
}

#' @rdname fm_CRS_sp
#' @param crsonly logical; if `TRUE`, remove any oblique` information
#' for `inla.CRS` class objects and return a pure `CRS` class object.
#' Default: `FALSE`.
#' @export
fm_CRS.inla.CRS <- function(x, ..., crsonly = FALSE) {
  y <- x[["crs"]]
  if (crsonly) {
    y <- y[["crs"]]
  }
  y
}

#' @rdname fm_CRS_sp
#' @export
fm_CRS.sf <- function(x, ..., crsonly = FALSE) {
  fm_CRS(sf::st_crs(x, ...))
}

#' @rdname fm_CRS_sp
#' @export
fm_CRS.sfc <- function(x, ..., crsonly = FALSE) {
  fm_CRS(sf::st_crs(x, ...))
}

#' @rdname fm_CRS_sp
#' @export
fm_CRS.sfg <- function(x, ..., crsonly = FALSE) {
  fm_CRS(sf::st_crs(x, ...))
}

#' @rdname fm_CRS_sp
#' @export
fm_CRS.inla.mesh <- function(x, ..., crsonly = FALSE) {
  fm_CRS(x[["crs"]], ..., crsonly = crsonly)
}

#' @rdname fm_CRS_sp
#' @export
fm_CRS.inla.mesh.lattice <- function(x, ..., crsonly = FALSE) {
  fm_CRS(x[["crs"]], ..., crsonly = crsonly)
}

#' @rdname fm_CRS_sp
#' @export
fm_CRS.inla.mesh.segment <- function(x, ..., crsonly = FALSE) {
  fm_CRS(x[["crs"]], ..., crsonly = crsonly)
}

#' @export
#' @rdname fm_CRS_sp
fm_CRS.CRS <- function(x, oblique = NULL,
                       ...) {
  if (!is.null(oblique)) {
    stopifnot(is.vector(oblique))
    if (length(oblique) < 4) {
      oblique <- c(oblique, rep(0, 4 - length(oblique)))
    }
    x <- list(crs = x, oblique = oblique)
    class(x) <- "inla.CRS"
  }
  x
}

#' @export
#' @rdname fm_CRS_sp
fm_CRS.default <- function(projargs = NULL, doCheckCRSArgs = NULL,
                           args = NULL, oblique = NULL,
                           SRS_string = NULL,
                           ...) {
  if (identical(projargs, "")) {
    projargs <- NULL
  }
  if (is.null(SRS_string) &&
    !is.null(projargs) &&
    !is.na(projargs) &&
    is.character(projargs)) {
    if (projargs %in% c("hammer", "lambert", "longlat", "mollweide")) {
      warning(paste0(
        "Use of old predefined projection '",
        projargs,
        "' is deprecated. Converting to '",
        projargs,
        "_norm'"
      ))
      projargs <- paste0(projargs, "_norm")
    }
    predef <- fm_wkt_predef()
    if (projargs %in% names(predef)) {
      SRS_string <- predef[[projargs]]
      projargs <- NULL
    } else {
      warning(
        paste0(
          "'fm_CRS' should be given a SRS_string for PROJ6 or a known keyword for a\n",
          "  predefined string given in projargs. Using 'fm_crs()' workaround."
        )
      )
      if (!is.null(args)) {
        x <- fm_crs(projargs)
        if (typeof(args) != "list") {
          stop("'args' must be NULL or a list of name=value pairs.")
        }
        xargs <- fm_CRS_as_list(x)
        for (name in names(args)) {
          xargs[[name]] <- args[[name]]
        }
        projargs <- fm_list_as_CRSargs(xargs)
      }
      SRS_string <- fm_crs(projargs)$wkt
      projargs <- NULL
    }
  }

  if (!is.null(SRS_string)) {
    if (!is.null(projargs)) {
      warning("SRS_string specified. Ignoring non-null projargs.")
    }
    x <- as(sf::st_crs(SRS_string), "CRS")
  } else if (inherits(projargs, "CRS")) {
    x <- projargs
  } else {
    x <- sp::CRS(NA_character_)
  }

  if (!is.null(oblique)) {
    stopifnot(is.vector(oblique))
    if (length(oblique) < 4) {
      oblique <- c(oblique, rep(0, 4 - length(oblique)))
    }
    x <- list(crs = x, oblique = oblique)
    class(x) <- "inla.CRS"
  }
  x
}

#' @return `fm_wkt_predef` returns a WKT2 string defining a projection
#' @examples
#' \dontrun{
#' names(fm_wkt_predef())
#' }
#' @export
#' @rdname fm_crs

fm_wkt_predef <- function() {
  list(
    hammer_norm = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["unknown",ELLIPSOID["unknown",0.707106781186548,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Reference meridian",0,ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]],CONVERSION["unknown",METHOD["PROJ hammer"]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]',
    lambert_norm = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["unknown",ELLIPSOID["unknown",1,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Reference meridian",0,ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]],CONVERSION["unknown",METHOD["Lambert Cylindrical Equal Area (Spherical)",ID["EPSG",9834]],PARAMETER["Latitude of 1st standard parallel",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8823]],PARAMETER["Longitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["False easting",0,LENGTHUNIT["metre",1],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]',
    longlat_norm = 'GEOGCRS["unknown",DATUM["unknown",ELLIPSOID["unknown",1,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Reference meridian",0,ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]],CS[ellipsoidal,2],AXIS["longitude",east,ORDER[1],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]],AXIS["latitude",north,ORDER[2],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]]',
    mollweide_norm = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["unknown",ELLIPSOID["unknown",0.707106781186548,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Reference meridian",0,ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]],CONVERSION["unknown",METHOD["Mollweide"],PARAMETER["Longitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["False easting",0,LENGTHUNIT["metre",1],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]',
    hammer_globe = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["Unknown based on Normal Sphere (r=6370997) ellipsoid",ELLIPSOID["Normal Sphere (r=6370997)",6370997,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]]],CONVERSION["unknown",METHOD["PROJ hammer"]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]]]',
    lambert_globe = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["Unknown based on Normal Sphere (r=6370997) ellipsoid",ELLIPSOID["Normal Sphere (r=6370997)",6370997,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]]],CONVERSION["unknown",METHOD["Lambert Cylindrical Equal Area (Spherical)",ID["EPSG",9834]],PARAMETER["Latitude of 1st standard parallel",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8823]],PARAMETER["Longitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["False easting",0,LENGTHUNIT["kilometre",1000],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["kilometre",1000],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]]]',
    longlat_globe = 'GEOGCRS["unknown",DATUM["Unknown based on Normal Sphere (r=6370997) ellipsoid",ELLIPSOID["Normal Sphere (r=6370997)",6370997,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]],CS[ellipsoidal,2],AXIS["longitude",east,ORDER[1],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]],AXIS["latitude",north,ORDER[2],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]]',
    mollweide_globe = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["Unknown based on Normal Sphere (r=6370997) ellipsoid",ELLIPSOID["Normal Sphere (r=6370997)",6370997,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]]],CONVERSION["unknown",METHOD["Mollweide"],PARAMETER["Longitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["False easting",0,LENGTHUNIT["kilometre",1000],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["kilometre",1000],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]]]',
    sphere = 'GEODCRS["unknown",DATUM["unknown",ELLIPSOID["unknown",1,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Reference meridian",0,ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]],CS[Cartesian,3],AXIS["(X)",geocentricX,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(Y)",geocentricY,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(Z)",geocentricZ,ORDER[3],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]',
    globe = 'GEODCRS["unknown",DATUM["Unknown based on Normal Sphere (r=6370997) ellipsoid",ELLIPSOID["Normal Sphere (r=6370997)",6370997,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]],CS[Cartesian,3],AXIS["(X)",geocentricX,ORDER[1],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]],AXIS["(Y)",geocentricY,ORDER[2],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]],AXIS["(Z)",geocentricZ,ORDER[3],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]]]'
  )
}





#' Internal WKT handling
#'
#' Conversion between WKT and a tree representation
#'
#' @param x A WKT2 string, or a `wkt_tree` list structure
#' @param \dots Unused
#'
#' @rdname wkt_tree
#' @export

fm_wkt_as_wkt_tree <- function(x, ...) {
  # Basic parsing of WKT string
  # ITEM[Param1, Param2, ...]
  # Param can be a constant or an ITEM[...]
  parse_item <- function(x) {
    # Parse item label
    n <- regexpr("\\[", x)
    item_label <- substr(x, 1, n - 1)
    x <- substr(x, n + 1, nchar(x))
    params <- list()
    # Parse parameters
    done <- FALSE
    while (!done) {
      # If [ comes before , or ], it's an Item, otherwise a constant
      n <- regexpr("[],[]", x)
      if (n < 1) {
        # Nothing found
        done <- TRUE
        break
      }
      if (substr(x, n, n) == "[") {
        n <- regexpr("[^ ]", x)
        x <- substr(x, n, nchar(x))
        item_info <- parse_item(x)
        params[[length(params) + 1]] <- item_info$item
        x <- item_info$x
      } else {
        const <- substr(x, 1, n - 1)
        params[[length(params) + 1]] <- const
        x <- substr(x, n, nchar(x))
      }
      n <- regexpr("[],]", x)
      done <- (substr(x, n, n) == "]")
      x <- substr(x, n + 1, nchar(x))
    }
    list(item = list(label = item_label, params = params), x = x)
  }

  x <- gsub("\n", "", x)
  item_info <- parse_item(x)
  item <- item_info$item
  item
}

#' @param pretty logical; If TRUE, use pretty formatting. Default: FALSE
#' @rdname wkt_tree
#' @export

fm_wkt_tree_as_wkt <- function(x, pretty = FALSE, ...) {
  construct_item <- function(x, level) {
    paste0(
      if (pretty) {
        paste0(rep("    ", level), collapse = "")
      } else {
        ""
      },
      x[["label"]],
      "[",
      paste0(
        vapply(
          x[["params"]],
          function(param) {
            if (!is.list(param)) {
              paste0(param)
            } else {
              paste0(
                if (pretty) {
                  "\n"
                } else {
                  ""
                },
                construct_item(param,
                  level = level + 1
                )
              )
            }
          },
          ""
        ),
        collapse = ","
      ),
      "]"
    )
  }
  construct_item(x, 0)
}

#' @param item character vector with item labels identifying a parameter item
#' entry.
#' @param duplicate For items that have more than one match, `duplicate`
#' indicates the index number of the desired version. Default: 1
#' @rdname wkt_tree
#' @export

fm_wkt_tree_get_item <- function(x, item, duplicate = 1) {
  for (k in seq_along(x[["params"]])) {
    if (is.list(x[["params"]][[k]]) &&
      (!is.null(x[["params"]][[k]][["label"]])) &&
      (x[["params"]][[k]][["label"]] %in% item)) {
      if (duplicate == 1) {
        return(x[["params"]][[k]])
      }
      duplicate <- duplicate - 1
    }
  }
  NULL
}

#' @param item_tree An item tree identifying a parameter item entry
#' @rdname wkt_tree
#' @export

fm_wkt_tree_set_item <- function(x, item_tree, duplicate = 1) {
  success <- FALSE
  for (k in seq_along(x[["params"]])) {
    if (is.list(x[["params"]][[k]]) && (x[["params"]][[k]][["label"]] == item_tree[["label"]])) {
      if (duplicate == 1) {
        x[["params"]][[k]] <- item_tree
      }
      duplicate <- duplicate - 1
      success <- TRUE
      break
    }
  }
  if (!success) {
    x[["params"]] <- c(x[["params"]], list(item_tree))
  }
  x
}




#' @export
#' @rdname fm_CRSargs
fm_CRS_as_list <- function(x, ...) {
  fm_CRSargs_as_list(fm_proj4string(x))
}


#' @export
#' @rdname fm_CRSargs
fm_list_as_CRS <- function(x, ...) {
  fm_CRS(fm_list_as_CRSargs(x))
}

#' Show expanded CRS arguments
#'
#' `r lifecycle::badge("deprecated")`
#' Wrappers for `sp::CRS` and `inla.CRS` objects to handle the
#' coordinate reference system argument string.
#' These methods should no longer be used with PROJ6/rgdal3;
#' see [fm_wkt()] and [fm_proj4string()] for a new approach.
#'
#' @aliases fm_CRSargs fm_CRS_as_list fm_CRSargs_as_list fm_list_as_CRS
#' fm_list_as_CRSargs
#' @param x An `sp::CRS` or `inla.CRS` object (for
#' `fm_CRSargs` and `fm_CRS_as_list`), a character string (for
#' `fm_CRSargs_as_list`), or a list (for `fm_list_as_CRS` and
#' `fm_list_as_CRSargs`).
#' @param \dots Additional arguments passed on to other methods.
#' @return For `fm_CRSargs` and `fm_list_as_CRSargs`, a character
#' string with PROJ.4 arguments.
#'
#' For `fm_CRS_as_list` and `fm_CRSargs_as_list`, a list of
#' name/value pairs.
#'
#' For `fm_list_as_CRS`, a `CRS` or `inla.CRS` object.
#' @author Finn Lindgren <finn.lindgren@@gmail.com>
#' @seealso [fm_CRS()]
#' @export
#' @keywords internal
#' @examples
#'
#' crs0 <- fm_CRS("longlat")
#' p4s <- fm_proj4string(crs0)
#' lst <- fm_CRSargs_as_list(p4s)
#' crs1 <- fm_list_as_CRS(lst)
#' lst$a <- 2
#' crs2 <- fm_CRS(p4s, args = lst)
#' print(fm_proj4string(crs0))
#' print(fm_proj4string(crs1))
#' print(fm_proj4string(crs2))
fm_CRSargs <- function(x, ...) {
  lifecycle::deprecate_warn("2.8.0", "fm_CRSargs()", "fm_proj4string()")

  fm_proj4string(x)
}


#' @return For `fm_list_as_CRSargs()`, a CRS proj4 string for name=value pair list
#' @rdname fm_CRSargs
fm_list_as_CRSargs <- function(x, ...) {
  paste(
    lapply(
      names(x),
      function(xx) {
        if (is.na(x[[xx]])) {
          paste("+", xx, sep = "")
        } else {
          paste("+", xx, "=", x[[xx]], sep = "")
        }
      }
    ),
    collapse = " "
  )
}

#' @return For `fm_CRSargs_as_list()`, a list of name=value pairs from CRS proj4 string
#' @rdname fm_CRSargs
#' @export
fm_CRSargs_as_list <- function(x, ...) {
  if (is.na(x)) {
    return(list())
  }
  if (!is.character(x)) {
    stop("proj4string must be of class character")
  }
  do.call(c, lapply(
    strsplit(
      x = strsplit(
        x = paste(" ", x, sep = ""),
        split = " \\+"
      )[[1]][-1],
      split = "="
    ),
    function(x) {
      xx <- list(x[2])
      names(xx) <- x[1]
      xx
    }
  ))
}






# fm_wkt ----

#' @describeIn fm_crs_wkt Returns a WKT2 string, for any input supported by [fm_crs()].
#' @export

fm_wkt <- function(crs) {
  fm_crs(crs, crsonly = TRUE)$wkt
}

#' @describeIn fm_crs_wkt Returns a proj4 string, for any input supported by [fm_crs()].
#' @export
fm_proj4string <- function(crs) {
  fm_crs(crs, crsonly = TRUE)$proj4string
}

#' @export
#' @describeIn fm_crs_wkt `r lifecycle::badge("deprecated")` Use [fm_wkt()]
#' instead.

fm_crs_get_wkt <- function(crs) {
  lifecycle::deprecate_warn(
    "2.8.0",
    "fm_crs_get_wkt()",
    "fm_wkt()"
  )
  fm_wkt(crs)
}


fm_rotmat3213 <- function(rot) {
  cs <- cos(rot[1])
  sn <- sin(rot[1])
  R <- matrix(c(
    cs, -sn, 0,
    sn, cs, 0,
    0, 0, 1
  ), 3, 3)
  cs <- cos(rot[2])
  sn <- sin(rot[2])
  R <- R %*% matrix(c(
    cs, 0, sn,
    0, 1, 0,
    -sn, 0, cs
  ), 3, 3)
  cs <- cos(rot[3])
  sn <- sin(rot[3])
  R <- R %*% matrix(c(
    1, 0, 0,
    0, cs, -sn,
    0, sn, cs
  ), 3, 3)
  cs <- cos(rot[4])
  sn <- sin(rot[4])
  R <- R %*% matrix(c(
    cs, -sn, 0,
    sn, cs, 0,
    0, 0, 1
  ), 3, 3)
  R
}

fm_rotmat3123 <- function(rot) {
  cs <- cos(rot[4])
  sn <- sin(rot[4])
  R <- matrix(c(
    cs, -sn, 0,
    sn, cs, 0,
    0, 0, 1
  ), 3, 3)
  cs <- cos(rot[3])
  sn <- sin(rot[3])
  R <- R %*% matrix(c(
    1, 0, 0,
    0, cs, -sn,
    0, sn, cs
  ), 3, 3)
  cs <- cos(rot[2])
  sn <- sin(rot[2])
  R <- R %*% matrix(c(
    cs, 0, sn,
    0, 1, 0,
    -sn, 0, cs
  ), 3, 3)
  cs <- cos(rot[1])
  sn <- sin(rot[1])
  R <- R %*% matrix(c(
    cs, -sn, 0,
    sn, cs, 0,
    0, 0, 1
  ), 3, 3)
  R
}


fm_crs_transform_oblique <- function(x, oblique, to.oblique = TRUE) {
  if (to.oblique) {
    ## Transform to oblique orientation
    ## 1) Rotate -oblique[1] around (0,0,1)
    ## 2) Rotate +oblique[2] around (0,1,0)
    ## 3) Rotate -oblique[3] around (1,0,0)
    ## 3) Rotate -oblique[4] around (0,0,1)
    x %*% fm_rotmat3213(c(-1, 1, -1, -1) * oblique * pi / 180)
  } else {
    ## Transform back from oblique orientation
    ## 1) Rotate +oblique[4] around (0,0,1)
    ## 2) Rotate +oblique[3] around (1,0,0)
    ## 3) Rotate -oblique[2] around (0,1,0)
    ## 4) Rotate +oblique[1] around (0,0,1)
    x %*% fm_rotmat3123(c(1, -1, 1, 1) * oblique * pi / 180)
  }
}



fm_wkt_tree_projection_type <- function(wt) {
  axis1 <- fm_wkt_tree_get_item(wt, "AXIS", 1)
  axis2 <- fm_wkt_tree_get_item(wt, "AXIS", 2)
  if (identical(axis1[["params"]][[1]], '"longitude"') &&
    identical(axis2[["params"]][[1]], '"latitude"')) {
    return("longlat")
  }
  conversion <- fm_wkt_tree_get_item(wt, "CONVERSION")
  if (!is.null(conversion)) {
    method <- fm_wkt_tree_get_item(conversion, "METHOD")
    if (identical(method[["params"]][[1]], '"Lambert Cylindrical Equal Area (Sherical)"')) {
      return("lambert")
    }
    if (identical(method[["params"]][[1]], '"Mollweide"')) {
      return("mollweide")
    }
    if (identical(method[["params"]][[1]], '"PROJ hammer"')) {
      return("hammer")
    }
    if (identical(method[["params"]][[1]], '"tmerc"')) {
      return("tmerc")
    }
  }
  NULL
}

fm_wkt_projection_type <- function(wkt) {
  wt <- fm_wkt_as_wkt_tree(wkt)
  fm_wkt_tree_projection_type(wt)
}

fm_crs_projection_type <- function(crs) {
  wkt <- fm_wkt(crs)
  fm_wkt_projection_type(wkt)
}

## +proj=longlat in (-180,180)x(-90,90)
## +proj=moll in (-2,2)x(-1,1) scaled by +a and +b, and +units
## +proj=lambert in (-pi,pi)x(-1,1) scaled by +a and +b, and +units
fm_crs_bounds <- function(crs, warn.unknown = FALSE) {
  wkt <- fm_wkt(crs)
  wt <- fm_wkt_as_wkt_tree(wkt)
  type <- fm_wkt_tree_projection_type(wt)

  if (is.null(type)) {
    if (fm_wkt_is_geocent(wkt)) {
      bounds <- list(type = "rectangle", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf))
    } else {
      if (warn.unknown) {
        warning("Could not determine shape of transformation bounds. Using infinite rectangle.")
      }
      bounds <- list(type = "rectangle", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf))
    }
  } else if (type == "longlat") {
    bounds <- list(type = "rectangle", xlim = c(-180, 180), ylim = c(-90, 90))
  } else if (type == "lambert") {
    axis <- c(pi, 1)
    radius <- fm_wkt_get_ellipsoid_radius(wkt)
    axis[1] <- axis[1] * radius
    # TODO: handle eccentricity
    axis[2] <- axis[2] * sqrt(radius) * sqrt(radius)
    # TODO: Handle units"
    bounds <- list(
      type = "rectangle",
      xlim = c(-1, 1) * axis[1],
      ylim = c(-1, 1) * axis[2]
    )
  } else if (type %in% c("mollweide", "hammer")) {
    axis <- c(2, 1)
    center <- c(0, 0)
    radius <- fm_wkt_get_ellipsoid_radius(wkt)
    axis[1] <- axis[1] * radius / sqrt(1 / 2)
    axis[2] <- axis[2] * radius / sqrt(1 / 2)
    # TODO: Handle "units"
    bounds <- list(
      type = "ellipse", axis = axis, center = center,
      xlim = center[1] + c(-1, 1) * axis[1],
      ylim = center[2] + c(-1, 1) * axis[2]
    )
  } else if (type == "tmerc") {
    bounds <- list(type = "rectangle", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf))
  } else {
    stop("'fm_crs_bounds' internal error: transformation detected but not handled.")
  }

  if (bounds$type == "rectangle") {
    bounds$polygon <- cbind(
      bounds$xlim[c(1, 2, 2, 1, 1)],
      bounds$ylim[c(1, 1, 2, 2, 1)]
    )
  } else if (bounds$type == "ellipse") {
    theta <- seq(0, 2 * pi, length = 1000)
    bounds$polygon <- cbind(
      bounds$center[1] + bounds$axis[1] * cos(theta),
      bounds$center[2] + bounds$axis[2] * sin(theta)
    )
  } else {
    stop("Unknown transformation type. This should not happen.")
  }
  bounds
}

## TRUE/FALSE for points inside/outside projection domain.

# ! sp code here
fm_crs_bounds_check <- function(x, bounds) {
  stopifnot(inherits(x, "matrix"))
  if (all(is.finite(bounds$xlim)) && all(is.finite(bounds$ylim))) {
    (sp::point.in.polygon(
      x[, 1], x[, 2],
      bounds$polygon[, 1], bounds$polygon[, 2]
    )
    > 0)
  } else {
    rep(TRUE, nrow(x))
  }
}


# ! show() for "CRS" and "crs" class appear identical
fm_internal_update_crs <- function(crs, newcrs, mismatch.allowed) {
  if (is.null(crs)) {
    newcrs
  } else {
    if (!mismatch.allowed && !identical(crs, newcrs)) {
      methods::show(crs)
      methods::show(newcrs)
      stop("CRS information mismatch.")
    }
    crs
  }
}


#' @title Check if two CRS objects are identical
#' @param crs0,crs1 Two `sf::crs`, `sp::CRS`, `fm_crs` or `inla.CRS` objects to be compared.
#' @param crsonly logical. If `TRUE` and any of `crs0` and `crs1` are `fm_crs` or `inla.CRS`
#' objects, extract and compare only the `sf::crs` or `sp::CRS` aspects. Default: `FALSE`
#' @export
#' @keywords internal
#' @seealso [fm_crs()], [fm_CRS()]
#' @examples
#'
#' crs0 <- crs1 <- fm_crs("longlat_globe")
#' fm_crs_oblique(crs1) <- c(0, 90)
#' print(c(
#'   fm_identical_CRS(crs0, crs0),
#'   fm_identical_CRS(crs0, crs1),
#'   fm_identical_CRS(crs0, crs1, crsonly = TRUE)
#' ))
fm_identical_CRS <- function(crs0, crs1, crsonly = FALSE) {
  if (crsonly) {
    crs0 <- fm_crs(crs0, crsonly = TRUE)
    crs1 <- fm_crs(crs1, crsonly = TRUE)
  }
  identical(crs0, crs1)
}





# fm_transform ----

#' @title Object coordinate transformation
#'
#' @description
#' Handle transformation of various inla objects according to coordinate
#' reference systems of `crs` (from `sf::st_crs()`), `fm_crs`,¬`sp::CRS` or
#' `INLA::inla.CRS` class.
#'
#' @param x
#' The object that should be transformed from it's current CRS to a new CRS
#' @param crs
#' The target crs object
#' @param crs0
#' The source crs object for spatial classes without crs information
#' @param passthrough
#' Default is FALSE.
#' Setting to TRUE allows objects with no CRS information to be passed
#' through without transformation. Use with care!
#' @param \dots
#' Potential additional arguments
#' @seealso [fm_CRS()]
#' @export
fm_transform <- function(x, crs = fm_crs(x), ...) {
  UseMethod("fm_transform")
}

#' @rdname fm_transform
#' @export
fm_transform.default <- function(x, crs = fm_crs(x), ..., crs0 = NULL) {
  stop(paste0(
    "fm_transform() for '",
    paste0(class(x), sep = ", "),
    "' not yet implemented."
  ))
}


fm_transform_raw <- function(x, from, to) {
  adjust_input <- function(x, crs) {
    if (fm_crs_is_geocent(crs) &&
      ncol(x) == 2) {
      if (nrow(x) > 0) {
        x <- cbind(x, 0)
      } else {
        x <- matrix(0, 0, ncol(x) + 1)
      }
    }
    x
  }

  adjust_output <- function(x, crs) {
    if (!fm_crs_is_geocent(crs) &&
      ncol(x) == 3) {
      if (nrow(x) > 0) {
        x <- x[, 1:2, drop = FALSE]
      } else {
        x <- matrix(0, 0, 2)
      }
    }
    x
  }

  x <- adjust_input(x, crs = to)
  if (nrow(x) == 0) {
    x
  } else {
    y <- sf::sf_project(
      pts = x,
      from = from,
      to = to,
      keep = TRUE
    )
    adjust_output(y, to)
  }
}


#' @rdname fm_transform
#' @export
fm_transform.matrix <- function(x, crs = NULL, ..., passthrough = FALSE, crs0 = NULL) {
  crs1 <- fm_crs(crs)
  crs0 <- fm_crs(crs0)
  if (fm_crs_is_null(crs0) || fm_crs_is_null(crs1)) {
    if (!passthrough) {
      if (fm_crs_is_null(crs0)) {
        stop("'crs0' is an invalid coordinate reference object.")
      }
      if (fm_crs_is_null(crs1)) {
        stop("'crs' is an invalid coordinate reference object.")
      }
    }
    return(x)
  }
  #  if (ncol(x) == 2) {
  #    x <- cbind(x, 0)
  #  }
  sphere_radius_0 <- fm_crs_get_ellipsoid_radius(crs0)
  sphere_radius_1 <- fm_crs_get_ellipsoid_radius(crs1)
  different_radii <- (sphere_radius_0 != sphere_radius_1)
  longlat_norm <- fm_crs("longlat_norm")
  longlat_0 <- fm_crs_set_ellipsoid_radius(longlat_norm, sphere_radius_0)
  longlat_1 <- fm_crs_set_ellipsoid_radius(longlat_norm, sphere_radius_1)

  crs_sphere <- fm_crs("sphere")
  onsphere_0 <- fm_identical_CRS(crs0, crs_sphere, crsonly = TRUE)
  onsphere_1 <- fm_identical_CRS(crs1, crs_sphere, crsonly = TRUE)
  is_geocentric_0 <- fm_crs_is_geocent(crs0)
  is_geocentric_1 <- fm_crs_is_geocent(crs1)
  if (is_geocentric_0) {
    ok <- rep(TRUE, nrow(x))
  } else {
    bounds <- fm_crs_bounds(crs0)
    if (identical(fm_crs_projection_type(crs0), "longlat")) {
      ## Wrap longitudes to [-180,180]
      needswrap <- (x[, 1] < -180) | (x[, 1] > 180)
      if (any(needswrap)) {
        x[needswrap, 1] <- ((x[needswrap, 1] + 180) %% 360) - 180
      }
    }
    ok <- fm_crs_bounds_check(x, bounds)
    if (!all(ok)) {
      xx <- x
    }
  }
  do_work_on_sphere <-
    inherits(crs0, "fm_crs") ||
      inherits(crs1, "fm_crs") ||
      different_radii
  x <- x[ok, , drop = FALSE]
  current_crs <- fm_crs(crs0, crsonly = TRUE)
  if (do_work_on_sphere) {
    if (!onsphere_0) {
      if (sphere_radius_0 != 1) {
        x <- fm_transform_raw(
          x,
          from = current_crs,
          to = longlat_0
        )
        # x can now be treated as being in longlat_norm coordinates
        current_crs <- longlat_norm
      }
      x <- fm_transform_raw(
        x,
        from = current_crs,
        to = crs_sphere
      )
      current_crs <- crs_sphere
    }
    if (!is.null(crs0$oblique)) {
      x <- fm_crs_transform_oblique(
        x,
        fm_crs_oblique(crs0),
        to.oblique = FALSE
      )
    }

    if (!is.null(crs1$oblique)) {
      x <- fm_crs_transform_oblique(
        x,
        fm_crs_oblique(crs1),
        to.oblique = TRUE
      )
    }
    if (sphere_radius_1 != 1) {
      x <- fm_transform_raw(
        x,
        from = current_crs,
        to = longlat_norm
      )
      # x can now be treated as being in longlat_1
      current_crs <- longlat_1
    }
  }

  x <- fm_transform_raw(
    x,
    from = current_crs,
    to = crs1
  )

  if (!all(ok)) {
    xx[ok, ] <- x
    xx[!ok, ] <- NA
    x <- xx
  }

  x
}

#' @rdname fm_transform
#' @export
fm_transform.list <- function(x, crs = fm_crs(x), ...) {
  if (!fm_crs_is_null(crs)) {
    x <- lapply(x, function(xx) fm_transform(xx, crs = crs, ...))
  }
  x
}


# Internal helper function for fm_transform.sf/sfc/sfg
fm_transform_sf <- function(x, crs, ..., passthrough) {
  crs1 <- fm_crs(crs)
  crs0 <- fm_crs(x)
  if (fm_crs_is_null(crs0) || fm_crs_is_null(crs1)) {
    if (!passthrough) {
      if (fm_crs_is_null(crs0)) {
        stop("'crs0' is an invalid coordinate reference object.")
      }
      if (fm_crs_is_null(crs1)) {
        stop("'crs' is an invalid coordinate reference object.")
      }
    }
    return(x)
  }

  if (inherits(x, "sfc_POINT")) {
    adjust_input <- function(x, crs) {
      if (fm_crs_is_geocent(crs) &&
        length(x) &&
        inherits(x[[1]], "XY")) {
        x <- sf::st_zm(x = x, drop = FALSE, what = "Z")
      }
      x
    }

    x <- adjust_input(x, crs1)
    coord <- sf::st_coordinates(x)
    M <- if ("M" %in% colnames(coord)) coord[, "M"] else NULL
    coord <- coord[, intersect(colnames(coord), c("X", "Y", "Z")), drop = FALSE]
    coord <- fm_transform(coord, crs = crs1, crs0 = crs0, passthrough = passthrough)
    if (is.null(M)) {
      the_dim <- c("X", "XY", "XYZ")[ncol(coord)]
    } else {
      the_dim <- c("XM", "XYM", "XYZM")[ncol(coord)]
      coord <- cbind(coord, M)
    }
    x <-
      do.call(sf::st_sfc, c(
        lapply(
          seq_len(nrow(coord)),
          function(k) sf::st_point(coord[k, , drop = FALSE], dim = the_dim)
        ),
        list(crs = sf::st_crs(crs1))
      ))
  } else {
    x <- sf::st_transform(x, sf::st_crs(crs1))
  }
  x
}

#' @export
#' @rdname fm_transform
fm_transform.sf <- function(x, crs = fm_crs(x), ..., passthrough = FALSE) {
  geo <- fm_transform(sf::st_geometry(x), crs = crs, ..., passthrough = passthrough)
  sf::st_geometry(x) <- geo
  x
}
#' @export
#' @rdname fm_transform
fm_transform.sfc <- function(x, crs = fm_crs(x), ..., passthrough = FALSE) {
  fm_transform_sf(x, crs = crs, ..., passthrough = passthrough)
}
#' @export
#' @rdname fm_transform
fm_transform.sfg <- function(x, crs = fm_crs(x), ..., passthrough = FALSE) {
  fm_transform_sf(x, crs = crs, ..., passthrough = passthrough)
}

#' @export
#' @rdname fm_transform
fm_transform.Spatial <- function(x, crs = fm_crs(x), ..., passthrough = FALSE) {
  orig_class <- class(x)
  as(
    as(
      fm_transform(
        sf::st_as_sf(x),
        crs = crs,
        ...,
        passthrough = passthrough
      ),
      "Spatial"
    ),
    orig_class[1]
  )
}



fm_crs_detect_manifold <- function(crs) {
  if (fm_crs_is_geocent(crs)) {
    manifold <- "S2"
  } else {
    manifold <- "R2"
  }
  manifold
}

#' @export
#' @rdname fm_transform
fm_transform.inla.mesh <- function(x,
                                   crs = fm_crs(x),
                                   ...,
                                   passthrough = FALSE,
                                   crs0 = fm_crs(x)) {
  x$loc <- fm_transform(x$loc, crs = crs, ..., crs0 = x$crs, passthrough = passthrough)
  x$manifold <- fm_crs_detect_manifold(crs)
  x$crs <- fm_CRS(crs)
  x
}

#' @export
#' @rdname fm_transform
fm_transform.inla.mesh.lattice <- function(x,
                                           crs = fm_crs(x),
                                           ...,
                                           passthrough = FALSE,
                                           crs0 = fm_crs(x)) {
  x$segm <- fm_transform(x$segm, crs = crs, crs0 = x$crs, ..., passthrough = passthrough)
  x$loc <- fm_transform(x$loc, crs = crs, crs0 = x$crs, ..., passthrough = passthrough)
  x$crs <- fm_CRS(crs)
  invisible(x)
}

#' @export
#' @rdname fm_transform
fm_transform.inla.mesh.segment <- function(x, crs = fm_crs(x), ..., passthrough = FALSE) {
  x$loc <- fm_transform(x$loc, crs = crs, crs0 = x$crs, ..., passthrough = passthrough)
  x$crs <- fm_CRS(crs)
  invisible(x)
}



# fm_spTransform ----

#' @describeIn inlabru-deprecated Handle transformation of various inla objects according to coordinate
#' reference systems of sp::CRS or INLA::inla.CRS class.
#'
#' @param x
#' The object that should be transformed from it's current CRS to a new CRS
#' @param crs0
#' The source sp::CRS or inla.CRS object
#' @param crs1
#' The target sp::CRS or inla.CRS object
#' @param CRSobj
#' The target sp::CRS or inla.CRS object
#' @param passthrough
#' Default is FALSE.
#' Setting to TRUE allows objects with no CRS information to be passed
#' through without transformation.
#' @param \dots
#' Potential additional arguments
#' @seealso [fm_transform()]
#' @export
fm_spTransform <- function(x, ...) {
  lifecycle::deprecate_soft("2.8.0", "fm_spTransform()", "fm_transform()")
  UseMethod("fm_spTransform")
}

#' @describeIn inlabru-deprecated The default method handles low level transformation of raw
#' coordinates.
#' @export
fm_spTransform.default <- function(x, crs0 = NULL, crs1 = NULL, passthrough = FALSE, ...) {
  fm_transform(x, crs = crs1, crs0 = crs0, passthrough = passthrough)
}

#' @export
#' @rdname inlabru-deprecated
fm_spTransform.SpatialPoints <- function(x, CRSobj, passthrough = FALSE, ...) {
  fm_transform(x, crs = CRSobj, passthrough = passthrough)
}

#' @export
#' @rdname inlabru-deprecated
fm_spTransform.SpatialPointsDataFrame <- function(x,
                                                  CRSobj,
                                                  passthrough = FALSE,
                                                  ...) {
  fm_transform(x, crs = CRSobj, passthrough = passthrough)
}

#' @export
#' @rdname inlabru-deprecated
fm_spTransform.inla.mesh.lattice <- function(x, CRSobj, passthrough = FALSE, ...) {
  fm_transform(x, crs = CRSobj, passthrough = passthrough)
}

#' @export
#' @rdname inlabru-deprecated
fm_spTransform.inla.mesh.segment <- function(x, CRSobj, passthrough = FALSE, ...) {
  fm_transform(x, crs = CRSobj, passthrough = passthrough)
}

#' @export
#' @rdname inlabru-deprecated
fm_spTransform.inla.mesh <- function(x, CRSobj, passthrough = FALSE, ...) {
  fm_transform(x, crs = CRSobj, passthrough = passthrough)
}
