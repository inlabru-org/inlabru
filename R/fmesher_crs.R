#' @importFrom sp coordinates proj4string `proj4string<-`

#' @title PROJ6 detection
#' @description Detect whether PROJ6 is available
#'
#' @return For `fm_has_PROJ6`, logical;
#' `TRUE` if PROJ6 is available, `FALSE` otherwise
#' @examples
#' fm_has_PROJ6()
#' @export
#' @rdname fm_has_PROJ6

fm_has_PROJ6 <- function() {
  stopifnot(requireNamespace("rgdal", quietly = TRUE))
  result <- tryCatch(rgdal::new_proj_and_gdal(),
    error = function(e) FALSE
  )
  result
}

#' @details `fm_not_for_PROJ6` is called to warn about using old PROJ4
#' features even though PROJ6 is available
#' @rdname fm_has_PROJ6

fm_not_for_PROJ6 <- function(fun = NULL) {
  if (fm_has_PROJ6()) {
    fun <- fm_caller_name(1, fun)
    msg <- paste0(
      "Call stack:\n",
      paste0(fm_call_stack(end = 1), collapse = "\n")
    )
    warning(paste0(
      "'",
      fun,
      "()' should not be used with PROJ6 and rgdal v3\n",
      msg
    ))
  }
}

#' @details `fm_not_for_PROJ4` is called to give an error when
#' calling methods that are only available for PROJ6
#' @rdname fm_has_PROJ6

fm_not_for_PROJ4 <- function(fun = NULL) {
  if (!fm_has_PROJ6()) {
    fun <- fm_caller_name(1, fun)
    stop(paste0(
      "'",
      fun,
      "()' is not supported for PROJ4"
    ))
  }
}

#' @details `fm_fallback_PROJ6` is called to warn about falling back
#' to using old PROJ4 methods when a PROJ6 method hasn't been implemented
#' @rdname fm_has_PROJ6

fm_fallback_PROJ6 <- function(fun = NULL) {
  if (fm_has_PROJ6()) {
    fun <- fm_caller_name(1, fun)
    warning(paste0(
      "'",
      fun,
      "()' method for PROJ6 not implemented. Falling back to PROJ4."
    ))
  }
}


#' @param fun The name of the function that requires PROJ6. Default: NULL,
#' which uses the name of the calling function.
#' @details `fm_requires_PROJ6` is called to give an error when PROJ6
#' is required but not available
#' @rdname fm_has_PROJ6

fm_requires_PROJ6 <- function(fun = NULL) {
  if (!fm_has_PROJ6()) {
    fun <- fm_caller_name(which = 1, override = fun)
    stop(paste0(
      "'",
      fun,
      "' requires PROJ6/RGDAL3"
    ))
  }
}



#' @title Extract CRS information
#' @description Wrapper for CRS(projargs) (PROJ4) and CRS(wkt) for
#' `sp::Spatial` objects.
#' @param x A `sp::Spatial` object
#' @return A `CRS` object, or NULL if no valid CRS identified
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @details This function is a convenience method to workaround PROJ4/PROJ6
#' differences, and the lack of a crs extraction method for Spatial objects.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   s <- sp::SpatialPoints(matrix(1:6, 3, 2), proj4string = fm_CRS("sphere"))
#'   fm_sp_get_crs(s)
#' }
#' }
#' @export
#' @rdname fm_sp_get_crs

fm_sp_get_crs <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (fm_has_PROJ6()) {
    suppressWarnings(crs <- sp::CRS(SRS_string = sp::wkt(x)))
  } else {
    crs <- sp::CRS(proj4string(x))
  }
  crs
}

fm_crs_is_null <- function(crs) {
  if (is.null(crs)) {
    TRUE
  } else if (fm_has_PROJ6()) {
    wkt <- fm_crs_get_wkt(crs)
    is.null(wkt)
  } else {
    #    is.na(crs)
    is.na(crs@projargs) && is.null(comment(crs))
  }
}

fm_ensure_crs <- function(crs) {
  if (fm_crs_is_null(crs)) {
    crs <- sp::CRS(NA_character_)
  }
  crs
}




#' @title Handling CRS/WKT
#' @description Get and set CRS object or WKT string properties.
#' @export
#' @rdname fm_crs_wkt

fm_wkt_is_geocent <- function(wkt) {
  if (is.null(wkt) || identical(wkt, "")) {
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
  if (fm_has_PROJ6()) {
    wkt <- fm_crs_get_wkt(crs)
    result <- fm_wkt_is_geocent(wkt)
  } else {
    args <- fm_CRS_as_list(crs)
    result <- identical(args[["proj"]], "geocent")
  }
  result
}


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
    stop("Ellipsoid settings not found")
  }

  datum <- fm_wkt_tree_get_item(wt, "DATUM")
  if (is.null(datum)) {
    stop("Ellipsoid settings not found")
  }
  ellipsoid <- fm_wkt_tree_get_item(datum, "ELLIPSOID")
  if (is.null(ellipsoid)) {
    stop("Ellipsoid settings not found")
  }
  as.numeric(ellipsoid[["params"]][[2]])
}

#' @rdname fm_crs_wkt
#' @export

fm_crs_get_ellipsoid_radius <- function(crs) {
  fm_requires_PROJ6()

  fm_wkt_get_ellipsoid_radius(fm_crs_get_wkt(crs))
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
      datum <- fm_wkt_tree_get_item(wt, "DATUM")
      null_datum <- is.null(datum)
      if (is.null(datum)) {
        stop("Ellipsoid settings not found")
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
#' @export

fm_crs_set_ellipsoid_radius <- function(crs, radius) {
  fm_requires_PROJ6()

  wkt <- fm_crs_get_wkt(crs)
  wkt <- fm_wkt_set_ellipsoid_radius(wkt, radius)
  new_crs <- fm_CRS(SRS_string = wkt)
  if (inherits(crs, "inla.CRS")) {
    crs$crs <- new_crs
    crs
  } else {
    new_crs
  }
}











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
#' @param doCheckCRSArgs default TRUE, must be set to FALSE by package
#' developers including `CRS` in an S4 class definition to avoid
#' uncontrollable loading of the `rgdal` namespace.
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
#'
#' if (require(rgdal)) {
#'   if (fm_has_PROJ6()) {
#'     crs1 <- fm_CRS("longlat_globe")
#'     crs2 <- fm_CRS("lambert_globe")
#'     crs3 <- fm_CRS("mollweide_norm")
#'     crs4 <- fm_CRS("hammer_globe")
#'     crs5 <- fm_CRS("sphere")
#'     crs6 <- fm_CRS("globe")
#'   } else {
#'     # Old definitions for pre-PROJ6:
#'     # Old radius-1 projections have a added "_norm" in the PROJ6 version of
#'     # the fm_CRS() predefined projections. They are detected and converted
#'     # to the new versions when RPOJ6 is available.
#'     crs1 <- fm_CRS("longlat") # PROJ6: longlat_norm
#'     crs2 <- fm_CRS("lambert") # PROJ6: lambert_norm
#'     crs3 <- fm_CRS("mollweide") # PROJ6: mollweide_norm
#'     crs4 <- fm_CRS("hammer") # PROJ6: hammer_norm
#'     crs5 <- fm_CRS("sphere")
#'     crs6 <- fm_CRS("globe")
#'   }
#' }
#' @export
#' @rdname fm_CRS

fm_CRS <- function(projargs = NULL, doCheckCRSArgs = TRUE,
                   args = NULL, oblique = NULL,
                   SRS_string = NULL,
                   ...) {
  if (fm_has_PROJ6()) {
    # PROJ6
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
        warning("'fm_CRS' should be given a SRS_string for PROJ6 or a known keyword for a predefined string given in projargs. Using fallback PROJ4 method.")
        if (!is.null(args)) {
          x <- sp::CRS(projargs, doCheckCRSArgs = doCheckCRSArgs)
          if (typeof(args) != "list") {
            stop("'args' must be NULL or a list of name=value pairs.")
          }
          xargs <- fm_CRS_as_list(x)
          for (name in names(args)) {
            xargs[[name]] <- args[[name]]
          }
          projargs <- fm_list_as_CRSargs(xargs)
        }
        SRS_string <- rgdal::showSRID(projargs, multiline = "NO")
        projargs <- NULL
      }
    }

    if (!is.null(SRS_string)) {
      if (!is.null(projargs)) {
        warning("SRS_string specified. Ignoring non-null projargs.")
      }
      x <- sp::CRS(SRS_string = SRS_string, doCheckCRSArgs = doCheckCRSArgs)
    } else if (inherits(projargs, "CRS")) {
      x <- projargs
    } else {
      x <- sp::CRS(NA_character_, doCheckCRSArgs = doCheckCRSArgs)
    }
  } else {
    # PROJ4
    if (is.null(projargs)) {
      projargs <- NA_character_
    }
    halfroot <- "+a=0.7071067811865476 +b=0.7071067811865476"
    predef <- list(
      hammer = paste("+proj=hammer +ellps=sphere +units=m", halfroot),
      lambert = "+proj=cea +ellps=sphere +lat_ts=0 +units=m +a=1 +b=1",
      longlat = "+proj=longlat +ellps=sphere +a=1 +b=1",
      mollweide = paste("+proj=moll +ellps=sphere +units=m", halfroot),
      sphere = "+proj=geocent +ellps=sphere +a=1 +b=1 +units=m",
      globe = "+proj=geocent +ellps=sphere +units=m"
    )
    if (is.character(projargs)) {
      if (projargs %in% names(predef)) {
        projargs <- predef[[projargs]]
      }
      x <- sp::CRS(projargs, doCheckCRSArgs = doCheckCRSArgs)
    } else if (inherits(projargs, "CRS")) {
      x <- projargs
    } else {
      stop(paste(
        "Unsupported projargs input class",
        paste(class(projargs), collapse = ",")
      ))
    }
    if (!is.null(args)) {
      if (typeof(args) != "list") {
        stop("'args' must be NULL or a list of name=value pairs.")
      }
      xargs <- fm_CRS_as_list(x)
      for (name in names(args)) {
        xargs[[name]] <- args[[name]]
      }
      x <- sp::CRS(fm_list_as_CRSargs(xargs), doCheckCRSArgs = doCheckCRSArgs)
    }
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
#' @rdname fm_CRS

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
  fm_requires_PROJ6()

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
  fm_requires_PROJ6()

  construct_item <- function(x, level) {
    paste0(
      if (pretty) {
        paste0(rep("    ", level), collapse = "")
      } else {
        ""
      },
      x[["label"]],
      "[",
      paste0(vapply(
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
  fm_not_for_PROJ6()
  fm_CRSargs_as_list(fm_CRSargs(x))
}


#' @export
#' @rdname fm_CRSargs
fm_list_as_CRS <- function(x, ...) {
  fm_not_for_PROJ6()
  fm_CRS(args = x)
}

#' Show expanded CRS arguments
#'
#' Wrappers for `sp::CRS` and `inla.CRS` objects to handle the
#' coordinate reference system argument string.
#' These methods should no longer be used with PROJ6/rgdal3;
#' see [fm_crs_get_wkt()] for a new approach.
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
#' @seealso [rgdal::CRSargs()], [fm_CRS()]
#' @export
#' @keywords internal
#' @examples
#'
#' if (require(rgdal) && !fm_has_PROJ6()) {
#'   crs0 <- fm_CRS("longlat")
#'   p4s <- fm_CRSargs(crs0)
#'   lst <- fm_CRSargs_as_list(p4s)
#'   crs1 <- fm_list_as_CRS(lst)
#'   lst$a <- 2
#'   crs2 <- fm_CRS(p4s, args = lst)
#'   print(fm_CRSargs(crs0))
#'   print(fm_CRSargs(crs1))
#'   print(fm_CRSargs(crs2))
#' }
fm_CRSargs <- function(x, ...) {
  fm_not_for_PROJ6()

  if (inherits(x, "inla.CRS")) {
    x <- x[["crs"]]
  }
  if (is.null(x)) {
    as.character(NA)
  } else {
    stopifnot(requireNamespace("rgdal", quietly = TRUE))
    rgdal::CRSargs(x)
  }
}


#' @return For `fm_list_as_CRSargs()`, a CRS proj4 string for name=value pair list
#' @rdname fm_CRSargs
fm_list_as_CRSargs <- function(x, ...) {
  paste(lapply(
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
  fm_not_for_PROJ6()

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




#' @param crs A `sp::CRS` or `inla.CRS` object
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
#' if (fm_has_PROJ6()) {
#'   c1 <- fm_CRS("globe")
#'   fm_crs_get_lengthunit(c1)
#'   c2 <- fm_crs_set_lengthunit(c1, "km")
#'   fm_crs_get_lengthunit(c2)
#' }
#' }
#' @export
#' @seealso [fm_sp_get_crs()]
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

#' @return For `fm_crs_get_wkt`, WKT2 string.
#' @export
#' @rdname fm_crs_wkt

fm_crs_get_wkt <- function(crs) {
  fm_requires_PROJ6()

  if (inherits(crs, "inla.CRS")) {
    crs <- crs[["crs"]]
  }

  if (is.null(crs)) {
    return(NULL)
  }

  comment(crs)
}

#' @return For `fm_crs_get_lengthunit`, a
#' list of length units used in the wkt string, excluding the ellipsoid radius
#' unit. (For legacy PROJ4 code, the raw units from the proj4string are
#' returned, if present.)
#' @export
#' @rdname fm_crs_wkt

fm_crs_get_lengthunit <- function(crs) {
  if (fm_has_PROJ6()) {
    x <- fm_wkt_get_lengthunit(fm_crs_get_wkt(crs))
  } else {
    if (inherits(crs, "inla.CRS")) {
      crs_ <- crs
      crs <- crs[["crs"]]
    } else {
      crs_ <- NULL
    }

    crs_args <- fm_CRS_as_list(crs)
    x <- crs_args[["units"]]
  }
  x
}

#' @return For `fm_crs_set_lengthunit`, a `sp::CRS` object with
#' altered length units.
#' Note that the length unit for the ellipsoid radius is unchanged.
#' @export
#' @rdname fm_crs_wkt

fm_crs_set_lengthunit <- function(crs, unit, params = NULL) {
  if (inherits(crs, "inla.CRS")) {
    crs_ <- crs
    crs <- crs[["crs"]]
  } else {
    crs_ <- NULL
  }
  if (fm_has_PROJ6()) {
    x <- sp::CRS(SRS_string = fm_wkt_set_lengthunit(fm_crs_get_wkt(crs),
      unit,
      params = params
    ))
  } else {
    crs_args <- fm_CRS_as_list(crs)
    crs_args[["units"]] <- unit
    x <- fm_list_as_CRS(crs_args)
  }
  if (!is.null(crs_)) {
    crs_[["crs"]] <- x
    x <- crs_
  }
  x
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
  fm_requires_PROJ6()
  wt <- fm_wkt_as_wkt_tree(wkt)
  fm_wkt_tree_projection_type(wt)
}

fm_crs_projection_type <- function(crs) {
  fm_requires_PROJ6()
  wkt <- fm_crs_get_wkt(crs)
  fm_wkt_projection_type(wkt)
}

## +proj=longlat in (-180,180)x(-90,90)
## +proj=moll in (-2,2)x(-1,1) scaled by +a and +b, and +units
## +proj=lambert in (-pi,pi)x(-1,1) scaled by +a and +b, and +units
fm_crs_bounds <- function(crs, warn.unknown = FALSE) {
  if (fm_has_PROJ6()) {
    # PROJ6

    wkt <- fm_crs_get_wkt(crs)
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
  } else {
    # PROJ4
    args <- fm_CRS_as_list(crs)
    if (args[["proj"]] == "longlat") {
      bounds <- list(type = "rectangle", xlim = c(-180, 180), ylim = c(-90, 90))
    } else if (args[["proj"]] == "cea") {
      axis <- c(pi, 1)
      if (!is.null(args[["a"]])) {
        axis[1] <- axis[1] * as.numeric(args$a)
      }
      if (!is.null(args[["b"]])) {
        axis[2] <- axis[2] * as.numeric(args$a)^0.5 * as.numeric(args$b)^0.5
      }
      ## TODO: Handle "lat_ts" and "units"
      bounds <- list(
        type = "rectangle",
        xlim = c(-1, 1) * axis[1], ylim = c(-1, 1) * axis[2]
      )
    } else if (args[["proj"]] %in% c("moll", "hammer")) {
      axis <- c(2, 1)
      center <- c(0, 0)
      if (!is.null(args[["a"]])) {
        axis[1] <- axis[1] * as.numeric(args$a) / sqrt(1 / 2)
        axis[2] <- axis[2] * as.numeric(args$a) / sqrt(1 / 2)
      }
      ## TODO: Handle "units"
      bounds <- list(
        type = "ellipse", axis = axis, center = center,
        xlim = center[1] + c(-1, 1) * axis[1],
        ylim = center[2] + c(-1, 1) * axis[2]
      )
    } else if (args[["proj"]] == "tmerc") {
      bounds <- list(type = "rectangle", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf))
    } else if (fm_crs_is_geocent(crs)) {
      bounds <- list(type = "rectangle", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf))
    } else {
      if (warn.unknown) {
        warning("Could not determine shape of transformation bounds. Using infinite rectangle.")
      }
      bounds <- list(type = "rectangle", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf))
    }
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


#' @title Chack if two CRS objects are identical
#' @param crs0,crs1 Two `sp::CRS` or `inla.CRS` objects to be compared.
#' @param crsonly logical. If `TRUE` and any of `crs0` and `crs1` are `inla.CRS`
#' objects, extract and compare only the `sp::CRS` objects. Default: `FALSE`
#' @export
#' @keywords internal

fm_identical_CRS <- function(crs0, crs1, crsonly = FALSE) {
  if (!crsonly) {
    identical(crs0, crs1)
  } else {
    if (inherits(crs0, "inla.CRS")) {
      crs0 <- crs0$crs
    }
    if (inherits(crs1, "inla.CRS")) {
      crs1 <- crs1$crs
    }
    identical(crs0, crs1)
  }
}

#' Handle transformation of various inla objects according to coordinate
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
#' @seealso [fm_CRS()]
#' @export
fm_spTransform <- function(x, ...) {
  UseMethod("fm_spTransform")
}

#' @details The default method handles low level transformation of raw
#' coordinates.
#' @export
#' @rdname fm_spTransform
fm_spTransform.default <- function(x, crs0, crs1, passthrough = FALSE, ...) {
  if (fm_has_PROJ6()) {
    # PROJ6
    ok0 <- (!is.null(crs0) &&
      ((inherits(crs0, "CRS") && !is.null(fm_crs_get_wkt(crs0))) ||
        (inherits(crs0, "inla.CRS"))))
    ok1 <- (!is.null(crs1) &&
      ((inherits(crs1, "CRS") && !is.null(fm_crs_get_wkt(crs1))) ||
        (inherits(crs1, "inla.CRS"))))
    if (ok0 && ok1) {
      if (ncol(x) == 2) {
        x <- cbind(x, 0)
      }
      sphere_radius_0 <- fm_crs_get_ellipsoid_radius(crs0)
      sphere_radius_1 <- fm_crs_get_ellipsoid_radius(crs1)
      different_radii <- (sphere_radius_0 != sphere_radius_1)
      longlat_norm <- fm_CRS("longlat_norm")
      longlat_0 <- fm_crs_set_ellipsoid_radius(longlat_norm, sphere_radius_0)
      longlat_1 <- fm_crs_set_ellipsoid_radius(longlat_norm, sphere_radius_1)

      crs_sphere <- fm_CRS("sphere")
      onsphere_0 <- fm_identical_CRS(crs0, crs_sphere, crsonly = TRUE)
      onsphere_1 <- fm_identical_CRS(crs1, crs_sphere, crsonly = TRUE)
      is_geocentric_0 <- fm_crs_is_geocent(crs0)
      is_geocentric_1 <- fm_crs_is_geocent(crs1)
      if (is_geocentric_0) {
        ok <- TRUE
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
        inherits(crs0, "inla.CRS") ||
          inherits(crs1, "inla.CRS") ||
          different_radii
      if (inherits(crs0, "inla.CRS")) {
        crs0crs <- crs0$crs
        crs0oblique <- crs0$oblique
      } else {
        crs0crs <- crs0
        crs0oblique <- NULL
      }
      if (inherits(crs1, "inla.CRS")) {
        crs1crs <- crs1$crs
        crs1oblique <- crs1$oblique
      } else {
        crs1crs <- crs1
        crs1oblique <- NULL
      }
      x <- sp::SpatialPoints(x[ok, , drop = FALSE], proj4string = crs0crs)
      if (do_work_on_sphere) {
        if (!onsphere_0) {
          if (sphere_radius_0 != 1) {
            x <- sp::spTransform(x, longlat_0)
            proj4string(x) <- sp::CRS(NA_character_) # Reset CRS to avoid warning
            proj4string(x) <- longlat_norm
          }
          x <- sp::spTransform(x, crs_sphere)
        }
        if (!is.null(crs0oblique)) {
          x <- sp::SpatialPoints(fm_crs_transform_oblique(coordinates(x),
            crs0oblique,
            to.oblique = FALSE
          ),
          proj4string = crs_sphere
          )
        }

        if (!is.null(crs1oblique)) {
          x <- sp::SpatialPoints(fm_crs_transform_oblique(coordinates(x),
            crs1oblique,
            to.oblique = TRUE
          ),
          proj4string = crs_sphere
          )
        }
        if (sphere_radius_1 != 1) {
          x <- sp::spTransform(x, longlat_norm)
          proj4string(x) <- sp::CRS(NA_character_) # Reset CRS to avoid warning
          proj4string(x) <- longlat_1
        }
      }

      x <- sp::spTransform(x, crs1crs)

      if (!all(ok)) {
        xx[ok, ] <- coordinates(x)
        xx[!ok, ] <- NA
        x <- xx
      }
    } else if (!passthrough) {
      if (!ok0) {
        stop("'crs0' is an invalid coordinate reference object.")
      }
      if (!ok1) {
        stop("'crs1' is an invalid coordinate reference object.")
      }
    }
    if (is.matrix(x)) {
      invisible(x)
    } else {
      invisible(coordinates(x))
    }
  } else {
    # PROJ4
    ok0 <- (!is.null(crs0) &&
      ((inherits(crs0, "CRS") && !is.na(fm_CRSargs(crs0))) ||
        (inherits(crs0, "inla.CRS"))))
    ok1 <- (!is.null(crs1) &&
      ((inherits(crs1, "CRS") && !is.na(fm_CRSargs(crs1))) ||
        (inherits(crs1, "inla.CRS"))))
    if (ok0 && ok1) {
      if (ncol(x) == 2) {
        x <- cbind(x, 0)
      }
      onsphere <- fm_identical_CRS(crs0, fm_CRS("sphere"), crsonly = TRUE)
      isgeocentric <- identical(fm_CRS_as_list(crs0)[["proj"]], "geocent")
      if (isgeocentric) {
        ok <- TRUE
      } else {
        bounds <- fm_crs_bounds(crs0)
        if (identical(fm_CRS_as_list(crs0)[["proj"]], "longlat")) {
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
      if (inherits(crs0, "inla.CRS")) {
        if (!onsphere) {
          x <- sp::spTransform(
            sp::SpatialPoints(x[ok, , drop = FALSE], proj4string = crs0$crs),
            fm_CRS("sphere")
          )
        }
        if (!is.null(crs0$oblique)) {
          x <- sp::SpatialPoints(fm_crs_transform_oblique(coordinates(x),
            crs0$oblique,
            to.oblique = FALSE
          ),
          proj4string = fm_CRS("sphere")
          )
        }
        onshpere <- TRUE
      } else {
        x <- sp::SpatialPoints(x[ok, , drop = FALSE], proj4string = crs0)
      }
      if (inherits(crs1, "inla.CRS")) {
        if (!onsphere) {
          x <- sp::spTransform(x, fm_CRS("sphere"))
        }
        if (!is.null(crs1$oblique)) {
          x <- sp::SpatialPoints(fm_crs_transform_oblique(coordinates(x),
            crs1$oblique,
            to.oblique = TRUE
          ),
          proj4string = fm_CRS("sphere")
          )
        }
        x <- sp::spTransform(x, crs1$crs)
      } else {
        x <- sp::spTransform(x, crs1)
      }
      if (!all(ok)) {
        xx[ok, ] <- coordinates(x)
        xx[!ok, ] <- NA
        x <- xx
      }
    } else if (!passthrough) {
      if (!ok0) {
        stop("'crs0' is an invalid coordinate reference object.")
      }
      if (!ok1) {
        stop("'crs1' is an invalid coordinate reference object.")
      }
    }
    if (is.matrix(x)) {
      invisible(x)
    } else {
      invisible(coordinates(x))
    }
  }
}

#' @export
#' @rdname fm_spTransform
fm_spTransform.SpatialPoints <- function(x, CRSobj, passthrough = FALSE, ...) {
  if (fm_has_PROJ6()) {
    crs_x <- fm_sp_get_crs(x)
    ok0 <- !is.null(fm_crs_get_wkt(crs_x))
    ok1 <- (!missing(CRSobj) && !is.null(CRSobj) &&
      (inherits(CRSobj, "CRS") && !is.null(fm_crs_get_wkt(CRSobj))))
    if (ok0 && ok1) {
      invisible(sp::SpatialPoints(fm_spTransform(
        coordinates(x),
        crs_x,
        CRSobj
      ),
      proj4string = CRSobj
      ))
    } else if (ok1) { ## Know: !ok0 && ok1
      if (!passthrough) {
        stop("Invalid origin CRS for sp::SpatialPoints")
      }
      invisible(sp::SpatialPoints(coordinates(x), proj4string = CRSobj))
    } else { ## Know: (ok0 || !ok0) && !ok1
      if (!passthrough) {
        stop("Invalid target CRS for sp::SpatialPoints")
      }
      invisible(sp::SpatialPoints(coordinates(x), proj4string = fm_CRS()))
    }
  } else {
    # PROJ4
    ok0 <- !is.na(proj4string(x))
    ok1 <- (!missing(CRSobj) && !is.null(CRSobj) &&
      (inherits(CRSobj, "CRS") && !is.na(fm_CRSargs(CRSobj))))
    if (ok0 && ok1) {
      invisible(sp::SpatialPoints(fm_spTransform(
        coordinates(x),
        sp::CRS(sp::proj4string(x)),
        CRSobj
      ),
      proj4string = CRSobj
      ))
    } else if (ok1) { ## Know: !ok0 && ok1
      if (!passthrough) {
        stop("Invalid origin CRS for sp::SpatialPoints")
      }
      invisible(sp::SpatialPoints(coordinates(x), proj4string = CRSobj))
    } else { ## Know: (ok0 || !ok0) && !ok1
      if (!passthrough) {
        stop("Invalid target CRS for sp::SpatialPoints")
      }
      invisible(sp::SpatialPoints(coordinates(x), proj4string = fm_CRS()))
    }
  }
}

#' @export
#' @rdname fm_spTransform
fm_spTransform.SpatialPointsDataFrame <- function(x,
                                                  CRSobj,
                                                  passthrough = FALSE,
                                                  ...) {
  fm_requires_PROJ6()

  ok1 <- (!missing(CRSobj) && !is.null(CRSobj) &&
    (inherits(CRSobj, "CRS") && !is.null(fm_crs_get_wkt(CRSobj))))
  if (!ok1 && !passthrough) {
    stop("Invalid target CRS for sp::SpatialPointsDataFrame")
  }

  x_no_df <- sp::SpatialPoints(coordinates(x),
    proj4string = fm_sp_get_crs(x)
  )
  x_no_df <- fm_spTransform(x_no_df,
    CRSobj = CRSobj,
    passthrough = passthrough
  )
  if (ok1) {
    invisible(sp::SpatialPointsDataFrame(coordinates(x_no_df),
      proj4string = CRSobj,
      data = x@data
    ))
  } else {
    invisible(sp::SpatialPointsDataFrame(coordinates(x_no_df),
      proj4string = fm_CRS(),
      data = x@data
    ))
  }
}

#' @export
#' @rdname fm_spTransform
fm_spTransform.inla.mesh.lattice <- function(x, CRSobj, passthrough = FALSE, ...) {
  x$segm <- fm_spTransform(x$segm, CRSobj, passthrough = passthrough)
  x$loc <- fm_spTransform(x$loc, x$crs, CRSobj, passthrough = passthrough)
  x$crs <- CRSobj
  invisible(x)
}

#' @export
#' @rdname fm_spTransform
fm_spTransform.inla.mesh.segment <- function(x, CRSobj, passthrough = FALSE, ...) {
  x$loc <- fm_spTransform(x$loc, x$crs, CRSobj, passthrough = passthrough)
  x$crs <- CRSobj
  invisible(x)
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
#' @rdname fm_spTransform
fm_spTransform.inla.mesh <- function(x, CRSobj, passthrough = FALSE, ...) {
  x$loc <- fm_spTransform(x$loc, x$crs, CRSobj, passthrough = passthrough)
  x$manifold <- fm_crs_detect_manifold(CRSobj)
  x$crs <- CRSobj
  invisible(x)
}


## Input: list of segments, all closed polygons.
fm_internal_sp2segment_join <- function(inp, grp = NULL, closed = TRUE) {
  crs <- NULL
  if (length(inp) > 0) {
    out.loc <- matrix(0, 0, ncol(inp[[1]]$loc))
    for (k in seq_along(inp)) {
      crs <- fm_internal_update_crs(crs, inp[[k]]$crs, mismatch.allowed = FALSE)
    }
  } else {
    out.loc <- matrix(0, 0, 2)
  }
  out.idx <- matrix(0L, 0, 2)
  if (is.null(grp)) {
    out.grp <- NULL
  } else {
    out.grp <- integer(0)
  }
  for (k in seq_along(inp)) {
    inp.loc <- inp[[k]]$loc
    inp.idx <- inp[[k]]$idx
    inp.grp <- inp[[k]]$grp
    offset <- nrow(out.loc)
    n <- nrow(as.matrix(inp.idx))
    if (closed) {
      if (!is.null(grp) && is.null(inp.grp)) {
        inp.grp <- rep(grp[k], n)
      }
      if (ncol(as.matrix(inp.idx)) == 1) {
        inp.idx <- cbind(inp.idx, inp.idx[c(2:n, 1)])
      }
    } else {
      if (!is.null(grp) && is.null(inp.grp)) {
        inp.grp <- rep(grp[k], n - 1)
      }
      if (ncol(as.matrix(inp.idx)) == 1) {
        inp.idx <- cbind(inp.idx[-n], inp.idx[-1])
      }
    }
    out.loc <- rbind(out.loc, inp.loc)
    out.idx <- rbind(out.idx, inp.idx + offset)
    if (!is.null(grp)) {
      out.grp <- c(out.grp, inp.grp)
    }
  }
  INLA::inla.mesh.segment(
    loc = out.loc, idx = out.idx, grp = out.grp, is.bnd = FALSE,
    crs = crs
  )
}


fm_as_inla_mesh_segment <-
  function(sp, ...) {
    UseMethod("fm_as_inla_mesh_segment")
  }

fm_sp2segment <-
  function(sp, ...) {
    UseMethod("fm_as_inla_mesh_segment")
  }



fm_as_inla_mesh_segment.SpatialPoints <-
  function(sp, reverse = FALSE, grp = NULL, is.bnd = TRUE, ...) {
    crs <- fm_sp_get_crs(sp)
    loc <- coordinates(sp)

    n <- dim(loc)[1L]
    if (reverse) {
      idx <- seq(n, 1L, length = n)
    } else {
      idx <- seq_len(n)
    }
    INLA::inla.mesh.segment(
      loc = loc, idx = idx, grp = grp, is.bnd = is.bnd,
      crs = crs
    )
  }

fm_as_inla_mesh_segment.SpatialPointsDataFrame <-
  function(sp, ...) {
    fm_as_inla_mesh_segment.SpatialLines(sp, ...)
  }



fm_as_inla_mesh_segment.Line <-
  function(sp, reverse = FALSE, crs = NULL, ...) {
    loc <- sp@coords
    n <- dim(loc)[1L]
    if (reverse) {
      idx <- seq(n, 1L, length = n)
    } else {
      idx <- seq_len(n)
    }
    INLA::inla.mesh.segment(loc = loc, idx = idx, is.bnd = FALSE, crs = crs)
  }

fm_as_inla_mesh_segment.Lines <-
  function(sp, join = TRUE, crs = NULL, ...) {
    segm <- as.list(lapply(
      sp@Lines,
      function(x) fm_as_inla_mesh_segment(x, crs = crs, ...)
    ))
    if (join) {
      segm <- fm_internal_sp2segment_join(segm, grp = NULL, closed = FALSE)
    }
    segm
  }

fm_as_inla_mesh_segment.SpatialLines <-
  function(sp, join = TRUE, grp = NULL, ...) {
    crs <- fm_sp_get_crs(sp)
    segm <- list()
    for (k in seq_len(length(sp@lines))) {
      segm[[k]] <- fm_as_inla_mesh_segment(sp@lines[[k]],
        join = TRUE,
        crs = crs, ...
      )
    }
    if (join) {
      if (missing(grp)) {
        grp <- seq_len(length(segm))
      }
      segm <- fm_internal_sp2segment_join(segm, grp = grp, closed = FALSE)
    }
    segm
  }

fm_as_inla_mesh_segment.SpatialLinesDataFrame <-
  function(sp, ...) {
    fm_as_inla_mesh_segment.SpatialLines(sp, ...)
  }

fm_as_inla_mesh_segment.SpatialPolygons <-
  function(sp, join = TRUE, grp = NULL, ...) {
    crs <- fm_sp_get_crs(sp)
    segm <- list()
    for (k in seq_len(length(sp@polygons))) {
      segm[[k]] <- fm_as_inla_mesh_segment(sp@polygons[[k]], join = TRUE, crs = crs)
    }
    if (join) {
      if (missing(grp)) {
        grp <- seq_len(length(segm))
      }
      segm <- fm_internal_sp2segment_join(segm, grp = grp)
    }
    segm
  }

fm_as_inla_mesh_segment.SpatialPolygonsDataFrame <-
  function(sp, ...) {
    fm_as_inla_mesh_segment.SpatialPolygons(sp, ...)
  }

fm_as_inla_mesh_segment.Polygons <-
  function(sp, join = TRUE, crs = NULL, ...) {
    segm <- as.list(lapply(
      sp@Polygons,
      function(x) fm_as_inla_mesh_segment(x, crs = crs)
    ))
    if (join) {
      segm <- fm_internal_sp2segment_join(segm, grp = NULL)
    }
    segm
  }

fm_as_inla_mesh_segment.Polygon <-
  function(sp, crs = NULL, ...) {
    loc <- sp@coords[-dim(sp@coords)[1L], , drop = FALSE]
    n <- dim(loc)[1L]
    if (sp@hole) {
      if (sp@ringDir == 1) {
        idx <- c(1L:n, 1L)
      } else {
        idx <- c(1L, seq(n, 1L, length.out = n))
      }
    } else
    if (sp@ringDir == 1) {
      idx <- c(1L, seq(n, 1L, length.out = n))
    } else {
      idx <- c(1L:n, 1L)
    }
    INLA::inla.mesh.segment(loc = loc, idx = idx, is.bnd = TRUE, crs = crs)
  }
