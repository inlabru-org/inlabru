#' @title Generate integration points
#'
#' @description
#' `r lifecycle::badge("deprecated")` in favour of [fmesher::fm_int()]
#'
#' This function generates points in one or two dimensions with a weight attached to each point.
#' The weighted sum of a function evaluated at these points is the integral of that function approximated
#' by linear basis functions. The parameter `samplers` describes the area(s) integrated over.
#'
#' In case of a single dimension `samplers` is supposed to be a two-column `matrix` where
#' each row describes the start and end points of the interval to integrate over. In the two-dimensional
#' case `samplers` can be either a `SpatialPolygon`, an `inla.mesh` or a
#' `SpatialLinesDataFrame` describing the area to integrate over. If a `SpatialLineDataFrame`
#' is provided, it has to have a column called 'weight' in order to indicate the width of the line.
#'
#' The domain parameter is an `inla.mesh.1d` or `inla.mesh` object that can be employed to
#' project the integration points to the vertices of the mesh. This reduces the final number of
#' integration points and reduces the computational cost of the integration. The projection can also
#' prevent numerical issues in spatial LGCP models where each observed point is ideally surrounded
#' by three integration point sitting at the corresponding mesh vertices. This is controlled
#' by `int.args$method="stable"` (default) or `"direct"`, where the latter uses the integration
#' points directly, without aggregating to the mesh vertices.
#'
#' For convenience, the
#' `domain` parameter can also be a single integer setting the number of equally spaced integration
#' points in the one-dimensional case.
#'
#' @aliases ipoints
#' @export
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com} and
#' \email{finn.lindgren@@gmail.com}
#'
#' @param samplers Description of the integration region boundary.
#' In 1D, a length 2 vector or two-column matrix where each row describes an interval,
#' or `NULL`
#' In 2D either a `SpatialPolygon` or a `SpatialLinesDataFrame` with a `weight` column
#' defining the width of the a transect line, and optionally further columns used by the
#' `group` argument, or `NULL`.  When `domain` is `NULL`, `samplers` may also
#' be an `inla.mesh.1d` or `inla.mesh` object, that is then treated as a `domain`
#' argument instead.
#' @param domain Either
#' * when `samplers` is a 1D interval(s) definition only, `domain` can be
#'   a single integer for the number of integration points to place in each 1D
#'   interval, overriding `int.args[["nsub1"]]`, and otherwise
#' * when `samplers` is `NULL`, `domain` can be a numeric vector of points,
#'   each given integration weight 1 (and no additional points are added
#'   in between),
#' * an `inla.mesh.1d` object for continuous 1D integration, or
#' * an `inla.mesh.2d` object for continuous 2D integration.
#' @param name Character array stating the name of the domains dimension(s).
#' If `NULL`, the names are taken from coordinate names from `samplers` for
#' `Spatial*` objects, otherwise "x", "y", "z" for 2D regions and
#' `"x"` for 1D regions
#' @param group Column names of the `samplers` object (if applicable) for which
#' the integration points are calculated independently and not merged when
#' aggregating to mesh nodes.
#' @param int.args List of arguments passed to `bru_int_polygon`.
#' * `method`: "stable" (to aggregate integration weights onto mesh nodes)
#'   or "direct" (to construct a within triangle/segment integration scheme
#'   without aggregating onto mesh nodes)
#' * `nsub1`, `nsub2`: integers controlling the number of internal integration
#'   points before aggregation. Points per triangle: `(nsub2+1)^2`.
#'   Points per knot segment: `nsub1`
#' * `poly_method`: if set to "legacy", selects an old polygon integration method
#'   that doesn't handle holes. No longer supported, and will generate an error.
#' @param project `r lifecycle::badge("deprecated")` Deprecated in favour of `int.args(method=...)`.
#' If TRUE, aggregate the integration points to mesh vertices. Default:
#' `project = (identical(int.args$method, "stable"))`
#'
#' @return A `data.frame`, `tibble`, `sf`, or `SpatialPointsDataFrame` of 1D and
#' 2D integration points, including a `weight` column and `.block` column.
#'
#' @examples
#' \donttest{
#' if (require("INLA", quietly = TRUE) &&
#'   require("ggplot2", quietly = TRUE) &&
#'   bru_safe_sp() &&
#'   require("sp") &&
#'   require("fmesher")) {
#'   # Create 50 integration points covering the dimension 'myDim' between 0 and 10.
#'
#'   ips <- ipoints(c(0, 10), 50, name = "myDim")
#'   plot(ips)
#'
#'
#'   # Create integration points for the two intervals [0,3] and [5,10]
#'
#'   ips <- ipoints(matrix(c(0, 3, 5, 10), nrow = 2, byrow = TRUE), 50)
#'   plot(ips)
#'
#'
#'   # Convert a 1D mesh into integration points
#'   mesh <- fm_mesh_1d(seq(0, 10, by = 1))
#'   ips <- ipoints(mesh, name = "time")
#'   plot(ips)
#'
#'
#'   # Obtain 2D integration points from a SpatialPolygon
#'
#'   data(gorillas, package = "inlabru")
#'   ips <- ipoints(gorillas$boundary)
#'   ggplot() +
#'     gg(gorillas$boundary) +
#'     gg(ips, aes(size = weight))
#'
#'
#'   #' Project integration points to mesh vertices
#'
#'   ips <- ipoints(gorillas$boundary, domain = gorillas$mesh)
#'   ggplot() +
#'     gg(gorillas$mesh) +
#'     gg(gorillas$boundary) +
#'     gg(ips, aes(size = weight))
#'
#'
#'   # Turn a 2D mesh into integration points
#'
#'   ips <- ipoints(gorillas$mesh)
#'   ggplot() +
#'     gg(gorillas$boundary) +
#'     gg(ips, aes(size = weight))
#' }
#' }
#'
#' @importFrom sp coordnames coordinates
ipoints <- function(samplers = NULL, domain = NULL, name = NULL, group = NULL,
                    int.args = NULL,
                    project = deprecated()) {
  lifecycle::deprecate_warn(
    "2.8.0.9004",
    "ipoints()",
    "fmesher::fm_int()",
    details = c(
      "`ipoints(samplers, domain)` has been replaced by more versatile `fm_int(domain, samplers, ...)` methods."
    )
  )

  if (!is.null(group)) {
    if (is.null(name) && inherits(domain, c("fm_mesh_2d", "inla.mesh"))) {
      if (inherits(samplers, "sf")) {
        name <- "geometry"
      } else {
        name <- "coordinates"
      }
    }
    domain <- c(list(domain), as.list(as.data.frame(samplers[group]))[group])
    names(domain) <- c(name, group)
    ips <- fm_int(
      domain = domain,
      samplers = samplers,
      int.args = int.args
    )
    return(ips)
  }


  int.args.default <- list(method = "stable", nsub1 = 30, nsub2 = 9)
  if (is.null(int.args)) {
    int.args <- list()
  }
  missing.args <- setdiff(names(int.args.default), names(int.args))
  int.args[missing.args] <- int.args.default[missing.args]
  if (!is.null(int.args[["nsub"]])) {
    int.args[["nsub1"]] <- int.args[["nsub"]]
  }
  if (!is.null(int.args[["nsub"]])) {
    int.args[["nsub2"]] <- int.args[["nsub"]]
  }

  if (lifecycle::is_present(project)) {
    lifecycle::deprecate_warn(
      "2.7.0",
      "ipoints(project)",
      details = paste0(
        "For project=", project,
        ", use int.args$method = '",
        list(
          "TRUE" = "stable",
          "FALSE" = "direct"
        )[as.character(project)],
        "' instead."
      )
    )
    project <- NULL
  }

  if (is.null(domain) &&
    inherits(samplers, c("fm_mesh_1d", "fm_mesh_2d", "inla.mesh.1d", "inla.mesh"))) {
    domain <- samplers
    samplers <- NULL
  }

  samplers_is_sf <- inherits(samplers, c("sf", "sfc"))
  if (samplers_is_sf) {
    samplers <- as(samplers, "Spatial")
  }

  is_2d <-
    (
      !is.null(samplers) &&
        inherits(
          samplers,
          c(
            "SpatialPoints",
            "SpatialPointsDataFrame",
            "SpatialPolygons",
            "SpatialPolygonsDataFrame",
            "SpatialLines",
            "SpatialLinesDataFrame"
          )
        )
    ) ||
      inherits(domain, c("fm_mesh_2d", "inla.mesh"))
  is_1d <- !is_2d &&
    (
      (!is.null(samplers) &&
        (is.numeric(samplers) ||
          is.character(samplers) ||
          is.factor(samplers)) ||
        (!is.null(domain) &&
          (is.numeric(domain) ||
            is.character(domain) ||
            is.factor(domain)) ||
          inherits(domain, c("fm_mesh_1d", "inla.mesh.1d"))))
    )
  if (!is_1d && !is_2d) {
    stop("Unable to determine integration domain definition")
  }

  if (is_1d && !is.null(samplers) && !is.null(domain) &&
    is.numeric(domain) && length(domain) == 1) {
    int.args[["nsub1"]] <- domain
    domain <- NULL
    int.args[["method"]] <- "direct"
  }
  if (is_2d && !is.null(samplers) && !is.null(domain) &&
    is.numeric(domain) && length(domain) == 1) {
    int.args[["nsub2"]] <- domain
    domain <- NULL
    int.args[["method"]] <- "direct"
  }

  # Do this check and transfer once more
  if (is.null(domain) &&
    inherits(
      samplers,
      c("fm_mesh_1d", "fm_mesh_2d", "inla.mesh.1d", "inla.mesh")
    )) {
    domain <- samplers
    samplers <- NULL
  }

  if (is_1d && is.null(name)) {
    name <- "x"
  }

  if (is.data.frame(samplers)) {
    if (!("weight" %in% names(samplers))) {
      samplers$weight <- 1
    }
    ips <- samplers
  } else if (is_1d && is.null(samplers) && (is.numeric(domain) ||
    is.character(domain) ||
    is.factor(domain))) {
    ips <- data.frame(
      x = as.vector(domain),
      weight = 1
    )
    colnames(ips) <- c(name, "weight")
    # TODO samplers numeric and you do not have domain, it is not safe to forget to supply domain, if stored in integer
  } else if (is_1d && is.null(domain) && (is.integer(samplers) ||
    is.character(samplers) ||
    is.factor(samplers))) {
    ips <- data.frame(
      x = as.vector(samplers),
      weight = 1
    )
    colnames(ips) <- c(name, "weight")
  } else if (is_1d &&
    is.null(samplers) &&
    inherits(domain, c("fm_mesh_1d", "inla.mesh.1d")) &&
    identical(int.args[["method"]], "stable")) {
    if (isTRUE(domain$cyclic)) {
      loc_trap <- c(domain$loc, domain$interval[2])
    } else {
      loc_trap <- domain$loc
    }
    loc_mid <- (loc_trap[-1] + loc_trap[-length(loc_trap)]) / 2
    weight_mid <- diff(loc_trap)
    weight_trap <- c(weight_mid / 2, 0) + c(0, weight_mid / 2)
    loc_simpson <- c(loc_trap, loc_mid)
    weight_simpson <- c(weight_trap / 3, weight_mid * 2 / 3)

    ips <- data.frame(
      x = loc_simpson,
      weight = weight_simpson,
      .block = 1
    )
    colnames(ips) <- c(name, "weight", ".block")
  } else if (is_1d) {
    domain_range <-
      if (inherits(domain, c("fm_mesh_1d", "inla.mesh.1d"))) {
        domain$interval
      } else {
        NULL
      }

    if (is.null(samplers)) {
      samplers <- matrix(domain_range, 1, 2)
    } else {
      if (is.null(dim(samplers))) {
        samplers <- matrix(samplers, nrow = 1)
      }
      if (ncol(samplers) != 2) {
        stop("Interval description matrix must have 2 elements or be a 2-column matrix.")
      }
      if (is.null(domain)) {
        domain <- fm_mesh_1d(sort(unique(as.vector(samplers))))
      }
    }
    domain <- fm_as_mesh_1d(domain)
    # Now samplers is a 2-column matrix, and domain is an `fm_mesh_1d` object.

    ips <- fmesher::fm_int(
      domain = domain,
      samplers = samplers,
      name = name,
      int.args = int.args
    )
  } else if (inherits(domain, c("fm_mesh_2d", "inla.mesh")) &&
    is.null(samplers) &&
    identical(int.args[["method"]], "stable")) {
    coord_names <- c("x", "y", "z")
    if (!is.null(name)) {
      coord_names[seq_along(name)] <- name
    }

    # transform to equal area projection
    if (!fm_crs_is_null(fm_crs(domain))) {
      crs <- fm_crs(domain)
      samplers <- fm_transform(domain, crs = fm_crs("+proj=cea +units=km"))
    }

    ips <- fm_vertices(domain, format = "sp")
    ips$weight <- fm_fem(domain, order = 1)$va

    # backtransform
    if (!fm_crs_is_null(fm_crs(domain))) {
      ips <- fm_transform(ips, crs = crs)
    }
    sp::coordnames(ips) <- coord_names[seq_len(NCOL(sp::coordinates(ips)))]
  } else if (inherits(samplers, "SpatialPointsDataFrame")) {
    if (!("weight" %in% names(samplers))) {
      warning("The integration points provided have no weight column. Setting weights to 1.")
      samplers$weight <- 1
    }

    ips <- samplers
  } else if (inherits(samplers, "SpatialPoints")) {
    ips <- samplers
    ips$weight <- 1
  } else if (inherits(samplers, "SpatialLines") ||
    inherits(samplers, "SpatialLinesDataFrame")) {
    # If SpatialLines are provided convert into SpatialLinesDataFrame and attach weight = 1
    if (inherits(samplers, "SpatialLines") &&
      !inherits(samplers, "SpatialLinesDataFrame")) {
      samplers <- SpatialLinesDataFrame(
        samplers,
        data = data.frame(weight = rep(1, length(samplers))),
        match.ID = FALSE
      )
    }

    # Set weight to 1 if not provided
    if (!("weight" %in% names(samplers))) {
      samplers$weight <- 1
    }

    the_domain <- list(coordinates = domain)
    if (!is.null(group)) {
      for (grp in group) {
        the_domain[[grp]] <- sort(unique(samplers[[grp]]))
      }
    }
    ips <- fm_int(domain = the_domain, samplers = samplers, int.args = int.args)

    coord_names <- c("x", "y", "z")
    #    if (!is.null(coordnames(samplers))) {
    #      coord_names[seq_along(coordnames(samplers))] <- coordnames(samplers)
    #    } else if (!is.null(name)) {
    #      coord_names[seq_along(name)] <- name
    #    }
    sp::coordnames(ips) <- coord_names[seq_len(NCOL(sp::coordinates(ips)))]
  } else if (is_2d &&
    (inherits(samplers, c(
      "SpatialPolygons",
      "SpatialPolygonsDataFrame"
    )) ||
      is.null(samplers))) {
    if (is.null(samplers)) {
      stop("Direct integration scheme for mesh domain with no samplers is not yet implemented.")
    }

    # If SpatialPolygons are provided convert into SpatialPolygonsDataFrame and attach weight = 1
    if (class(samplers)[1] == "SpatialPolygons") {
      samplers <- SpatialPolygonsDataFrame(samplers,
        data = data.frame(weight = rep(1, length(samplers))),
        match.ID = FALSE
      )
    } else if (is.null(samplers@data[["weight"]])) {
      samplers@data[["weight"]] <- 1
    }

    cnames <- coordnames(samplers)
    samplers_crs <- fm_CRS(samplers)

    # Convert samplers and domain to equal area CRS
    if (!fm_crs_is_null(fm_crs(domain))) {
      samplers <- fm_transform(samplers, crs = fm_crs("+proj=cea +units=km"))
    }

    # This old code doesn't handle holes properly.
    polyloc <- do.call(rbind, lapply(
      seq_len(length(samplers)),
      function(k) {
        cbind(
          x = rev(sp::coordinates(samplers@polygons[[k]]@Polygons[[1]])[, 1]),
          y = rev(sp::coordinates(samplers@polygons[[k]]@Polygons[[1]])[, 2]),
          group = k
        )
      }
    ))

    # This handles holes
    poly_segm <- fm_as_segm(samplers, join = FALSE)
    poly_segm <- lapply(
      seq_along(poly_segm),
      function(k) {
        segm <- poly_segm[[k]]
        segm[["grp"]] <- rep(k, NROW(segm[["idx"]]))
        segm[["is.bnd"]] <- TRUE
        segm
      }
    )

    # If domain is NULL, make a mesh with the polygons as boundary
    if (is.null(domain)) {
      warning("Computing integration points from polygon; specify a mesh for better numerical control.")
      max.edge <- max(diff(range(polyloc[, 1])), diff(range(polyloc[, 2]))) / 20
      domain <- fm_mesh_2d_inla(boundary = list(samplers), max.edge = max.edge)
      domain$crs <- fm_CRS(samplers)
      domain_crs <- fm_CRS(domain)
    } else {
      domain_crs <- fm_CRS(domain)
      if (!fm_crs_is_null(fm_crs(domain))) {
        domain <- fm_transform(domain, crs = fm_crs("+proj=cea +units=km"))
      }
    }


    if (identical(int.args[["poly_method"]], "legacy")) {
      stop("the legacy integration method method is no longer supported")
    } else {
      ips <- bru_int_polygon(
        domain,
        method = int.args$method,
        nsub = int.args$nsub2,
        samplers = samplers
      )
    }


    if (is.null(group)) {
      df <- data.frame(
        weight = ips[, "weight"] * samplers@data[ips$.block, "weight"],
        .block = ips$.block
      )
    } else {
      df <- data.frame(
        samplers@data[ips$.block, group, drop = FALSE],
        weight = ips[, "weight"] * samplers@data[ips$.block, "weight"],
        .block = ips$.block
      )
    }
    if (is.null(ips$z)) {
      ips <- sp::SpatialPointsDataFrame(ips[, c("x", "y")],
        data = df,
        match.ID = FALSE, proj4string = fm_CRS(domain)
      )
    } else {
      ips <- sp::SpatialPointsDataFrame(ips[, c("x", "y", "z")],
        data = df,
        match.ID = FALSE, proj4string = fm_CRS(domain)
      )
    }

    if (!fm_crs_is_null(domain_crs) && !fm_crs_is_null(samplers_crs)) {
      ips <- fm_transform(ips, crs = domain_crs)
    }

    coord_names <- c("x", "y", "z")
    #    if (!is.null(coordnames(samplers))) {
    #      coord_names[seq_along(coordnames(samplers))] <- coordnames(samplers)
    #    } else if (!is.null(name)) {
    #      coord_names[seq_along(name)] <- name
    #    }
    sp::coordnames(ips) <- coord_names[seq_len(NCOL(sp::coordinates(ips)))]
  } else {
    stop("No integration handling code reached; please notify the package developer.")
  }

  if (samplers_is_sf && inherits(ips, "Spatial")) {
    ips <- sf::st_as_sf(ips)
  }

  ips
}





#' @title (Blockwise) cross product of integration points
#'
#' @description
#' Calculates the groupwise cross product of integration points in different
#' dimensions and multiplies their weights accordingly.
#' If the object defining points in a particular dimension has no
#' weights attached to it all weights are assumed to be 1.
#'
#' Legacy wrapper for [fm_cprod()]
#'
#' @seealso [fm_cprod()]
#' @export
#' @keywords internal
#'
#'
#' @param ... `data.frame`, `sf`, or `SpatialPointsDataFrame` objects, each one
#' usually obtained by a call to an [fm_int()] method.
#' @param na.rm logical; if `TRUE`, the rows with weight `NA` from the
#' non-overlapping full_join will be removed; if `FALSE`, set the undefined weights to `NA`.
#' If `NULL` (default), act as `TRUE`, but warn if any elements needed removing.
#' @param .blockwise logical; if `FALSE`, computes full tensor product integration.
#' If `TRUE`, computes within-block tensor product integration (used internally
#' by [fm_int()]).
#' Default `FALSE`
#' @return A `data.frame`, `sf`, or `SpatialPointsDataFrame` of multidimensional
#' integration points and their weights
#'
#' @examples
#' \donttest{
#' # ipoints needs INLA
#' if (bru_safe_inla()) {
#'   # Create integration points in dimension 'myDim' and 'myDiscreteDim'
#'   ips1 <- fm_int(fm_mesh_1d(0:20),
#'     rbind(c(0, 3), c(3, 8)),
#'     name = "myDim"
#'   )
#'   ips2 <- fm_int(domain = c(1, 2, 4), name = "myDiscreteDim")
#'
#'   # Calculate the cross product
#'   ips <- cprod(ips1, ips2)
#'
#'   # Plot the integration points
#'   plot(ips$myDim, ips$myDiscreteDim, cex = 10 * ips$weight)
#' }
#' }
#'
cprod <- function(..., na.rm = NULL, .blockwise = FALSE) {
  lifecycle::deprecate_soft(
    "2.8.0",
    "cprod()",
    "fmesher::fm_cprod()"
  )
  fm_cprod(..., na.rm = na.rm, .blockwise = .blockwise)
}
