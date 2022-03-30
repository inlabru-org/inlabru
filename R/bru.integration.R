#' @title Generate integration points
#'
#' @description
#' This function generates points in one or two dimensions with a weight attached to each point.
#' The weighted sum of a function evaluated at these points is the integral of that function approximated
#' by linear basis functions. The parameter `samplers` describes the area(s) integrated over.
#'
#' In case of a single dimension `samplers` is supposed to be a two-column `matrix` where
#' each row describes the start and end point of the interval to integrate over. In the two-dimensional
#' case `samplers` can be either a `SpatialPolygon`, an `inla.mesh` or a
#' `SpatialLinesDataFrame` describing the area to integrate over. If a `SpatialLineDataFrame`
#' is provided it has to have a column called 'weight' in order to indicate the width of the line.
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
#' `Spatial*` objects, otherwise "x", "y", "coordinateZ" for 2D regions and
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
#'   that doesn't handle holes. Currently only used for debugging purposes.
#' @param project Deprecated in favour of `int.args(method=...)`.
#' If TRUE, aggregate the integration points to mesh vertices. Default:
#' `project = (identical(int.args$method, "stable"))`
#'
#' @return A `data.frame` or `SpatialPointsDataFrame` of 1D and 2D integration points, respectively.
#'
#' @examples
#' \donttest{
#' if (require("INLA", quietly = TRUE) &&
#'   require("ggplot2", quietly = TRUE)) {
#'
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
#'   mesh <- inla.mesh.1d(seq(0, 10, by = 1))
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
                    project = NULL) {
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

  if (!is.null(project)) {
    if (project && !identical(int.args$method, "stable")) {
      stop("ipoints(project=TRUE) is deprecated, and int.args$methods != 'stable'")
    } else if (!project && identical(int.args$method, "stable")) {
      stop("ipoints(project=FALSE) is deprecated, and int.args$methods == 'stable'")
    }
    warning(
      "ipoints(project=", ifelse(project, "TRUE", "FALSE"),
      ") is deprecated. Will use int.args$method = '",
      int.args[["method"]],
      "' instead."
    )
  }

  if (is.null(domain) &&
    inherits(samplers, c("inla.mesh.1d", "inla.mesh"))) {
    domain <- samplers
    samplers <- NULL
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
      inherits(domain, "inla.mesh")
  is_1d <- !is_2d &&
    (
      (!is.null(samplers) &&
        is.numeric(samplers)) ||
        (!is.null(domain) &&
          (is.numeric(domain) ||
            inherits(domain, "inla.mesh.1d")))
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
    inherits(samplers, c("inla.mesh.1d", "inla.mesh"))) {
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
  } else if (is_1d && is.null(samplers) && is.numeric(domain)) {
    ips <- data.frame(
      x = as.vector(domain),
      weight = 1
    )
    colnames(ips) <- c(name, "weight")
  } else if (is_1d && is.null(domain) && is.integer(samplers)) {
    ips <- data.frame(
      x = as.vector(samplers),
      weight = 1
    )
    colnames(ips) <- c(name, "weight")
  } else if (is_1d && is.null(samplers) && inherits(domain, "inla.mesh.1d") &&
    identical(int.args[["method"]], "stable")) {
    ips <- data.frame(
      x = domain$loc,
      weight = Matrix::diag(INLA::inla.mesh.fem(domain)$c0)
    )
    colnames(ips) <- c(name, "weight")
  } else if (is_1d) {
    domain_range <-
      if (inherits(domain, "inla.mesh.1d")) {
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
        domain <- INLA::inla.mesh.1d(sort(unique(as.vector(samplers))))
      }
    }
    # Now samplers is a 2-column matrix, and domain is an `inla.mesh.1d` object.

    ips <- list()

    # Integrate over each samplers subinterval, with nsub integration points
    # per domain knot interval
    if (domain$degree >= 2) {
      warning("Integration points projected onto knots may lead to instability for degree >= 2 basis functions.")
    }
    nsub <- int.args[["nsub1"]]
    u <- rep(
      (seq_len(nsub) - 0.5) / nsub,
      domain$n - 1
    )
    int_loc <-
      domain$loc[rep(seq_len(domain$n - 1), each = nsub)] * (1 - u) +
      domain$loc[rep(seq_len(domain$n - 1) + 1, each = nsub)] * u
    int_w <-
      (domain$loc[rep(seq_len(domain$n - 1) + 1, each = nsub)] -
        domain$loc[rep(seq_len(domain$n - 1), each = nsub)]) /
        nsub

    for (j in seq_len(nrow(samplers))) {
      subsamplers <- samplers[j, ]

      if (identical(int.args[["method"]], "stable")) {
        A_w <- INLA::inla.spde.make.A(domain,
          int_loc,
          weights =
            int_w *
              (int_loc >= min(subsamplers)) *
              (int_loc <= max(subsamplers))
        )
        ips[[j]] <- data.frame(
          loc = domain$loc,
          weight = Matrix::colSums(A_w)
        )
      } else {
        inside <-
          (int_loc >= min(subsamplers)) &
            (int_loc <= max(subsamplers))

        ips[[j]] <- data.frame(
          loc = int_loc[inside],
          weight = int_w[inside]
        )
      }
      colnames(ips[[j]]) <- c(name, "weight")
    }

    ips <- do.call(rbind, ips)
  } else if (inherits(domain, "inla.mesh") &&
    is.null(samplers) &&
    identical(int.args[["method"]], "stable")) {
    coord_names <- c("x", "y", "coordinateZ")
    if (!is.null(name)) {
      coord_names[seq_along(name)] <- name
    }

    # transform to equal area projection
    if (!fm_crs_is_null(domain$crs)) {
      crs <- domain$crs
      samplers <- stransform(domain, crs = CRS("+proj=cea +units=km"))
    }

    ips <- vertices(domain)
    ips$weight <- INLA::inla.mesh.fem(domain, order = 1)$va

    # backtransform
    if (!fm_crs_is_null(domain$crs)) {
      ips <- stransform(ips, crs = crs)
    }
    coordnames(ips) <- coord_names[seq_len(NCOL(coordinates(ips)))]
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
        data = data.frame(weight = rep(1, length(samplers)))
      )
    }

    # Set weight to 1 if not provided
    if (!("weight" %in% names(samplers))) {
      samplers$weight <- 1
    }

    ips <- int.slines(
      samplers,
      domain,
      group = group,
      project = identical(int.args[["method"]], "stable")
    )

    coord_names <- c("x", "y", "coordinateZ")
    if (!is.null(coordnames(samplers))) {
      coord_names[seq_along(coordnames(samplers))] <- coordnames(samplers)
    } else if (!is.null(name)) {
      coord_names[seq_along(name)] <- name
    }
    coordnames(ips) <- coord_names[seq_len(NCOL(coordinates(ips)))]
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
    samplers_crs <- fm_sp_get_crs(samplers)

    # Convert samplers and domain to equal area CRS
    if (!fm_crs_is_null(domain$crs)) {
      samplers <- stransform(samplers, crs = sp::CRS("+proj=cea +units=km"))
    }

    # This old code doesn't handle holes properly.
    polyloc <- do.call(rbind, lapply(
      seq_len(length(samplers)),
      function(k) {
        cbind(
          x = rev(coordinates(samplers@polygons[[k]]@Polygons[[1]])[, 1]),
          y = rev(coordinates(samplers@polygons[[k]]@Polygons[[1]])[, 2]),
          group = k
        )
      }
    ))

    # This handles holes
    poly_segm <- INLA::inla.sp2segment(samplers, join = FALSE)
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
      domain <- INLA::inla.mesh.2d(boundary = samplers, max.edge = max.edge)
      domain$crs <- fm_sp_get_crs(samplers)
    } else {
      if (!fm_crs_is_null(domain$crs)) {
        domain <- stransform(domain, crs = CRS("+proj=cea +units=km"))
      }
    }
    domain_crs <- fm_ensure_crs(domain$crs)

    if (identical(int.args[["poly_method"]], "legacy")) {
      ips <- int.polygon(domain,
        loc = polyloc[, 1:2], group = polyloc[, 3],
        method = int.args$method, nsub = int.args$nsub2
      )
    } else {
      if (!is.null(int.args$use_new) && !int.args$use_new) {
        ips <- bru_int_polygon(
          domain,
          polylist = poly_segm,
          method = int.args$method,
          nsub = int.args$nsub2
        )
      } else {
        ips <- bru_int_polygon(
          domain,
          polylist = poly_segm,
          method = int.args$method,
          nsub = int.args$nsub2,
          samplers = samplers
        )
      }
    }


    df <- data.frame(
      samplers@data[ips$group, group, drop = FALSE],
      weight = ips[, "weight"] * samplers@data[ips$group, "weight"]
    )
    if (is.null(ips$coordinateZ)) {
      ips <- SpatialPointsDataFrame(ips[, c("x", "y")],
        data = df,
        match.ID = FALSE, proj4string = domain_crs
      )
    } else {
      ips <- SpatialPointsDataFrame(ips[, c("x", "y", "coordinateZ")],
        data = df,
        match.ID = FALSE, proj4string = domain_crs
      )
    }

    if (!fm_crs_is_null(domain_crs) && !fm_crs_is_null(samplers_crs)) {
      ips <- stransform(ips, crs = samplers_crs)
    }

    coord_names <- c("x", "y", "coordinateZ")
    if (!is.null(coordnames(samplers))) {
      coord_names[seq_along(coordnames(samplers))] <- coordnames(samplers)
    } else if (!is.null(name)) {
      coord_names[seq_along(name)] <- name
    }
    coordnames(ips) <- coord_names[seq_len(NCOL(coordinates(ips)))]
  } else {
    stop("No integration handling code reached; please notify the package developer.")
  }

  ips
}

#' @title Cross product of integration points
#'
#' @description
#' Calculates the cross product of integration points in different dimensions
#' and multiplies their weights accordingly. If the object defining points in a particular
#' dimension has no weights attached to it all weights are assumend to be 1.
#'
#' @aliases cprod
#' @export
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'
#' @param ... `data.frame` or `SpatialPointsDataFrame` objects, each one usually obtained by a call to the [ipoints] function.
#' @return A `data.frame` or `SpatialPointsDataFrame` of multidimensional integration points and their weights
#'
#' @examples
#' \donttest{
#' # ipoints needs INLA
#' if (bru_safe_inla()) {
#'   # Create integration points in dimension 'myDim' and 'myDiscreteDim'
#'   ips1 <- ipoints(rbind(c(0, 3), c(3, 8)), 17, name = "myDim")
#'   ips2 <- ipoints(domain = c(1, 2, 4), name = "myDiscreteDim")
#'
#'   # Calculate the cross product
#'   ips <- cprod(ips1, ips2)
#'
#'   # Plot the integration points
#'   plot(ips$myDim, ips$myDiscreteDim, cex = 10 * ips$weight)
#' }
#' }
#'
cprod <- function(...) {
  ipl <- list(...)
  ipl <- ipl[!vapply(ipl, is.null, TRUE)]
  if (length(ipl) == 0) {
    return(NULL)
  }

  if (length(ipl) == 1) {
    ips <- ipl[[1]]
  } else {
    ips1 <- ipl[[1]]
    if (length(ipl) > 2) {
      ips2 <- do.call(cprod, ipl[2:length(ipl)])
    } else {
      ips2 <- ipl[[2]]
    }
    if (!"weight" %in% names(ips1)) {
      ips1$weight <- 1
    }
    if (!"weight" %in% names(ips2)) {
      ips2$weight <- 1
    }

    by <- setdiff(intersect(names(ips1), names(ips2)), "weight")
    if (inherits(ips1, "Spatial")) {
      ips <- sp::merge(ips1, ips2, by = by, duplicateGeoms = TRUE)
    } else if (inherits(ips2, "Spatial")) {
      ips <- sp::merge(ips2, ips1, by = by, duplicateGeoms = TRUE)
    } else {
      ips <- base::merge(ips1, ips2, by = by)
    }
    ips$weight <- ips$weight.x * ips$weight.y
    ips[["weight.x"]] <- NULL
    ips[["weight.y"]] <- NULL
    row.names(ips) <- as.character(seq_len(NROW(ips)))
  }
  ips
}

# Integration points for log Gaussian Cox process models using INLA
#
# prerequisits:
#
# - List of integration dimension names, extend and quadrature
# - Samplers: These may live in a subset of the dimensions, usually space and time
#             ("Where and wehen did a have a look at the point process")
# - Actually this is a simplified view. Samplers should have start and end time !
#
# Procedure:
# - Select integration strategy by type of samplers:
#       1) SpatialPointsDataFrame: Assume these are already integration points
#       2) SpatialLinesDataFrame: Use simplified integration along line with (width provided by samplers)
#       3) SpatialPolygonDataFrame: Use full integration over polygons
#
# - Create integration points from samplers. Do NOT perform simplification projection here!
# - Simplify integration points.
#   1) Group by non-mesh dimensions, e.g. time, weather
#   2) For each group simplify with respect to mesh-dimensions, e.g. space
#   3) Merge
#
#   Local dependencies:
#     `int.polygon()`
#
# @aliases ipoints
# @export
# @param samplers A `Spatial[Points/Lines/Polygons]DataFrame` object
# @param domain A list of named integration definitions, each either a numeric
# vector of points given integration weight 1, an `inla.mesh.1d` object, or an
# `inla.mesh.2d` object. Only those domains that are not given in the `samplers`
# data.frame are used, plus the coordinates object, used for the spatial aspect
# of the `samplers` object.
# @param int.args List of arguments passed on to \code{ipoints}
# @return Integration points


ipmaker <- function(samplers, domain, dnames,
                    int.args = list(method = "stable", nsub = NULL)) {
  if ("coordinates" %in% dnames) {
    spatial <- TRUE
  } else {
    spatial <- FALSE
  }

  # Dimensions provided via samplers (except "coordinates")
  samp.dim <- intersect(names(samplers), dnames)

  # Dimensions provided via domain but not via samplers
  nosamp.dim <- setdiff(names(domain), c(samp.dim, "coordinates"))

  # Check if a domain definition is missing
  missing.dims <- setdiff(dnames, c(names(domain), samp.dim))
  if (length(missing.dims > 0)) {
    stop(paste0(
      "Domain definitions missing for dimensions: ",
      paste0(missing.dims, collapse = ", ")
    ))
  }
  extra.dims <- setdiff(names(domain), c(samp.dim, nosamp.dim))
  if (length(missing.dims > 0)) {
    warning(paste0(
      "Unexpected extra domain defintions: ",
      paste0(extra.dims, collapse = ", ")
    ))
  }

  if (spatial) {
    ips <- ipoints(samplers, domain$coordinates,
      group = samp.dim, int.args = int.args
    )
  } else {
    ips <- NULL
  }

  lips <- lapply(nosamp.dim, function(nm) ipoints(NULL, domain[[nm]], name = nm, int.args = int.args))
  ips <- do.call(cprod, c(list(ips), lips))
}




# Project data to mesh vertices under the assumption of lineariity
#
#
# @aliases vertex.projection
# @export
# @param points A SpatialPointsDataFrame object
# @param mesh An inla.mesh object
# @param columns A character array of the points columns which whall be projected
# @param group Character array identifying columns in \code{points}. These
# colouns are interpreted as factors and the projection is performed
# independently for each combination of factor levels.
# @return SpatialPointsDataFrame of mesh vertices with projected data attached

vertex.projection <- function(points, mesh, columns = names(points), group = NULL, fill = NULL) {
  if (is.null(group) | (length(group) == 0)) {
    res <- INLA::inla.fmesher.smorg(mesh$loc, mesh$graph$tv, points2mesh = coordinates(points))
    tri <- res$p2m.t

    data <- list()
    for (k in seq_along(columns)) {
      cn <- columns[k]
      nw <- points@data[, columns] * res$p2m.b
      w.by <- by(as.vector(nw), as.vector(mesh$graph$tv[tri, ]), sum, simplify = TRUE)
      data[[cn]] <- as.vector(w.by)
    }

    data <- data.frame(data)
    coords <- mesh$loc[as.numeric(names(w.by)), , drop = FALSE]
    data$vertex <- as.numeric(names(w.by))

    ret <- SpatialPointsDataFrame(coords,
      proj4string = fm_sp_get_crs(points),
      data = data,
      match.ID = FALSE
    )
    coordnames(ret) <- coordnames(points)

    # If fill is not NULL, add vertices to which no data was projected
    # and set their projected data according to `fill`

    if (!is.null(fill)) {
      vrt <- vertices(mesh)
      vrt <- vrt[setdiff(vrt$vertex, data$vertex), ]
      if (nrow(vrt) > 0) {
        for (nm in setdiff(names(data), "vertex")) vrt[[nm]] <- fill
        ret <- rbind(ret, vrt)
      }
      ret <- ret[match(1:mesh$n, ret$vertex), ]
    }
  } else {
    fn <- function(X) {
      ret <- vertex.projection(X, mesh, columns = columns)
      for (g in group) {
        ret[[g]] <- X[[g]][1]
      }
      ret
    }
    ret <- lapply(
      unique(points[[group]]),
      function(x) {
        fn(points[points[[group]] == x, , drop = FALSE])
      }
    )
    ret <- do.call(rbind, ret)
  }
  ret
}


# Project data to mesh vertices under the assumption of linearity
#
#
# @aliases vertex.projection
# @export
# @param points A SpatialPointsDataFrame object
# @param mesh An inla.mesh object
# @param columns A character array of the points columns which whall be projected
# @param group Character array identifying columns in \code{points}. These
# columns are interpreted as factors and the projection is performed
# independently for each combination of factor levels.
# @return SpatialPointsDataFrame of mesh vertices with projected data attached
# @example
#
# pts = data.frame(x = 50 * runif(10), weight = abs(rnorm(100)))
# msh = inla.mesh.1d(seq(0,50,by=1))
# pts$year = c(rep(1,5), rep(2,5))
# ip =  vertex.projection.1d(pts, msh)
# ggplot(ip) + geom_point(aes(x=x, y=weight))
#
# ip =  vertex.projection.1d(pts, msh, group = "year", fill = 0, column = "weight")
# head(ip)
# ggplot(ip) + geom_point(aes(x=x, y=weight, color = year))

vertex.projection.1d <- function(points, mesh, group = NULL, column = "weight", simplify = TRUE, fill = NULL) {
  dname <- setdiff(names(points), c(column, group))
  if (length(dname) > 1) {
    dname <- dname[1]
  }

  xx <- points[, dname]
  ww <- points[, column]
  iv <- findInterval(xx, mesh$loc)

  # Left and right vertex location
  left <- mesh$loc[iv]
  right <- mesh$loc[iv + 1]

  # Relative location within the two neighboring vertices
  w.right <- (xx - left) / (right - left)
  w.left <- 1 - w.right

  # Projected integration points
  ips <- rbind(
    data.frame(x = left, vertex = iv),
    data.frame(x = right, vertex = iv + 1)
  )
  ips[column] <- c(ww * w.left, ww * w.right)


  # Simplify
  if (simplify) {
    bygroup <- list(vertex = ips$vertex)
    if (!is.null(group)) {
      bygroup <- c(bygroup, as.list(rbind(points[, group, drop = FALSE], points[, group, drop = FALSE])))
    }
    ips <- aggregate(ips[, column, drop = FALSE], by = bygroup, FUN = sum)
  }

  # Add x-coordinate
  ips[dname] <- mesh$loc[ips$vertex]

  # Fill
  if (!is.null(fill)) {
    miss <- setdiff(seq_len(length(mesh$loc)), ips$vertex)
    mips <- data.frame(vertex = miss, x = mesh$loc[miss])
    mips[, column] <- fill
    ips <- rbind(ips, merge(mips, ips[, group, drop = FALSE]))
  }

  ips
}


#' Weighted summation (integration) of data frame subsets
#'
#' A typical task in statistical inference to integrate a (multivariate) function along one or
#' more dimensions of its domain. For this purpose, the function is evaluated at some points
#' in the domain and the values are summed up using weights that depend on the area being
#' integrated over. This function performs the weighting and summation conditional for each level
#' of the dimensions that are not integrated over. The parameter `dims` states the the
#' dimensions to integrate over. The set of dimensions that are held fixed is the set difference
#' of all column names in `data` and the dimensions stated by `dims`.
#'
#' @aliases int
#' @export
#' @param data A `data.frame` or `Spatial` object. Has to have a `weight` column with numeric values.
#' @param values Numerical values to be summed up, usually the result of function evaluations.
#' @param dims Column names (dimension names) of the `data` object to integrate over.
#' @return A `data.frame` of integrals, one for each level of the cross product of all dimensions not being integrated over.
#'
#' @examples
#' \donttest{
#' # ipoints needs INLA
#' if (bru_safe_inla(quietly = TRUE)) {
#'   # Create integration points in two dimensions, x and y
#'
#'   ips <- cprod(
#'     ipoints(c(0, 10), 10, name = "x"),
#'     ipoints(c(1, 5), 10, name = "y")
#'   )
#'
#'   # The sizes of the domains are 10 and 4 for x and y, respectively.
#'   # Integrating f(x,y) = 1 along x and y should result in the total
#'   # domain size 40
#'
#'   int(ips, rep(1, nrow(ips)), c("x", "y"))
#' }
#' }
#'
int <- function(data, values, dims = NULL) {
  keep <- setdiff(names(data), c(dims, "weight"))
  if (length(keep) > 0 & !is.null(dims)) {
    agg <- aggregate(values * data$weight, by = as.list(data[, keep, drop = FALSE]), FUN = sum)
    names(agg)[ncol(agg)] <- "integral" # paste0("integral_{",dims,"}(",deparse(values),")")
  } else {
    agg <- sum(data$weight * values)
  }

  agg
}
