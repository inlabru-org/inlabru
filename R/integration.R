# Gaussian quadrature and other integration point constructors
#
# Construct integration points for each of lines defined by the start and end
# points provided.
# The following schemes are available:
# "equidistant" : Equidistant integration points without boundary. All weights
# are identical and sum uf to the length of a line.
# "gaussian": Points and weight according to the Gaussian quadrature rule.
# Currently only n=1 and n=2 are supported (Exact integration for linear and
# quadratic functions).
# "twosided-gaussian": Experimental
#
# @export
# @param sp Start points of lines
# @param ep End points of lines
# @param scheme Integration scheme (gaussian or equidistant)
# @param n Number of integration points
# @return List with integration poins (ips), weights (w) and weights including
#   line length (wl)
# @author Fabian E. Bachl \email{f.e.bachl@@bath.ac.uk}

int.quadrature <- function(sp = NULL,
                           ep = NULL,
                           scheme = "gaussian",
                           n.points = 1,
                           geometry = "euc",
                           coords = NULL) {
  if (is.null(colnames(sp)) && !is.null(coords)) {
    sp <- data.frame(sp)
    colnames(sp) <- coords
  }
  if (is.null(colnames(ep)) && !is.null(coords) && !(length(ep) == 0)) {
    ep <- data.frame(ep)
    colnames(ep) <- coords
  }

  if (is.vector(sp)) {
    n.lines <- 1
    len <- function(x) {
      abs(x)
    }
  } else {
    n.lines <- dim(sp)[1]
    len <- function(x) {
      apply(x^2, MARGIN = 1, function(y) {
        return(sqrt(sum(y)))
      })
    }
  }
  if (scheme == "gaussian") {
    if (n.points == 1) {
      # points
      if (geometry == "euc") {
        ips <- sp + (ep - sp) / 2
        # weights
        w <- rep(1, n.lines)
        wl <- w * len(ep - sp)
      } else if (geometry == "geo") {
        stop("Geometry geo not supported")
      }

      # PD vector
      # pdv = cbind(-(ep[,2]-sp[,2]),ep[,1]-sp[,1])
      # Index of line a point comes from
      line.idx <- 1:n.lines
    } else if (n.points == 2) {
      if (geometry == "geo") {
        stop("Geometry geo not supported")
      }
      # points
      ips1 <- sp + (-0.5 * sqrt(1 / 3) + 1 / 2) * (ep - sp)
      ips2 <- sp + (0.5 * sqrt(1 / 3) + 1 / 2) * (ep - sp)
      ips <- rbind(ips1, ips2)

      # weights
      w <- rep(1, dim(ips)[1]) / 2
      wl <- w * len(ep - sp)

      # PD vector
      # dvec = cbind(-(ep[,2]-sp[,2]),ep[,1]-sp[,1])
      # dvec = rbind(dvec,dvec)

      # Index of line a point comes from
      line.idx <- c(1:n.lines, 1:n.lines)
    } else if (n.points == 3) {
      if (geometry == "geo") {
        stop("Geometry geo not supported")
      }
      # points
      ips1 <- sp + (-0.5 * sqrt(3 / 5) + 1 / 2) * (ep - sp)
      ips2 <- sp + 0.5 * (ep - sp)
      ips3 <- sp + (0.5 * sqrt(3 / 5) + 1 / 2) * (ep - sp)
      ips <- rbind(ips1, ips2, ips3)

      # weights
      w <- c(rep(5 / 9, n.lines), rep(8 / 9, n.lines), rep(5 / 9), n.lines) / 2
      wl <- w * len(ep - sp)

      # Index of line a point comes from
      line.idx <- c(1:n.lines, 1:n.lines)
    } else {
      stop("Gaussian quadrature with n>3 not implemented")
    }
  } else if (scheme == "twosided-gaussian") {
    if (geometry == "geo") {
      stop("Geometry geo not supported")
    }
    ips1 <- int.quadrature(sp,
      sp + 0.5 * (ep - sp),
      scheme = "gaussian",
      n.points = n.points,
      geometry,
      coords
    )
    ips2 <- int.quadrature(sp + 0.5 * (ep - sp),
      ep,
      scheme = "gaussian",
      n.points = n.points,
      geometry,
      coords
    )
    ips <- rbind(ips1$ips, ips2$ips)
    w <- 0.5 * rbind(ips1$w, ips2$w)
    wl <- 0.5 * rbind(ips1$wl, ips2$wl)
    line.idx <- rbind(ips1$line.idx, ips2$line.idx)
  } else if (scheme == "equidistant") {
    if (geometry == "geo") {
      stop("Geometry geo not supported")
    }
    # points
    mult <- seq(0 - 1 / (2 * n.points), 1 + 1 / (2 * n.points),
      length.out = n.points + 2
    )[2:(n.points + 1)]
    nips <- list()
    for (k in 1:n.points) {
      nips[[k]] <- sp + mult[k] * (ep - sp)
    }
    ips <- do.call(rbind, nips)
    # weights
    w <- rep(1, dim(ips)[1]) / n.points
    wl <- w * (len(ep - sp)[rep(1:n.lines, n.points)])
    # Index of line a point comes from
    line.idx <- rep(1:n.lines, n.points)
  } else if (scheme == "trapezoid") {
    if (geometry == "geo") {
      stop("Geometry geo not supported")
    }
    # points
    mult <- seq(0, 1, length.out = n.points)
    nips <- list()
    for (k in 1:n.points) {
      nips[[k]] <- sp + mult[k] * (ep - sp)
    }
    ips <- do.call(rbind, nips)
    # weights
    w <- rep(1, dim(ips)[1]) / (n.points - 1)
    w[1] <- 0.5 * w[1]
    w[length(w)] <- 0.5 * w[length(w)]
    wl <- w * (len(ep - sp)[rep(1:n.lines, n.points)])
    # Index of line a point comes from
    line.idx <- rep(1:n.lines, n.points)
  } else if (scheme == "fixed") {
    ips <- data.frame(tmp = sp)
    # colnames(ips) = coords
    w <- rep(1, length(sp))
    wl <- rep(1, length(sp))
    line.idx <- rep(NaN, length(sp))
  }
  return(list(ips = ips, w = w, wl = wl, line.idx = line.idx))
}









# New integration methods ----





#' @describeIn inlabru-deprecated
#' Aggregate integration weights onto mesh nodes
#'
#' `r lifecycle::badge("deprecated")`
#' Use [fmesher::fm_vertex_projection()] instead.
#'
#' @param mesh Mesh on which to integrate
#' @param integ `list` of `loc`, integration points,
#'   and `weight`, integration weights,
#'   or a `SpatialPointsDataFrame`. Only the coordinates and `weight` column
#'   are handled.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @keywords internal
integration_weight_aggregation <- function(mesh, integ) {
  lifecycle::deprecate_stop(
    "2.8.0",
    "integration_weight_aggregation()",
    "fmesher::fm_vertex_projection()"
  )

  fm_vertex_projection(points = integ, mesh = mesh)
}

#' @describeIn inlabru-deprecated
#' Integration scheme for mesh triangle interiors
#'
#' `r lifecycle::badge("deprecated")`
#' Use [fmesher::fm_int_mesh_2d_core()] instead.
#'
#' @param mesh Mesh on which to integrate
#' @param tri_subset Optional triangle index vector for integration on a subset
#' of the mesh triangles (Default `NULL`)
#' @param nsub number of subdivision points along each triangle edge, giving
#'    `(nsub + 1)^2` proto-integration points used to compute
#'   the vertex weights
#'   (default `NULL=9`, giving 100 integration points for each triangle)
#' @return * `mesh_triangle_integration` returns a `list` with elements `loc`
#' and `weight` with integration points for the mesh
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @keywords internal
mesh_triangle_integration <- function(mesh, tri_subset = NULL, nsub = NULL) {
  lifecycle::deprecate_stop(
    "2.8.0",
    "mesh_triangle_integration()",
    "fmesher::fm_int_mesh_2d_core()"
  )

  fmesher::fm_int_mesh_2d_core(
    mesh = mesh,
    tri_subset = tri_subset,
    nsub = nsub
  )
}



#' @describeIn inlabru-deprecated
#' Integration points for polygons inside an inla.mesh
#'
#' `r lifecycle::badge("deprecated")` Use [fmesher::fm_int()] instead.
#'
#' @keywords internal

bru_int_polygon <- function(...) {
  lifecycle::deprecate_stop(
    "2.8.0",
    "bru_int_polygon()",
    "fmesher::fm_int()"
  )
}
