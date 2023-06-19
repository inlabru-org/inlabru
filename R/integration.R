# Gaussian quadrature and other integration point constructors
#
# Construct integration points for each of lines defined by the start and end points provided.
# The following schemes are available:
# "equidistant" : Equidistant integration points without boundary. All weights are identical and sum uf to the length of a line.
# "gaussian": Points and weight according to the Gaussian quadrature rule. Currently only n=1 and n=2 are supported (Exact integration for linear and quadratic functions).
# "twosided-gaussian": Experimental
#
# @aliases int.quadrature
# @export
# @param sp Start points of lines
# @param ep End points of lines
# @param scheme Integration scheme (gaussian or equidistant)
# @param n Number of integration points
# @return List with integration poins (ips), weights (w) and weights including line length (wl)
# @author Fabian E. Bachl \email{f.e.bachl@@bath.ac.uk}

int.quadrature <- function(sp = NULL, ep = NULL, scheme = "gaussian", n.points = 1, geometry = "euc", coords = NULL) {
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
    ips1 <- int.quadrature(sp, sp + 0.5 * (ep - sp), scheme = "gaussian", n.points = n.points, geometry, coords)
    ips2 <- int.quadrature(sp + 0.5 * (ep - sp), ep, scheme = "gaussian", n.points = n.points, geometry, coords)
    ips <- rbind(ips1$ips, ips2$ips)
    w <- 0.5 * rbind(ips1$w, ips2$w)
    wl <- 0.5 * rbind(ips1$wl, ips2$wl)
    line.idx <- rbind(ips1$line.idx, ips2$line.idx)
  } else if (scheme == "equidistant") {
    if (geometry == "geo") {
      stop("Geometry geo not supported")
    }
    # points
    mult <- seq(0 - 1 / (2 * n.points), 1 + 1 / (2 * n.points), length.out = n.points + 2)[2:(n.points + 1)]
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


# Arc length from Cartesian coordinates
sphere_geodesic_length <- function(sp, ep) {
  # Needs to get coordinates on S2 with arbitrary radius; the caller
  # can then scale by the appropriate ellipsoid radius.
  #   len = atan2(crossproduct(sp, ep), dotproduct(sp, ep))
  # but check fmesher code for potentially more stable method when sp \approx ep.
  stop("Not implemented")
}

# 2022-11-29: Currently assumes the data is Spatial
int.slines <- function(data, mesh, .block = NULL, project = TRUE) {
  # Extract start and end coordinates
  qq <- coordinates(data)
  sp <- do.call(
    rbind,
    lapply(
      qq,
      function(k) {
        do.call(
          rbind,
          lapply(k, function(x) x[1:(nrow(x) - 1), , drop = FALSE])
        )
      }
    )
  )
  ep <- do.call(
    rbind,
    lapply(
      qq,
      function(k) {
        do.call(
          rbind,
          lapply(k, function(x) x[2:(nrow(x)), , drop = FALSE])
        )
      }
    )
  )

  idx <- do.call(
    rbind,
    lapply(
      seq_along(qq),
      function(k) {
        do.call(
          cbind,
          lapply(qq[[k]], function(x) rep(k, nrow(x) - 1))
        )
      }
    )
  )
  idx <- cbind(idx, idx)

  sampler_crs <- fm_crs(data)
  target_crs <- fm_crs(mesh)
  if (!fm_crs_is_null(sampler_crs) &&
    fm_crs_is_null(target_crs)) {
    target_crs <- sampler_crs
  }

  if (!is.null(mesh)) {
    # Filter out points outside the mesh...
    sp <- fm_transform(sp, crs = target_crs, crs0 = sampler_crs, passthrough = TRUE)
    ep <- fm_transform(ep, crs = target_crs, crs0 = sampler_crs, passthrough = TRUE)
    proj1 <- fm_evaluator(mesh, loc = sp, crs = target_crs)
    proj2 <- fm_evaluator(mesh, loc = ep, crs = target_crs)
    ok <- (proj1$proj$ok & proj2$proj$ok)
    if (!all(ok)) {
      warning("Found spatial lines with start or end point ouside of the mesh. Omitting.")
    }
    sp <- sp[ok, , drop = FALSE]
    ep <- ep[ok, , drop = FALSE]
    idx <- idx[ok, , drop = FALSE]

    # Split at mesh edges
    line.spl <- split_lines(mesh, sp, ep, TRUE)
    sp <- line.spl$sp
    ep <- line.spl$ep
    idx <- idx[line.spl$split.origin, ]
  }

  # At this point, sp and ep are in the target_crs

  # Determine integration points along lines

  if (fm_crs_is_null(sampler_crs)) {
    ips <- (sp + ep) / 2
    w <- rowSums((ep - sp)^2)^0.5
  } else {
    # Has CRS
    longlat.crs <- fm_crs("longlat_globe")
    geocentric.crs <- fm_crs("sphere")
    sp3d <- fm_transform(sp, crs = geocentric.crs, crs0 = target_crs)
    ep3d <- fm_transform(ep, crs = geocentric.crs, crs0 = target_crs)
    mp3d <- (sp3d + ep3d) / rowSums((sp3d + ep3d)^2)^0.5

    ips <- fm_transform(mp3d, crs = target_crs, crs0 = geocentric.crs)
    w <- sp::spDists(
      fm_transform(sp3d, crs = longlat.crs, crs0 = geocentric.crs)[, 1:2, drop = FALSE],
      fm_transform(ep3d, crs = longlat.crs, crs0 = geocentric.crs)[, 1:2, drop = FALSE],
      diagonal = TRUE, longlat = TRUE
    )
  }

  # Wrap everything up and perform projection according to distance and given .block argument
  ips <- data.frame(ips)
  d_ips <- ncol(ips)
  # Temporary names
  colnames(ips) <- c("x", "y", "z")[seq_len(d_ips)]

  # Weights
  ips <- cbind(ips, weight = w)
  if ("weight" %in% names(data)) {
    ips$weight <- ips$weight * data$weight[idx[, 1]]
  }

  if (!is.null(.block)) {
    ips <- cbind(ips, as.data.frame(data)[idx[, 1], .block, drop = FALSE])
  }

  ips <- sp::SpatialPointsDataFrame(
    ips[, 1:d_ips, drop = FALSE],
    data = ips[, -(1:d_ips), drop = FALSE],
    proj4string = fm_CRS(target_crs)
  )
  if (!is.null(coordnames(data))) {
    name <- coordnames(data)
    if (length(name) < d_ips) {
      name <- c(name, "z")
    }
    coordnames(ips) <- name
  }

  # Project to mesh vertices
  if (project && !is.null(mesh)) {
    ips <- fm_vertex_projection(ips, mesh)
  }

  ips
}



join_segm <- function(...) {
  segm_list <- list(...)
  loc <- matrix(0, 0, 3)
  idx <- matrix(0, 0, 2)
  for (k in seq_along(segm_list)) {
    idx <- rbind(idx, segm_list[[k]]$idx + nrow(loc))
    loc <- rbind(loc, segm_list[[k]]$loc)
  }

  # Collapse duplicate points
  new_loc <- loc
  new_idx <- seq_len(nrow(loc))
  prev_idx <- 0
  for (k in seq_len(nrow(loc))) {
    if (any(is.na(new_loc[k, ]))) {
      new_idx[k] <- NA
    } else {
      if (prev_idx == 0) {
        prev_dist <- 1
      } else {
        prev_dist <- ((new_loc[seq_len(prev_idx), 1] - new_loc[k, 1])^2 +
          (new_loc[seq_len(prev_idx), 2] - new_loc[k, 2])^2 +
          (new_loc[seq_len(prev_idx), 3] - new_loc[k, 3])^2)^0.5
      }
      if (all(prev_dist > 0)) {
        prev_idx <- prev_idx + 1
        new_idx[k] <- prev_idx
        new_loc[prev_idx, ] <- new_loc[k, ]
      } else {
        new_idx[k] <- which.min(prev_dist)
      }
    }
  }
  idx <- matrix(new_idx[idx], nrow(idx), 2)
  # Remove NA and atomic lines
  ok <-
    !is.na(idx[, 1]) &
      !is.na(idx[, 2]) &
      idx[, 1] != idx[, 2]
  idx <- idx[ok, , drop = FALSE]
  # Set locations
  loc <- new_loc[seq_len(prev_idx), , drop = FALSE]

  INLA::inla.mesh.segment(
    loc = loc,
    idx = idx,
    is.bnd = FALSE
  )
}



#' Construct the intersection mesh of a mesh and a polygon
#'
#' @param mesh `inla.mesh` object to be intersected
#' @param poly `inla.mesh.segment` object with a closed polygon
#'   to intersect with the mesh
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @keywords internal
intersection_mesh <- function(mesh, poly) {
  if (ncol(poly$loc) < 3) {
    poly$loc <- cbind(poly$loc, 0)
  }

  all_edges <- INLA::inla.mesh.segment(
    loc = mesh$loc,
    idx = cbind(
      as.vector(t(mesh$graph$tv)),
      as.vector(t(mesh$graph$tv[, c(2, 3, 1), drop = FALSE]))
    ),
    is.bnd = FALSE
  )

  mesh_cover <- INLA::inla.mesh.create(
    loc = rbind(mesh$loc, poly$loc),
    interior = c(list(all_edges))
  )

  split <- INLA::inla.fmesher.smorg(mesh_cover$loc,
    mesh_cover$graph$tv,
    splitlines = list(
      loc = poly$loc,
      idx = poly$idx
    )
  )
  split_segm <- INLA::inla.mesh.segment(
    loc = split$split.loc,
    idx = split$split.idx,
    is.bnd = FALSE
  )

  joint_segm <- join_segm(split_segm, all_edges)

  mesh_joint_cover <- INLA::inla.mesh.create(
    interior = list(joint_segm),
    extend = TRUE
  )

  mesh_poly <- INLA::inla.mesh.create(boundary = poly)

  loc_tri <-
    (mesh_joint_cover$loc[mesh_joint_cover$graph$tv[, 1], , drop = FALSE] +
      mesh_joint_cover$loc[mesh_joint_cover$graph$tv[, 2], , drop = FALSE] +
      mesh_joint_cover$loc[mesh_joint_cover$graph$tv[, 3], , drop = FALSE]) / 3
  ok_tri <-
    fm_evaluator(mesh, loc = loc_tri)$proj$ok &
      fm_evaluator(mesh_poly, loc = loc_tri)$proj$ok
  if (any(ok_tri)) {
    loc_subset <- unique(sort(as.vector(mesh_joint_cover$graph$tv[ok_tri, , drop = FALSE])))
    new_idx <- integer(mesh$n)
    new_idx[loc_subset] <- seq_along(loc_subset)
    tv_subset <- matrix(new_idx[mesh_joint_cover$graph$tv[ok_tri, , drop = FALSE]],
      ncol = 3
    )
    loc_subset <- mesh_joint_cover$loc[loc_subset, , drop = FALSE]
    mesh_subset <- INLA::inla.mesh.create(
      loc = loc_subset,
      tv = tv_subset,
      extend = FALSE
    )
  } else {
    mesh_subset <- NULL
  }

  mesh_subset
}



# Basic robust integration weights for mesh/polygon intersections
#
# @param mesh Mesh on which to integrate
# @param bnd `inla.mesh.segment` defining the integration domain
# @param nsub number of subdivision points along each triangle edge, giving
#    `(nsub + 1)^2` proto-integration points used to compute
#   the vertex weights
#   (default `nsub=9`, giving 100 integration points for each triangle)
# @return `list` with elements `loc` and `weight` with
#   integration points for the intersection of the mesh and polygon
# @author Finn Lindgren \email{finn.lindgren@@gmail.com}
# @keywords internal
make_stable_integration_points <- function(mesh, bnd, nsub = NULL) {
  # Construct a barycentric grid of subdivision triangle midpoints
  if (is.null(nsub)) {
    nsub <- 9
  }
  stopifnot(nsub >= 0)
  nB <- (nsub + 1)^2

  b <- seq(1 / 3, 1 / 3 + nsub, length = nsub + 1) / (nsub + 1)
  bb <- as.matrix(expand.grid(b, b))
  # Points above the diagonal should be reflected into the lower triangle:
  refl <- rowSums(bb) > 1
  if (any(refl)) {
    bb[refl, ] <- cbind(1 - bb[refl, 2], 1 - bb[refl, 1])
  }
  # Construct complete barycentric coordinates:
  barycentric_grid <- cbind(1 - rowSums(bb), bb)

  # Construct integration points
  nT <- nrow(mesh$graph$tv)
  loc <- matrix(0.0, nT * nB, 3)
  idx_end <- 0
  for (tri in seq_len(nT)) {
    idx_start <- idx_end + 1
    idx_end <- idx_start + nB - 1
    loc[seq(idx_start, idx_end, length = nB), ] <-
      as.matrix(barycentric_grid %*%
        mesh$loc[mesh$graph$tv[tri, ], , drop = FALSE])
  }

  # Construct integration weights
  weight <- rep(INLA::inla.mesh.fem(mesh, order = 1)$ta / nB, each = nB)

  # Filter away points outside integration domain boundary:
  mesh_bnd <- INLA::inla.mesh.create(boundary = bnd)
  ok <- fm_evaluator(mesh_bnd, loc = loc)$proj$ok

  list(
    loc = loc[ok, , drop = FALSE],
    weight = weight[ok]
  )
}

# Integration points for polygons inside an inla.mesh
#
# This method doesn't handle polygons with holes. Use [bru_int_polygon()]
# instead.
#
# @param mesh An inla.mesh object
# @param loc Locations defining the polygons
# @param group If loc defines multiple polygons then this is the ID of the group for each location in loc
# @param method Which integration method to use ("stable", with aggregation to mesh vertices, or "direct")
# @param ... Arguments passed to the low level integration method (`make_stable_integration_points`)
# @author Fabian E. Bachl \email{f.e.bachl@@bath.ac.uk} and Finn Lindgren \email{finn.lindgren@@gmail.com}
# @keywords internal

int.polygon <- function(mesh, loc, group = NULL, method = NULL, ...) {
  if (is.null(group)) {
    group <- rep(1, nrow(loc))
  }
  method <- match.arg(method, c("stable", "direct"))

  ipsl <- list()
  # print(paste0("Number of polygons to integrate over: ", length(unique(group)) ))
  for (g in unique(group)) {
    gloc <- loc[group == g, , drop = FALSE]

    # Combine polygon with mesh boundary to get mesh covering the intersection.
    bnd <- INLA::inla.mesh.segment(loc = gloc, is.bnd = TRUE)
    integ <- make_stable_integration_points(mesh, bnd, ...)

    if (method %in% c("stable")) {
      # Project integration points and weights to mesh nodes
      integ <- fm_vertex_projection(integ, mesh)
    }

    # Keep points inside the mesh with positive weights
    ok <-
      fm_evaluator(mesh, integ$loc)$proj$ok &
        (integ$weight > 0)

    ips <- data.frame(integ$loc[ok, 1:2, drop = FALSE])
    colnames(ips) <- c("x", "y")
    ips$weight <- integ$weight[ok]

    ips$.block <- rep(g, nrow(ips))
    ipsl <- c(ipsl, list(ips))
  }

  do.call(rbind, ipsl)
}



# Integration points for polygons inside an inla.mesh
#
# @param mesh An inla.mesh object
# @param polylist A list of `inla.mesh.segment` objects
# @param method Which integration method to use ("stable",
#   with aggregation to mesh vertices, or "direct")
# @param ... Arguments passed to the low level integration method (`make_stable_integration_points`)
# @author Finn Lindgren \email{finn.lindgren@@gmail.com}
# @keywords internal

bru_int_polygon_old <- function(mesh, polylist, method = NULL, ...) {
  method <- match.arg(method, c("stable", "direct"))

  ipsl <- list()
  # print(paste0("Number of polygons to integrate over: ", length(polylist) ))
  for (g in seq_along(polylist)) {
    poly <- polylist[[g]]

    # Combine polygon with mesh boundary to get mesh covering the intersection.
    integ <- make_stable_integration_points(mesh, poly, ...)

    # Keep points inside the mesh with positive weights
    ok <-
      fm_evaluator(mesh, integ$loc)$proj$ok &
        (integ$weight > 0)

    if (any(ok)) {
      integ <- list(
        loc = integ$loc[ok, 1:2, drop = FALSE],
        weight = integ$weight[ok]
      )

      if (method %in% c("stable")) {
        # Project integration points and weights to mesh nodes
        integ <- integration_weight_aggregation(mesh, integ)
      }

      ips <- data.frame(
        x = integ$loc[, 1],
        y = integ$loc[, 2],
        weight = integ$weight,
        .block = g
      )

      ipsl <- c(ipsl, list(ips))
    }
  }

  do.call(rbind, ipsl)
}


# New integration methods ----





#' Aggregate integration weights onto mesh nodes
#'
#' `r lifecycle::badge("deprecated")` Use [fm_vertex_projection()] instead.
#'
#' @param mesh Mesh on which to integrate
#' @param integ `list` of `loc`, integration points,
#'   and `weight`, integration weights,
#'   or a `SpatialPointsDataFrame`. Only the coordinates and `weight` column
#'   are handled.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @keywords internal
#' @export
integration_weight_aggregation <- function(mesh, integ) {
  lifecycle::deprecate_warn(
    "2.8.0",
    "integration_weight_aggregation()",
    "fm_vertex_projection()"
  )

  fm_vertex_projection(points = integ, mesh = mesh)
}

#' Integration scheme for mesh triangle interiors
#'
#' `r lifecycle::badge("deprecated")` Use [fm_int_inla_mesh_core()] instead.
#'
#' @param mesh Mesh on which to integrate
#' @param tri_subset Optional triangle index vector for integration on a subset
#' of the mesh triangles (Default `NULL`)
#' @param nsub number of subdivision points along each triangle edge, giving
#'    `(nsub + 1)^2` proto-integration points used to compute
#'   the vertex weights
#'   (default `NULL=9`, giving 100 integration points for each triangle)
#' @return `list` with elements `loc` and `weight` with
#'   integration points for the mesh
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @keywords internal
mesh_triangle_integration <- function(mesh, tri_subset = NULL, nsub = NULL) {
  lifecycle::deprecate_warn(
    "2.8.0",
    "mesh_triangle_integration()",
    "fm_int_inla_mesh_core()"
  )

  fm_int_inla_mesh_core(mesh = mesh, tri_subset = tri_subset, nsub = NULL)
}



#' Integration points for polygons inside an inla.mesh
#'
#' @export
#' @param mesh An inla.mesh object
#' @param method Which integration method to use ("stable",
#'   with aggregation to mesh vertices, or "direct")
#' @param samplers If non-NULL, a SpatialPolygons* object
#' @param ... Arguments passed to the low level integration method (`make_triangle_integration`)
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @keywords internal

bru_int_polygon <- function(mesh,
                            method = NULL,
                            samplers = NULL,
                            ...) {
  method <- match.arg(method, c("stable", "direct"))

  ipsl <- list()

  # Compute direct integration points
  # TODO: Allow blockwise construction to avoid
  # overly large temporary coordinate matrices (via tri_subset)
  integ <- fm_int_inla_mesh_core(mesh, ...)

  # Keep points with positive weights (This should be all,
  # but if there's a degenerate triangle, this gets rid of it)
  ok <- (integ$weight > 0)
  integ$loc <- integ$loc[ok, , drop = FALSE]
  integ$weight <- integ$weight[ok]

  if (!is.null(samplers)) {
    mesh_crs <- fm_CRS(mesh)
    samplers_crs <- fm_CRS(samplers)
    integ_sp <- sp::SpatialPoints(integ$loc, proj4string = mesh_crs)
    if (!identical(mesh_crs, samplers_crs) &&
      !fm_crs_is_null(mesh_crs) &&
      !fm_crs_is_null(samplers_crs)) {
      integ_sp <- fm_transform(integ_sp,
        crs = samplers_crs,
        passthrough = TRUE
      )
    }

    #    idx_ <- sp::over(samplers, integ_sp, returnList = TRUE)
    samplers_sf <- sf::st_as_sf(samplers)
    integ_sp_sf <- sf::st_as_sf(integ_sp)
    idx <- sf::st_contains(samplers_sf, integ_sp_sf, sparse = TRUE)

    for (g in seq_along(idx)) {
      if (length(idx[[g]]) > 0) {
        integ_ <- list(
          loc = integ$loc[idx[[g]], , drop = FALSE],
          weight = integ$weight[idx[[g]]]
        )

        if (method %in% c("stable")) {
          # Project integration points and weights to mesh nodes
          integ_ <- fm_vertex_projection(integ_, mesh)
        }

        if (ncol(integ_$loc) > 2) {
          ips <- data.frame(
            x = integ_$loc[, 1],
            y = integ_$loc[, 2],
            z = integ_$loc[, 3],
            weight = integ_$weight,
            .block = g
          )
        } else {
          ips <- data.frame(
            x = integ_$loc[, 1],
            y = integ_$loc[, 2],
            weight = integ_$weight,
            .block = g
          )
        }

        ipsl <- c(ipsl, list(ips))
      }
    }
  } else {
    if (method %in% c("stable")) {
      # Project integration points and weights to mesh nodes
      integ <- fm_vertex_projection(integ, mesh)
    }

    if (ncol(integ$loc) > 2) {
      ipsl <- list(data.frame(
        x = integ$loc[, 1],
        y = integ$loc[, 2],
        z = integ$loc[, 3],
        weight = integ$weight,
        .block = 1
      ))
    } else {
      ipsl <- list(data.frame(
        x = integ$loc[, 1],
        y = integ$loc[, 2],
        weight = integ$weight,
        .block = 1
      ))
    }
  }

  do.call(rbind, ipsl)
}







# From plotsample test, comparing before and after code refactor to avoid
# recreating the full integration scheme for every polygon
#
# print(bench::mark(
#   A=bru_int_polygon_old(domain,
#     polylist=poly_segm,method=int.args$method,
#     nsub=int.args$nsub2),
#   B=bru_int_polygon(
#     domain,
#     polylist = poly_segm,
#     method = int.args$method,
#     nsub = int.args$nsub2),
#   check = FALSE))
## A tibble: 2 × 13
## expression      min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result memory time  gc
## <bch:expr> <bch:tm> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list> <list> <lis> <lis>
##  1 A            14.03s 14.03s    0.0713        NA     2.28     1    32     14.03s <NULL> <NULL> <ben… <tib…
##  2 B             8.05s  8.05s    0.124         NA     2.61     1    21      8.05s <NULL> <NULL> <ben… <tib…

# With a SPDF samplers object:
#
# bench::mark(
#   ips_old = ipoints(
#     samplers = gorillas$plotsample$plots,
#     domain = gorillas$mesh,
#     int.args = list(use_new = FALSE)
#   ),
#   ips_new = ipoints(
#     samplers = gorillas$plotsample$plots,
#     domain = gorillas$mesh,
#     int.args = list(use_new = TRUE)
#   ),
#   check = FALSE)
# expression      min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result memory time    gc
# <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list> <list> <list>  <list>
#   1 ips_old      10.04s   10.04s    0.0996        NA    1.10      1    11     10.04s <NULL> <NULL> <bench… <tibbl…
# 2 ips_new       2.91s    2.91s    0.343         NA    0.343     1     1      2.91s <NULL> <NULL> <bench… <tibbl…
