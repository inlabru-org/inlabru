#' @title Coercion methods to and from meshes
#' @rdname fm_as
#' @param ... Arguments passed on to other methods
#' @export
fm_as_sfc <- function(x, ...) {
  UseMethod("fm_as_sfc")
}

#' @rdname fm_as
#' @export
fm_as_inla_mesh_segment <-
  function(...) {
    UseMethod("fm_as_inla_mesh_segment")
  }

#' @rdname fm_as
#' @export
fm_as_inla_mesh <- function(...) {
  UseMethod("fm_as_inla_mesh")
}


#' @rdname fm_as
#' @aliases fm_as_sfc fm_as_sfc.inla.mesh
#'
#' @param x An object to be coerced/transformed/converted into another class
#' @param multi logical; if `TRUE`, attempt to a `sfc_MULTIPOLYGON`, otherwise
#' a set of `sfc_POLYGON`. Default `FALSE`
#' @returns * `fm_as_sfc`: An `sfc_MULTIPOLYGON` or `sfc_POLYGON` object
#' @exportS3Method fm_as_sfc inla.mesh
#' @export
fm_as_sfc.inla.mesh <- function(x, ..., multi = FALSE) {
  stopifnot(inherits(x, "inla.mesh"))
  if (multi) {
    geom <- sf::st_sfc(
      sf::st_multipolygon(
        lapply(
          seq_len(nrow(x$graph$tv)),
          function(k) {
            list(x$loc[x$graph$tv[k, c(1, 2, 3, 1)], , drop = FALSE])
          }
        ),
        dim = "XYZ"
      ),
      check_ring_dir = TRUE
    )
  } else {
    geom <- sf::st_sfc(
      lapply(
        seq_len(nrow(x$graph$tv)),
        function(k) {
          sf::st_polygon(
            list(x$loc[x$graph$tv[k, c(1, 2, 3, 1)], , drop = FALSE]),
            dim = "XYZ"
          )
        }
      )
    )
  }
  sf::st_crs(geom) <- fm_crs(x$crs)
  geom
}


#' @rdname fm_as
#' @aliases fm_as_inla_mesh fm_as_inla_mesh.sfc_MULTIPOLYGON
#'
#' @returns * `fm_as_inla_mesh`: An `inla.mesh` mesh object
#' @export
fm_as_inla_mesh.sfc_MULTIPOLYGON <- function(x, ...) {
  if (length(x) > 1) {
    warning("More than one MULTIPOLYGON detected, but conversion method only uses one.",
      immediate. = TRUE
    )
  }
  # Ensure correct CCW ring orientation; sf doesn't take into account
  # that geos has CW as canonical orientation
  x <- sf::st_sfc(x, check_ring_dir = TRUE)
  tv <- matrix(seq_len(3 * length(x[[1]])), length(x[[1]]), 3, byrow = TRUE)
  loc <- do.call(
    rbind,
    lapply(
      x[[1]],
      function(xx) {
        if ((length(xx) > 1) ||
          (nrow(xx[[1]]) > 4)) {
          stop("Invalid geometry; non-triangle detected.")
        }
        xx[[1]][1:3, , drop = FALSE]
      }
    )
  )
  crs <- fm_CRS(sf::st_crs(x))
  mesh <- INLA::inla.mesh.create(
    loc = loc, tv = tv, ...,
    crs = crs
  )
  mesh
}

#' @rdname fm_as
#' @aliases fm_as_inla_mesh fm_as_inla_mesh.sfc_POLYGON
#'
#' @export
fm_as_inla_mesh.sfc_POLYGON <- function(x, ...) {
  # Ensure correct CCW ring orientation; sf doesn't take into account
  # that geos has CW as canonical orientation
  sfc <- sf::st_sfc(x, check_ring_dir = TRUE)
  tv <- matrix(seq_len(3 * NROW(x)), NROW(x), 3, byrow = TRUE)
  loc <- do.call(
    rbind,
    lapply(
      x,
      function(xx) {
        if ((length(xx) > 1) ||
          (nrow(xx[[1]]) > 4)) {
          stop("Invalid geometry; non-triangle detected.")
        }
        xx[[1]][1:3, , drop = FALSE]
      }
    )
  )
  crs <- fm_CRS(sf::st_crs(x))
  mesh <- INLA::inla.mesh.create(
    loc = loc, tv = tv, ...,
    crs = crs
  )
  mesh
}





#' @rdname fm_as
#' @param reverse logical; When TRUE, reverse the order of the input points.
#'   Default `FALSE`
#' @param grp if non-null, should be an integer vector of grouping labels for
#'   one for each segment.
#'    Default `NULL`
#' @param is.bnd logical; if `TRUE`, set the boundary flag for the segments.
#'   Default `TRUE`
#' @export
fm_as_inla_mesh_segment.sfc_POINT <-
  function(x, reverse = FALSE, grp = NULL, is.bnd = TRUE, ...) {
    sfc <- x
    crs <- sf::st_crs(sfc)

    #    if (st_check_dim(sfc)) {
    #      warning(
    #        "XYZ, XYM and XYZM sfg classes are not fully supported. In general the Z and M coordinates will be ignored"
    #      )
    #    }

    crs <- fm_CRS(crs) # required for INLA::inla.mesh.segment

    loc <- sf::st_coordinates(sfc)
    coord_names <- intersect(c("X", "Y", "Z"), colnames(loc))
    loc <- unname(loc[, coord_names, drop = FALSE])

    n <- dim(loc)[1L]
    if (is.bnd) {
      idx <- c(seq_len(n), 1L)
    } else {
      idx <- seq_len(n)
    }
    if (reverse) {
      idx <- rev(idx)
      if (!is.null(grp)) {
        grp <- rev(grp)
      }
    }
    INLA::inla.mesh.segment(
      loc = loc, idx = idx, grp = grp, is.bnd = is.bnd,
      crs = crs
    )
  }

#' @rdname fm_as
#' @param join logical; if `TRUE`, join input segments with common vertices.
#'    Default `TRUE`
#' @export
fm_as_inla_mesh_segment.sfc_LINESTRING <-
  function(x, join = TRUE, grp = NULL, reverse = FALSE, ...) {
    sfc <- x
    # Note: Z should be fully supported in what we do with 3D coordinates ourselves.
    # It's when applying st_ methods that the check needs to be done, not when crating
    # objects, as we _do_ support 3D meshes in inla.mesh, and _should_ support those
    # in inlabru.
    #    if (st_check_dim(sfc)) {
    #      warning(
    #        "XYZ, XYM and XYZM sfg classes are not fully supported. In general the Z and M coordinates will be ignored"
    #      )
    #    }

    crs <- sf::st_crs(sfc)
    crs <- fm_CRS(crs) # required for INLA::inla.mesh.segment

    segm <- list()
    if (is.null(grp)) {
      grp <- seq_len(length(sfc))
    }
    for (k in seq_len(length(sfc))) {
      loc <- sf::st_coordinates(sfc[k])
      coord_names <- intersect(c("X", "Y", "Z"), colnames(loc))
      loc <- unname(loc[, coord_names, drop = FALSE])

      n <- dim(loc)[1L]
      if (reverse) {
        idx <- seq(n, 1L, length = n)
      } else {
        idx <- seq_len(n)
      }
      segm[[k]] <- INLA::inla.mesh.segment(
        loc = loc,
        idx = idx,
        grp = grp[k],
        is.bnd = FALSE,
        crs = crs
      )
    }

    if (join) {
      segm <- fm_internal_sp2segment_join(segm, grp = grp)
    }
    segm
  }

#' @rdname fm_as
#' @export
fm_as_inla_mesh_segment.sfc_POLYGON <-
  function(x, join = TRUE, grp = NULL, ...) {
    # Ensure correct CCW ring orientation; sf doesn't take into account
    # that geos has CW as canonical orientation
    sfc <- sf::st_sfc(x, check_ring_dir = TRUE)
    crs <- sf::st_crs(sfc)
    crs <- fm_CRS(crs) # required for INLA::inla.mesh.segment

    segm <- list()
    if (is.null(grp)) {
      grp <- seq_len(length(sfc))
    }
    for (k in seq_len(length(sfc))) {
      loc <- sf::st_coordinates(sfc[k])
      coord_names <- intersect(c("X", "Y", "Z"), colnames(loc))
      L1info <- loc[, "L1", drop = TRUE]
      L2info <- loc[, "L2", drop = TRUE]
      loc <- unname(loc[, coord_names, drop = FALSE])
      # If winding directions are correct, all info is already available
      # For 3D, cannot check winding, so must assume correct.
      segm_k <-
        lapply(
          unique(L1info),
          function(i) {
            subset <- which(L1info == i)
            # sfc_POLYGON repeats the initial point within each L1
            n <- length(subset) - 1
            subset <- subset[-(n + 1)]
            idx <- c(seq_len(n), 1L)
            INLA::inla.mesh.segment(
              loc = loc[subset, , drop = FALSE],
              idx = idx,
              grp = grp[k],
              is.bnd = TRUE,
              crs = crs
            )
          }
        )
      segm[[k]] <- fm_internal_sp2segment_join(segm_k)
    }

    if (join) {
      segm <- fm_internal_sp2segment_join(segm, grp = grp)
    }
    segm
  }

#' @rdname fm_as
#' @export
fm_as_inla_mesh_segment.sfc_MULTIPOLYGON <-
  function(x, join = TRUE, grp = NULL, ...) {
    # Ensure correct CCW ring orientation; sf doesn't take into account
    # that geos has CW as canonical orientation
    sfc <- sf::st_sfc(x, check_ring_dir = TRUE)
    crs <- sf::st_crs(sfc)
    crs <- fm_CRS(crs) # required for INLA::inla.mesh.segment

    segm <- list()
    if (is.null(grp)) {
      grp <- seq_len(length(sfc))
    }
    for (k in seq_len(length(sfc))) {
      loc <- sf::st_coordinates(sfc[k])
      coord_names <- intersect(c("X", "Y", "Z"), colnames(loc))
      Linfo <- loc[, c("L1", "L2"), drop = FALSE]
      loc <- unname(loc[, coord_names, drop = FALSE])
      # If winding directions are correct, all info is already available
      # For 3D, cannot check winding, so must assume correct.
      uniqueLinfo <- unique(Linfo)
      segm_k <-
        lapply(
          seq_len(nrow(uniqueLinfo)),
          function(i) {
            subset <- which((Linfo[, 1] == uniqueLinfo[i, 1]) &
              (Linfo[, 2] == uniqueLinfo[i, 2]))
            # sfc_POLYGON repeats the initial point
            n <- length(subset) - 1
            subset <- subset[-(n + 1)]
            idx <- c(seq_len(n), 1L)
            INLA::inla.mesh.segment(
              loc = loc[subset, , drop = FALSE],
              idx = idx,
              grp = grp[k],
              is.bnd = TRUE,
              crs = crs
            )
          }
        )
      segm[[k]] <- fm_internal_sp2segment_join(segm_k)
    }

    if (join) {
      segm <- fm_internal_sp2segment_join(segm, grp = grp)
    }
    segm
  }

#' @rdname fm_as
#' @export
fm_as_inla_mesh_segment.sfc_GEOMETRY <-
  function(x, grp = NULL, join = TRUE, ...) {
    if (is.null(grp)) {
      grp <- seq_len(length(x))
    }
    segm <-
      lapply(
        seq_along(x),
        function(k) {
          fm_as_inla_mesh_segment(x[k], grp = grp[[k]], join = join, ...)
        }
      )
    if (join) {
      segm <- fm_internal_sp2segment_join(segm, grp = grp)
    }
    segm
  }

#' @rdname fm_as
#' @export
fm_as_inla_mesh_segment.sf <-
  function(x, ...) {
    sfc <- sf::st_geometry(x)
    fm_as_inla_mesh_segment(sfc, ...)
  }

#' @rdname fm_as
#' @export
fm_as_inla_mesh.sf <-
  function(x, ...) {
    sfc <- sf::st_geometry(x)
    fm_as_inla_mesh(sfc, ...)
  }
