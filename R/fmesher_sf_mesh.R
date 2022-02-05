#' @title Coercion methods to and from meshes
#' @rdname fm_as
#' @param ... Arguments passed on to other methods
#' @export
fm_as_sfc <- function(x) {
  UseMethod("fm_as_sfc")
}


#' @rdname fm_as
#' @export
fm_as_inla_mesh <- function(...) {
  UseMethod("fm_as_inla_mesh")
}


#' @rdname fm_as
#' @export
fm_as_sf_crs <- function(x) {
  if (inherits(x, "crs")) {
    x
  } else if (inherits(x, "CRS")) {
    sf::st_crs(x)
  } else if (is.null(x)) {
    sf::st_crs(x)
  } else {
    warning(paste0("Unsupported source crs class ",
                   paste(class(x), sep = ",")),
            immediate. = TRUE)
    x
  }
}
#' @rdname fm_as
#' @export
fm_as_sp_crs <- function(x) {
  if (inherits(x, "CRS")) {
    x
  } else if (inherits(x, "crs")) {
    if (is.na(x)) {
      NULL
    } else {
      fm_CRS(SRS_string = x$wkt)
    }
  } else if (is.null(x)) {
    NULL
  } else {
    warning(paste0("Unsupported source crs class ",
                 paste(class(x), sep = ",")),
          immediate. = TRUE)
  x
}
}

#' @rdname fm_as
#' @aliases fm_as_sfc fm_as_sfc.inla.mesh
#'
#' @param x An `inla.mesh` mesh object, or an `sfg` `MULTIPOLYGON` object
#' @returns * `fm_as_sfc`: An `sfc_MULTIPOLYGON` object
#' @exportS3Method fm_as_sfc inla.mesh
#' @export
fm_as_sfc.inla.mesh <- function(x, ...) {
  stopifnot(inherits(x, "inla.mesh"))
  geom <- sf::st_geometry(
    sf::st_multipolygon(
      lapply(
        seq_len(nrow(x$graph$tv)),
        function(k) {
          list(x$loc[x$graph$tv[k, c(1, 2, 3, 1)], , drop = FALSE])
        }
      ),
      dim = "XYZ"
    )
  )
  sf::st_crs(geom) <- fm_as_sf_crs(x$crs)
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
            immediate. = TRUE)
  }
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
  crs <- fm_as_sp_crs(sf::st_crs(x))
  mesh <- INLA::inla.mesh.create(loc = loc, tv = tv, ...,
                                 crs = crs)
  mesh
}





#' @rdname fm_as
#' @export
fm_as_inla_mesh_segment.sfc_POINT <-
  function(sfc, reverse = FALSE, grp = NULL, is.bnd = TRUE, ...) {
    crs <- sf::st_crs(sfc)

#    if (st_check_dim(sfc)) {
#      warning(
#        "XYZ, XYM and XYZM sfg classes are not fully supported. In general the Z and M coordinates will be ignored"
#      )
#    }

    crs <- fm_as_sp_crs(crs) # required for INLA::inla.mesh.segment

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
#' @export
fm_as_inla_mesh_segment.sfc_LINESTRING <-
  function(sfc, join = TRUE, grp = NULL, reverse = FALSE, ...) {
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
    crs <- fm_as_sp_crs(crs) # required for INLA::inla.mesh.segment

    segm <- list()
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
        is.bnd = FALSE,
        crs = crs
      )
    }

    # Think this can stay the same?
    if (join) {
      if (missing(grp)) {
        grp <- seq_len(length(segm))
      }
      segm <- fm_internal_sp2segment_join(segm, grp = grp, closed = FALSE)
    }
    segm
  }

#' @rdname fm_as
#' @export
fm_as_inla_mesh_segment.sfc_POLYGON <-
  function(sp, join = TRUE, grp = NULL, ...) {
    crs <- sf::st_crs(sfc)

    if (st_check_dim(sfc)) {
      warning(
        "XYZ, XYM and XYZM are not fully supported. In general the Z and M coordinates will be ignored"
      )
    }

    crs <- fm_as_sp_crs(crs) # required for INLA::inla.mesh.segment

    segm <- list()
    for (k in seq_len(length(sfc))) {
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
      # segm[[k]] <- fm_as_inla_mesh_segment(sp@polygons[[k]], join = TRUE, crs = crs)
    }

    if (join) {
      if (missing(grp)) {
        grp <- seq_len(length(segm))
      }
      segm <- fm_internal_sp2segment_join(segm, grp = grp)
    }
    segm
  }

#' @rdname fm_as
#' @export
fm_as_inla_mesh_segment.sf <-
  function(sf, reverse = FALSE, grp = NULL, is.bnd = TRUE, ...) {
    sfc <- sf::st_geometry(sf)
    fm_as_inla_mesh_segment(sfc)
  }

#' @rdname fm_as
#' @export
fm_as_inla_mesh.sf <-
  function(sf, ...) {
    sfc <- sf::st_geometry(sf)
    fm_as_inla_mesh(sfc)
  }
