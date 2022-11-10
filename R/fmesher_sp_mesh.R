#' @importFrom sp coordinates proj4string `proj4string<-`

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
  is.bnd <- all(vapply(inp, function(x) x$is.bnd, TRUE))
  INLA::inla.mesh.segment(
    loc = out.loc, idx = out.idx, grp = out.grp, is.bnd = is.bnd,
    crs = crs
  )
}




#' @export
#' @rdname fm_as
#' @param sp An `sp` style S4 object to be converted
fm_sp2segment <-
  function(sp, ...) {
    UseMethod("fm_as_inla_mesh_segment")
  }




#' @export
#' @rdname fm_as
#' @param crs A crs object
#' @param closed logical; whether to treat a point sequence as a closed polygon.
#' Default: `FALSE`
fm_as_inla_mesh_segment.matrix <-
  function(sp, reverse = FALSE, grp = NULL, is.bnd = FALSE, crs = NULL, closed = FALSE, ...) {
    loc <- sp
    n <- dim(loc)[1L]
    if (closed) {
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
      loc = loc, idx = idx, grp = grp, is.bnd = is.bnd, crs = crs
    )
  }


#' @export
#' @rdname fm_as
fm_as_inla_mesh_segment.SpatialPoints <-
  function(sp, reverse = FALSE, grp = NULL, is.bnd = TRUE, closed = FALSE, ...) {
    crs <- fm_sp_get_crs(sp)
    loc <- coordinates(sp)

    n <- dim(loc)[1L]
    if (closed) {
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
      loc = loc, idx = idx, grp = grp, is.bnd = is.bnd, crs = crs
    )
  }

#' @export
#' @rdname fm_as
fm_as_inla_mesh_segment.SpatialPointsDataFrame <-
  function(sp, ...) {
    fm_as_inla_mesh_segment.SpatialPoints(sp, ...)
  }



#' @export
#' @rdname fm_as
fm_as_inla_mesh_segment.Line <-
  function(sp, reverse = FALSE, grp = NULL, crs = NULL, ...) {
    loc <- sp@coords
    n <- dim(loc)[1L]
    if (reverse) {
      idx <- seq(n, 1L, length = n)
      if (!is.null(grp)) {
        grp <- rev(grp)
      }
    } else {
      idx <- seq_len(n)
    }
    INLA::inla.mesh.segment(
      loc = loc, idx = idx, grp = grp, is.bnd = FALSE, crs = crs
    )
  }

#' @export
#' @rdname fm_as
fm_as_inla_mesh_segment.Lines <-
  function(sp, join = TRUE, grp = NULL, crs = NULL, ...) {
    segm <- as.list(lapply(
      seq_len(length(sp@Lines)),
      function(k) {
        x <- sp@Lines[[k]]
        if (!is.null(grp)) {
          grp_ <- grp[k]
        } else {
          grp_ <- NULL
        }
        fm_as_inla_mesh_segment(x, grp = grp_, crs = crs, ...)
      }
    ))
    if (join) {
      segm <- fm_internal_sp2segment_join(segm, grp = grp, closed = FALSE)
    }
    segm
  }

#' @export
#' @rdname fm_as
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

#' @export
#' @rdname fm_as
fm_as_inla_mesh_segment.SpatialLinesDataFrame <-
  function(sp, ...) {
    fm_as_inla_mesh_segment.SpatialLines(sp, ...)
  }

#' @export
#' @rdname fm_as
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

#' @export
#' @rdname fm_as
fm_as_inla_mesh_segment.SpatialPolygonsDataFrame <-
  function(sp, ...) {
    fm_as_inla_mesh_segment.SpatialPolygons(sp, ...)
  }

#' @export
#' @rdname fm_as
fm_as_inla_mesh_segment.Polygons <-
  function(sp, join = TRUE, crs = NULL, grp = NULL, ...) {
    segm <- as.list(lapply(
      sp@Polygons,
      function(x) fm_as_inla_mesh_segment(x, crs = crs, ...)
    ))
    if (join) {
      if (missing(grp)) {
        grp <- seq_len(length(segm))
      }
      segm <- fm_internal_sp2segment_join(segm, grp = grp)
    }
    segm
  }

#' @export
#' @rdname fm_as
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
