#' @importFrom sp coordinates proj4string `proj4string<-`

## Input: list of segments
fm_internal_sp2segment_join <- function(inp, grp = NULL) {
  crs <- NULL
  if (length(inp) > 0) {
    out.loc <- matrix(0, 0, ncol(inp[[1]]$loc))
    is.bnd <- inp[[1]]$is.bnd
    for (k in seq_along(inp)) {
      if (inp[[k]]$is.bnd != is.bnd) {
        stop("Cannot join a mix of boundary and interior segment")
      }
      crs <- fm_internal_update_crs(crs, inp[[k]]$crs, mismatch.allowed = FALSE)
    }
  } else {
    is.bnd <- FALSE
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
    if (ncol(as.matrix(inp.idx)) == 1) {
      stop("inla.mesh.segment idx should be a matrix but a vector has been detected.")
    }
    out.loc <- rbind(out.loc, inp.loc)
    out.idx <- rbind(out.idx, inp.idx + offset)
    if (!is.null(grp)) {
      out.grp <- c(out.grp, rep(grp[k], nrow(inp.idx)))
    } else {
      out.grp <- c(out.grp, inp.grp)
    }
  }
  INLA::inla.mesh.segment(
    loc = out.loc, idx = out.idx, grp = out.grp, is.bnd = is.bnd,
    crs = crs
  )
}


fm_as_inla_mesh_segment <-
  function(...) {
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

fm_as_inla_mesh_segment.SpatialPointsDataFrame <-
  function(sp, ...) {
    fm_as_inla_mesh_segment.SpatialPoints(sp, ...)
  }



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
    INLA::inla.mesh.segment(loc = loc, idx = idx, grp = grp, is.bnd = FALSE, crs = crs)
  }

fm_as_inla_mesh_segment.Lines <-
  function(sp, join = TRUE, crs = NULL, ...) {
    segm <- as.list(lapply(
      sp@Lines,
      function(x) fm_as_inla_mesh_segment(x, crs = crs, ...)
    ))
    if (join) {
      segm <- fm_internal_sp2segment_join(segm, grp = NULL)
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
      segm <- fm_internal_sp2segment_join(segm, grp = grp)
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
      function(x) fm_as_inla_mesh_segment(x, crs = crs, ...)
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
