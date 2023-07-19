#' @include fmesher_mesh.R
#' @include deprecated.R

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
  fm_segm(
    loc = out.loc, idx = out.idx, grp = out.grp, is.bnd = is.bnd,
    crs = crs
  )
}




#' @export
#' @describeIn inlabru-deprecated `r lifecycle::badge("deprecated")` in favour of
#' [fm_as_segm()]
fm_sp2segment <-
  function(...) {
    fm_as_segm(...)
  }




#' @export
#' @rdname fm_as_segm
#' @param crs A crs object
#' @param closed logical; whether to treat a point sequence as a closed polygon.
#' Default: `FALSE`
fm_as_segm.matrix <-
  function(x, reverse = FALSE, grp = NULL, is.bnd = FALSE, crs = NULL, closed = FALSE, ...) {
    loc <- x
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
    fm_segm(
      loc = loc, idx = idx, grp = grp, is.bnd = is.bnd, crs = fm_CRS(crs)
    )
  }


#' @export
#' @rdname fm_as_segm
fm_as_segm.SpatialPoints <-
  function(x, reverse = FALSE, grp = NULL, is.bnd = TRUE, closed = FALSE, ...) {
    crs <- fm_CRS(x)
    loc <- sp::coordinates(x)

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
    fm_segm(
      loc = loc, idx = idx, grp = grp, is.bnd = is.bnd, crs = crs
    )
  }

#' @export
#' @rdname fm_as_segm
fm_as_segm.SpatialPointsDataFrame <-
  function(x, ...) {
    fm_as_segm.SpatialPoints(x, ...)
  }



#' @export
#' @rdname fm_as_segm
fm_as_segm.Line <-
  function(x, reverse = FALSE, grp = NULL, crs = NULL, ...) {
    loc <- x@coords
    n <- dim(loc)[1L]
    if (reverse) {
      idx <- seq(n, 1L, length.out = n)
      if (!is.null(grp)) {
        grp <- rev(grp)
      }
    } else {
      idx <- seq_len(n)
    }
    fm_segm(
      loc = loc, idx = idx, grp = grp, is.bnd = FALSE, crs = fm_CRS(crs)
    )
  }

#' @export
#' @rdname fm_as_segm
fm_as_segm.Lines <-
  function(x, join = TRUE, grp = NULL, crs = NULL, ...) {
    segm <- as.list(lapply(
      seq_len(length(x@Lines)),
      function(k) {
        x <- x@Lines[[k]]
        if (!is.null(grp)) {
          grp_ <- grp[k]
        } else {
          grp_ <- NULL
        }
        fm_as_segm(x, grp = grp_, crs = crs, ...)
      }
    ))
    if (join) {
      segm <- fm_internal_sp2segment_join(segm, grp = grp, closed = FALSE)
    }
    segm
  }

#' @export
#' @rdname fm_as_segm
fm_as_segm.SpatialLines <-
  function(x, join = TRUE, grp = NULL, ...) {
    crs <- fm_CRS(x)
    segm <- list()
    for (k in seq_len(length(x@lines))) {
      segm[[k]] <- fm_as_segm(x@lines[[k]],
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
#' @rdname fm_as_segm
fm_as_segm.SpatialLinesDataFrame <-
  function(x, ...) {
    fm_as_segm.SpatialLines(x, ...)
  }

#' @export
#' @rdname fm_as_segm
fm_as_segm.SpatialPolygons <-
  function(x, join = TRUE, grp = NULL, ...) {
    crs <- fm_CRS(x)
    segm <- list()
    for (k in seq_len(length(x@polygons))) {
      segm[[k]] <-
        fm_as_segm(x@polygons[[k]], join = TRUE, crs = crs)
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
#' @rdname fm_as_segm
fm_as_segm.SpatialPolygonsDataFrame <-
  function(x, ...) {
    fm_as_segm.SpatialPolygons(x, ...)
  }

#' @export
#' @rdname fm_as_segm
fm_as_segm.Polygons <-
  function(x, join = TRUE, crs = NULL, grp = NULL, ...) {
    segm <- as.list(lapply(
      x@Polygons,
      function(x) fm_as_segm(x, crs = fm_CRS(crs), ...)
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
#' @rdname fm_as_segm
fm_as_segm.Polygon <-
  function(x, crs = NULL, ...) {
    loc <- x@coords[-dim(x@coords)[1L], , drop = FALSE]
    n <- dim(loc)[1L]
    if (x@hole) {
      if (x@ringDir == 1) {
        idx <- c(1L:n, 1L)
      } else {
        idx <- c(1L, seq(n, 1L, length.out = n))
      }
    } else if (x@ringDir == 1) {
      idx <- c(1L, seq(n, 1L, length.out = n))
    } else {
      idx <- c(1L:n, 1L)
    }
    fm_segm(
      loc = loc,
      idx = idx,
      is.bnd = TRUE,
      crs = fm_CRS(crs)
    )
  }
