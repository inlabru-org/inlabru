# GENERICS
# refine = function(...) {UseMethod("refine")}
# tsplit = function(...) {UseMethod("tsplit")}


#' Query if a point is inside the mesh boundary
#'
#'
#' @aliases is.inside
#' @export
#' @keywords internal
#' @param mesh an inla.mesh object.
#' @param loc Points in space stored either as data.frame, a two-column matrix of x and y coordinates or a SpatialPoints object.
#' @param mesh.coords Coordinate names of the mesh. Use only if loc is a data.frame with respective column names.
#' @return Single column matrix of Boolean values indicating if a point is inside the mesh.
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'
#' @examples
#' \dontrun{
#' if (bru_safe_inla(quietly = TRUE)) {
#'   # Load Gorilla data
#'
#'   data("gorillas", package = "inlabru")
#'
#'   # Check if all Gorilla nests are inside the mesh
#'
#'   all(is.inside(gorillas$mesh, gorillas$nests))
#'
#'   # Also works for locations not stored as SpatialPoints object
#'
#'   loc <- coordinates(gorillas$nests)
#'   all(is.inside(gorillas$mesh, loc))
#' }
#' }
#'
is.inside <- function(mesh, loc, mesh.coords = NULL) {
  if (!is.null(mesh.coords)) {
    loc <- as.matrix(loc[, mesh.coords])
  }
  proj <- fm_evaluator(mesh, loc = loc)

  return(proj$proj$ok)
}

# Query if a point is inside a polygon AND inside the mesh;
#
#
# @aliases is.in.polygon
# @export
# @param mesh an inla.mesh object
# @param ploc Points defining a polygon
# @param loc Points to quer
# @param mask.mesh Mask points outside mesh, default: TRUE
# @param mesh.coords Coordinate names of the mesh. Use only if loc is a data.frame with respective column names.
# @return inside Boolean, TRUE if inside polygon
# @author Fabian E. Bachl \email{f.e.bachl@@bath.ac.uk}

is.inside.polygon <- function(mesh, ploc, loc, mesh.coords = NULL, mask.mesh = TRUE) {
  if (!is.null(mesh.coords) && is.data.frame(loc)) {
    loc <- as.matrix(loc[, mesh.coords, drop = FALSE])
  }
  if (!is.null(mesh.coords) && is.data.frame(ploc)) {
    ploc <- as.matrix(ploc[, mesh.coords, drop = FALSE])
  }

  mask <- sp::point.in.polygon(loc[, 1], loc[, 2], ploc[, 1], ploc[, 2]) > 0
  if (mask.mesh) {
    mask2 <- is.inside(mesh, loc)
    return(mask & mask2)
  } else {
    return(mask)
  }
}

#' @title Extract vertex locations from an `inla.mesh`
#'
#' @description Converts the vertices of an `inla.mesh` object into a `SpatialPointsDataFrame`.
#'
#' @aliases vertices.inla.mesh
#' @export
#' @param object An `inla.mesh` object.
#' @return A SpatialPointsDataFrame of mesh vertex locations. The `vrt` column indicates the internal vertex id.
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'
#' @examples
#' \donttest{
#' if (require(ggplot2, quietly = TRUE)) {
#'   data("mrsea", package = "inlabru")
#'   vrt <- vertices.inla.mesh(mrsea$mesh)
#'   ggplot() +
#'     gg(mrsea$mesh) +
#'     gg(vrt, color = "red")
#' }
#' }
#'
vertices.inla.mesh <- function(object) {
  object$crs <- fm_crs(object$crs)
  if (is.na(object$crs)) {
    object$crs <- sp::CRS(NA_character_)
  } else {
    object$crs <- fm_CRS(object$crs)
  }

  vrt <- data.frame(object$loc)
  if (!is.null(colnames(vrt))) {
    colnames(vrt) <- c("x", "y", "z")
  }
  if (all(vrt[, 3] == 0)) {
    vrt <- vrt[, 1:2]
  }
  vrt <- SpatialPoints(vrt, proj4string = object$crs)
  vrt <- SpatialPointsDataFrame(
    vrt,
    data = data.frame(vertex = seq_len(nrow(object$loc)))
  )

  vrt # return
}


#' @title Generate `SpatialPixels` covering an `inla.mesh`
#'
#' @description Generate `SpatialPixels` covering an `inla.mesh`
#'
#' @aliases pixels
#' @export
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'
#' @param mesh An `inla.mesh` object
#' @param nx Number of pixels in x direction
#' @param ny Number of pixels in y direction
#' @param mask If logical and TRUE, remove pixels that are outside the mesh.
#' If `mask` is a `Spatial` object, only return pixels covered by this object.
#' @return `SpatialPixelsDataFrame` covering the mesh
#'
#' @examples
#' \donttest{
#' if (require(ggplot2, quietly = TRUE)) {
#'   data("mrsea", package = "inlabru")
#'   pxl <- pixels(mrsea$mesh, nx = 50, ny = 50, mask = mrsea$boundary)
#'   ggplot() +
#'     gg(pxl, fill = "grey", alpha = 0.5) +
#'     gg(mrsea$mesh)
#' }
#' }
#'
pixels <- function(mesh, nx = 150, ny = 150, mask = TRUE) {
  lifecycle::deprecate_soft(
    "2.8.0",
    "pixels()",
    "fm_pixels(format = 'sp')",
    details = "The fm_pixels() function can generate sf, terra, and sp output."
  )
  fm_pixels(mesh, nx = nx, ny = ny, mask = mask, format = "sp")
}



#' @title Generate lattice points covering a mesh
#'
#' @description Generate `terra`, `sf`, or `sp` lattice locations
#'
#' @aliases pixels
#' @export
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com} and
#'  Finn Lindgren \email{finn.lindgren@@gmail.com}
#'
#' @param mesh An `inla.mesh` object
#' @param nx Number of pixels in x direction
#' @param ny Number of pixels in y direction
#' @param mask If logical and TRUE, remove pixels that are outside the mesh.
#' If `mask` is an `sf` or `Spatial` object, only return pixels covered by this object.
#' @param format character; "sf", "terra" or "sp"
#' @return `sf`, `SpatRaster`, or `SpatialPixelsDataFrame` covering the mesh
#'
#' @examples
#' \donttest{
#' if (require(ggplot2, quietly = TRUE)) {
#'   data("mrsea", package = "inlabru")
#'   pxl <- fm_pixels(mrsea$mesh, nx = 50, ny = 50, mask = mrsea$boundary,
#'                    format = "terra")
#'   ggplot() +
#'     gg(pxl, fill = "grey", alpha = 0.5) +
#'     gg(mrsea$mesh)
#' }
#' }
#'
fm_pixels <- function(mesh, nx = 150, ny = 150, mask = TRUE,

                      format = "sf") {
  format <- match.arg(format, c("sf", "terra", "sp"))
  if (!identical(mesh$manifold, "R2")) {
    stop("inlabru::fm_pixels() currently works for R2 meshes only.")
  }

  if (length(nx) == 1) {
    x <- seq(min(mesh$loc[, 1]), max(mesh$loc[, 1]), length = nx)
  } else {
    x <- nx
  }
  if (length(ny) == 1) {
    y <- seq(min(mesh$loc[, 2]), max(mesh$loc[, 2]), length = ny)
  } else {
    y <- ny
  }

  pixels <- expand.grid(x = x, y = y)
  pixels <- sf::st_as_sf(pixels, coords = c("x", "y"), crs = fm_crs(mesh))

  pixels_within <- rep(TRUE, NROW(pixels))
  if (is.logical(mask)) {
    if (mask) {
      pixels_within <- is.inside(mesh, pixels)
      pixels <- pixels[pixels_within, , drop = FALSE]
    }
  } else {
    if (inherits(mask, "SpatialPolygonsDataFrame")) {
      mask <- as(mask, "SpatialPolygons")
    }
    mask <- sf::st_as_sf(mask)
    pixels_within <- sf::st_within(pixels, mask)
    pixels <- pixels[lengths(pixels_within) > 0, , drop = FALSE]
  }

  if (identical(format, "sp")) {
    pixels <- as(pixels, "Spatial")
    pixels <- as(pixels, "SpatialPixelsDataFrame")
  } else if (identical(format, "terra")) {
    pixels <- as(pixels, "Spatial")
    pixels <- as(pixels, "SpatialPixels")
    pixels <- terra::rast(pixels)
    pixels$.mask <- as.vector(as.matrix(pixels_within))
  }

  pixels
}



# Triangle indices of points given a mesh
#
# @aliases triangle
# @export
# @param mesh A inla.mesh
# @param loc Locations using the coordinate system of the mesh
# @return tri Triangle indices
# @examples \\dontrun{}
# @author Fabian E. Bachl \email{f.e.bachl@@bath.ac.uk}

triangle <- function(mesh, loc) {
  warning("'triangle()' is an internal unused function. To be replaced by method using inla.mesh.smorg, and later fmesher")
  mcross <- function(a, b) {
    return((a[, 1]) * (b[2]) - (a[, 2]) * (b[1]))
  }
  tri <- numeric(length = dim(loc)[1])
  tv <- mesh$graph$tv
  for (j in seq_len(nrow(tv))) {
    v <- mesh$loc[tv[j, ], c(1, 2)]

    a <- v[1, ]
    b <- v[2, ]
    c <- v[3, ]

    ap <- loc - cbind(rep(a[1], nrow(loc)), rep(a[2], nrow(loc)))
    bp <- loc - cbind(rep(b[1], nrow(loc)), rep(b[2], nrow(loc)))
    cp <- loc - cbind(rep(c[1], nrow(loc)), rep(c[2], nrow(loc)))

    ab <- a - b
    bc <- b - c
    ca <- c - a

    c1 <- mcross(ap, ab)
    c2 <- mcross(bp, bc)
    c3 <- mcross(cp, ca)

    # AP x AB, BP x BC, and CP x CA must have same sign
    inside <- which(((sign(c1) == -1) & (sign(c2) == -1) & (sign(c3) == -1)) | ((sign(c1) == 1) & (sign(c2) == 1) & (sign(c3) == 1)))

    tri[inside] <- j
  }
  return(tri)
}




#' Refine an inla.mesh object
#'
#'
#' @aliases refine.inla.mesh
#' @keywords internal
#'
#' @param mesh an inla.mesh object
#' @param refine A list of refinement options passed on to
#' `INLA::inla.mesh.create`
#' @return mesh A refined inla.mesh object
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}

refine.inla.mesh <- function(mesh, refine = list(max.edge = 1)) {
  warning("refine.inla.mesh() is experimental and will be replaced by a new method")
  rmesh <- INLA::inla.mesh.create(
    loc = mesh$loc,
    interior = INLA::inla.mesh.interior(mesh),
    boundary = INLA::inla.mesh.boundary(mesh),
    refine = refine
  )
  return(rmesh)
}

#' Split triangles of a mesh into four triangles
#'
#' * Warning: The original triangle edges are not always preserved.
#' * Warning: Works in euclidean coordinates. Not suitable for sphere.
#' * Warning: Experimental; will be replaced by a new method
#'
#' @aliases tsplit.inla.mesh
#' @keywords internal
#' @param mesh an inla.mesh object
#' @param n number of splitting recursions
#' @return mesh A refined inla.mesh object
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}

tsplit.inla.mesh <- function(mesh, n = 1) {
#  lifecycle::deprecate_soft(
#    "2.8.0",
#    "tsplit.inla.mesh()",
#    "fm_subdivide()",
#    details = "tsplit.inla.mesh(mesh, n) has been replaced by fm_subdivide(mesh, n = 2^n - 1)"
#  )

  p1 <- mesh$loc[mesh$graph$tv[, 1], ]
  p2 <- mesh$loc[mesh$graph$tv[, 2], ]
  p3 <- mesh$loc[mesh$graph$tv[, 3], ]

  m1 <- (p2 + p3) / 2
  m2 <- (p1 + p3) / 2
  m3 <- (p1 + p2) / 2
  tri.loc <- rbind(mesh$loc, m1, m2, m3)

  split.edges <- function(segm) {
    if (is.null(segm) || (nrow(segm$idx) == 0)) {
      return(segm)
    }
    n.loc <- nrow(segm$loc)
    n.idx <- nrow(segm$idx)
    loc <- rbind(segm$loc,
                 (segm$loc[segm$idx[, 1], ] +
                    segm$loc[segm$idx[, 2], ]) / 2)
    idx <- rbind(
      cbind(segm$idx[, 1], n.loc + seq_len(n.idx)),
      cbind(n.loc + seq_len(n.idx), segm$idx[, 2])
    )

    segm2 <-
      INLA::inla.mesh.segment(loc = loc,
                              idx = idx,
                              grp = c(segm$grp, segm$grp),
                              is.bnd = segm$is.bnd,
                              crs = fm_CRS(segm))

    segm2
  }

  interior2 <- split.edges(INLA::inla.mesh.interior(mesh)[[1]])
  boundary2 <- split.edges(INLA::inla.mesh.boundary(mesh)[[1]])

  mesh2 <- INLA::inla.mesh.create(loc = tri.loc,
                                  interior = interior2,
                                  boundary = boundary2,
                                  crs = fm_CRS(mesh))

  if (n <= 1) {
    return(mesh2)
  } else {
    return(tsplit.inla.mesh(mesh2, n - 1))
  }
}


# Split triangles of a mesh into subtriangles
#
# @param mesh an inla.mesh object
# @param n number of added points along each edge
# @return A refined inla.mesh object
# @author Finn Lindgren \email{finn.lindgren@@gmail.com}
# @export

fm_subdivide <- function(mesh, n = 1) {
  if (n < 1) {
    return(mesh)
  }

  split.edges <- function(segm, n) {
    if (is.null(segm) || (nrow(segm$idx) == 0)) {
      return(segm)
    }
    n.loc <- nrow(segm$loc)
    n.idx <- nrow(segm$idx)
    loc <- do.call(rbind,
                   c(list(segm$loc),
                     lapply(seq_len(n),
                            function(k) {
                              (segm$loc[segm$idx[, 1], ] * k / (n + 1) +
                                 segm$loc[segm$idx[, 2], ] * (n - k + 1) / (n +
                                                                              1))
                            })))
    idx <- do.call(rbind,
                   c(list(cbind(
                     segm$idx[, 1], n.loc + seq_len(n.idx)
                   )),
                   lapply(
                     seq_len(n - 1),
                     function(k) {
                     cbind(n.loc * k + seq_len(n.idx),
                           n.loc * (k + 1) + seq_len(n.idx))
                     }
                   ),
                   list(cbind(
                     n.loc * n + seq_len(n.idx), segm$idx[, 2]
                   ))))

    segm2 <-
      INLA::inla.mesh.segment(loc = loc,
                              idx = idx,
                              grp = rep(segm$grp, n + 1),
                              is.bnd = segm$is.bnd)

    segm2
  }

  p1 <- mesh$loc[mesh$graph$tv[, 1], ]
  p2 <- mesh$loc[mesh$graph$tv[, 2], ]
  p3 <- mesh$loc[mesh$graph$tv[, 3], ]

  tri.inner.loc <-
    do.call(rbind,
            lapply(seq_len(n + 2) - 1,
                   function(k2) {
                     do.call(rbind,
                             lapply(seq_len(n + 2- k2) - 1,
                                    function(k3) {
                                      w1 <- (n + 1 - k2 - k3) / (n + 1)
                                      w2 <- k2 / (n + 1)
                                      w3 <- k3 / (n + 1)
                                      p1 * w1 + p2 * w2 + p3 * w3
                                    }))
                   }))

  n.tri <- nrow(p1)
  tri.edges <- INLA::inla.mesh.segment(
    loc = rbind(p1, p2, p3),
    idx = rbind(cbind(seq_len(n.tri), seq_len(n.tri) + n.tri),
                cbind(seq_len(n.tri) + n.tri, seq_len(n.tri) + 2 * n.tri),
                cbind(seq_len(n.tri) + 2 * n.tri, seq_len(n.tri))),
    is.bnd = FALSE
  )

  boundary2 <- split.edges(INLA::inla.mesh.boundary(mesh)[[1]], n = n)

  if (identical(mesh$manifold, "S2")) {
    radius <- mean(rowSums(mesh$loc^2)^0.5)
    renorm <- function(loc) {
      loc * (radius / rowSums(loc^2)^0.5)
    }
    tri.inner.loc <- renorm(tri.inner.loc)
    tri.edges2$loc <- renorm(tri.edges2$loc)
    boundary2$loc <- renorm(boundary2$loc)
  }

  mesh2 <- INLA::inla.mesh.create(loc = tri.inner.loc,
                                  interior = tri.edges,
                                  boundary = boundary2,
                                  refine = list(min.angle = 0,
                                                max.edge = Inf),
                                  crs = fm_CRS(mesh))

  mesh2
}
