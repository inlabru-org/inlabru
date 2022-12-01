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
  if (inherits(loc, "Spatial")) {
    loc <- coordinates(loc)
  }
  if (!is.null(mesh.coords) && is.data.frame(loc)) {
    loc <- as.matrix(loc[, mesh.coords, drop = FALSE])
  }
  p2m <- INLA::inla.fmesher.smorg(loc = mesh$loc, tv = mesh$graph$tv, points2mesh = loc)
  return(!(p2m$p2m.t == 0))
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
    object$crs <- fm_as_sp_crs(object$crs)
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
  if (!identical(mesh$manifold, "R2")) {
    stop("inlabru::pixels() currently works for R2 meshes only.")
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
  coordinates(pixels) <- c("x", "y")
  if (!fm_crs_is_null(mesh$crs)) {
    proj4string(pixels) <- mesh$crs
  }

  if (is.logical(mask)) {
    if (mask) {
      pixels <- pixels[is.inside(mesh, coordinates(pixels))]
    }
  } else {
    if (inherits(mask, "SpatialPolygonsDataFrame")) {
      mask <- as(mask, "SpatialPolygons")
    }
    pixels <- pixels[!is.na(sp::over(pixels, mask))]
  }

  pixels <- sp::SpatialPixelsDataFrame(
    pixels,
    data = data.frame(matrix(nrow = NROW(pixels), ncol = 0))
  )

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
#' * Warning: does not reconstruct interior boundary
#' * Warning2: Works in euclidean coordinates. Not suitable for sphere.
#' * Warning3: Experimental; will be replaced by a new method
#'
#' @aliases tsplit.inla.mesh
#' @keywords internal
#' @param mesh an inla.mesh object
#' @param n number of splitting recursions
#' @return mesh A refined inla.mesh object
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}

tsplit.inla.mesh <- function(mesh, n = 1) {
  warning("tsplit.inla.mesh() is experimental and will be replaced by a new method")

  n <- 1

  p1 <- mesh$loc[mesh$graph$tv[, 1], ]
  p2 <- mesh$loc[mesh$graph$tv[, 2], ]
  p3 <- mesh$loc[mesh$graph$tv[, 3], ]

  m1 <- p1 + 0.5 * (p2 - p1)
  m2 <- p1 + 0.5 * (p3 - p1)
  m3 <- p2 + 0.5 * (p3 - p2)
  all.loc <- rbind(mesh$loc, m1, m2, m3)

  # TODO: Make this compliant with the inla.mesh boundary specifications;
  #   - Order is not guaranteed
  #   - Multiple disconnected components may occur
  bnd.mid <- mesh$loc[mesh$segm$bnd$idx[, 1], ] + 0.5 * (mesh$loc[mesh$segm$bnd$idx[, 2], ] - mesh$loc[mesh$segm$bnd$idx[, 1], ])
  all.bnd <- rbind(mesh$segm$bnd$loc, bnd.mid)


  mesh2 <- INLA::inla.mesh.create(loc = all.loc, boundary = all.bnd)

  if (n == 1) {
    return(mesh2)
  } else {
    return(tsplit.inla.mesh(mesh2, n - 1))
  }
}
