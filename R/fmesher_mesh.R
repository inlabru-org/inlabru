#' @include deprecated.R

#' @title Generate lattice points covering a mesh
#'
#' @description Generate `terra`, `sf`, or `sp` lattice locations
#'
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
#'   pxl <- fm_pixels(mrsea$mesh,
#'     nx = 50, ny = 50, mask = mrsea$boundary,
#'     format = "terra"
#'   )
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
    x <- seq(min(mesh$loc[, 1]), max(mesh$loc[, 1]), length.out = nx)
  } else {
    x <- nx
  }
  if (length(ny) == 1) {
    y <- seq(min(mesh$loc[, 2]), max(mesh$loc[, 2]), length.out = ny)
  } else {
    y <- ny
  }

  pixels <- expand.grid(x = x, y = y)
  pixels <- sf::st_as_sf(pixels, coords = c("x", "y"), crs = fm_crs(mesh))

  pixels_within <- rep(TRUE, NROW(pixels))
  if (is.logical(mask)) {
    if (mask) {
      pixels_within <- fm_is_within(pixels, mesh)
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
    loc <- do.call(
      rbind,
      c(
        list(segm$loc),
        lapply(
          seq_len(n),
          function(k) {
            (segm$loc[segm$idx[, 1], ] * k / (n + 1) +
              segm$loc[segm$idx[, 2], ] * (n - k + 1) / (n +
                1))
          }
        )
      )
    )
    idx <- do.call(
      rbind,
      c(
        list(cbind(
          segm$idx[, 1], n.loc + seq_len(n.idx)
        )),
        lapply(
          seq_len(n - 1),
          function(k) {
            cbind(
              n.loc * k + seq_len(n.idx),
              n.loc * (k + 1) + seq_len(n.idx)
            )
          }
        ),
        list(cbind(
          n.loc * n + seq_len(n.idx), segm$idx[, 2]
        ))
      )
    )

    segm2 <-
      INLA::inla.mesh.segment(
        loc = loc,
        idx = idx,
        grp = rep(segm$grp, n + 1),
        is.bnd = segm$is.bnd
      )

    segm2
  }

  p1 <- mesh$loc[mesh$graph$tv[, 1], ]
  p2 <- mesh$loc[mesh$graph$tv[, 2], ]
  p3 <- mesh$loc[mesh$graph$tv[, 3], ]

  tri.inner.loc <-
    do.call(
      rbind,
      lapply(
        seq_len(n + 2) - 1,
        function(k2) {
          do.call(
            rbind,
            lapply(
              seq_len(n + 2 - k2) - 1,
              function(k3) {
                w1 <- (n + 1 - k2 - k3) / (n + 1)
                w2 <- k2 / (n + 1)
                w3 <- k3 / (n + 1)
                p1 * w1 + p2 * w2 + p3 * w3
              }
            )
          )
        }
      )
    )

  n.tri <- nrow(p1)
  tri.edges <- INLA::inla.mesh.segment(
    loc = rbind(p1, p2, p3),
    idx = rbind(
      cbind(seq_len(n.tri), seq_len(n.tri) + n.tri),
      cbind(seq_len(n.tri) + n.tri, seq_len(n.tri) + 2 * n.tri),
      cbind(seq_len(n.tri) + 2 * n.tri, seq_len(n.tri))
    ),
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

  mesh2 <- INLA::inla.mesh.create(
    loc = tri.inner.loc,
    interior = tri.edges,
    boundary = boundary2,
    refine = list(
      min.angle = 0,
      max.edge = Inf
    ),
    crs = fm_CRS(mesh)
  )

  mesh2
}



#' @title store points in different formats
#'
#' @description Convert a matrix of points into different formats.
#'
#' @param format character; `"sf"`, `"df"`, `"sp"`
#' @return
#' An `sf`, `data.frame`, or `SpatialPointsDataFrame` object, with
#' optional added information.
#' @export
#' @keywords internal
fm_store_points <- function(loc, crs = NULL, info = NULL, format = NULL) {
  format <- match.arg(
    format,
    c("sf", "df", "sp")
  )

  crs <- fm_crs(crs)

  points <- as.data.frame(loc)
  colnames(points) <- c("x", "y", "z")[seq_len(ncol(points))]
  if (!fm_crs_is_null(crs) && !fm_crs_is_geocent(crs)) {
    points <- points[, 1:2, drop = FALSE]
  }

  if (identical(format, "df")) {
    points <- cbind(points, info)
  } else if (identical(format, "sp")) {
    points <- sp::SpatialPointsDataFrame(
      points,
      data = info,
      proj4string = fm_CRS(crs)
    )
  } else if (identical(format, "sf")) {
    points <- sf::st_as_sf(
      cbind(points, info),
      coords = seq_len(ncol(points)),
      crs = crs
    )
  }

  points # return
}


#' @title Extract vertex locations from an `inla.mesh`
#'
#' @description Extracts the vertices of an `inla.mesh` object.
#'
#' @export
#' @param x An `inla.mesh` object.
#' @param format character; `"sf"`, `"df"`, `"sp"`
#' @return
#' An `sf`, `data.frame`, or `SpatialPointsDataFrame` object, with the vertex
#' coordinates, and a `.vertex` column with the vertex indices.
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com},
#' Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [fm_centroids()]
#'
#' @examples
#' \donttest{
#' if (require(ggplot2, quietly = TRUE)) {
#'   data("mrsea", package = "inlabru")
#'   vrt <- fm_vertices(mrsea$mesh, format = "sp")
#'   ggplot() +
#'     gg(mrsea$mesh) +
#'     gg(vrt, color = "red")
#' }
#' }
#'
fm_vertices <- function(x, format = NULL) {
  fm_store_points(
    loc = x$loc,
    info = data.frame(.vertex = seq_len(nrow(x$loc))),
    crs = fm_crs(x),
    format = format
  )
}

#' @title Extract triangle centroids from an `inla.mesh`
#'
#' @description Computes the centroids of the triangles of an `inla.mesh`
#' object.
#'
#' @export
#' @param x An `inla.mesh` object.
#' @param format character; `"sf"`, `"df"`, `"sp"`
#' @return
#' An `sf`, `data.frame`, or `SpatialPointsDataFrame` object, with the vertex
#' coordinates, and a `.triangle` column with the triangle indices.
#'
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [fm_vertices()]
#'
#' @examples
#' \donttest{
#' if (require(ggplot2, quietly = TRUE)) {
#'   data("mrsea", package = "inlabru")
#'   vrt <- fm_centroids(mrsea$mesh, format = "sp")
#'   ggplot() +
#'     gg(mrsea$mesh) +
#'     gg(vrt, color = "red")
#' }
#' }
#'
fm_centroids <- function(x, format = NULL) {
  ## Extract triangle centroids
  loc <- (x$loc[x$graph$tv[, 1], , drop = FALSE] +
    x$loc[x$graph$tv[, 2], , drop = FALSE] +
    x$loc[x$graph$tv[, 3], , drop = FALSE]) / 3

  if (identical(x$manifold, "S2")) {
    loc <- loc / rowSums(loc^2)^0.5 * sum(x$loc[1, ]^2)^0.5
  }

  fm_store_points(
    loc = loc,
    info = data.frame(.triangle = seq_len(nrow(loc))),
    crs = fm_crs(x),
    format = format
  )
}



# Convert loc information to raw matrix coordinates for the mesh
fm_onto_mesh <- function(mesh, loc, crs = NULL) {
  if (!is.matrix(loc) && !fm_crs_is_null(crs)) {
    warning("loc is non-matrix but crs specified; will be ignored")
  }
  if (inherits(loc, c(
    "SpatialPoints", "SpatialPointsDataFrame",
    "sf", "sfc", "sfg"
  ))) {
    crs <- fm_crs(loc)
  }
  mesh_crs <- fm_crs(mesh)

  loc_needs_normalisation <- FALSE
  if (!fm_crs_is_null(crs) && !fm_crs_is_null(mesh_crs)) {
    if (fm_crs_is_geocent(mesh_crs)) {
      if (!is.matrix(loc)) {
        if (!fm_identical_CRS(crs, mesh_crs)) {
          loc <- fm_transform(loc, crs = mesh_crs, crs0 = crs)
        }
      }
    } else if (!fm_identical_CRS(crs, mesh_crs)) {
      loc <- fm_transform(loc, crs = mesh_crs, crs0 = crs, passthrough = TRUE)
    }
  } else if (identical(mesh$manifold, "S2")) {
    loc_needs_normalisation <- TRUE
  }

  if (inherits(loc, c("SpatialPoints", "SpatialPointsDataFrame"))) {
    loc <- sp::coordinates(loc)
  } else if (inherits(loc, c("sf", "sfc", "sfg"))) {
    loc <- sf::st_coordinates(loc)
    c_names <- colnames(loc)
    c_names <- intersect(c_names, c("X", "Y", "Z"))
    loc <- loc[, c_names, drop = FALSE]
  } else if (!is.matrix(loc)) {
    warning(
      paste0(
        "Unclear if the 'loc' class ('",
        paste0(class(loc), collapse = "', '"),
        "') is of a type we know how to handle."
      ),
      immediate. = TRUE
    )
  }
  if (loc_needs_normalisation) {
    loc <- loc / rowSums(loc^2)^0.5 * mean(rowSums(mesh$loc^2)^0.5)
  }

  loc
}






#' @title Compute an extension of a spatial object
#'
#' @description
#' Constructs a potentially nonconvex extension of a spatial object by
#' performing dilation by `convex + concave` followed by
#' erosion by `concave`. This is equivalent to dilation by `convex` followed
#' by closing (dilation + erosion) by `concave`.
#'
#' @param x A spatial object
#' @param ... Arguments passed on to the sub-methods
#' @param convex How much to extend
#' @param concave The minimum allowed reentrant curvature. Default equal to `convex`
#' @param preserveTopology logical; argument to `sf::st_simplify()`
#' @param dTolerance If not null, the `dTolerance` argument to `sf::st_simplify()`
#' @returns An extended object
#' @references Gonzalez and Woods (1992), Digital Image Processing
#' @export
#' @examples
#' inp <- sf::st_as_sf(as.data.frame(matrix(1:6, 3, 2)), coords = 1:2)
#' out <- fm_nonconvex_hull(inp, convex = 1)
#' plot(out)
fm_nonconvex_hull <- function(x, ...) {
  UseMethod("fm_nonconvex_hull")
}

#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.matrix <- function(x, ...) {
  fm_nonconvex_hull.sfc(sf::st_multipoint(x), ...)
}

#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.sf <- function(x, ...) {
  fm_nonconvex_hull.sfc(sf::st_geometry(x), ...)
}

#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.Spatial <- function(x, ...) {
  fm_nonconvex_hull.sfc(sf::st_as_sfc(x), ...)
}

#' @rdname fm_nonconvex_hull
#' @export
fm_nonconvex_hull.sfg <- function(x, ...) {
  fm_nonconvex_hull.sfc(sf::st_sfc(x), ...)
}

#' @rdname fm_nonconvex_hull
#' @details Differs from `sf::st_buffer(x, convex)` followed by
#' `sf::st_concave_hull()` (available from GEOS 3.11)
#' in how the amount of allowed concavity is controlled.
#' @export
fm_nonconvex_hull.sfc <- function(x,
                                  convex = -0.15,
                                  concave = convex,
                                  preserveTopology = TRUE,
                                  dTolerance = NULL,
                                  ...) {
  if ((convex < 0) || (concave < 0)) {
    approx_diameter <- fm_diameter(x)
    if (convex < 0) {
      convex <- approx_diameter * abs(convex)
    }
    if (concave < 0) {
      concave <- approx_diameter * abs(concave)
    }
  }

  nQuadSegs <- 64
  y <- sf::st_buffer(x, dist = convex + concave, nQuadSegs = nQuadSegs)
  y <- sf::st_union(y)
  y <- sf::st_buffer(y, dist = -concave, nQuadSegs = nQuadSegs)
  if (!is.null(dTolerance)) {
    y <- sf::st_simplify(y,
      preserveTopology = preserveTopology,
      dTolerance = dTolerance
    )
  }
  y <- sf::st_union(y)
  y
}

#' @title Compute approximate spatial object diameter
#' @param x A spatial object
#' @param ... Currently unused
#' @export
fm_diameter <- function(x, ...) {
  UseMethod("fm_diameter")
}

#' @rdname fm_diameter
#' @export
fm_diameter.sf <- function(x, ...) {
  fm_diameter.sfc(sf::st_geometry(x))
}

#' @rdname fm_diameter
#' @export
fm_diameter.sfg <- function(x, ...) {
  fm_diameter.sfc(sf::st_sfc(x))
}

#' @rdname fm_diameter
#' @export
fm_diameter.sfc <- function(x, ...) {
  fm_diameter.matrix(sf::st_coordinates(x))
}

#' @rdname fm_diameter
#' @export
fm_diameter.matrix <- function(x, ...) {
  if (ncol(x) == 1) {
    lim <- range(x)
    approx_diameter <- diff(lim)
  } else if (ncol(x) == 2) {
    lim <- rbind(range(x[, 1]), range(x[, 2]))
    approx_diameter <- max(diff(lim[1, ]), diff(lim[2, ]))
  } else if (ncol(x) >= 3) {
    lim <- rbind(range(x[, 1]), range(x[, 2]), range(x[, 3]))
    approx_diameter <- max(diff(lim[1, ]), diff(lim[2, ]), diff(lim[3, ]))
  }
  approx_diameter
}


#' @title Compute extensions of a spatial object
#'
#' @description
#' Constructs a potentially nonconvex extension of a spatial object by
#' performing dilation by `convex + concave` followed by
#' erosion by `concave`. This is equivalent to dilation by `convex` followed
#' by closing (dilation + erosion) by `concave`.
#' @seealso [fm_nonconvex_hull()]
#'
#' @param x A spatial object
#' @param convex numeric vector; How much to extend
#' @param concave numeric vector; The minimum allowed reentrant curvature. Default equal to `convex`
#' @param dTolerance If not null, the `dTolerance` argument to `sf::st_simplify()`,
#' passed on to [fm_nonconvex_hull()].
#' The default is `pmin(convex, concave) / 40`, chosen to
#' give approximately 4 or more subsegments per circular quadrant.
#' @param ... Optional further arguments to pass on to [fm_nonconvex_hull()].
#' @returns A list of `sfg` objects.
#' @export
#' @examples
#' if (bru_safe_inla()) {
#'   inp <- sf::st_as_sf(as.data.frame(matrix(1:6, 3, 2)), coords = 1:2)
#'   out <- fm_extensions(inp, convex = c(0.75, 2))
#'   bnd <- lapply(out, fm_as_segm)
#'   plot(INLA::inla.mesh.2d(boundary = bnd, max.edge = c(0.25, 1)), asp = 1)
#' }
#'
fm_extensions <- function(x,
                          convex = -0.15,
                          concave = convex,
                          dTolerance = NULL,
                          ...) {
  if (any(convex < 0) || any(concave < 0) || any(dTolerance < 0)) {
    approx_diameter <- fm_diameter(x)
  }
  len <- max(length(convex), length(concave), length(dTolerance))
  scale_fun <- function(val) {
    if (any(val < 0)) {
      val[val < 0] <- approx_diameter * abs(val[val < 0])
    }
    if (length(val) < len) {
      val <- c(val, rep(val[length(val)], len - length(val)))
    }
    val
  }
  convex <- scale_fun(convex)
  concave <- scale_fun(concave)
  if (is.null(dTolerance)) {
    dTolerance <- pmin(convex, concave) / 40
  } else {
    dTolerance <- scale_fun(dTolerance)
  }

  y <- lapply(
    seq_along(convex),
    function(k) {
      fm_nonconvex_hull(
        x,
        convex = convex[k],
        concave = concave[k],
        dTolerance = dTolerance[k],
        ...
      )
    }
  )
  y
}



# fm_segm ####

#' @title Make a spatial segment object
#' @describeIn fm_segm Create a new `fm_segm` object.
#' @export
#' @param ... Currently passed on to `inla.mesh.segment`
#' @family object creation and conversion
fm_segm <- function(...) {
  UseMethod("fm_segm")
}

#' @rdname fm_segm
#' @export
fm_segm.default <- function(...) {
  fm_as_segm(INLA::inla.mesh.segment(...))
}

#' @describeIn fm_segm Join multiple `fm_segm` objects into a single `fm_segm`
#' object.
#' @param grp.default When joining segments, use this group label for segments
#' that have `grp == NULL`.
#' @export
fm_segm.fm_segm <- function(..., grp.default = 0) {
  fm_as_segm(INLA::inla.mesh.segment(..., grp.default = grp.default))
}
#' @describeIn fm_segm Join multiple `inla.mesh.segment` objects into a single `fm_segm`
#' object.
#' @param grp.default When joining segments, use this group label for segments
#' that have `grp == NULL`.
#' @export
#' @method fm_segm inla.mesh.segment
fm_segm.inla.mesh.segment <- function(..., grp.default = 0) {
  fm_as_segm(INLA::inla.mesh.segment(..., grp.default = grp.default))
}


#' @title Convert objects to `fm_segm`
#' @describeIn fm_as_segm Convert an object to `fm_segm`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @export
#' @family object creation and conversion
fm_as_segm <- function(x, ...) {
  UseMethod("fm_as_segm")
}
#' @rdname fm_as_segm
#' @export
#' @method fm_as_segm inla.mesh.segment
fm_as_segm.inla.mesh.segment <- function(x, ...) {
  class(x) <- c("fm_segm", class(x))
  x
}




# fm_mesh ####

#' @title Convert objects to fmesher objects
#' @description
#' Used for conversion from general objects
#' (usually `inla.mesh` and other INLA specific classes)
#' to `fmesher` classes.
#'
#' @param x Object to be converted
#' @param ... Arguments forwarded to submethods
#' @rdname fm_as_fm
#' @export
#' @family object creation and conversion
fm_as_fm <- function(x, ...) {
  UseMethod("fm_as_fm")
}
#' @rdname fm_as_fm
#' @export
fm_as_fm.fm_mesh_1d <- function(x, ...) {
  #  class(x) <- c("fm_mesh_1d", setdiff(class(x), "fm_mesh_1d"))
  x
}
#' @rdname fm_as_fm
#' @export
fm_as_fm.fm_mesh_2d <- function(x, ...) {
  #  class(x) <- c("fm_mesh_2d", setdiff(class(x), "fm_mesh_2d"))
  x
}
#' @rdname fm_as_fm
#' @export
fm_as_fm.fm_segm <- function(x, ...) {
  #  class(x) <- c("fm_segm", setdiff(class(x), "fm_segm"))
  x
}
#' @rdname fm_as_fm
#' @export
fm_as_fm.fm_lattice_2d <- function(x, ...) {
  #  class(x) <- c("fm_lattice_2d", setdiff(class(x), "fm_lattice_2d"))
  x
}
#' @rdname fm_as_fm
#' @export
#' @method fm_as_fm inla.mesh.1d
fm_as_fm.inla.mesh.1d <- function(x, ...) {
  fm_as_mesh_1d(x, ...)
}
#' @rdname fm_as_fm
#' @export
#' @method fm_as_fm inla.mesh
fm_as_fm.inla.mesh <- function(x, ...) {
  fm_as_mesh_2d(x, ...)
}

#' @rdname fm_as_fm
#' @export
#' @method fm_as_fm inla.mesh.segment
fm_as_fm.inla.mesh.segment <- function(x, ...) {
  fm_as_segm(x, ...)
}

#' @rdname fm_as_fm
#' @export
#' @method fm_as_fm inla.mesh.lattice
fm_as_fm.inla.mesh.lattice <- function(x, ...) {
  fm_as_lattice_2d(x, ...)
}

# fm_mesh_1d ####

#' @title Make a 1D mesh object
#' @export
#' @param ... Currently passed on to `inla.mesh.1d`
#' @family object creation and conversion
fm_mesh_1d <- function(...) {
  UseMethod("fm_mesh_1d")
}

#' @rdname fm_mesh_1d
#' @export
fm_mesh_1d.default <- function(...) {
  fm_as_mesh_1d(INLA::inla.mesh.1d(...))
}

#' @title Convert objects to `fm_segm`
#' @describeIn fm_as_mesh_1d Convert an object to `fm_mesh_1d`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @export
#' @family object creation and conversion
#' @export
fm_as_mesh_1d <- function(...) {
  UseMethod("fm_as_mesh_1d")
}
#' @rdname fm_as_mesh_1d
#' @param x Object to be converted
#' @export
fm_as_mesh_1d.fm_mesh_1d <- function(x, ...) {
  #  class(x) <- c("fm_mesh_1d", setdiff(class(x), "fm_mesh_1d"))
  x
}
#' @rdname fm_as_mesh_1d
#' @param x Object to be converted
#' @export
#' @method fm_as_mesh_1d inla.mesh.1d
fm_as_mesh_1d.inla.mesh.1d <- function(x, ...) {
  class(x) <- c("fm_mesh_1d", class(x))
  x
}


# fm_mesh_2d ####

#' @title Make a 2D mesh object
#' @export
#' @param ... Currently passed on to `inla.mesh.2d`
#' @family object creation and conversion
fm_mesh_2d <- function(...) {
  UseMethod("fm_mesh_2d")
}

#' @rdname fm_mesh_2d
#' @export
fm_mesh_2d.default <- function(...) {
  fm_as_mesh_2d(INLA::inla.mesh.2d(...))
}

#' @title Convert objects to `fm_mesh_2d`
#' @describeIn fm_as_mesh_2d Convert an object to `fm_mesh_2d`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @export
#' @family object creation and conversion
#' @export
fm_as_mesh_2d <- function(...) {
  UseMethod("fm_as_mesh_2d")
}
#' @rdname fm_as_mesh_2d
#' @param x Object to be converted
#' @export
fm_as_mesh_2d.fm_mesh_2d <- function(x, ...) {
#  class(x) <- c("fm_mesh_2d", setdiff(class(x), "fm_mesh_2d"))
  x
}
#' @rdname fm_as_mesh_2d
#' @export
#' @method fm_as_mesh_2d inla.mesh
fm_as_mesh_2d.inla.mesh <- function(x, ...) {
  class(x) <- c("fm_mesh_2d", class(x))
  x
}




# fm_lattice_2d ####

#' @title Make a lattice object
#' @export
#' @param ... Currently passed on to `inla.mesh.lattice`
#' @family object creation and conversion
fm_lattice_2d <- function(...) {
  UseMethod("fm_lattice_2d")
}

#' @rdname fm_lattice_2d
#' @export
fm_lattice_2d.default <- function(...) {
  fm_as_lattice_2d(INLA::inla.mesh.lattice(...))
}

#' @title Convert objects to `fm_lattice_2d`
#' @describeIn fm_as_lattice_2d Convert an object to `fm_lattice_2d`.
#' @param x Object to be converted.
#' @param ... Arguments passed on to submethods
#' @export
#' @family object creation and conversion
#' @export
fm_as_lattice_2d <- function(...) {
  UseMethod("fm_as_lattice_2d")
}
#' @rdname fm_as_lattice_2d
#' @param x Object to be converted
#' @export
fm_as_lattice_2d.fm_lattice_2d <- function(x, ...) {
  #  class(x) <- c("fm_lattice_2d", setdiff(class(x), "fm_lattice_2d"))
  x
}
#' @rdname fm_as_lattice_2d
#' @param x Object to be converted
#' @export
#' @method fm_as_lattice_2d inla.mesh.lattice
fm_as_lattice_2d.inla.mesh.lattice <- function(x, ...) {
  class(x) <- c("fm_lattice_2d", class(x))
  x
}



# Deprecated ####

#' @describeIn inlabru-deprecated Conversion to inla.mesh.segment
#' `r lifecycle::badge("deprecated")` in favour of [fm_as_segm()].
#' @returns An `fm_segm` object
#' @export
fm_as_inla_mesh_segment <-
  function(...) {
    lifecycle::deprecate_soft("2.8.0.9008",
                              "fm_as_inla_mesh_segment()",
                              "fm_as_segm()")
    fm_as_segm(...)
  }

#' @describeIn inlabru-deprecated Conversion to inla.mesh.
#' `r lifecycle::badge("deprecated")` in favour of [fm_as_mesh_2d()].
#' @returns An `fm_mesh_2d` object
#' @export
fm_as_inla_mesh <- function(...) {
  lifecycle::deprecate_soft("2.8.0.9008",
                            "fm_as_inla_mesh()",
                            "fm_as_mesh_2d()")
  fm_as_mesh_2d(...)
}
