#' Visualize a globe using RGL
#'
#' @description
#' Creates a textured sphere and lon/lat coordinate annotations.
#' This function requires the `rgl` and `sphereplot` packages.
#'
#' @export
#' @param R Radius of the globe
#' @param R.grid Radius of the annotation sphere.
#' @param specular Light color of specular effect.
#' @param axes If TRUE, plot x, y and z axes.
#' @param box If TRUE, plot a box around the globe.
#' @param xlab,ylab,zlab Axes labels
#'
#' @return No value, used for plotting side effect.
#'
#' @family inlabru RGL tools
#'
#' @example inst/examples/rgl.R

globe <- function(R = 1,
                  R.grid = 1.05,
                  specular = "black",
                  axes = FALSE,
                  box = FALSE,
                  xlab = "", ylab = "", zlab = "") {
  # coordinates for texture
  n.smp <- 50
  lat <- matrix(-asin(seq(-1, 1, length.out = n.smp)),
    n.smp,
    n.smp,
    byrow = TRUE
  )
  long <- matrix(seq(-180, 180, length.out = n.smp) * pi / 180, n.smp, n.smp)
  x <- R * cos(lat) * cos(long)
  y <- R * cos(lat) * sin(long)
  z <- R * sin(lat)

  # globe and texture
  requireNamespace("rgl")
  rgl::persp3d(x, y, z,
    col = "white",
    texture = system.file("misc/Lambert_ocean.png", package = "inlabru"),
    specular = "black",
    axes = axes,
    box = box,
    xlab = xlab,
    ylab = ylab,
    zlab = zlab,
    normal_x = x,
    normal_y = y,
    normal_z = z
  )

  # spheric grid
  requireNamespace("sphereplot")
  sphereplot::rgl.sphgrid(longtype = "D", add = TRUE, radius = R.grid)
}

#' Render objects using RGL
#'
#' `glplot()` is a generic function for renders various kinds of spatial
#' objects, i.e. `Spatial*` data and `fm_mesh_2d` objects. The function invokes
#' particular methods which depend on the class of the first argument.
#'
#' @name glplot
#' @export
#' @param object an object used to select a method.
#' @param ... further arguments passed to or from other methods.
#'
#' @family inlabru RGL tools
#'
#' @example inst/examples/rgl.R

glplot <- function(object, ...) {
  UseMethod("glplot")
}

#' @describeIn glplot This function will calculate the cartesian coordinates of
#'   the points provided and use points3d() in order to render them.
#'
#' @export
#'
#' @param add If TRUE, add the points to an existing plot. If FALSE, create new
#'   plot.
#' @param color vector of R color characters. See material3d() for details.
#'
#' @family inlabru RGL tools


glplot.SpatialPoints <- function(object, add = TRUE, color = "red", ...) {
  if (length(sp::coordnames(object)) < 3) {
    ll <- data.frame(object)
    ll$TMP.ZCOORD <- 0
    sp::coordinates(ll) <- c(sp::coordnames(object), "TMP.ZCOORD")
    sp::proj4string(ll) <- fm_CRS(object)
    object <- ll
  }

  object <- fm_transform(object, crs = fm_crs("sphere"))
  cc <- sp::coordinates(object)
  requireNamespace("rgl")
  rgl::points3d(
    x = cc[, 1],
    y = cc[, 2],
    z = cc[, 3],
    add = add,
    color = color,
    ...
  )
}

#' @describeIn glplot This function will calculate a cartesian representation of
#'   the lines provided and use `lines3d()` in order to render them.
#'
#' @export
#'
#' @family inlabru RGL tools

glplot.SpatialLines <- function(object, add = TRUE, ...) {
  qq <- sp::coordinates(object)
  sp <- do.call(rbind, lapply(qq, function(k) {
    do.call(rbind, lapply(k, function(x) {
      x[1:(nrow(x) - 1), ]
    }))
  }))
  ep <- do.call(rbind, lapply(qq, function(k) {
    do.call(rbind, lapply(k, function(x) {
      x[2:(nrow(x)), ]
    }))
  }))
  sp <- data.frame(x = sp[, 1], y = sp[, 2], z = 0)
  ep <- data.frame(x = ep[, 1], y = ep[, 2], z = 0)

  sp::coordinates(sp) <- c("x", "y", "z")
  sp::coordinates(ep) <- c("x", "y", "z")
  sp::proj4string(sp) <- fm_CRS(object)
  sp::proj4string(ep) <- fm_CRS(object)

  sp <- fm_transform(sp, crs = fm_crs("sphere"))
  ep <- fm_transform(ep, crs = fm_crs("sphere"))

  cs <- sp::coordinates(sp)
  ce <- sp::coordinates(ep)
  na <- matrix(NA, ncol = 3, nrow = nrow(cs))

  mm <- matrix(t(cbind(cs, ce, na)),
    ncol = 3,
    nrow = 3 * nrow(ce),
    byrow = TRUE
  )

  requireNamespace("rgl")

  rgl::lines3d(mm, add = add, ...)
}


#' @describeIn glplot This function transforms the mesh to 3D cartesian
#'   coordinates and uses `inla.plot.mesh()` with `rgl=TRUE` to plot the result.
#'
#' @export
#'
#' @param col Color specification. A single named color, a vector of scalar
#'   values, or a matrix of RGB values.
#' @param ... Parameters passed on to plot_rgl.fm_mesh_2d()
#'
#' @family inlabru RGL tools

glplot.fm_mesh_2d <- function(object, add = TRUE, col = NULL, ...) {
  if (object$manifold == "S2") {
    # mesh$loc = mesh$loc
  } else {
    object <- fm_transform(object, crs = fm_crs("sphere"), passthrough = TRUE)
  }

  if (is.null(col)) {
    fmesher::plot_rgl(object, add = add, ...)
  } else {
    fmesher::plot_rgl(object, add = add, col = col, ...)
  }
}
#' @rdname glplot
#' @export
glplot.inla.mesh <- function(object, add = TRUE, col = NULL, ...) {
  glplot(fm_as_fm(object), add = add, col = col, ...)
}
