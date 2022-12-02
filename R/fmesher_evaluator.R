# Point/mesh connection methods ####

#' @title Methods for projecting to/from an inla.mesh
#'
#' @description Calculate a lattice projection to/from an `inla.mesh`
#'
#' @param mesh An `inla.mesh` or `inla.mesh.1d` object.
#' @param loc Projection locations.  Can be a matrix, `SpatialPoints`,
#' `SpatialPointsDataFrame`, `sf`, `sfc`, or `sfg` object.
#' @param lattice An `inla.mesh.lattice()` object.
#' @param xlim X-axis limits for a lattice. For R2 meshes, defaults to covering
#' the domain.
#' @param ylim Y-axis limits for a lattice. For R2 meshes, defaults to covering
#' the domain.
#' @param dims Lattice dimensions.
#' @param projector An `fm_evaluator` object.
#' @param field Basis function weights, one per mesh basis function, describing
#' the function to be evaluated at the projection locations
#' @param projection One of `c("default", "longlat", "longsinlat",
#' "mollweide")`.
#' @param crs An optional CRS or inla.CRS object associated with `loc`
#' and/or `lattice`.
#' @param \dots Additional arguments passed on to methods.
#' @return For `fm_evaluate(mesh, ...)`, a list with projection
#' information.  For `fm_evaluator(mesh, ...)`, an
#' `fm_evaluator` object.  For `fm_evaluate(projector,
#' field, ...)`, a field projected from the mesh onto the locations given by
#' the projector object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso `inla.mesh`, `inla.mesh.1d`,
#' `inla.mesh.lattice`
#' @examples
#' if (bru_safe_inla()) {
#'   n <- 20
#'   loc <- matrix(runif(n * 2), n, 2)
#'   mesh <- INLA::inla.mesh.create(loc, refine = list(max.edge = 0.05))
#'   proj <- fm_evaluator(mesh)
#'   field <- cos(mesh$loc[, 1] * 2 * pi * 3) * sin(mesh$loc[, 2] * 2 * pi * 7)
#'   image(proj$x, proj$y, fm_evaluate(proj, field))
#' }
#' \donttest{
#' if (bru_safe_inla() &&
#'   require(rgl)) {
#'   plot(mesh, rgl = TRUE, col = field, draw.edges = FALSE, draw.vertices = FALSE)
#' }
#' }
#' @export fm_evaluate
fm_evaluate <- function(...) {
  UseMethod("fm_evaluate")
}

#' @export
#' @rdname fm_evaluate
fm_evaluate.inla.mesh <- function(mesh, field, ...) {
  stopifnot(!missing(field) && !is.null(field))

  proj <- fm_evaluator(mesh, ...)
  fm_evaluate(proj, field = field)
}


#' @export
#' @rdname fm_evaluate
fm_evaluate.inla.mesh.1d <- function(mesh, field, ...) {
  stopifnot(!missing(field) && !is.null(field))

  proj <- fm_evaluator(mesh, ...)
  fm_evaluate(proj, field)
}


#' @export
#' @rdname fm_evaluate
fm_evaluate.fm_evaluator <-
  function(projector, field, ...) {
    if (is.data.frame(field)) {
      field <- as.matrix(field)
    }

    if (is.null(dim(field))) {
      if (is.null(projector$lattice)) {
        data <- as.vector(projector$proj$A %*% as.vector(field))
        data[!projector$proj$ok] <- NA
        return(data)
      } else {
        data <- as.vector(projector$proj$A %*% as.vector(field))
        data[!projector$proj$ok] <- NA
        return(matrix(
          data,
          projector$lattice$dims[1],
          projector$lattice$dims[2]
        ))
      }
    } else {
      data <- projector$proj$A %*% field
      data[!projector$proj$ok, ] <- NA
      return(data)
    }
  }


#' @export
#' @rdname fm_evaluate
fm_evaluator <- function(...) {
  UseMethod("fm_evaluator")
}


#' @export
#' @rdname fm_evaluate
fm_evaluator_inla_mesh <- function(mesh, loc = NULL, crs = NULL, ...) {
  stopifnot(inherits(mesh, "inla.mesh"))

  # Support INLA <= 22.11.27 by converting globes to spheres
  # TODO: Handle the > 22.11.27 more efficiently
  if (!fm_crs_is_null(mesh$crs)) {
    if (fm_crs_is_geocent(mesh$crs)) {
      crs.sphere <- fm_CRS("sphere")
      if (!fm_identical_CRS(mesh$crs, crs.sphere)) {
        ## Convert the mesh to a perfect sphere.
        mesh <- fm_transform(mesh, crs = crs.sphere)
      }
      if (!is.matrix(loc)) {
        if (!fm_identical_CRS(crs, crs.sphere)) {
          loc <- fm_transform(loc, crs = crs.sphere, crs0 = crs)
        }
      } else {
        if (fm_crs_is_null(crs)) {
          loc <- loc / rowSums(loc^2)^0.5
        } else {
          loc <- fm_transform(loc, crs = crs.sphere, crs0 = crs)
        }
      }
    } else if (!fm_identical_CRS(crs, mesh$crs)) {
      mesh <- fm_transform(mesh, crs = mesh$crs, crs0 = crs, passthrough = TRUE)
    }
  } else if (identical(mesh$manifold, "S2")) {
    mesh$loc <- mesh$loc / rowSums(mesh$loc^2)^0.5
    loc <- loc / rowSums(loc^2)^0.5
  }

  if (inherits(loc, c("SpatialPoints", "SpatialPointsDataFrame"))) {
    loc <- sp::coordinates(loc)
  } else if (inherits(loc, c("sf", "sfc", "sfg"))) {
    loc <- sf::st_coordinates(loc)
    c_names <- colnames(loc)
    c_names <- intersect(c_names, c("X", "Y", "Z"))
    loc <- loc[, c_names, drop = FALSE]
  }

  jj <-
    which(rowSums(matrix(is.na(as.vector(loc)),
      nrow = nrow(loc),
      ncol = ncol(loc)
    )) == 0)
  smorg <- (INLA::inla.fmesher.smorg(mesh$loc,
    mesh$graph$tv,
    points2mesh = loc[jj, , drop = FALSE]
  ))
  ti <- matrix(0L, nrow(loc), 1)
  b <- matrix(0, nrow(loc), 3)
  ti[jj, 1L] <- smorg$p2m.t
  b[jj, ] <- smorg$p2m.b

  ok <- (ti[, 1L] > 0L)
  ti[ti[, 1L] == 0L, 1L] <- NA

  ii <- which(ok)
  A <- (Matrix::sparseMatrix(
    dims = c(nrow(loc), mesh$n),
    i = rep(ii, 3),
    j = as.vector(mesh$graph$tv[ti[ii, 1L], ]),
    x = as.vector(b[ii, ])
  ))

  list(t = ti, bary = b, A = A, ok = ok)
}


#' @export
#' @rdname fm_evaluate
fm_evaluator_inla_mesh_1d <- function(mesh, loc, ...) {
  stopifnot(inherits(mesh, "inla.mesh.1d"))

  A <- INLA::inla.mesh.1d.A(mesh, loc = loc)

  return(list(
    A = A,
    ok = (loc >= mesh$interval[1]) & (loc <= mesh$interval[2])
  ))
}




#' @export
#' @rdname fm_evaluate
fm_evaluator_lattice <- function(mesh,
                                 xlim = NULL,
                                 ylim = NULL,
                                 dims = c(100, 100),
                                 projection = NULL,
                                 crs = NULL,
                                 ...) {
  stopifnot(inherits(mesh, "inla.mesh"))
  if (identical(mesh$manifold, "R2") &&
    (is.null(mesh$crs) || is.null(crs))) {
    units <- "default"
    lim <- list(
      xlim = if (is.null(xlim)) range(mesh$loc[, 1]) else xlim,
      ylim = if (is.null(ylim)) range(mesh$loc[, 2]) else ylim
    )
  } else if (identical(mesh$manifold, "S2") &&
    (is.null(mesh$crs) || is.null(crs))) {
    projection <-
      match.arg(projection, c(
        "longlat", "longsinlat",
        "mollweide"
      ))
    units <- projection
    lim <- INLA::inla.mesh.map.lim(loc = mesh$loc, projection = projection)
  } else {
    lim <- fm_crs_bounds(crs)
    if (identical(mesh$manifold, "R2")) {
      lim0 <- list(
        xlim = if (is.null(xlim)) range(mesh$loc[, 1]) else xlim,
        ylim = if (is.null(ylim)) range(mesh$loc[, 2]) else ylim
      )
      lim$xlim[1] <- max(lim$xlim[1], lim0$xlim[1])
      lim$xlim[2] <- min(lim$xlim[2], lim0$xlim[2])
      lim$ylim[1] <- max(lim$ylim[1], lim0$ylim[1])
      lim$ylim[2] <- min(lim$ylim[2], lim0$ylim[2])
    }
  }
  if (missing(xlim) && is.null(xlim)) {
    xlim <- lim$xlim
  }
  if (missing(ylim) && is.null(ylim)) {
    ylim <- lim$ylim
  }
  x <- seq(xlim[1], xlim[2], length.out = dims[1])
  y <- seq(ylim[1], ylim[2], length.out = dims[2])
  if (is.null(mesh$crs) || is.null(crs)) {
    lattice <- INLA::inla.mesh.lattice(x = x, y = y, units = units)
  } else {
    lattice <- INLA::inla.mesh.lattice(x = x, y = y, crs = crs)
  }
  lattice
}

#' @export
#' @describeIn fm_evaluate The `...` arguments are passed on to `fm_evaluator_lattice()`
#' if no `loc` or `lattice` is provided.
fm_evaluator.inla.mesh <- function(mesh,
                                   loc = NULL,
                                   lattice = NULL,
                                   crs = NULL,
                                   ...) {
  if (missing(loc) || is.null(loc)) {
    if (missing(lattice) || is.null(lattice)) {
      lattice <- fm_evaluator_lattice(mesh,
        crs = crs,
        ...
      )
    }
    dims <- lattice$dims
    x <- lattice$x
    y <- lattice$y
    crs <- lattice$crs

    if (is.null(mesh$crs) || is.null(lattice$crs)) {
      proj <- fm_evaluator_inla_mesh(mesh, lattice$loc)
    } else {
      proj <- fm_evaluator_inla_mesh(mesh,
        loc = lattice$loc,
        crs = lattice$crs
      )
    }
    projector <- list(x = x, y = y, lattice = lattice, loc = NULL, proj = proj, crs = crs)
    class(projector) <- "fm_evaluator"
  } else {
    proj <- fm_evaluator_inla_mesh(mesh, loc = loc, crs = crs)
    projector <- list(x = NULL, y = NULL, lattice = NULL, loc = loc, proj = proj, crs = crs)
    class(projector) <- "fm_evaluator"
  }

  return(projector)
}


#' @export
#' @rdname fm_evaluate
fm_evaluator.inla.mesh.1d <- function(mesh,
                                      loc = NULL,
                                      xlim = mesh$interval,
                                      dims = 100,
                                      ...) {
  stopifnot(inherits(mesh, "inla.mesh.1d"))

  if (missing(loc) || is.null(loc)) {
    loc <- seq(xlim[1], xlim[2], length.out = dims[1])
  }

  proj <- fm_evaluator_inla_mesh_1d(mesh, loc)
  projector <- list(x = loc, lattice = NULL, loc = loc, proj = proj)
  class(projector) <- "fm_evaluator"

  return(projector)
}
