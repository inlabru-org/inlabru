#' Sample from an inhomogeneous Poisson process
#'
#' This function provides point samples from one- and two-dimensional inhomogeneous Poisson processes. The
#' log intensity has to be provided via its values at the nodes of an \code{inla.mesh.1d} or
#' \code{inla.mesh} object. In between mesh nodes the log intensity is assumed to be linear.
#'
#' For 2D processes on a sphere the \code{R} parameter can be used to adjust to sphere's radius implied by
#' the mesh. If the intensity is very high the standard \code{strategy} "spherical" can cause memory issues.
#' Using the "sliced-spherical" strategy can help in this case.
#'
#' @aliases sample.lgcp
#' @export
#'
#' @param mesh An \link[INLA]{inla.mesh} object
#' @param loglambda vector or matrix; A vector of log intensities at the mesh vertices
#'   (for higher order basis functions, e.g.
#'   for \code{inla.mesh.1d} meshes, \code{loglambda} should be given as \code{mesh$m} basis
#'   function weights rather than the values at the \code{mesh$n} vertices)
#'   A single scalar is expanded to a vector of the appropriate length.
#'   If a matrix is supplied, one process sample for each column is produced.
#' @param strategy Only relevant for 2D meshes. One of \code{'triangulated'}, \code{'rectangle'},
#'   \code{'sliced-spherical'}, \code{'spherical'}. The \code{'rectangle'} method is only valid for
#'   CRS-less flat 2D meshes.
#'   If \code{NULL} or \code{'auto'}, the the likely fastest method is chosen;
#'   \code{'rectangle'} for flat 2D meshes with no CRS,
#'   \code{'sliced-spherical'} for CRS \code{'longlat'} meshes, and
#'   \code{'triangulated'} for all other meshes.
#' @param R Numerical value only applicable to spherical and geographical meshes. It is interpreted as
#'   \code{R} is the equivalent Earth radius, in km, used to scale the lambda intensity.
#'     For CRS enabled meshes, the default is 6371. For CRS-less spherical meshes, the default is 1.
#' @param samplers A `SpatialPolygonsDataFrame` or `inla.mesh` object.
#'   Simulated points that fall outside these polygons are discarded.
#' @param ignore.CRS logical; if \code{TRUE}, ignore any CRS information in the mesh. Default \code{FALSE}.
#'   This affects \code{R} and the permitted values for \code{strategy}.
#'
#' @return A \code{data.frame} (1D case),
#'   SpatialPoints (2D flat and 3D spherical surface cases)
#'   SpatialPointsDataFrame (2D/3D surface cases with multiple samples).
#'   For multiple samples, the \code{data.frame} output has a
#'   column \code{'sample'} giving the index for each sample.
#' object of point locations.
#'
#' @details
#' \itemize{
#' \item For crs-less meshes on R2: Lambda is interpreted in the raw coordinate system. Output has an NA CRS.
#' \item For crs-less meshes on S2: Lambda with raw units, after scaling the mesh to radius \code{R}, if specified.
#'   Output is given on the same domain as the mesh, with an NA CRS.
#' \item For crs meshes on R2: Lambda is interpreted as per km^2, after scaling the globe to the Earth radius 6371 km,
#'   or \code{R}, if specified. Output given in the same CRS as the mesh.
#' \item For crs meshes on S2: Lambda is interpreted as per km^2, after scaling the globe to the Earth radius 6371 km,
#'   or \code{R}, if specified. Output given in the same CRS as the mesh.
#' }
#'
#' @author Daniel Simpson <\email{dp.simpson@@gmail.com}> (base rectangle and spherical algorithms),
#' Fabian E. Bachl <\email{bachlfab@@gmail.com}> (inclusion in inlabru, sliced spherical sampling),
#' Finn Lindgren <\email{finn.lindgren@@gmail.com}> (extended CRS support, triangulated sampling)
#'
#' @examples
#' \donttest{
#' # The INLA package is required
#' if (require("INLA", quietly = TRUE)) {
#'
#' vertices = seq(0, 3, by = 0.1)
#' mesh = inla.mesh.1d(vertices)
#' loglambda = 5-0.5*vertices
#' pts = sample.lgcp(mesh, loglambda)
#' pts$y = 0
#' plot(vertices, exp(loglambda), type = "l", ylim = c(0,150))
#' points(pts, pch = "|" )
#'
#' }
#' }
#'
#' @examples
#' \donttest{
#' # The INLA package is required
#' if (require("INLA", quietly = TRUE)) {
#'
#' data("gorillas", package = "inlabru")
#' pts = sample.lgcp(gorillas$mesh,
#'                   loglambda = 1.5,
#'                   samplers = gorillas$boundary)
#' ggplot() + gg(gorillas$mesh) + gg(pts)
#'
#' }
#' }

sample.lgcp <- function(mesh, loglambda, strategy = NULL, R = NULL, samplers = NULL,
                        ignore.CRS = FALSE) {
  requireINLA()
  if (inherits(mesh, "inla.mesh.1d")) {
    xmin <- mesh$interval[1]
    xmax <- mesh$interval[2]
    area <- xmax - xmin

    multi.samples <- is.matrix(loglambda)
    if (multi.samples) {
      result <- list()
      n.samples <- ncol(loglambda)
      loglambda.matrix <- loglambda
    } else {
      n.samples <- 1
    }

    multi.sample <- list()
    for (sample in seq_len(n.samples)) {
      if (multi.samples) {
        loglambda <- as.vector(loglambda.matrix[, sample, drop = TRUE])
      }

      # The B-spline basis expansion has a convex hull property
      # which implies that the maximum weight is not greater or equal to the maximum function value.
      wmax <- max(loglambda)

      Npoints <- rpois(1, lambda = area * exp(wmax))
      if (Npoints > 0) {
        points <- runif(n = Npoints, min = xmin, max = xmax)
        proj <- INLA::inla.mesh.project(mesh, points)
        if (length(loglambda) == 1) {
          lambda_ratio <- exp(as.vector(rowSums(proj$A) * loglambda) - wmax)
        } else {
          lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - wmax)
        }
        keep <- (runif(Npoints) <= lambda_ratio)
        waste_ratio <- sum(keep) / length(keep)
        if (multi.samples) {
          multi.sample[[sample]] <- data.frame(x = points[keep], sample = sample)
        } else {
          multi.sample[[sample]] <- data.frame(x = points[keep])
        }
      } else {
        waste_ratio <- 0
        if (multi.samples) {
          multi.sample[[sample]] <- data.frame(x = numeric(0), sample = integer(0))
        } else {
          multi.sample[[sample]] <- data.frame(x = numeric(0))
        }
      }
    }
    ret <- do.call(rbind, multi.sample)
  } else if (inherits(mesh, "inla.mesh")) {
    multi.samples <- is.matrix(loglambda)
    if (multi.samples) {
      result <- list()
      n.samples <- ncol(loglambda)
      if (nrow(loglambda) == 1) {
        loglambda.matrix <- matrix(loglambda, mesh$n, n.samples, byrow = TRUE)
      } else {
        loglambda.matrix <- loglambda
      }
    } else {
      n.samples <- 1
      loglambda <- as.vector(loglambda)
      if (length(loglambda) == 1) {
        loglambda <- rep(loglambda, mesh$n)
      }
    }

    if (ignore.CRS) {
      mesh$crs <- NULL
    }
    input.crs <- INLA::inla.CRS(INLA::inla.CRSargs(mesh$crs))
    input.crs.list <- INLA::inla.as.list.CRSargs(INLA::inla.CRSargs(input.crs))
    use.crs <- !is.null(input.crs.list$proj) && !ignore.CRS
    is.geocent <- (mesh$manifold == "S2")

    if (is.geocent || use.crs) {
      strategies <- c("triangulated", "sliced-spherical", "spherical")
    } else {
      strategies <- c("triangulated", "rectangle")
    }
    if (!is.null(strategy) && (strategy == "auto")) {
      strategy <- NULL
    }
    if (is.null(strategy)) {
      if (is.geocent) {
        strategy <- "triangulated"
      } else if (!use.crs) {
        strategy <- "rectangle"
      } else if (identical(input.crs.list$proj, "longlat")) {
        strategy <- "sliced-spherical"
      }
    }
    strategy <- match.arg(strategy, strategies)

    if (is.geocent) {
      space.R <- mean(rowSums(mesh$loc ^ 2) ^ 0.5)
      if (is.null(input.crs.list$units)) {
        space.units <- "m"
      } else {
        space.units <- input.crs.list$units
      }
      internal.crs <- INLA::inla.CRS("sphere", args = list(a = 1, b = 1, units = "m"))
      mesh$loc <- mesh$loc / space.R
      mesh$crs <- internal.crs
    } else {
      if (use.crs) {
        internal.crs <- INLA::inla.CRS("sphere", args = list(a = 1, b = 1, units = "m"))
      } else {
        internal.crs <- CRS(as.character(NA))
        mesh$crs <- NULL
      }
    }

    area.R <- 1
    if (!missing(R) && !is.null(R)) {
      area.R <- R
    } else {
      if (use.crs) {
        area.R <- 6371
      } else if (is.geocent) {
        area.R <- space.R
      }
    }

    for (sample in seq_len(n.samples)) {
      if (multi.samples) {
        loglambda <- as.vector(loglambda.matrix[, sample, drop = TRUE])
      }

      if (strategy == "triangulated") {
        if (is.geocent) {
          area.mesh <- mesh
        } else if (use.crs) {
          area.mesh <- INLA::inla.spTransform(mesh, CRSobj = internal.crs)
        } else {
          area.mesh <- mesh
          area.R <- 1
        }
        areas <- INLA::inla.fmesher.smorg(
          area.mesh$loc,
          area.mesh$graph$tv,
          fem = 0,
          output = "ta"
        )$ta * area.R ^ 2

        loglambda_tri <- matrix(loglambda[mesh$graph$tv], nrow(mesh$graph$tv), ncol(mesh$graph$tv))
        loglambda_max <- apply(loglambda_tri, 1, max)

        Npoints <- rpois(length(areas), lambda = exp(loglambda_max) * areas)

        if (sum(Npoints) > 0) {
          triangle <- rep(seq_along(areas), Npoints)

          # Sample uniformly on base triangle (0,0),(1,0),(0,1)
          points <- matrix(runif(2 * sum(Npoints)), sum(Npoints), 2)
          flip <- rowSums(points) > 1.0
          points[flip, ] <- 1 - points[flip, 2:1, drop = FALSE]
          # Convert to in-triangle points
          points <-
            mesh$loc[mesh$graph$tv[triangle, 1], , drop = FALSE] +
            (mesh$loc[mesh$graph$tv[triangle, 2], , drop = FALSE] -
              mesh$loc[mesh$graph$tv[triangle, 1], , drop = FALSE]) * points[, 1] +
            (mesh$loc[mesh$graph$tv[triangle, 3], , drop = FALSE] -
              mesh$loc[mesh$graph$tv[triangle, 1], , drop = FALSE]) * points[, 2]
          if (is.geocent) {
            points <- points / rowSums(points ^ 2) ^ 0.5
          }
        }

        # Construct SpatialPoints
        if (is.geocent) {
          target.crs <- internal.crs
        } else if (use.crs) {
          target.crs <- input.crs
        } else {
          target.crs <- CRS(as.character(NA))
        }
        if (sum(Npoints) > 0) {
          points <- sp::SpatialPoints(points, proj4string = target.crs)

          A <- INLA::inla.mesh.project(mesh, points)$A
          lambda_ratio <- exp(as.vector(A %*% loglambda) - loglambda_max[triangle])
          keep <- (runif(sum(Npoints)) <= lambda_ratio)
          ret <- points[keep]
          waste_ratio <- sum(keep) / length(keep)
        } else {
          points <- sp::SpatialPoints(matrix(0, 1, 3), proj4string = target.crs)[-1]
          ret <- points
          waste_ratio <- 0
        }
      } else if (strategy == "rectangle") {
        # Construct bounding rectangle
        loc <- mesh$loc
        xmin <- min(loc[, 1])
        xmax <- max(loc[, 1])
        ymin <- min(loc[, 2])
        ymax <- max(loc[, 2])
        area <- (xmax - xmin) * (ymax - ymin)

        # Simulate number of points
        lambda_max <- max(loglambda)
        Npoints <- rpois(1, lambda = area * exp(lambda_max))

        if (Npoints == 0) {
          ret <- SpatialPoints(matrix(0, 1, 2), proj4string = internal.crs)[-1]
          waste_ratio <- 0
        } else {
          # Simulate uniform points on the bounding rectangle
          points <- sp::SpatialPoints(
            cbind(
              x = runif(n = Npoints, min = xmin, max = xmax),
              y = runif(n = Npoints, min = ymin, max = ymax)
            ),
            proj4string = internal.crs
          )

          # Do some thinning
          proj <- INLA::inla.mesh.project(mesh, points)
          lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - lambda_max)
          keep <- proj$ok & (runif(Npoints) <= lambda_ratio)
          ret <- points[keep]
          waste_ratio <- sum(keep) / length(keep)
        }
      } else if (strategy == "spherical") {
        # Simulate number of points
        area <- 4 * pi * area.R ^ 2
        lambda_max <- max(loglambda)
        Npoints <- rpois(n = 1, lambda = area * exp(lambda_max))
        if (Npoints > 5000000) {
          stop(paste0("Too many points!: ", Npoints))
        }

        if (Npoints == 0) {
          ret <- SpatialPoints(matrix(0, 1, 3), proj4string = internal.crs)[-1]
          waste_ratio <- 0
        } else {
          # Choose z uniformly distributed in [-1,1].
          z <- runif(n = Npoints, min = -1, max = 1)
          # Choose t uniformly distributed on [0, 2*pi).
          t <- runif(n = Npoints, min = 0, max = 2 * pi)
          r <- sqrt(1 - z ^ 2)
          x <- r * cos(t)
          y <- r * sin(t)

          points <- data.frame(x, y, z)
          coordinates(points) <- c("x", "y", "z")
          proj4string(points) <- internal.crs

          proj <- INLA::inla.mesh.project(mesh, points)
          lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - lambda_max)
          keep <- proj$ok & (runif(Npoints) <= lambda_ratio)
          ret <- points[keep]
          waste_ratio <- sum(keep) / length(keep)
        }
      } else if (strategy == "sliced-spherical") {
        if (identical(input.crs.list$proj, "longlat")) {
          # Simulate number of points
          lon.range <- range(mesh$loc[, 1])
          lat.range <- range(mesh$loc[, 2])
          lon.rel <- diff(lon.range) / 360
          lat.rel <- diff(lat.range) / 180
        } else {
          lon.range <- c(-180, 180)
          lat.range <- c(-90, 90)
          lon.rel <- 1
          lat.rel <- 1
        }
        slat.range <- pmax(-1, pmin(1, sin(lat.range * pi / 180)))
        radlon.range <- lon.range * pi / 180
        area <- (4 * pi * area.R ^ 2) * lon.rel * lat.rel

        lambda_max <- max(loglambda)
        Npoints <- rpois(n = 1, lambda = area * exp(lambda_max))

        sampled.points <- list()
        nblocks <- max(1, Npoints / 1000000)
        sp <- ceiling(seq(1, Npoints + 1, length.out = max(2, nblocks)))
        n.keep <- 0
        n.sampled <- 0

        for (k in seq_len(length(sp) - 1)) {
          n.points <- sp[k + 1] - sp[k]
          if (n.points == 0) {
            sampled.points[[k]] <- SpatialPoints(matrix(0, 1, 3), proj4string = internal.crs)[-1]
            break
          }

          # Choose z uniformly distributed in [-1,1].
          z <- runif(n = n.points, min = slat.range[1], max = slat.range[2]) # transforms into latitude
          # Choose t uniformly distributed on [0, 2*pi).
          angle <- runif(n = n.points, min = radlon.range[1], max = radlon.range[2])
          r <- sqrt(1 - z ^ 2)
          x <- r * cos(angle)
          y <- r * sin(angle)

          points <- data.frame(x, y, z)
          coordinates(points) <- c("x", "y", "z")
          proj4string(points) <- internal.crs

          proj <- INLA::inla.mesh.project(mesh, points)
          lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - lambda_max)
          keep <- proj$ok & (runif(Npoints) <= lambda_ratio)
          sampled.points[[k]] <- points[keep]

          n.keep <- n.keep + sum(keep)
          n.sampled <- n.sampled + length(keep)
        }

        # What to return
        if (sum(vapply(sampled.points, length, 1L)) > 0) {
          ret <- do.call(rbind, sampled.points)
        } else {
          ret <- sampled.points[[1]]
        }
        waste_ratio <- n.keep / n.sampled
      }

      if (multi.samples) {
        result[[sample]] <-
          sp::SpatialPointsDataFrame(
            ret,
            data = data.frame(sample = rep(sample, length(ret)))
          )
      }
    }
    if (multi.samples) {
      ret <- do.call(rbind, result[vapply(result, length, 1L) > 0])
    }

    if (is.geocent) {
      if (length(ret) > 0) {
        ret <- sp::spTransform(ret, INLA::inla.CRS("sphere", args = list(a = space.R, b = space.R, units = "m")))
      } else if (multi.samples) {
        ret <- SpatialPointsDataFrame(matrix(0, 1, 3), data = data.frame(sample = 1))[-1]
      } else {
        ret <- SpatialPoints(matrix(0, 1, 3))[-1]
      }
      if (use.crs) {
        proj4string(ret) <- input.crs
      } else {
        proj4string(ret) <- CRS(as.character(NA))
      }
    } else {
      if (use.crs) {
        if (length(ret) > 0) {
          ret <- spTransform(ret, input.crs)
        } else if (multi.samples) {
          ret <- SpatialPointsDataFrame(matrix(0, 1, 2), data = data.frame(sample = 1))[-1]
          proj4string(ret) <- input.crs
        } else {
          ret <- SpatialPoints(matrix(0, 1, 2))[-1]
          proj4string(ret) <- input.crs
        }
      } else {
        proj4string(ret) <- CRS(as.character(NA))
      }
    }

    # Only retain points within the samplers
    if (!is.null(samplers) && (length(ret) > 0)) {
      if (inherits(samplers, "inla.mesh")) {
        proj <- INLA::inla.mesh.project(samplers, points)
        ret <- ret[proj$ok]
      } else {
        ret <- ret[!is.na(over(ret, samplers))]
      }
    }
  } else {
    stop(paste0(
      "The `mesh` must inherit from `inla.mesh.1d` or `inla.mesh`.\n  class = c('",
      paste0(class(mesh), collapse = "', '"),
      "')"
    ))
  }

  ret
}