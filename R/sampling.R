#' Sample from an inhomogeneous Poisson process
#'
#' This function provides point samples from one- and two-dimensional
#' inhomogeneous Poisson processes. The log intensity has to be provided via its
#' values at the nodes of an `fm_mesh_1d` or `fm_mesh_2d` object. In between
#' mesh nodes the log intensity is assumed to be linear.
#'
#' For 2D processes on a sphere the `R` parameter can be used to adjust to
#' sphere's radius implied by the mesh. If the intensity is very high the
#' standard `strategy` "spherical" can cause memory issues. Using the
#' "sliced-spherical" strategy can help in this case.
#'
#' @export
#'
#' @param mesh An [fmesher::fm_mesh_1d] or [fmesher::fm_mesh_2d] object
#' @param loglambda vector or matrix; A vector of log intensities at the mesh
#'   vertices (for higher order basis functions, e.g. for `fm_mesh_1d` meshes,
#'   `loglambda` should be given as `mesh$m` basis function weights rather than
#'   the values at the `mesh$n` vertices) A single scalar is expanded to a
#'   vector of the appropriate length. If a matrix is supplied, one process
#'   sample for each column is produced.
#' @param strategy Only relevant for 2D meshes. One of `'triangulated'`,
#'   `'rectangle'`, `'sliced-spherical'`, `'spherical'`. The `'rectangle'`
#'   method is only valid for CRS-less flat 2D meshes.
#'   If `NULL` or `'auto'`, the the likely fastest method is chosen;
#'   `'rectangle'` for flat 2D meshes with no CRS,
#'   `'sliced-spherical'` for CRS `'longlat'` meshes, and
#'   `'triangulated'` for all other meshes.
#' @param R Numerical value only applicable to spherical and geographical
#'   meshes. It is interpreted as `R` is the equivalent Earth radius, in km,
#'   used to scale the lambda intensity. For CRS enabled meshes, the default is
#'   6371. For CRS-less spherical meshes, the default is 1.
#' @param samplers A `SpatialPolygonsDataFrame` or `fm_mesh_2d` object.
#'   Simulated points that fall outside these polygons are discarded.
#' @param ignore.CRS logical; if `TRUE`, ignore any CRS information in the mesh.
#'   Default `FALSE`. This affects `R` and the permitted values for `strategy`.
#'
#' @return A `data.frame` (1D case),
#'   or `sf` (2D flat and 3D spherical surface cases,
#'   or 2D/2.5D surface cases with multiple samples).
#'   For multiple samples, the `data.frame` output has a
#'   column `'sample'` giving the index for each sample.
#' object of point locations.
#'
#' For the old (pre version `2.13.0.9034`) `sp` output format, use `as(ret,
#' "Spatial")` or `sf::as_Spatial(ret)`.
#'
#' @details
#' * For crs-less meshes on R2: Lambda is interpreted in the raw coordinate
#'   system. Output has an NA CRS.
#' * For crs-less meshes on S2: Lambda with raw units, after scaling the mesh
#'   to radius `R`, if specified.
#'   Output is given on the same domain as the mesh, with an NA CRS.
#' * For crs meshes on R2: Lambda is interpreted as per km^2, after scaling the
#'   globe to the Earth radius 6371 km, or `R`, if specified. Output given in
#'   the same CRS as the mesh.
#' * For crs meshes on S2: Lambda is interpreted as per km^2, after scaling the
#'   globe to the Earth radius 6371 km, or `R`, if specified. Output given in
#'   the same CRS as the mesh.
#'
#' @author Daniel Simpson \email{dp.simpson@@gmail.com} (base rectangle and
#'   spherical algorithms), Fabian E. Bachl \email{bachlfab@@gmail.com}
#'   (inclusion in inlabru, sliced spherical sampling), Finn Lindgren
#'   \email{finn.lindgren@@gmail.com} (extended CRS support, triangulated
#'   sampling)
#'
#' @examples
#' \donttest{
#' # The INLA package is required
#' if (bru_safe_inla()) {
#'   vertices <- seq(0, 3, by = 0.1)
#'   mesh <- fm_mesh_1d(vertices)
#'   loglambda <- 5 - 0.5 * vertices
#'   pts <- sample.lgcp(mesh, loglambda)
#'   pts$y <- 0
#'   plot(vertices, exp(loglambda), type = "l", ylim = c(0, 150))
#'   points(pts, pch = "|")
#' }
#' }
#'
#' \donttest{
#' # The INLA package is required
#' if (bru_safe_inla() &&
#'   require(ggplot2, quietly = TRUE) &&
#'   bru_safe_terra(quietly = TRUE) &&
#'   require("sf", quietly = TRUE)) {
#'   gorillas <- gorillas_sf
#'   pts <- sample.lgcp(gorillas$mesh,
#'     loglambda = 1.5,
#'     samplers = gorillas$boundary
#'   )
#'   ggplot() +
#'     gg(gorillas$mesh) +
#'     gg(pts)
#'
#'   pts <- sample.lgcp(gorillas$mesh,
#'     loglambda = base::scale(gorillas$mesh$loc) * 2,
#'     samplers = gorillas$boundary,
#'     ignore.CRS = TRUE
#'   )
#'   ggplot() +
#'     gg(gorillas$mesh) +
#'     gg(pts, aes(color = as.factor(sample))) +
#'     labs(color = "Sample")
#' }
#' }
#'
sample.lgcp <- function(mesh,
                        loglambda,
                        strategy = NULL,
                        R = NULL,
                        samplers = NULL,
                        ignore.CRS = FALSE) {
  mesh <- fm_as_fm(mesh)
  if (inherits(mesh, "fm_mesh_1d")) {
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

      # The B-spline basis expansion has a convex hull property which implies
      # that the maximum weight is not greater or equal to the maximum function
      # value.
      wmax <- max(loglambda)

      Npoints <- rpois(1, lambda = area * exp(wmax))
      if (Npoints > 0) {
        points <- runif(n = Npoints, min = xmin, max = xmax)
        proj <- fm_basis(mesh, points, full = TRUE)
        if (length(loglambda) == 1) {
          lambda_ratio <- exp(as.vector(
            Matrix::rowSums(proj$A) * loglambda
          ) - wmax)
        } else {
          lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - wmax)
        }
        keep <- (runif(Npoints) <= lambda_ratio)
        waste_ratio <- sum(keep) / length(keep)
        if (multi.samples) {
          multi.sample[[sample]] <- data.frame(
            x = points[keep],
            sample = sample
          )
        } else {
          multi.sample[[sample]] <- data.frame(x = points[keep])
        }
      } else {
        waste_ratio <- 0
        if (multi.samples) {
          multi.sample[[sample]] <- data.frame(
            x = numeric(0),
            sample = integer(0)
          )
        } else {
          multi.sample[[sample]] <- data.frame(x = numeric(0))
        }
      }
    }
    ret <- do.call(rbind, multi.sample)
  } else if (inherits(mesh, "fm_mesh_2d")) {
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
    input.crs <- fm_crs(mesh)
    use.crs <- !is.na(input.crs) && !ignore.CRS
    is.geocent <- fmesher::fm_crs_is_geocent(mesh)

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
      } else if (identical(sf::st_crs(input.crs)$proj, "longlat")) {
        strategy <- "sliced-spherical"
      }
    }
    strategy <- match.arg(strategy, strategies)

    if (is.geocent) {
      space.R <- mean(rowSums(mesh$loc^2)^0.5)
      internal.crs <- fm_crs("sphere")
      mesh$loc <- mesh$loc / space.R
      mesh$crs <- internal.crs
    } else {
      if (use.crs) {
        internal.crs <- fm_crs("sphere")
      } else {
        internal.crs <- fm_crs(NA_character_)
        mesh$crs <- NULL
      }
    }

    area.R <- 1
    if (!missing(R) && !is.null(R)) {
      area.R <- R
    } else {
      if (use.crs) {
        #        if (is.na(input.crs)) {
        #          space.units <- fm_length_unit(input.crs)
        #        }
        area.R <- 6371
        if (is.geocent) {
          if (abs(1 - space.R / area.R) > 1e-2) {
            warning(
              "The mesh has radius '", space.R,
              "', but crs information is available. ",
              "Using radius 6371 for area calculations."
            )
          }
        }
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
          area.mesh <- fm_transform(mesh, crs = internal.crs)
        } else {
          area.mesh <- mesh
          area.R <- 1
        }
        areas <- fm_fem(area.mesh, order = 0)$ta * area.R^2

        loglambda_tri <- matrix(
          loglambda[mesh$graph$tv],
          nrow(mesh$graph$tv),
          ncol(mesh$graph$tv)
        )
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
              mesh$loc[mesh$graph$tv[triangle, 1], , drop = FALSE]) *
              points[, 1] +
            (mesh$loc[mesh$graph$tv[triangle, 3], , drop = FALSE] -
              mesh$loc[mesh$graph$tv[triangle, 1], , drop = FALSE]) *
              points[, 2]
          if (is.geocent) {
            points <- points / rowSums(points^2)^0.5
          }
        }

        # Construct SpatialPoints
        if (is.geocent) {
          target.crs <- internal.crs
        } else if (use.crs) {
          target.crs <- fm_crs(input.crs)
        } else {
          target.crs <- fm_crs(NA_character_)
        }
        if (sum(Npoints) > 0) {
          points <- sf::st_as_sf(
            as.data.frame(points),
            coords = seq_len(ncol(points)),
            crs = target.crs
          )

          A <- fm_basis(mesh, points)
          lambda_ratio <- exp(as.vector(A %*% loglambda) -
            loglambda_max[triangle])
          keep <- (runif(sum(Npoints)) <= lambda_ratio)
          ret <- points[keep, , drop = FALSE]
          waste_ratio <- sum(keep) / length(keep)
        } else {
          points <- sf::st_as_sf(
            as.data.frame(matrix(0, 1, 3)),
            coords = seq_len(3),
            crs = target.crs
          )[-1, , drop = FALSE]
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
          ret <- sf::st_as_sf(
            as.data.frame(matrix(0, 1, 2)),
            coords = seq_len(2),
            crs = internal.crs
          )[-1, , drop = FALSE]
          waste_ratio <- 0
        } else {
          # Simulate uniform points on the bounding rectangle
          points <- sf::st_as_sf(
            data.frame(
              x = runif(n = Npoints, min = xmin, max = xmax),
              y = runif(n = Npoints, min = ymin, max = ymax)
            ),
            coords = seq_len(2),
            crs = internal.crs
          )

          # Do some thinning
          proj <- fm_basis(mesh, points, full = TRUE)
          lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - lambda_max)
          keep <- proj$ok & (runif(Npoints) <= lambda_ratio)
          ret <- points[keep, , drop = FALSE]
          waste_ratio <- sum(keep) / length(keep)
        }
      } else if (strategy == "spherical") {
        # Simulate number of points
        area <- 4 * pi * area.R^2
        lambda_max <- max(loglambda)
        Npoints <- rpois(n = 1, lambda = area * exp(lambda_max))
        if (Npoints > 5000000) {
          stop(paste0("Too many points!: ", Npoints))
        }

        if (Npoints == 0) {
          ret <- sf::st_as_sf(
            as.data.frame(matrix(0, 1, 3)),
            coords = seq_len(3),
            crs = internal.crs
          )[-1, , drop = FALSE]
          waste_ratio <- 0
        } else {
          # Choose z uniformly distributed in [-1,1].
          z <- runif(n = Npoints, min = -1, max = 1)
          # Choose t uniformly distributed on [0, 2*pi).
          t <- runif(n = Npoints, min = 0, max = 2 * pi)
          r <- sqrt(1 - z^2)
          x <- r * cos(t)
          y <- r * sin(t)

          points <- sf::st_as_sf(
            as.data.frame(x, y, z),
            coords = seq_len(3),
            crs = internal.crs
          )

          proj <- fm_basis(mesh, points, full = TRUE)
          lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - lambda_max)
          keep <- proj$ok & (runif(Npoints) <= lambda_ratio)
          ret <- points[keep, , drop = FALSE]
          waste_ratio <- sum(keep) / length(keep)
        }
      } else if (strategy == "sliced-spherical") {
        if (identical(input.crs$proj, "longlat")) {
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
        area <- (4 * pi * area.R^2) * lon.rel * lat.rel

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
            sampled.points[[k]] <-
              sf::st_as_sf(
                as.data.frame(matrix(0, 1, 3)),
                coords = seq_len(3),
                crs = internal.crs
              )[-1, , drop = FALSE]
            break
          }

          # Choose z uniformly distributed in [-1,1].
          # Transforms into latitude
          z <- runif(n = n.points, min = slat.range[1], max = slat.range[2])
          # Choose t uniformly distributed on [0, 2*pi).
          angle <- runif(
            n = n.points,
            min = radlon.range[1],
            max = radlon.range[2]
          )
          r <- sqrt(1 - z^2)
          x <- r * cos(angle)
          y <- r * sin(angle)

          points <- sf::st_as_sf(
            as.data.frame(x, y, z),
            coords = seq_len(3),
            crs = internal.crs
          )

          proj <- fm_basis(mesh, points, full = TRUE)
          lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - lambda_max)
          keep <- proj$ok & (runif(Npoints) <= lambda_ratio)
          sampled.points[[k]] <- points[keep, , drop = FALSE]

          n.keep <- n.keep + sum(keep)
          n.sampled <- n.sampled + length(keep)
        }

        # What to return
        if (sum(lengths(sampled.points)) > 0) {
          ret <- do.call(rbind, sampled.points)
        } else {
          ret <- sampled.points[[1]]
        }
        waste_ratio <- n.keep / n.sampled
      }

      if (multi.samples) {
        result[[sample]] <- cbind(ret, sample = rep(sample, NROW(ret)))
      }
    }
    if (multi.samples) {
      ret <- do.call(rbind, result)
    }

    if (is.geocent) {
      if (NROW(ret) > 0) {
        ret <- fm_transform(ret, crs = fm_crs("globe"))
      } else if (multi.samples) {
        if (use.crs) {
          ret_crs <- fm_crs(input.crs)
        } else {
          ret_crs <- fm_crs(NA_character_)
        }
        ret <- sf::st_as_sf(
          cbind(as.data.frame(matrix(0, 1, 3)), sample = 1),
          coords = seq_len(3),
          crs = ret_crs
        )[-1, , drop = FALSE]
      } else {
        ret <- sf::st_as_sf(
          as.data.frame(matrix(0, 1, 3)),
          coords = seq_len(3),
          crs = ret_crs
        )[-1, , drop = FALSE]
      }
    } else {
      if (use.crs) {
        if (NROW(ret) > 0) {
          ret <- fm_transform(ret, input.crs)
        } else if (multi.samples) {
          ret <- sf::st_as_sf(
            cbind(as.data.frame(matrix(0, 1, 2)), sample = 1),
            coords = seq_len(2),
            crs = fm_crs(input.crs)
          )[-1, , drop = FALSE]
        } else {
          ret <- sf::st_as_sf(
            as.data.frame(matrix(0, 1, 2)),
            coords = seq_len(2),
            crs = fm_crs(input.crs)
          )[-1, , drop = FALSE]
        }
      } else {
        sf::st_crs(ret) <- fm_crs(NA_character_)
      }
    }

    # Only retain points within the samplers
    if (!is.null(samplers) && (NROW(ret) > 0)) {
      if (inherits(samplers, "fm_mesh_2d")) {
        proj <- fm_basis(samplers, ret, full = TRUE)
        ret <- ret[proj$ok, , drop = FALSE]
      } else if (inherits(samplers, "Spatial")) {
        ret_ <- sf::as_Spatial(ret)
        if (use.crs) {
          ret_ <- fm_transform(ret_, fm_crs(samplers))
        } else {
          fm_crs(ret_) <- fm_crs(NA_character_)
          fm_crs(samplers) <- fm_crs(NA_character_)
        }
        ok <- !is.na(sp::over(ret_, samplers))
        ret <- ret[ok, , drop = FALSE]
      } else {
        ret_ <- sf::st_as_sf(ret)
        if (use.crs) {
          ret_ <- fm_transform(ret_, fm_crs(samplers))
        } else {
          fm_crs(ret_) <- fm_crs(NA_character_)
          fm_crs(samplers) <- fm_crs(NA_character_)
        }
        idx <- sf::st_within(ret_, samplers)
        ok <- vapply(idx, function(x) length(x) > 0, TRUE)
        ret <- ret[ok, , drop = FALSE]
      }
    }
  } else {
    stop(paste0(
      "The `mesh` must be convertible by `fm_as_fm()` to `fm_mesh_1d` ",
      "or `fm_mesh_2d`.\n",
      "  class = c('",
      paste0(class(mesh), collapse = "', '"),
      "')"
    ))
  }

  ret
}
