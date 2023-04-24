#' @title Multi-domain integration
#'
#' @description Construct integration points on tensor product spaces
#'
#' @export
fm_int <- function(domain, samplers = NULL, ...) {
  UseMethod("fm_int")
}

#' @export
#' @describeIn fm_int multi-domain integration
fm_int.list <- function(domain, samplers = NULL, ...) {
  stop("TODO: convert apmaker")
}

#' @export
#' @describeIn fm_int Discrete double or integer space integration
fm_int.numeric <- function(domain, samplers = NULL, name = "x", ...) {
  if (is.null(samplers)) {
    ips <- data.frame(
      x = as.vector(domain),
      weight = 1
    )
    colnames(ips)[1] <- name
    return(ips)
  }

  stopifnot(is.numeric(samplers))

  ok <- samplers %in% domain
  ips <- data.frame(
    x = as.vector(samplers),
    weight = 1
  )
  ips[["x"]][!ok] <- if (is.integer(samplers)) NA_integer_ else NA_real_
  storage.mode(ips[["x"]]) <- storage.mode(domain)
  colnames(ips)[1] <- name
  ips
}

#' @export
#' @describeIn fm_int Discrete character space integration
fm_int.character <- function(domain, samplers = NULL, name = "x", ...) {
  if (is.null(samplers)) {
    ips <- data.frame(
      x = as.vector(domain),
      weight = 1
    )
    colnames(ips)[1] <- name
    return(ips)
  }

  stopifnot(is.character(samplers))

  ok <- samplers %in% domain
  ips <- data.frame(
    x = as.vector(samplers),
    weight = 1
  )
  ips[["x"]][!ok] <- NA_character_
  colnames(ips)[1] <- name
  ips
}

#' @export
#' @describeIn fm_int Discrete factor space integration
fm_int.factor <- function(domain, samplers = NULL, name = "x", ...) {
  if (is.null(samplers)) {
    ips <- data.frame(
      x = factor(as.vector(domain), levels = levels(domain)),
      weight = 1
    )
    colnames(ips)[1] <- name
    return(ips)
  }

  stopifnot(is.factor(samplers) || is.character(samplers))

  ips <- data.frame(
    x = factor(as.vector(samplers), levels = levels(domain)),
    weight = 1
  )
  colnames(ips)[1] <- name
  ips
}


#' @export
#' @describeIn fm_int `inla.mesh.1d` integration
fm_int.inla.mesh.1d <- function(domain, samplers = NULL, name = "x", int.args = NULL, ...) {
  int.args.default <- list(method = "stable", nsub1 = 30, nsub2 = 9)
  if (is.null(int.args)) {
    int.args <- list()
  }
  missing.args <- setdiff(names(int.args.default), names(int.args))
  int.args[missing.args] <- int.args.default[missing.args]
  if (!is.null(int.args[["nsub"]])) {
    int.args[["nsub1"]] <- int.args[["nsub"]]
  }

  if (is.null(samplers)) {
    samplers <- rbind(domain$interval)
  } else if (is.null(dim(samplers))) {
    samplers <- matrix(samplers, 1, 2)
  }
  if (ncol(samplers) != 2) {
    stop("Interval description matrix must have 2 elements or be a 2-column matrix.")
  }

  ips <- list()
  for (j in seq_len(nrow(samplers))) {
    subsampler <- samplers[j, ]

    if (identical(int.args[["method"]], "stable")) {
      if (isTRUE(domain$cyclic)) {
        loc_trap <-
          sort(unique(c(domain$loc, as.vector(subsampler), domain$interval)))
      } else {
        loc_trap <- sort(unique(c(domain$loc, as.vector(subsampler))))
      }

      # Simpson's rule integration
      loc_mid <- (loc_trap[-1] + loc_trap[-length(loc_trap)]) / 2
      # Detect mid-points inside the samplers
      if (isTRUE(domain$cyclic) && (subsampler[1] > subsampler[2])) {
        inside <- (loc_mid < min(subsampler)) |
          (loc_mid > max(subsampler))
      } else {
        inside <- (loc_mid >= min(subsampler)) &
          (loc_mid <= max(subsampler))
      }
      weight_mid <- diff(loc_trap)
      weight_mid[!inside] <- 0.0

      weight_trap <- c(weight_mid / 2, 0) + c(0, weight_mid / 2)
      loc_simpson <- c(loc_trap, loc_mid)
      weight_simpson <- c(weight_trap / 3, weight_mid * 2 / 3)

      ok <- Matrix::rowSums(fm_evaluator(domain, loc_simpson)$proj$A) > 0

      ips[[j]] <- data.frame(
        x = loc_simpson[ok & (weight_simpson > 0)],
        weight = weight_simpson[ok & (weight_simpson > 0)],
        group = j
      )
      colnames(ips[[j]])[1] <- name
    } else {
      nsub <- int.args[["nsub1"]]
      u <- rep(
        (seq_len(nsub) - 0.5) / nsub,
        domain$n - 1
      )
      int_loc <-
        domain$loc[rep(seq_len(domain$n - 1), each = nsub)] * (1 - u) +
        domain$loc[rep(seq_len(domain$n - 1) + 1, each = nsub)] * u
      int_w <-
        (domain$loc[rep(seq_len(domain$n - 1) + 1, each = nsub)] -
           domain$loc[rep(seq_len(domain$n - 1), each = nsub)]) /
        nsub

      inside <-
        (int_loc >= min(subsampler)) &
        (int_loc <= max(subsampler))

      ips[[j]] <- data.frame(
        loc = int_loc[inside],
        weight = int_w[inside],
        group = j
      )
    }
    colnames(ips[[j]])[1] <- name
  }

  ips <- do.call(rbind, ips)

  ips
}







#' @export
#' @describeIn fm_int `inla.mesh.2d` integration
fm_int.inla.mesh <- function(domain, samplers = NULL, name = NULL, int.args = NULL, ...) {
  int.args.default <- list(method = "stable", nsub1 = 30, nsub2 = 9)
  if (is.null(int.args)) {
    int.args <- list()
  }
  missing.args <- setdiff(names(int.args.default), names(int.args))
  int.args[missing.args] <- int.args.default[missing.args]
  if (!is.null(int.args[["nsub"]])) {
    int.args[["nsub2"]] <- int.args[["nsub"]]
  }

  samplers_is_sp <- inherits(samplers, "Spatial")
  if (!is.null(samplers) && !samplers_is_sp) {
    samplers <- as(samplers, "Spatial")
  }

  if (inherits(domain, "inla.mesh") &&
      is.null(samplers) &&
      identical(int.args[["method"]], "stable")) {
    coord_names <- c("x", "y", "z")
    if (!is.null(name)) {
      coord_names[seq_along(name)] <- name
    }

    # transform to equal area projection
    if (!fm_crs_is_null(domain$crs)) {
      crs <- domain$crs
      samplers <- fm_transform(domain, crs = fm_crs("+proj=cea +units=km"))
    }

    ips <- vertices.inla.mesh(domain)
    ips$weight <- INLA::inla.mesh.fem(domain, order = 1)$va

    ips$group <- 1

    # backtransform
    if (!fm_crs_is_null(domain$crs)) {
      ips <- fm_transform(ips, crs = crs)
    }
    coordnames(ips) <- coord_names[seq_len(NCOL(coordinates(ips)))]
  } else if (inherits(samplers, "SpatialPointsDataFrame")) {
    if (!("weight" %in% names(samplers))) {
      warning("The integration points provided have no weight column. Setting weights to 1.")
      samplers$weight <- 1
    }

    ips <- samplers
  } else if (inherits(samplers, "SpatialPoints")) {
    ips <- samplers
    ips$weight <- 1
  } else if (inherits(samplers, "SpatialLines") ||
             inherits(samplers, "SpatialLinesDataFrame")) {
    # If SpatialLines are provided convert into SpatialLinesDataFrame and attach weight = 1
    if (inherits(samplers, "SpatialLines") &&
        !inherits(samplers, "SpatialLinesDataFrame")) {
      samplers <- SpatialLinesDataFrame(
        samplers,
        data = data.frame(weight = rep(1, length(samplers)))
      )
    }

    # Set weight to 1 if not provided
    if (!("weight" %in% names(samplers))) {
      samplers$weight <- 1
    }

    ips <- int.slines(
      samplers,
      domain,
      group = group,
      project = identical(int.args[["method"]], "stable")
    )

    coord_names <- c("x", "y", "z")
    if (!is.null(coordnames(samplers))) {
      coord_names[seq_along(coordnames(samplers))] <- coordnames(samplers)
    } else if (!is.null(name)) {
      coord_names[seq_along(name)] <- name
    }
    coordnames(ips) <- coord_names[seq_len(NCOL(coordinates(ips)))]
  } else if ((inherits(samplers, c(
    "SpatialPolygons",
    "SpatialPolygonsDataFrame"
  )) ||
  is.null(samplers))) {

    if (!is.null(samplers)) {
      # If SpatialPolygons are provided convert into SpatialPolygonsDataFrame and attach weight = 1
      if (class(samplers)[1] == "SpatialPolygons") {
        samplers <- SpatialPolygonsDataFrame(samplers,
                                             data = data.frame(weight = rep(1, length(samplers))),
                                             match.ID = FALSE
        )
      } else if (is.null(samplers@data[["weight"]])) {
        samplers@data[["weight"]] <- 1
      }

      cnames <- coordnames(samplers)
      samplers_crs <- fm_CRS(samplers)

      # Convert samplers and domain to equal area CRS
      if (!fm_crs_is_null(domain$crs)) {
        samplers <- fm_transform(samplers, crs = fm_crs("+proj=cea +units=km"))
      }
    }

    if (!fm_crs_is_null(domain$crs)) {
      domain <- fm_transform(domain, crs = fm_crs("+proj=cea +units=km"))
    }
    domain_crs <- fm_crs(domain$crs)
    domain_crs <- fm_CRS(domain_crs)

    if (identical(int.args[["poly_method"]], "legacy")) {
      lifecycle::deprecate_stop("2.8.0", I("int.args$poly_method == 'legacy'"))
    } else {
      if (!is.null(int.args$use_new) && !int.args$use_new) {

        # This handles holes
        poly_segm <- fm_as_inla_mesh_segment(samplers, join = FALSE)
        poly_segm <- lapply(
          seq_along(poly_segm),
          function(k) {
            segm <- poly_segm[[k]]
            segm[["grp"]] <- rep(k, NROW(segm[["idx"]]))
            segm[["is.bnd"]] <- TRUE
            segm
          }
        )

        ips <- bru_int_polygon_old(
          domain,
          polylist = poly_segm,
          method = int.args$method,
          nsub = int.args$nsub2
        )
      } else {
        ips <- bru_int_polygon(
          domain,
          method = int.args$method,
          nsub = int.args$nsub2,
          samplers = samplers
        )
      }
    }

    if (!is.null(samplers)) {
      df <- data.frame(
        samplers@data[ips$group, group, drop = FALSE],
        weight = ips[, "weight"] * samplers@data[ips$group, "weight"]
      )
    } else {
      df <- as.data.frame(ips)[, !(names(ips) %in% c("x", "y", "z")), drop = FALSE]
    }
    if (is.null(ips$z)) {
      ips <- sp::SpatialPointsDataFrame(ips[, c("x", "y")],
                                        data = df,
                                        match.ID = FALSE, proj4string = domain_crs
      )
    } else {
      ips <- sp::SpatialPointsDataFrame(ips[, c("x", "y", "z")],
                                        data = df,
                                        match.ID = FALSE, proj4string = domain_crs
      )
    }

    if (!fm_crs_is_null(domain_crs) && !fm_crs_is_null(samplers_crs)) {
      ips <- fm_transform(ips, crs = samplers_crs)
    }

    coord_names <- c("x", "y", "z")
    if (!is.null(samplers) && !is.null(coordnames(samplers))) {
      coord_names[seq_along(coordnames(samplers))] <- coordnames(samplers)
    } else if (!is.null(name)) {
      coord_names[seq_along(name)] <- name
    }
    coordnames(ips) <- coord_names[seq_len(NCOL(coordinates(ips)))]
  } else {
    stop("No integration handling code reached; please notify the package developer.")
  }

  if (!samplers_is_sp && inherits(ips, "Spatial")) {
    ips <- sf::st_as_sf(ips)
  }

  ips
}
