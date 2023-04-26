#' @title Multi-domain integration
#'
#' @description Construct integration points on tensor product spaces
#' @param domain Functional space specification; single domain or a named list
#' of domains
#' @param samplers For single domain `fm_int` methods, an object specifying one or more
#' subsets of the domain, and optional weighting and group information, as variables `weight` and `group`
#' in a data frame.
#' For `fm_int.list`, a list of sampling definitions.
#' @param name For single-domain methods, the variable name to use for the
#' integration points. Default 'x'
#' @param \dots Additional arguments passed on to other methods
#'
#' @export
#' @examples
#' if (bru_safe_inla()) {
#'   # Integration on the interval (2, 3.5) with Simpson's rule
#'   fm_int(INLA::inla.mesh.1d(0:4), samplers = cbind(2, 3.5))
#' }
#'
fm_int <- function(domain, samplers = NULL, ...) {
  UseMethod("fm_int")
}

#' @export
#' @describeIn fm_int Multi-domain integration
fm_int.list <- function(domain, samplers = NULL, ...) {
  if (is.null(names(domain))) {
    stop("For 'fm_int.list', the domain must be a named list.")
  }

  if (!is.null(samplers) && !inherits(samplers, "list")) {
    samplers <- list(samplers)
  }

  #####################################################################################
  # 20221213 This function to figure out which samplers the domain integration are
  # class check is not needed here but to be done in S3 method
  # 20221212 We have to check the matching here, that is
  # 1) certain domain classes have to match certain samplers classes.
  # TODO 2) in which 1) has to check the domain one by one if it is a list
  # if (any(unlist(lapply(domain, function(x) {
  #   inherits(x, c("character", "factor", "numeric", "inla.mesh.1d"))
  # })))) {
  #   stopifnot_class(samplers, c("character", "factor", "numeric", "data.frame"))
  # }
  # if (any(unlist(lapply(domain, function(x) {
  #   inherits(x, cc(
  #     "inla.mesh",
  #     "inla.mesh.lattice",
  #     "raster",
  #     "SpatRaster"
  #   ))
  # })))) {
  #   stopifnot_class(samplers, c("sf", "sfc", "Spatial"))
  # }

  # Change a mix of sp and sf objects to sf
  sf_samplers <- unlist(lapply(samplers, function(x) inherits(x, c("sf", "sfc"))))
  sp_samplers <- unlist(lapply(samplers, function(x) inherits(x, "Spatial")))
  if (any(sp_samplers)) {
    if (any(sf_samplers)) {
      warning("Both `sf` and `sp` objects in the samplers are detected. Output will be `sf`.")
    }
    samplers[sp_samplers] <- lapply(samplers[sp_samplers], sf::st_as_sf)
    if (!("coordinates" %in% names(domain))) {
      stop("`sp` input detected but no `coordinates` domain present.")
    }
    names(domain)[names(domain) %in% "coordinates"] <- "geometry"
  }

  # TODO 20220126 lapply to extract the names, the current one is not sufficient
  # TODO Sort multidomain samplers, single domain samplers
  # TODO Multidomain samplers happens when a sampler across several domains. How to detect that?
  # TODO we should then do the name check for domain and samplers here
  # TODO remove sampler domains and full domain samplers
  # TODO some thoughts for S3 methods, there should be an extra layer ie function to sort samplers and domain arguments and difine multidomain, singledomain and full domain s3 class
  #######################
  names_domain <- names(domain)
  names_lsamplers <- names(samplers)
  if (is.null(names_lsamplers)) {
    names_lsamplers <- rep("", length(samplers))
  }
  index_single_samplers <- which(names_lsamplers != "")
  index_multi_samplers <- which(names_lsamplers == "")
  names_samplers <- as.list(names_lsamplers)
  names_samplers[index_multi_samplers] <- names_list(samplers[index_multi_samplers])
  names_reserved <- c("weight") # coordinate and geometry is not required here

  if (length(intersect(names_domain, names_reserved)) > 0) {
    stop(paste0("The reserved names ",
                paste0(intersect(names_domain, names_reserved), collapse = ", "),
                " cannot be used as domain names."))
  }

  lips_samplers <- list()

  #######################
  # multidomain samplers, ie unnamed element(s) in samplers, for each sampler and then for each domain(lapply)
  # TODO still have to deal with secondary geometry
  for (i in index_multi_samplers) {
    if (is.null(names_samplers[[i]])) {
      stop(paste0("The unnamed sampler #", i, " in the samplers is NULL"))
    }
    names_intersect <- intersect(names_samplers[[i]], names_domain)
    lips_multidomainsampler <- lapply(
      names_intersect,
      function(nm) fm_int(
        domain = domain[[nm]],
        samplers = samplers[[i]][[nm]],
        name = nm,
        #        group = names_intersect, # block=group should be the grouping, say season,
        ...
      ))
    lips_samplers[[i]] <- do.call(cprod, lips_multidomainsampler)
  }


  #######################
  # singledomain samplers, ie named element(s) in samplers
  for (i in index_single_samplers) {
    nm <- intersect(names_samplers[[i]], names_domain)
    stopifnot(length(nm) == 1)
    lips_samplers[[i]] <-
      fm_int(
        domain = domain[[nm]],
        samplers = samplers[[i]],
        name = nm,
        #      group = names_intersect, # block=group should be the grouping, say season,
        ...
      )
  }

  # Full domain samplers
  names_full_domain_samplers <- setdiff(names_domain, unlist(names_samplers))
  lips_full_domain_samplers <-
    lapply(
      names_full_domain_samplers,
      function(nm) {
        fm_int(
          domain = domain[[nm]],
          name = nm,
          #      group = names_intersect, # block=group should be the grouping, say season,
          ...
        )
      })

  ips <- do.call(cprod, c(lips_samplers, lips_full_domain_samplers))

  if (any(sp_samplers) && !any(sf_samplers)) {
    ips <- sf::as_Spatial(ips)
  }

  ips
}


#' @export
#' @describeIn fm_int Discrete double or integer space integration
fm_int.numeric <- function(domain, samplers = NULL, name = "x", ...) {
  if (is.null(samplers)) {
    ips <- data.frame(
      x = as.vector(domain),
      weight = 1,
      group = 1L
    )
    colnames(ips)[1] <- name
    return(ips)
  }

  if (!is.data.frame(samplers)) {
    samplers <- data.frame(
      x = samplers,
      weight = 1,
      group = 1L
    )
    colnames(samplers)[1] <- name
  } else {
    if (is.null(samplers[["weight"]])) {
      samplers[["weight"]] <- 1
    }
    if (is.null(samplers[["group"]])) {
      samplers[["group"]] <- 1L
    }
  }

  storage.mode(samplers[[name]]) <- storage.mode(domain)

  ok <- samplers[[name]] %in% domain
  ips <- samplers[ok, , drop = FALSE]
  ips
}

#' @export
#' @describeIn fm_int Discrete character space integration
fm_int.character <- function(domain, samplers = NULL, name = "x", ...) {
  if (is.null(samplers)) {
    ips <- data.frame(
      x = as.vector(domain),
      weight = 1,
      group = 1L
    )
    colnames(ips)[1] <- name
    return(ips)
  }

  if (!is.data.frame(samplers)) {
    samplers <- data.frame(
      x = samplers,
      weight = 1,
      group = 1L
    )
    colnames(samplers)[1] <- name
  } else {
    if (is.null(samplers[["weight"]])) {
      samplers[["weight"]] <- 1
    }
    if (is.null(samplers[["group"]])) {
      samplers[["group"]] <- 1L
    }
  }

  storage.mode(samplers[[name]]) <- storage.mode(domain)

  ok <- samplers[[name]] %in% domain
  ips <- samplers[ok, , drop = FALSE]
  ips
}

#' @export
#' @describeIn fm_int Discrete factor space integration
fm_int.factor <- function(domain, samplers = NULL, name = "x", ...) {
  if (is.null(samplers)) {
    ips <- data.frame(
      x = as.vector(domain),
      weight = 1,
      group = 1L
    )
    colnames(ips)[1] <- name
    return(ips)
  }

  if (!is.data.frame(samplers)) {
    samplers <- data.frame(
      x = factor(as.vector(samplers), levels = levels(domain)),
      weight = 1,
      group = 1L
    )
    colnames(samplers)[1] <- name
  } else {
    if (is.null(samplers[["weight"]])) {
      samplers[["weight"]] <- 1
    }
    if (is.null(samplers[["group"]])) {
      samplers[["group"]] <- 1L
    }
  }

  ok <- samplers[[name]] %in% domain
  ips <- samplers[ok, , drop = FALSE]
  ips
}


#' @export
#' @describeIn fm_int `SpatRaster` integration. Not yet implemented.
fm_int.SpatRaster <- function(domain, samplers = NULL, name = "x", ...) {
  stop("'SpatRaster' integration is not yet implemented.")
}

#' @export
#' @describeIn fm_int `inla.mesh.lattice` integration. Not yet implemented.
fm_int.inla.mesh.lattice <- function(domain, samplers = NULL, name = "x", ...) {
  stop("'inla.mesh.lattice' integration is not yet implemented.")
}

#' @param int.args A list of integration options
#' @export
#' @describeIn fm_int `inla.mesh.1d` integration.
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
    samplers <- data.frame(
      x1 = domain$interval[1],
      x2 = domain$interval[2],
      weight = 1,
      group = 1L
    )
    colnames(samplers)[1:2] <- paste0(name, 1:2)
  } else if (is.null(dim(samplers))) {
    samplers <- data.frame(
      x1 = samplers[1],
      x2 = samplers[2],
      weight = 1,
      group = 1L
    )
    colnames(samplers)[1:2] <- paste0(name, 1:2)
  } else {
    samplers <- as.data.frame(samplers)
    colnames(samplers)[1:2] <- paste0(name, 1:2)
    if (!("weight" %in% colnames(samplers))) {
      samplers <- cbind(samplers, weight = 1)
    }
    if (!("group" %in% colnames(samplers))) {
      samplers <- cbind(samplers, group = seq_len(NROW(samplers)))
    }
  }

  ips <- list()
  for (j in seq_len(nrow(samplers))) {
    subsampler <- unlist(samplers[j, 1:2, drop = TRUE])
    theweight <- samplers[j, "weight"]
    thegroup <- samplers[j, "group"]

    if (isTRUE(domain$cyclic)) {
      if (diff(subsampler) >= diff(domain$interval)) {
        subsampler <- domain$interval
      } else {
        subsampler[1] <- domain$interval[1] +
          (subsampler[1] - domain$interval[1]) %% diff(domain$interval)
        subsampler[2] <- subsampler[1] + diff(subsampler) %% diff(domain$interval)
        if (diff(subsampler) == 0.0) {
          subsampler <- domain$interval
        } else if (subsampler[2] > domain$interval[2]) {
          subsampler[2] <- domain$interval[1] +
            (subsampler[2] - domain$interval[1]) %% diff(domain$interval)
        }
      }
    } else if (diff(subsampler) <= 0.0) {
      # Empty interval, skip to next subsampler
      next
    }

    if (identical(int.args[["method"]], "stable")) {
      if (isTRUE(domain$cyclic)) {
        loc_trap <- sort(unique(pmin(
          domain$interval[2],
          pmax(
            domain$interval[1],
            c(
              domain$loc,
              as.vector(subsampler),
              domain$interval
            )
          )
        )))
      } else {
        loc_trap <- sort(unique(pmin(
          domain$interval[2],
          pmax(
            domain$interval[1],
            c(
              domain$loc,
              as.vector(subsampler)
            )
          )
        )))
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
        weight = weight_simpson[ok & (weight_simpson > 0)] * theweight,
        group = thegroup
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

      if (isTRUE(domain$cyclic) && (subsampler[1] > subsampler[2])) {
        inside <- (int_loc < min(subsampler)) |
          (int_loc > max(subsampler))
      } else {
        inside <- (int_loc >= min(subsampler)) &
          (int_loc <= max(subsampler))
      }

      ips[[j]] <- data.frame(
        loc = int_loc[inside],
        weight = int_w[inside] * theweight,
        group = thegroup
      )
    }
    colnames(ips[[j]])[1] <- name
  }

  ips <- do.call(rbind, ips)

  if (NROW(ips) == 0) {
    ips <- data.frame(x = numeric(0), weight = numeric(0), group = integer(0))
    colnames(ips)[1] <- name
  }

  ips
}




#' Subset integration on a mesh
#'
#' Integration methods for spatial samplers on `inla.mesh` meshes.
#'
#' @inheritParams fm_int
#' @export
#' @keywords internal
fm_int_inla_mesh <- function(samplers,
                             domain,
                             name = NULL,
                             int.args = NULL,
                             ...) {
  stopifnot(inherits(domain, "inla.mesh"))
  UseMethod("fm_int_inla_mesh")
}

#' @export
#' @describeIn fm_int_inla_mesh Full domain integration
fm_int_inla_mesh.NULL <- function(samplers,
                                  domain,
                                  name = NULL,
                                  int.args = NULL,
                                  ...) {
  ipsl <- bru_int_polygon(domain,
                          samplers = NULL,
                          method = int.args[["method"]],
                          nsub = int.args$nsub2)

  ips <- sf::st_as_sf(ipsl,
                      coords = intersect(c("x", "y", "z"), names(ipsl)),
                      crs = fm_crs(domain))

  ips
}

#' @export
#' @describeIn fm_int_inla_mesh General integration
fm_int_inla_mesh.default <- function(samplers,
                                     domain,
                                     name = NULL,
                                     int.args = NULL,
                                     ...) {
  if (is.null(samplers) &&
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
      warning("The integration points provided have no weight column. Setting weight to 1.")
      samplers$weight <- 1
    }
    if (!("group" %in% names(samplers))) {
      samplers$group <- 1L
    }

    ips <- samplers
  } else if (inherits(samplers, "SpatialPoints")) {
    ips <- samplers
    ips$weight <- 1
    ips$group <- 1L
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
    if (!("group" %in% names(samplers))) {
      samplers$group <- 1
    }

    ips <- int.slines(
      samplers,
      domain,
      group = "group",
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
        samplers@data[ips$group, "group", drop = FALSE],
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

  ips <- sf::st_as_sf(ips)

  ips
}

#' @export
#' @describeIn fm_int_inla_mesh `Spatial` integration
fm_int_inla_mesh.Spatial <- function(samplers,
                                     domain,
                                     name = NULL,
                                     int.args = NULL,
                                     ...) {
  samplers <- sf::st_as_sf(samplers)

  ips <-
    fm_int_inla_mesh(samplers,
                     domain = domain,
                     name = name,
                     int.args = int.args,
                     ...)

  ips <- as(ips, "Spatial")

  ips
}

#' @export
#' @describeIn fm_int `inla.mesh` integration
fm_int.inla.mesh <- function(domain,
                             samplers = NULL,
                             name = NULL,
                             int.args = NULL,
                             ...) {
  int.args.default <- list(method = "stable", nsub1 = 30, nsub2 = 9)
  if (is.null(int.args)) {
    int.args <- list()
  }
  missing.args <- setdiff(names(int.args.default), names(int.args))
  int.args[missing.args] <- int.args.default[missing.args]
  if (!is.null(int.args[["nsub"]])) {
    int.args[["nsub2"]] <- int.args[["nsub"]]
  }

  ips <- fm_int_inla_mesh(samplers,
                          domain = domain,
                          name = name,
                          int.args = int.args,
                          ...)
  return(ips)

  samplers_is_sp <- inherits(samplers, "Spatial")
  if (!is.null(samplers) && !samplers_is_sp) {
    samplers <- as(samplers, "Spatial")
  }

  if (is.null(samplers) &&
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
      warning("The integration points provided have no weight column. Setting weight to 1.")
      samplers$weight <- 1
    }
    if (!("group" %in% names(samplers))) {
      samplers$group <- 1L
    }

    ips <- samplers
  } else if (inherits(samplers, "SpatialPoints")) {
    ips <- samplers
    ips$weight <- 1
    ips$group <- 1L
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
    if (!("group" %in% names(samplers))) {
      samplers$group <- 1
    }

    ips <- int.slines(
      samplers,
      domain,
      group = "group",
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
        samplers@data[ips$group, "group", drop = FALSE],
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
