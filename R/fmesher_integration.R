#' @title Multi-domain integration
#'
#' @description Construct integration points on tensor product spaces
#' @param domain Functional space specification; single domain or a named list
#' of domains
#' @param samplers For single domain `fm_int` methods, an object specifying one or more
#' subsets of the domain, and optional weighting in a `weight` variable.
#' For `fm_int.list`, a list of sampling definitions, where data frame elements
#' may contain information for multiple domains, in which case each row represent
#' a separate tensor product integration subspace.
#' @param name For single-domain methods, the variable name to use for the
#' integration points. Default 'x'
#' @param \dots Additional arguments passed on to other methods
#'
#' @returns A `data.frame`, `tibble`, `sf`, or `SpatialPointsDataFrame` of 1D and
#' 2D integration points, including a `weight` column and `.block` column.

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

#' Multi-domain sampler integration
#'
#' Combine integration over different domains
#'
#' @export
#' @keywords internal
fm_int_multi_sampler <- function(domain, samplers, ...) {
  if (is.null(names(domain))) {
    stop("For 'fm_int_multi_sampler', the domain must be a named list.")
  }

  names_domain <- names(domain)
  names_samplers <- names(samplers)
  names_reserved <- c("weight", ".block")

  if (length(intersect(names_domain, names_reserved)) > 0) {
    stop(paste0(
      "The reserved name(s) ",
      paste0("'",
        intersect(names_domain, names_reserved),
        "'",
        collapse = ", "
      ),
      " cannot be used as domain names."
    ))
  }

  names_intersect <- intersect(names(samplers), names(domain))
  ips_list <- lapply(
    names_intersect,
    function(nm) {
      fm_int(
        domain = domain[[nm]],
        samplers = samplers[[nm]],
        name = nm,
        ...
      )
    }
  )
  ips <- do.call(cprod, c(ips_list, list(.blockwise = TRUE)))

  if ("weight" %in% names(samplers)) {
    ips$weight <- ips$weight * samplers$weight[ips$.block]
  }

  ips
}


#' @export
#' @describeIn fm_int Multi-domain integration
fm_int.list <- function(domain, samplers = NULL, ...) {
  weight_name <- "weight"

  if (is.null(names(domain))) {
    stop("For 'fm_int.list', the domain must be a named list.")
  }

  if (!is.null(samplers) && !inherits(samplers, "list")) {
    samplers <- list(samplers)
  }

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
  names_samplers[index_multi_samplers] <- lapply(samplers[index_multi_samplers], names)
  names_reserved <- c("weight", ".block") # coordinate and geometry is not required here

  if (length(intersect(names_domain, names_reserved)) > 0) {
    stop(paste0(
      "The reserved names ",
      paste0(intersect(names_domain, names_reserved), collapse = ", "),
      " cannot be used as domain names."
    ))
  }

  lips_samplers <- list()

  #######################
  # multidomain samplers, ie unnamed element(s) in samplers, for each sampler and then for each domain(lapply)
  # TODO still have to deal with secondary geometry
  for (i in index_multi_samplers) {
    if (is.null(names_samplers[[i]])) {
      stop(paste0("The unnamed sampler #", i, " in the samplers is NULL"))
    }
    lips_samplers[[i]] <-
      fm_int_multi_sampler(
        domain = domain,
        samplers = samplers[[i]],
        ...
      )
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
          ...
        )
      }
    )


  ips <- do.call(cprod, c(
    lips_samplers,
    lips_full_domain_samplers,
    list(.blockwise = FALSE)
  ))

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
      .block = 1L
    )
    colnames(ips)[1] <- name
    return(ips)
  }

  if (!is.data.frame(samplers)) {
    samplers <- data.frame(
      x = samplers,
      weight = 1,
      .block = seq_len(NROW(samplers))
    )
    colnames(samplers)[1] <- name
  } else {
    if (is.null(samplers[["weight"]])) {
      samplers[["weight"]] <- 1
    }
    samplers[[".block"]] <- seq_len(NROW(samplers))
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
      .block = 1L
    )
    colnames(ips)[1] <- name
    return(ips)
  }

  if (!is.data.frame(samplers)) {
    samplers <- data.frame(
      x = samplers,
      weight = 1,
      .block = seq_len(NROW(samplers))
    )
    colnames(samplers)[1] <- name
  } else {
    if (is.null(samplers[["weight"]])) {
      samplers[["weight"]] <- 1
    }
    samplers[[".block"]] <- seq_len(NROW(samplers))
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
      .block = 1L
    )
    colnames(ips)[1] <- name
    return(ips)
  }

  if (!is.data.frame(samplers)) {
    samplers <- data.frame(
      x = factor(as.vector(samplers), levels = levels(domain)),
      weight = 1,
      .block = seq_len(NROW(samplers))
    )
    colnames(samplers)[1] <- name
  } else {
    if (is.null(samplers[["weight"]])) {
      samplers[["weight"]] <- 1
    }
    samplers[[".block"]] <- seq_len(NROW(samplers))
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



# inla.mesh.1d integration ####

#' @param int.args A list of integration options
#' @export
#' @describeIn fm_int `inla.mesh.1d` integration. Supported samplers:
#' * `NULL` for integration over the entire domain;
#' * A length 2 vector defining an interval;
#' * A 2-column matrix with a single interval in each row;
#' * A tibble with a named column containing a matrix, and optionally a
#'  `weight` column.
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
    samplers <- tibble::tibble(
      x = cbind(domain$interval[1], domain$interval[2]),
      weight = 1,
      .block = 1L
    )
    colnames(samplers)[1] <- name
  } else if (is.null(dim(samplers))) {
    samplers <- tibble::tibble(
      x = cbind(samplers[1], samplers[2]),
      weight = 1,
      .block = 1L
    )
    colnames(samplers)[1] <- name
  } else if (is.matrix(samplers)) {
    samplers <- tibble::tibble(
      x = samplers,
      weight = 1,
      .block = seq_len(NROW(samplers))
    )
    colnames(samplers)[1] <- name
  } else {
    samplers <- tibble::as_tibble(samplers)
    if (!(name %in% colnames(samplers))) {
      stop(paste0("Domain name '", name, "' missing from sampler."))
    }
    if (!("weight" %in% colnames(samplers))) {
      samplers$weight <- 1
    }
    samplers$.block <- seq_len(NROW(samplers))
  }

  ips <- list()
  for (j in seq_len(nrow(samplers))) {
    subsampler <- samplers[[name]][j, , drop = TRUE]
    theweight <- samplers[j, "weight", drop = TRUE]
    the.block <- samplers[j, ".block", drop = TRUE]

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
        .block = the.block
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
        .block = the.block
      )
    }
    colnames(ips[[j]])[1] <- name
  }

  ips <- do.call(rbind, ips)

  if (NROW(ips) == 0) {
    ips <- data.frame(x = numeric(0), weight = numeric(0), .block = integer(0))
    colnames(ips)[1] <- name
  }

  ips
}


# inla.mesh integration ####

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
    nsub = int.args$nsub2
  )

  ips <- sf::st_as_sf(ipsl,
    coords = intersect(c("x", "y", "z"), names(ipsl)),
    crs = fm_crs(domain)
  )

  if (!is.null(name)) {
    ips <- tibble::as_tibble(ips)
    names(ips)[names(ips) == "geometry"] <- name
    ips <- sf::st_as_sf(ips, sf_column_name = name)
  }

  # TODO: figure out a way to decide the output format
  ips <- sf::as_Spatial(ips)

  ips
}

#' @export
#' @describeIn fm_int_inla_mesh `sf` integration
fm_int_inla_mesh.sf <- function(samplers,
                                domain,
                                name = NULL,
                                int.args = NULL,
                                ...) {
  if (is.null(name)) {
    name <- attr(samplers, "sf_column")
  }
  if (!("weight" %in% names(samplers))) {
    weight <- rep(1, NROW(samplers))
  } else {
    weight <- samplers$weight
  }

  fm_int_inla_mesh(sf::st_geometry(samplers),
    domain,
    name = name,
    int.args = int.args,
    .weight = weight,
    ...
  )
}

#' @param .weight Optional weight vector for `sfc_*` integration
#' @param .block Optional block grouping vector for `sfc_*` integration
#' @export
#' @describeIn fm_int_inla_mesh `sfc_POINT` integration
fm_int_inla_mesh.sfc_POINT <- function(samplers,
                                       domain,
                                       name = NULL,
                                       int.args = NULL,
                                       .weight = rep(1, NROW(samplers)),
                                       ...) {
  ips <- tibble::tibble(
    geometry = samplers,
    weight = .weight,
    .block = seq_len(NROW(samplers))
  )
  names(ips)[names(ips) == "geometry"] <- name
  ips <- sf::st_as_sf(ips, sf_column_name = name)

  # TODO: remove points outside the domain

  ips
}

#' @export
#' @describeIn fm_int_inla_mesh `sfc_LINESTRING` integration
fm_int_inla_mesh.sfc_LINESTRING <- function(samplers,
                                            domain,
                                            name = NULL,
                                            int.args = NULL,
                                            .weight = rep(1, NROW(samplers)),
                                            ...) {
  samplers <- sf::as_Spatial(samplers)
  samplers$weight <- .weight
  ips <- fm_int_inla_mesh(samplers,
    domain = domain,
    name = name,
    int.args = int.args, ...
  )
  ips <- sf::st_as_sf(ips)
  ips <- tibble::as_tibble(ips)
  names(ips)[names(ips) == "geometry"] <- name
  ips <- sf::st_as_sf(ips, sf_column_name = name)
  ips
}

#' @export
#' @describeIn fm_int_inla_mesh `sfc_POLYGON` integration
fm_int_inla_mesh.sfc_POLYGON <- function(samplers,
                                         domain,
                                         name = NULL,
                                         int.args = NULL,
                                         .weight = rep(1, NROW(samplers)),
                                         ...) {
  weight = .weight
  .block = seq_len(NROW(samplers))

#  samplers_crs <- fm_crs(samplers)

#  # Convert samplers and domain to equal area CRS
#  if (!fm_crs_is_null(domain$crs)) {
#    samplers <- fm_transform(samplers, crs = fm_crs("+proj=cea +units=km"))
#  }

#  if (!fm_crs_is_null(domain$crs)) {
#    domain <- fm_transform(domain, crs = fm_crs("+proj=cea +units=km"))
#  }
#  domain_crs <- fm_crs(domain$crs)

  if (identical(int.args[["poly_method"]], "legacy")) {
    lifecycle::deprecate_stop("2.8.0", I("int.args$poly_method == 'legacy'"))
  } else if (!is.null(int.args$use_new) && !int.args$use_new) {
    lifecycle::deprecate_stop("2.8.0", I("int.args$use_new == FALSE"))
  }

  ips <- bru_int_polygon_sf(
    domain,
    method = int.args$method,
    nsub = int.args$nsub2,
    samplers = samplers
  )

  ips$weight <- ips$weight * .weight[ips$.block]
  ips$.block = .block[ips$.block]

#  if (!fm_crs_is_null(domain_crs) && !fm_crs_is_null(samplers_crs)) {
#    ips <- fm_transform(ips, crs = domain_crs)
#  }

  ips <- tibble::as_tibble(ips)
  names(ips)[names(ips) == "geometry"] <- name
  ips <- sf::st_as_sf(ips, sf_column_name = name)
  ips
}

#' @export
#' @describeIn fm_int_inla_mesh `sfc_MULTIPOLYGON` integration
fm_int_inla_mesh.sfc_MULTIPOLYGON <- function(samplers,
                                              domain,
                                              name = NULL,
                                              int.args = NULL,
                                              .weight = rep(1, NROW(samplers)),
                                              ...) {
  samplers <- sf::as_Spatial(samplers)
  samplers$weight <- .weight
  ips <- fm_int_inla_mesh(samplers,
    domain = domain,
    name = name,
    int.args = int.args, ...
  )
  ips <- sf::st_as_sf(ips)
  ips <- tibble::as_tibble(ips)
  names(ips)[names(ips) == "geometry"] <- name
  ips <- sf::st_as_sf(ips, sf_column_name = name)
  ips
}



#' @export
#' @describeIn fm_int_inla_mesh SpatialPoints integration
fm_int_inla_mesh.SpatialPoints <- function(samplers,
                                           domain,
                                           name = NULL,
                                           int.args = NULL,
                                           ...) {
  if (!inherits(samplers, "SpatialPointsDataFrame")) {
    # If SpatialPoints are provided convert into SpatialPointsDataFrame and attach weight = 1
    samplers <- SpatialLinesDataFrame(
      samplers,
      data = data.frame(
        weight = rep(1, length(samplers)),
        .block <- seq_len(NROW(samplers))
      )
    )
  }
  if (!("weight" %in% names(samplers))) {
    samplers$weight <- 1
  }
  if (!(".block" %in% names(samplers))) {
    samplers$.block <- seq_len(NROW(samplers))
  }

  # TODO: remove points outside the domain

  samplers
}

#' @export
#' @describeIn fm_int_inla_mesh SpatialLines integration
fm_int_inla_mesh.SpatialLines <- function(samplers,
                                          domain,
                                          name = NULL,
                                          int.args = NULL,
                                          ...) {
  if (!inherits(samplers, "SpatialLinesDataFrame")) {
    samplers <- SpatialLinesDataFrame(
      samplers,
      data = data.frame(
        weight = rep(1, length(samplers)),
        .block = seq_len(NROW(samplers))
      )
    )
  }

  # Set weight to 1 if not provided
  if (!("weight" %in% names(samplers))) {
    samplers$weight <- 1
  }
  if (!(".block" %in% names(samplers))) {
    samplers$.block <- seq_len(NROW(samplers))
  }

  ips <- int.slines(
    samplers,
    domain,
    .block = ".block",
    project = identical(int.args[["method"]], "stable")
  )

  coord_names <- c("x", "y", "z")
  if (!is.null(coordnames(samplers))) {
    coord_names[seq_along(coordnames(samplers))] <- coordnames(samplers)
  } else if (!is.null(name)) {
    coord_names[seq_along(name)] <- name
  }
  coordnames(ips) <- coord_names[seq_len(NCOL(coordinates(ips)))]
  ips
}

#' @export
#' @describeIn fm_int_inla_mesh SpatialPolygons integration
fm_int_inla_mesh.SpatialPolygons <- function(samplers,
                                             domain,
                                             name = NULL,
                                             int.args = NULL,
                                             ...) {
  if (!inherits(samplers, "SpatialPolygonsDataFrame")) {
    samplers <- SpatialPolygonsDataFrame(
      samplers,
      data = data.frame(
        weight = rep(1, length(samplers)),
        .block = seq_len(NROW(samplers))
      ),
      match.ID = FALSE
    )
  }

  # Set weight to 1 if not provided
  if (!("weight" %in% names(samplers))) {
    samplers$weight <- 1
  }
  if (!(".block" %in% names(samplers))) {
    samplers$.block <- seq_len(NROW(samplers))
  }

  cnames <- coordnames(samplers)
  samplers_crs <- fm_CRS(samplers)

  # Convert samplers and domain to equal area CRS
  if (!fm_crs_is_null(domain$crs)) {
    samplers <- fm_transform(samplers, crs = fm_crs("+proj=cea +units=km"))
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

  df <- data.frame(
    weight = ips[, "weight"] * samplers@data[ips$.block, "weight"],
    .block = samplers@data[ips$.block, ".block", drop = TRUE]
  )
  if (is.null(ips$z)) {
    ips <- sp::SpatialPointsDataFrame(ips[, c("x", "y")],
      data = df,
      match.ID = FALSE,
      proj4string = domain_crs
    )
  } else {
    ips <- sp::SpatialPointsDataFrame(ips[, c("x", "y", "z")],
      data = df,
      match.ID = FALSE,
      proj4string = domain_crs
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
  ips
}



## @export
## @describeIn fm_int_inla_mesh `Spatial` integration
# fm_int_inla_mesh.Spatial <- function(samplers,
#                                     domain,
#                                     name = NULL,
#                                     int.args = NULL,
#                                     ...) {
#  samplers <- sf::st_as_sf(samplers)
#
#  ips <-
#    fm_int_inla_mesh(samplers,
#                     domain = domain,
#                     name = name,
#                     int.args = int.args,
#                     ...)
#
#  ips <- as(ips, "Spatial")
#
#  ips
# }

#' @export
#' @describeIn fm_int `inla.mesh` integration. Any sampler class with an
#' associated [fm_int_inla_mesh()] method is supported.
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
    ...
  )

  return(ips)
}
