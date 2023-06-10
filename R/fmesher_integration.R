#' Split lines at mesh edges
#'
#' @aliases split_lines
#' @export
#' @param mesh An inla.mesh object
#' @param sp Start points of lines
#' @param ep End points of lines
#' @param filter.zero.length Filter out segments with zero length? (Bool)
#' @param ... argments to int.quadrature
#' @return List of start and end points resulting from splitting the given lines
#' @author Fabian E. Bachl \email{f.e.bachl@@bath.ac.uk}
#' @keywords internal
split_lines <- function(mesh, sp, ep, filter.zero.length = TRUE) {
  idx <- seq_len(NROW(sp))
  if (NROW(sp) > 0) {
    # Filter out segments not on the mesh
    t1 <- INLA::inla.fmesher.smorg(
      loc = mesh$loc, tv = mesh$graph$tv,
      points2mesh = as.matrix(data.frame(sp, z = 0))
    )$p2m.t
    t2 <- INLA::inla.fmesher.smorg(
      loc = mesh$loc, tv = mesh$graph$tv,
      points2mesh = as.matrix(data.frame(ep, z = 0))
    )$p2m.t
    # if (any(t1==0) | any(t2==0)) { warning("points outside boundary! filtering...")}
    sp <- sp[!((t1 == 0) | (t2 == 0)), , drop = FALSE]
    ep <- ep[!((t1 == 0) | (t2 == 0)), , drop = FALSE]
    idx <- idx[!((t1 == 0) | (t2 == 0))]
  }

  if (NROW(sp) == 0) {
    return(list(
      sp = sp, ep = ep,
      split.origin = NULL,
      idx = idx,
      split.loc = NULL
    ))
  }

  loc <- as.matrix(rbind(sp, ep))

  # Split the segments into parts
  if (NCOL(loc) == 2) {
    loc <- cbind(loc, rep(0, NROW(loc)))
  }
  np <- dim(sp)[1]
  sp.idx <- t(rbind(seq_len(np), np + seq_len(np)))
  splt <- INLA::inla.fmesher.smorg(
    mesh$loc, mesh$graph$tv,
    splitlines = list(loc = loc, idx = sp.idx)
  )
  # plot(data$mesh)
  # points(loc)
  # points(splt$split.loc,col="blue)

  # Start points of new segments
  sp <- splt$split.loc[splt$split.idx[, 1], seq_len(dim(sp)[2]), drop = FALSE]
  # End points of new segments
  ep <- splt$split.loc[splt$split.idx[, 2], seq_len(dim(ep)[2]), drop = FALSE]
  idx <- idx[splt$split.idx[, 1]]
  origin <- splt$split.origin

  # Filter out zero length segments
  if (filter.zero.length) {
    sl <- apply((ep - sp)^2, MARGIN = 1, sum)
    sp <- sp[!(sl == 0), , drop = FALSE]
    ep <- ep[!(sl == 0), , drop = FALSE]
    origin <- origin[!(sl == 0)]
    idx <- idx[!(sl == 0)]
  }

  return(list(
    sp = sp, ep = ep,
    split.origin = origin,
    idx = idx,
    split.loc = splt$split.loc
  ))
}


#' @title (Blockwise) cross product of integration points
#'
#' @description
#' Calculates the groupwise cross product of integration points in different
#' dimensions and multiplies their weights accordingly.
#' If the object defining points in a particular dimension has no
#' weights attached to it all weights are assumed to be 1.
#'
#'
#' @export
#' @keywords internal
#'
#' @param ... `data.frame`, `sf`, or `SpatialPointsDataFrame` objects, each one
#' usually obtained by a call to an [fm_int()] method.
#' @param na.rm logical; if `TRUE`, the rows with weight `NA` from the
#' non-overlapping full_join will be removed; if `FALSE`, set the undefined weights to `NA`.
#' If `NULL` (default), act as `TRUE`, but warn if any elements needed removing.
#' @param .blockwise logical; if `FALSE`, computes full tensor product integration.
#' If `TRUE`, computes within-block tensor product integration (used internally
#' by [fm_int()]).
#' Default `FALSE`
#' @return A `data.frame`, `sf`, or `SpatialPointsDataFrame` of multidimensional
#' integration points and their weights
#'
#' @examples
#' \donttest{
#' # fm_int needs INLA
#' if (bru_safe_inla()) {
#'   # Create integration points in dimension 'myDim' and 'myDiscreteDim'
#'   ips1 <- fm_int(INLA::inla.mesh.1d(1:20),
#'     rbind(c(0, 3), c(3, 8)),
#'     name = "myDim"
#'   )
#'   ips2 <- fm_int(domain = c(1, 2, 4), name = "myDiscreteDim")
#'
#'   # Calculate the cross product
#'   ips <- fm_cprod(ips1, ips2)
#'
#'   # Plot the integration points
#'   plot(ips$myDim, ips$myDiscreteDim, cex = 10 * ips$weight)
#' }
#' }
#'
#' @importFrom stats na.omit
fm_cprod <- function(..., na.rm = NULL, .blockwise = FALSE) {
  ipl <- list(...)

  # Transform sp to sf
  # TODO make a test. and give a warning for NA non-overlapping outcome?
  # check for each element, or on the subset, change only for sp anonymous function on lapply
  ipl_sp <- vapply(ipl, function(x) inherits(x, "Spatial"), TRUE)
  ipl[ipl_sp] <- lapply(ipl[ipl_sp], sf::st_as_sf)

  ipl <- ipl[!vapply(ipl, is.null, TRUE)]
  if (length(ipl) == 0) {
    return(NULL)
  }

  if (length(ipl) == 1) {
    ips <- ipl[[1]]
  } else {
    ips1 <- ipl[[1]]
    if (length(ipl) > 2) {
      ips2 <- do.call(fm_cprod, ipl[2:length(ipl)])
    } else {
      ips2 <- ipl[[2]]
    }
    if (!"weight" %in% names(ips1)) {
      ips1$weight <- 1
    }
    if (!"weight" %in% names(ips2)) {
      ips2$weight <- 1
    }
    if (!".block" %in% names(ips1)) {
      ips1$.block <- seq_len(NROW(ips1))
    }
    if (!".block" %in% names(ips2)) {
      ips2$.block <- seq_len(NROW(ips2))
    }

    by <- setdiff(intersect(names(ips1), names(ips2)), "weight")
    if (!.blockwise) {
      by <- setdiff(by, ".block")
    }

    # `sf::st_join` performs spatial join/filter; `dplyr::*_join` expects `x` of
    # class `sf` and `y` of class `data.frame`. The trick `as.tibble(sf_obj)` allows
    # `dplyr::full_join` and turn it back to `sf` with active geometry as the ips1.
    # Z <- full_join(as_tibble(X), as_tibble(Y), by = "group")
    # st_as_sf(Z)
    # https://stackoverflow.com/questions/64365792/dplyr-full-join-on-geometry-columns-of-sf-objects
    if (inherits(ips1, c("sf", "sfc")) ||
      inherits(ips2, c("sf", "sfc"))) {
      if (length(by) == 0) {
        ips <-
          sf::st_as_sf(dplyr::cross_join(
            tibble::as_tibble(ips1),
            tibble::as_tibble(ips2)
          ))
      } else {
        ips <-
          sf::st_as_sf(dplyr::full_join(tibble::as_tibble(ips1),
            tibble::as_tibble(ips2),
            by = by,
            relationship = "many-to-many"
          ))
      }
    } else {
      # equivalent to base::merge(ips1, ips2, by = by, all = TRUE)
      if (length(by) == 0) {
        ips <-
          dplyr::cross_join(ips1, ips2)
      } else {
        ips <-
          dplyr::full_join(ips1, ips2,
            by = by,
            relationship = "many-to-many"
          )
      }
    }

    ips$weight <- ips$weight.x * ips$weight.y
    ips[["weight.x"]] <- NULL
    ips[["weight.y"]] <- NULL
    tibble::remove_rownames(ips)

    if (!.blockwise) {
      ips$.block <- paste0(ips$.block.x, ",", ips$.block.y)
      ips[[".block.x"]] <- NULL
      ips[[".block.y"]] <- NULL
    }
  }
  if (any(is.na(ips$weight)) && !isFALSE(na.rm)) {
    if (is.null(na.rm)) {
      warning(
        paste0(
          "Block information mismatch resulting in NA weights, and 'na.rm' was not supplied.",
          " These rows will be removed."
        )
      )
    }
    ips <- na.omit(ips)
  }

  # TODO Transform back to sp only if they are required. ips is a tibble sf tbl data.frame.
  # It does not make sense to revert certain indices back after merging. Hence, I revert the entire object back to sp.
  if (any(ipl_sp)) {
    ips <- sf::as_Spatial(ips)
    # WARNING SHOULD BE HERE, MORE FRIENDLY ERROR, METHOD ST_AS_SF
    # TODO deprecated_soft warning for `sp` presence
    # lifecycle::deprecate_soft(
    #   when = "2.7.0",
    #   what = "inlabru::ipl_sp()", # ipl_sp() does not exist but it does not allow mentioning deprecated function in another package
    #   details = "Support for `sp` is gradually deprecated in favour of `sf` and `terra`"
    # )
  }
  ips
}


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
#' if (bru_safe_inla() && bru_safe_sp()) {
#'   # Integration on the interval (2, 3.5) with Simpson's rule
#'   ips <- fm_int(INLA::inla.mesh.1d(0:4), samplers = cbind(2, 3.5))
#'   plot(ips)
#'
#'   # Create integration points for the two intervals [0,3] and [5,10]
#'
#'   ips <- fm_int(
#'     INLA::inla.mesh.1d(0:10),
#'     matrix(c(0, 3, 5, 10), nrow = 2, byrow = TRUE)
#'   )
#'   plot(ips)
#'
#'   # Convert a 1D mesh into integration points
#'   mesh <- INLA::inla.mesh.1d(seq(0, 10, by = 1))
#'   ips <- fm_int(mesh, name = "time")
#'   plot(ips)
#'
#'
#'   if (require("ggplot2", quietly = TRUE)) {
#'     data("gorillas", package = "inlabru")
#'     #' Integrate on a 2D mesh with polygon boundary subset
#'     ips <- fm_int(gorillas$mesh, gorillas$boundary)
#'     ggplot() +
#'       gg(gorillas$mesh) +
#'       gg(gorillas$boundary) +
#'       gg(ips, aes(size = weight)) +
#'       scale_size_area()
#'   }
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
  ips <- do.call(fm_cprod, c(ips_list, list(.blockwise = TRUE)))

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

  ips <- do.call(fm_cprod, c(
    lips_samplers,
    lips_full_domain_samplers,
    list(.blockwise = FALSE)
  ))

  if (any(sp_samplers) && !any(sf_samplers)) {
    ips <- sf::as_Spatial(ips)
    cnames <- sp::coordnames(ips)
    sp::coordnames(ips) <- c("x", "y", "z")[seq_along(cnames)]
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

#' @param int.args List of arguments passed to line and integration methods.
#' * `method`: "stable" (to aggregate integration weights onto mesh nodes)
#'   or "direct" (to construct a within triangle/segment integration scheme
#'   without aggregating onto mesh nodes)
#' * `nsub1`, `nsub2`: integers controlling the number of internal integration
#'   points before aggregation. Points per triangle: `(nsub2+1)^2`.
#'   Points per knot segment: `nsub1`
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

#' @title Project integration points to mesh vertices
#'
#' @description
#' Compute information for assigning points to the vertices of the covering triangle
#'
#' @export
#' @param points A `SpatialPointsDataFrame`, `sf`, or `list` object
#' @param mesh An `inla.mesh` object
#' @return `SpatialPointsDataFrame`, `sf`, or `list` of mesh vertices with
#' projected data attached
#' @importFrom rlang .data
#' @keywords internal
#' @export

fm_vertex_projection <- function(points, mesh) {
  if (inherits(points, "sf") ||
    inherits(points, "Spatial")) {
    n_points <- NROW(points)
    res <- fm_evaluator(mesh, points)
  } else {
    n_points <- NROW(points$loc)
    res <- fm_evaluator(mesh, points$loc)
  }
  tri <- res$proj$t
  bary <- res$proj$bary

  if (is.null(points$weight)) {
    points$weight <- rep(1L, n_points)
  }
  if (is.null(points$.block)) {
    points$.block <- rep(1L, n_points)
  }

  ok <- !is.na(tri)
  ok[ok] <- (tri[ok] > 0)
  if (any(!ok)) {
    warning("Some integration points were outside the mesh; check your coordinate systems.")
  }

  data <-
    data.frame(
      .vertex = as.vector(mesh$graph$tv[tri[ok], ]),
      weight = as.vector(points$weight[ok] * bary[ok, ]),
      .block = rep(points$.block[ok], times = 3)
    )

  data <-
    dplyr::summarise(
      dplyr::group_by(data, .data$.vertex, .data$.block),
      weight = sum(.data$weight),
      .groups = "drop"
    )
  coords <- mesh$loc[data$.vertex, , drop = FALSE]
  data <- dplyr::select(data, c("weight", ".block", ".vertex"))

  if (inherits(points, "Spatial")) {
    ret <- sp::SpatialPointsDataFrame(
      coords[, seq_along(sp::coordnames(points)), drop = FALSE],
      proj4string = fm_CRS(mesh),
      data = data,
      match.ID = FALSE
    )
    sp::coordnames(ret) <- sp::coordnames(points)
  } else if (inherits(points, "sf")) {
    colnames(coords) <- c("X", "Y", "Z")[seq_len(ncol(coords))]
    d <- length(intersect(colnames(sf::st_coordinates(points)), c("X", "Y", "Z")))
    data <- cbind(
      tibble::as_tibble(coords[, seq_len(d), drop = FALSE]),
      tibble::as_tibble(data)
    )
    ret <- sf::st_as_sf(
      data,
      coords = seq_len(d),
      crs = fm_crs(mesh)
    )
  } else {
    colnames(coords) <- c("X", "Y", "Z")[seq_len(ncol(coords))]
    data$loc <- coords
    ret <- data
  }

  ret
}

#' Subset integration on a mesh
#'
#' Integration methods for spatial samplers on `inla.mesh` meshes.
#'
#' @returns An `sf` point object with columns `weight` and `.block`
#' @inheritParams fm_int
#' @export
#' @keywords internal
fm_int_inla_mesh <- function(samplers,
                             domain,
                             name = NULL,
                             int.args = NULL,
                             ...) {
  stopifnot(inherits(domain, "inla.mesh"))

  if (missing(samplers) || is.null(samplers)) {
    return(
      fm_int_inla_mesh_NULL(
        samplers = NULL,
        domain = domain,
        name = name,
        int.args = int.args,
        ...
      )
    )
  }

  UseMethod("fm_int_inla_mesh")
}

#' @describeIn fm_int_inla_mesh Full domain integration
fm_int_inla_mesh_NULL <- function(samplers,
                                  domain,
                                  name = NULL,
                                  int.args = NULL,
                                  ...) {
  stopifnot(is.null(samplers))

  ips <- fm_int_inla_mesh_polygon(
    domain = domain,
    samplers = NULL,
    int.args = int.args
  )

  if (!is.null(name) && (name != attr(ips, "sf_column"))) {
    ips <- dplyr::rename(ips, "{name}" := "geometry")
  }

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

  fm_int_inla_mesh(
    sf::st_geometry(samplers),
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
  if (is.null(name)) {
    name <- "geometry"
  }
  ips <- tibble::tibble(
    "{name}" := samplers,
    weight = .weight,
    .block = seq_len(NROW(samplers))
  )
  ips <- sf::st_as_sf(ips, sf_column_name = name)

  # TODO: remove points outside the domain

  ips
}


#' @export
#' @describeIn fm_int_inla_mesh `sfc_MULTIPOINT` integration
#' @importFrom rlang :=
fm_int_inla_mesh.sfc_MULTIPOINT <- function(samplers,
                                            domain,
                                            name = NULL,
                                            int.args = NULL,
                                            .weight = rep(1, NROW(samplers)),
                                            ...) {
  coords <- tibble::as_tibble(sf::st_coordinates(samplers))
  coords <- dplyr::rename(coords, .block = "L1")
  coords$weight <- .weight[coords$.block]
  ips <- sf::st_as_sf(
    coords,
    coords = intersect(names(coords), c("X", "Y", "Z", "M")),
    crs = fm_crs(samplers)
  )

  # TODO: remove points outside the domain

  if (!is.null(name) && (name != attr(ips, "sf_column"))) {
    ips <- dplyr::rename(ips, "{name}" := "geometry")
  }
  ips
}









fm_int_inla_mesh_lines <- function(samplers,
                                   domain,
                                   name = NULL,
                                   int.args = NULL,
                                   .weight = rep(1, NROW(samplers)),
                                   ...) {
  project <- identical(int.args$method, "stable")

  weight <- .weight
  .block <- seq_len(NROW(samplers))

  # Extract start and end coordinates
  coords <- sf::st_coordinates(samplers)
  if ("L2" %in% colnames(coords)) {
    # MULTILINESTRING
    feature <- coords[, "L2"]
    part <- coords[, "L1"]
    L_idx <- which(colnames(coords) %in% c("L1", "L2"))
  } else {
    feature <- coords[, "L1"]
    part <- rep(1, nrow(coords))
    L_idx <- which(colnames(coords) %in% "L1")
  }

  segment <- which((diff(part) > 0) | (diff(feature) > 0))
  coordnames <- intersect(colnames(coords), c("X", "Y", "Z", "M"))

  sp <- coords[-c(segment, nrow(coords)), coordnames, drop = FALSE]
  ep <- coords[-c(1L, 1L + segment), coordnames, drop = FALSE]
  idx <- feature[-c(segment, nrow(coords))]

  sampler_crs <- fm_crs(samplers)
  target_crs <- fm_crs(domain)
  if (!fm_crs_is_null(sampler_crs) &&
    fm_crs_is_null(target_crs)) {
    target_crs <- sampler_crs
  }

  # Filter out points outside the mesh...
  sp <- fm_transform(sp, crs = target_crs, crs0 = sampler_crs, passthrough = TRUE)
  ep <- fm_transform(ep, crs = target_crs, crs0 = sampler_crs, passthrough = TRUE)
  proj1 <- fm_evaluator(domain, loc = sp, crs = target_crs)
  proj2 <- fm_evaluator(domain, loc = ep, crs = target_crs)
  ok <- (proj1$proj$ok & proj2$proj$ok)
  if (!all(ok)) {
    warning("Found spatial lines with start or end point ouside of the mesh. Omitting.")
  }
  sp <- sp[ok, , drop = FALSE]
  ep <- ep[ok, , drop = FALSE]
  idx <- idx[ok]

  # Split at mesh edges
  line.spl <- split_lines(domain, sp, ep, TRUE)
  sp <- line.spl$sp
  ep <- line.spl$ep
  idx <- idx[line.spl$split.origin]

  # At this point, sp and ep are in the target_crs

  # Determine integration points along lines

  if (fm_crs_is_null(sampler_crs)) {
    ips <- (sp + ep) / 2
    w <- rowSums((ep - sp)^2)^0.5
  } else {
    # Has CRS
    longlat.crs <- fm_crs("longlat_globe")
    geocentric.crs <- fm_crs("sphere")
    sp3d <- fm_transform(sp, crs = geocentric.crs, crs0 = target_crs)
    ep3d <- fm_transform(ep, crs = geocentric.crs, crs0 = target_crs)
    mp3d <- (sp3d + ep3d) / rowSums((sp3d + ep3d)^2)^0.5

    ips <- fm_transform(mp3d, crs = target_crs, crs0 = geocentric.crs)
    w <- sp::spDists(
      fm_transform(sp3d, crs = longlat.crs, crs0 = geocentric.crs)[, 1:2, drop = FALSE],
      fm_transform(ep3d, crs = longlat.crs, crs0 = geocentric.crs)[, 1:2, drop = FALSE],
      diagonal = TRUE, longlat = TRUE
    )
  }

  # Wrap everything up and perform projection according to distance and given .block argument
  ips <- data.frame(ips)
  d_ips <- ncol(ips)
  # Temporary names
  colnames(ips) <- c("x", "y", "z")[seq_len(d_ips)]

  # Weights
  ips <- cbind(ips, weight = w)
  ips$weight <- ips$weight * weight[idx]
  ips$.block <- .block[idx]

  ips <- sf::st_as_sf(as.data.frame(ips),
    coords = seq_len(d_ips),
    crs = target_crs
  )

  # Project to mesh vertices
  if (project) {
    ips <- fm_vertex_projection(ips, domain)
  }

  if (!is.null(name) && (name != attr(ips, "sf_column"))) {
    ips <- dplyr::rename(ips, "{name}" := "geometry")
  }

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
  ips <- fm_int_inla_mesh_lines(samplers, domain, name, int.args, .weight, ...)

  if (!is.null(name) && (name != attr(ips, "sf_column"))) {
    ips <- dplyr::rename(ips, "{name}" := "geometry")
  }

  ips
}

#' @export
#' @describeIn fm_int_inla_mesh `sfc_MULTILINESTRING` integration
fm_int_inla_mesh.sfc_MULTILINESTRING <- function(samplers,
                                                 domain,
                                                 name = NULL,
                                                 int.args = NULL,
                                                 .weight = rep(1, NROW(samplers)),
                                                 ...) {
  ips <- fm_int_inla_mesh_lines(samplers, domain, name, int.args, .weight, ...)

  if (!is.null(name) && (name != attr(ips, "sf_column"))) {
    ips <- dplyr::rename(ips, "{name}" := "geometry")
  }

  ips
}



#' Integration scheme for mesh triangle interiors
#'
#' @param mesh Mesh on which to integrate
#' @param tri_subset Optional triangle index vector for integration on a subset
#' of the mesh triangles (Default `NULL`)
#' @param nsub number of subdivision points along each triangle edge, giving
#'    `(nsub + 1)^2` proto-integration points used to compute
#'   the vertex weights
#'   (default `NULL=9`, giving 100 integration points for each triangle)
#' @return `list` with elements `loc` and `weight` with
#'   integration points for the mesh
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @keywords internal
fm_int_inla_mesh_core <- function(mesh, tri_subset = NULL, nsub = NULL) {
  # Construct a barycentric grid of subdivision triangle midpoints
  if (is.null(nsub)) {
    nsub <- 9
  }
  stopifnot(nsub >= 0)
  nB <- (nsub + 1)^2

  nT <- nrow(mesh$graph$tv)
  if (is.null(tri_subset)) {
    tri_subset <- seq_len(nT)
  }

  is_spherical <- identical(mesh$manifold, "S2")

  # Barycentric integration coordinates
  b <- seq(1 / 3, 1 / 3 + nsub, length = nsub + 1) / (nsub + 1)
  bb <- as.matrix(expand.grid(b, b))
  # Points above the diagonal should be reflected into the lower triangle:
  refl <- rowSums(bb) > 1
  if (any(refl)) {
    bb[refl, ] <- cbind(1 - bb[refl, 2], 1 - bb[refl, 1])
  }
  # Construct complete barycentric coordinates:
  barycentric_grid <- cbind(1 - rowSums(bb), bb)

  # Construct integration points
  loc <- matrix(0.0, length(tri_subset) * nB, ncol(mesh$loc))
  idx_end <- 0
  for (tri in tri_subset) {
    idx_start <- idx_end + 1
    idx_end <- idx_start + nB - 1
    loc[seq(idx_start, idx_end, length = nB), ] <-
      as.matrix(barycentric_grid %*%
        mesh$loc[mesh$graph$tv[tri, ], , drop = FALSE])
  }

  if (is_spherical) {
    # Normalise
    radius <- sum(mesh$loc[1, ]^2)^0.5
    mesh$loc <- mesh$loc / radius
    loc <- loc / rowSums(loc^2)^0.5
  }

  # Construct integration weights
  tri_area <- INLA::inla.mesh.fem(mesh, order = 1)$ta[tri_subset]

  if (is_spherical) {
    tri_area <- tri_area * radius^2
    loc <- loc * radius
  }

  list(
    loc = loc,
    weight = rep(tri_area / nB, each = nB)
  )
}


fm_int_inla_mesh_polygon <- function(samplers,
                                     domain,
                                     name = NULL,
                                     int.args = NULL,
                                     .weight = rep(1, NROW(samplers)),
                                     ...) {
  method <- match.arg(int.args[["method"]], c("stable", "direct"))

  ipsl <- list()

  # Compute direct integration points
  # TODO: Allow blockwise construction to avoid
  # overly large temporary coordinate matrices (via tri_subset)
  integ <- fm_int_inla_mesh_core(domain, nsub = int.args[["nsub2"]])

  # Keep points with positive weights (This should be all,
  # but if there's a degenerate triangle, this gets rid of it)
  ok <- (integ$weight > 0)
  integ$loc <- integ$loc[ok, , drop = FALSE]
  integ$weight <- integ$weight[ok]

  domain_crs <- fm_crs(domain)

  if (!is.null(samplers)) {
    samplers_crs <- fm_crs(samplers)
    integ_sf <- sf::st_as_sf(as.data.frame(integ$loc),
      coords = seq_len(ncol(integ$loc)),
      crs = domain_crs
    )
    if (!identical(domain_crs, samplers_crs) &&
      !fm_crs_is_null(domain_crs) &&
      !fm_crs_is_null(samplers_crs)) {
      integ_sf <- fm_transform(integ_sf,
        crs = samplers_crs,
        passthrough = TRUE
      )
    }

    idx <- sf::st_contains(samplers, integ_sf, sparse = TRUE)

    for (g in seq_along(idx)) {
      if (length(idx[[g]]) > 0) {
        integ_ <- list(
          loc = integ$loc[idx[[g]], , drop = FALSE],
          weight = integ$weight[idx[[g]]]
        )

        if (method %in% c("stable")) {
          # Project integration points and weights to mesh nodes
          integ_ <- fm_vertex_projection(integ_, domain)
        }

        if (ncol(integ_$loc) > 2) {
          ips <- sf::st_as_sf(
            tibble::tibble(
              x = integ_$loc[, 1],
              y = integ_$loc[, 2],
              z = integ_$loc[, 3],
              weight = integ_$weight,
              .block = g
            ),
            coords = c("x", "y", "z"),
            crs = domain_crs
          )
        } else {
          ips <- sf::st_as_sf(
            tibble::tibble(
              x = integ_$loc[, 1],
              y = integ_$loc[, 2],
              weight = integ_$weight,
              .block = g
            ),
            coords = c("x", "y"),
            crs = domain_crs
          )
        }

        ipsl <- c(ipsl, list(ips))
      }
    }
  } else {
    if (method %in% c("stable")) {
      # Project integration points and weights to mesh nodes
      integ <- fm_vertex_projection(integ, domain)
    }

    if (ncol(integ$loc) > 2) {
      ipsl <- list(sf::st_as_sf(
        tibble::tibble(
          x = integ$loc[, 1],
          y = integ$loc[, 2],
          z = integ$loc[, 3],
          weight = integ$weight,
          .block = 1L
        ),
        coords = c("x", "y", "z"),
        crs = domain_crs
      ))
    } else {
      ipsl <- list(sf::st_as_sf(
        tibble::tibble(
          x = integ$loc[, 1],
          y = integ$loc[, 2],
          weight = integ$weight,
          .block = 1L
        ),
        coords = c("x", "y"),
        crs = domain_crs
      ))
    }
  }

  ips <- do.call(rbind, ipsl)

  if (!is.null(name) && (name != attr(ips, "sf_column"))) {
    ips <- dplyr::rename(ips, "{name}" := "geometry")
  }

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
  weight <- .weight
  .block <- seq_len(NROW(samplers))

  ips <- fm_int_inla_mesh_polygon(
    domain = domain,
    int.args = int.args,
    samplers = samplers
  )

  ips$weight <- ips$weight * .weight[ips$.block]
  ips$.block <- .block[ips$.block]

  if (!is.null(name) && (name != attr(ips, "sf_column"))) {
    ips <- dplyr::rename(ips, "{name}" := "geometry")
  }

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
  weight <- .weight
  .block <- seq_len(NROW(samplers))

  ips <- fm_int_inla_mesh_polygon(
    domain = domain,
    int.args = int.args,
    samplers = samplers
  )

  ips$weight <- ips$weight * .weight[ips$.block]
  ips$.block <- .block[ips$.block]

  if (!is.null(name) && (name != attr(ips, "sf_column"))) {
    ips <- dplyr::rename(ips, "{name}" := "geometry")
  }

  ips
}



#' @export
#' @describeIn fm_int_inla_mesh `sfc_GEOMERY` integration
fm_int_inla_mesh.sfc_GEOMETRY <- function(samplers,
                                          domain,
                                          name = NULL,
                                          int.args = NULL,
                                          .weight = rep(1, NROW(samplers)),
                                          ...) {
  geometry_class <- vapply(
    seq_along(samplers),
    function(x) {
      class(samplers[x])[1]
    },
    ""
  )
  .block <- seq_len(NROW(samplers))

  ips <- list()
  for (g_class in unique(geometry_class)) {
    subset <- geometry_class == g_class
    ips[[g_class]] <-
      fm_int_inla_mesh(samplers[subset],
        domain = domain,
        name = name,
        int.args = int.args,
        .weight = .weight[subset]
      )
    ips[[g_class]][[".block"]] <- .block[subset][ips[[g_class]][[".block"]]]
  }
  ips <- do.call(dplyr::bind_rows, ips)

  if (!is.null(name) && (name != attr(ips, "sf_column"))) {
    ips <- dplyr::rename(ips, "{name}" := "geometry")
  }

  ips
}









#' @export
#' @describeIn fm_int_inla_mesh `Spatial` integration
fm_int_inla_mesh.Spatial <- function(samplers,
                                     domain,
                                     name = NULL,
                                     int.args = NULL,
                                     format = NULL,
                                     ...) {
  samplers <- sf::st_as_sf(samplers)

  ips <-
    fm_int_inla_mesh(
      samplers,
      domain = domain,
      name = name,
      int.args = int.args,
      ...
    )

  ips
}

#' @export
#' @describeIn fm_int `inla.mesh` integration. Any sampler class with an
#' associated [fm_int_inla_mesh()] method is supported.
#' @param format character; determines the output format, as either "sf"
#'   (default when the sampler is `NULL`) or "sp". When `NULL`, determined by
#'   the sampler type.
fm_int.inla.mesh <- function(domain,
                             samplers = NULL,
                             name = NULL,
                             int.args = NULL,
                             format = NULL,
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

  if (is.null(format) && inherits(samplers, "Spatial")) {
    format <- "sp"
  }
  if (!is.null(format)) {
    if ((format == "sf") && !inherits(ips, "sf")) {
      ips <- sf::st_as_sf(ips)
      if (!is.null(name) && (name != attr(ips, "sf_column"))) {
        ips <- dplyr::rename(ips, "{name}" := attr(ips, "sf_column"))
      }
    } else if ((format == "sp") && !inherits(ips, "Spatial")) {
      ips <- as(ips, "Spatial")
      cnames <- sp::coordnames(ips)
      sp::coordnames(ips) <- c("x", "y", "z")[seq_along(cnames)]
    }
  }

  ips
}
