#' @title Groupwise Cross product of integration points
#'
#' @description
#' Calculates the groupwise cross product of integration points in different
#' dimensions and multiplies their weights accordingly, where the group is told
#' by polygon ID. If the object defining points in a particular dimension has no
#' weights attached to it all weights are assumed to be 1.
#'
#' @aliases group_cprod
#' @keywords internal
#'
#'
#' @param ... `data.frame` or `SpatialPointsDataFrame` objects, each one usually obtained by a call to the [ipoints] function.
#' @return A `data.frame` or `SpatialPointsDataFrame` of multidimensional integration points and their weights
#'
#' @examples
#' \donttest{
#' # ipoints needs INLA
#' if (bru_safe_inla()) {
#'   # Create integration points in dimension 'myDim' and 'myDiscreteDim'
#'   ips1 <- ipoints(rbind(c(0, 3), c(3, 8)), 17, name = "myDim")
#'   ips2 <- ipoints(domain = c(1, 2, 4), name = "myDiscreteDim")
#'
#'   # Calculate the cross product
#'   ips <- cprod(ips1, ips2)
#'
#'   # Plot the integration points
#'   plot(ips$myDim, ips$myDiscreteDim, cex = 10 * ips$weight)
#' }
#' }
#'

# TODO For sp, change sp to sf, do the full_join, change sp back. Onion structure.
group_cprod <- function(..., group = NULL) {
  ipl <- list(...)

  # Transfrom sp to sf
  if (any(lapply(ipl, function(x) inherits(x, "Spatial")))) {
    ipl_sp <- TRUE
    ipl <- lapply(ipl, sf::st_as_sf)
  }

  ipl <- ipl[!vapply(ipl, is.null, TRUE)]
  if (length(ipl) == 0) {
    return(NULL)
  }

  if (length(ipl) == 1) {
    ips <- ipl[[1]]
  } else {
    ips1 <- ipl[[1]]
    if (length(ipl) > 2) {
      ips2 <- do.call(group_cprod, ipl[2:length(ipl)])
    } else {
      ips2 <- ipl[[2]]
    }
    if (!"weight" %in% names(ips1)) {
      ips1$weight <- 1
    }
    if (!"weight" %in% names(ips2)) {
      ips2$weight <- 1
    }

    # This line from the cprod handles the grouping already.
    by <- setdiff(intersect(names(ips1), names(ips2)), "weight")

# `sf::st_join` performs spatial join/filter; `dplyr::*_join` expects `x` of
# class `sf` and `y` of class `data.frame`. The trick `as.tibble(sf_obj)` allows
# `dplyr::full_join` and turn it back to `sf` with active geometry as the ips1.
# Z <- full_join(as_tibble(X), as_tibble(Y), by = "group")
# st_as_sf(Z)
# https://stackoverflow.com/questions/64365792/dplyr-full-join-on-geometry-columns-of-sf-objects
    if (inherits(ips1, "sf", "sfc") ||
        inherits(ips2, "sf", "sfc")) {
      ips <-
        sf::st_as_sf(dplyr::full_join(tibble::as_tibble(ips1),
                                      tibble::as_tibble(ips2),
                                      by = by))
    } else {
      ips <-
        base::merge(ips1, ips2, by = by, all = TRUE) # all = TRUE for full_join
    }

    ips$weight <- ips$weight.x * ips$weight.y
    ips[["weight.x"]] <- NULL
    ips[["weight.y"]] <- NULL
    row.names(ips) <- as.character(seq_len(NROW(ips)))
  }

  # Transform back to sp
  if(ipl_sp){
    ips <-  sp::as.Spatial(ips)
  }
  ips
}
