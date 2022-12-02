#' @title Groupwise Cross product of integration points
#'
#' @description
#' Calculates the groupwise cross product of integration points in different
#' dimensions and multiplies their weights accordingly, where the group is told
#' by polygon ID. If the object defining points in a particular dimension has no
#' weights attached to it all weights are assumed to be 1.
#'
#' @aliases group_cprod
#' @export
#'
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
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
#' 20221201 This is the full_join trick I need here
#' Z <- full_join(as_tibble(X), as_tibble(Y), by = "grp")
#' st_as_sf(Z)
group_cprod <- function(..., group = NULL) {
  ipl <- list(...)
  ipl <- ipl[!vapply(ipl, is.null, TRUE)]
  if (length(ipl) == 0) {
    return(NULL)
  }

  if (length(ipl) == 1) {
    ips <- ipl[[1]]
  } else {
    ips1 <- ipl[[1]]
    if (length(ipl) > 2) {
      ips2 <- do.call(cprod, ipl[2:length(ipl)])
    } else {
      ips2 <- ipl[[2]]
    }
    if (!"weight" %in% names(ips1)) {
      ips1$weight <- 1
    }
    if (!"weight" %in% names(ips2)) {
      ips2$weight <- 1
    }

    by <- setdiff(intersect(names(ips1), names(ips2)), "weight")
    if (inherits(ips1, "Spatial")) {
      ips <- sp::merge(ips1, ips2, by = by, duplicateGeoms = TRUE)
    } else if (inherits(ips2, "Spatial")) {
      ips <- sp::merge(ips2, ips1, by = by, duplicateGeoms = TRUE)
    } else {
      ips <- base::merge(ips1, ips2, by = by)
    }
    ips$weight <- ips$weight.x * ips$weight.y
    ips[["weight.x"]] <- NULL
    ips[["weight.y"]] <- NULL
    row.names(ips) <- as.character(seq_len(NROW(ips)))
  }
  ips
}
