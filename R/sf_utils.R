# helper functions for working with sf objects

#' Calculate signed area for polygon
#'
#' @aliases st_signed_area
#' @export
#' @param sfg A POLYGON sfg object
#' @return Returns the signed area.  Negative values indicate
#' anti-clockwise winding direction.
#' @author Andrew Seaton \email{Andrew.Seaton.2@@glasgow.ac.uk}
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}

st_signed_area <- function(sfg) {
  if (!inherits(sfg, c("POLYGON", "sfg"))) {
    stop("Signed area only implemented for POLYGON sfg objects")
  }

  coords <- as.matrix(sfg)
  i <- seq_len(nrow(coords) - 1)
  edges <- cbind(coords[i, , drop = FALSE], coords[i + 1, , drop = FALSE])
  area <- sum((edges[, 3] - edges[, 1]) * (edges[, 2] + edges[, 4]) / 2)
  return(area)
}

#' Check for "XYZ", "XYM" and "XYZM" sfg classes
#'
#' @aliases st_check_dim
#' @param sfc An sfc object
#' @return LOGICAL indicating if any sfg element of the sfc object has class
#' "XYZ", "XYM" or "XYZM". Internal function used to check for 3 and 4
#' dimensional objects.

st_check_dim <- function(sfc) {
  check <- vapply(sfc,
    FUN = function(x) inherits(x, c("XYZ", "XYM", "XYZM")),
    TRUE
  )

  sum(check)
}

#' Check sfg polygon satisfies standards for POLYGON simple features
#'
#' @description It seems as though st_polygon does not check this.
#' For now only implements a basic check for disjoint regions using st_within()
#'
#' @aliases st_check_polygon
#' @export
#' @param sfg A POLYGON sfg object
#' @return LOGICAL; `TRUE` if the `sfg` holes are entirely inside the outer ring, and
#' are disjoint, otherwise `FALSE`. When `FALSE`, the attribute `Message` is set
#' to a character vector describing the detected reasons.

st_check_polygon <- function(sfg) {
  if (!inherits(sfg, c("POLYGON", "sfg"))) {
    stop(
      "Requires POLYGON sfg object"
    )
  }

  np <- length(sfg)

  # 1st is outer boundary
  main <- sf::st_polygon(list(as.matrix(sfg[[1]])))
  ok <- TRUE
  msg <- NULL

  # Rest should be holes
  if (length(sfg) > 1) {
    holes <- lapply(
      sfg[-1],
      FUN = as.matrix
    )
    holes <- lapply(
      holes,
      FUN = function(x) sf::st_geometry(sf::st_polygon(list(x)))
    )
    holes <- do.call(
      c, holes
    )

    check_within <- sf::st_within(
      holes, main,
      sparse = FALSE
    )
    if (!all(check_within)) {
      ok <- FALSE
      msg <- c(msg, "Holes not entirely within the outer ring")
    }

    if (length(holes) > 1) {
      check_overlap <- vapply(
        seq_along(holes),
        function(k) {
          !any(sf::st_intersects(
            holes[-k], holes[k],
            sparse = FALSE
          ))
        },
        TRUE
      )
      if (!all(check_overlap)) {
        ok <- FALSE
        msg <- c(msg, "Some holes overlap")
      }
    }
  }

  attr(ok, "Message") <- msg
  ok
}
