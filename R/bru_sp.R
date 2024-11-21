# Convert data frame to SpatialLinesDataFrame
#
#
# @keywords internal
# @export
# @param data A SpatialGridDataFrame or SpatialPixelDataFrame
# @param where Spatial*DataFrame of locations for which to fill in values from
# \code{data}. If NULL, use \code{data} to determine the locations.
# @return Spatial object

sfill <- function(data, where = NULL) {
  check_spatstat("spatstat.geom")
  bru_safe_sp(force = TRUE)

  if (is.null(where)) {
    where <- data
  }
  vallist <- list()
  for (k in seq_len(ncol(data@data))) {
    dpoints <- sp::SpatialPoints(data)
    vals <- data@data[, k]
    dpoints <- dpoints[!is.na(vals), ]
    vals <- vals[!is.na(vals)]

    data.ow <- spatstat.geom::owin(
      range(sp::coordinates(dpoints)[, 1]),
      range(sp::coordinates(dpoints)[, 2])
    )
    data.ppp <- spatstat.geom::as.ppp(sp::coordinates(dpoints), data.ow)
    where.ow <- spatstat.geom::owin(
      range(sp::coordinates(where)[, 1]),
      range(sp::coordinates(where)[, 2])
    )
    where.ppp <- spatstat.geom::as.ppp(sp::coordinates(where), where.ow)

    nn <- spatstat.geom::nncross(where.ppp, data.ppp)[, "which"]
    vallist[[k]] <- vals[nn]
  }
  ret <- data.frame(do.call(data.frame, vallist))
  colnames(ret) <- colnames(data@data)
  ret <- sp::SpatialPixelsDataFrame(where, data = ret)
  ret
}


#' Convert data frame to SpatialLinesDataFrame
#'
#' A line in 2D space is defined by a start and an end point, each associated
#' with 2D coordinates. This function takes a `data.frame` as input and assumes
#' that each row defines a line in space. In order to do so, the data frame must
#' have at least four columns and the `start.cols` and `end.cols` parameters
#' must be used to point out the names of the columns that define the start and
#' end coordinates of the line. The data is then converted to a
#' `SpatialLinesDataFrame` `DF`. If a coordinate reference system `crs` is
#' provided it is attached to `DF`. If also `to.crs` is provided, the coordinate
#' system of `DF` is transformed accordingly. Additional columns of the input
#' data, e.g. covariates, are retained and attached to `DF`.
#'
#' @export
#' @param data A data.frame
#' @param start.cols Character array pointing out the columns of `data` that
#'   hold the start points of the lines
#' @param end.cols Character array pointing out the columns of `data` that hold
#'   the end points of the lines
#' @param crs Coordinate reference system of the original `data`
#' @param to.crs Coordinate reference system for the output.
#' @param format Format of the output object. Either `"sp"` (default) or `"sf"`
#' @return `sp::SpatialLinesDataFrame` or [sf::sf]
#' @keywords internal
#'
#' @examples
#' \donttest{
#' # Create a data frame defining three lines
#' lns <- data.frame(
#'   xs = c(1, 2, 3), ys = c(1, 1, 1), # start points
#'   xe = c(2, 3, 4), ye = c(2, 2, 2)
#' ) # end points
#'
#'
#' # Conversion to sf without CRS
#' spl <- sline(lns,
#'   start.cols = c("xs", "ys"),
#'   end.cols = c("xe", "ye"),
#'   format = "sf"
#' )
#'
#' if (require(ggplot2, quietly = TRUE)) {
#'   # Plot the lines
#'   ggplot() +
#'     gg(spl)
#' }
#' }
#'
sline <- function(data, start.cols, end.cols, crs = fm_crs(), to.crs = NULL,
                  format = c("sp", "sf")) {
  format <- match.arg(format)

  sp <- as.data.frame(data[, start.cols])
  ep <- as.data.frame(data[, end.cols])

  colnames(sp) <- c("x", "y")
  colnames(ep) <- c("x", "y")

  if (identical(format, "sp")) {
    bru_safe_sp(force = TRUE)

    lilist <- lapply(seq_len(nrow(sp)), function(k) {
      sp::Lines(list(sp::Line(rbind(sp[k, ], ep[k, ]))), ID = k)
    })
    spl <- sp::SpatialLines(lilist, proj4string = fm_CRS(crs))

    df <- data[, setdiff(names(data), c(start.cols, end.cols))]
    rownames(df) <- seq_len(nrow(df))

    slines <- sp::SpatialLinesDataFrame(spl, data = df)
  } else {
    lin <- lapply(seq_len(nrow(sp)), function(k) {
      sf::st_linestring(as.matrix(rbind(
        sp[k, , drop = FALSE],
        ep[k, , drop = FALSE]
      )))
    })
    lin <- sf::st_sfc(lin, crs = fm_crs(crs))
    df <- data[, setdiff(names(data), c(start.cols, end.cols))]
    rownames(df) <- seq_len(nrow(df))
    slines <- sf::st_sf(geometry = lin, df)
  }

  # If requested, change CRS
  if (!is.null(to.crs)) slines <- fm_transform(slines, to.crs)

  slines
}


#' Convert a data.frame of boundary points into a SpatialPolgonsDataFrame
#'
#' A polygon can be described as a sequence of points defining the polygon's
#' boundary. When given such a sequence (anti clockwise!) this function creates
#' a `SpatialPolygonsDataFrame` or `sf` holding the polygon decribed. By
#' default, the first two columns of `data` are assumed to define the x and y
#' coordinates of the points. This behaviour can be changed using the `cols`
#' parameter, which points out the names of the columns holding the coordinates.
#' The coordinate reference system of the resulting spatial polygon can be set
#' via the `crs` parameter. Posterior conversion to a different CRS is supported
#' using the `to.crs` parameter.
#'
#' @keywords internal
#' @export
#' @param data A data.frame of points describing the boundary of the polygon
#' (unique points, no holes)
#' @param cols Column names of the x and y coordinates within the data
#' @param crs Coordinate reference system of the points
#' @param to.crs Coordinate reference system for the `SpatialLines`/`sf` ouput.
#' @param format Format of the output object. Either "sp" (default) or "sf"
#' @return `sp::SpatialPolygonsDataFrame` or [sf::sf]
#'
#' @examples
#' \donttest{
#' # Create data frame of boundary points (anti clockwise!)
#' pts <- data.frame(
#'   x = c(1, 2, 1.7, 1.3),
#'   y = c(1, 1, 2, 2)
#' )
#'
#' # Convert to sf
#' pol <- spoly(pts, format = "sf")
#'
#' if (require(ggplot2, quietly = TRUE)) {
#'   # Plot it!
#'   ggplot() +
#'     gg(pol)
#' }
#' }
#'
spoly <- function(data,
                  cols = colnames(data)[1:2],
                  crs = fm_crs(),
                  to.crs = NULL,
                  format = c("sp", "sf")) {
  format <- match.arg(format)
  if (identical(format, "sp")) {
    bru_safe_sp(force = TRUE)

    po <- sp::Polygon(data[, cols], hole = FALSE)
    pos <- sp::Polygons(list(po), ID = "tmp")
    predpoly <- sp::SpatialPolygons(list(pos), proj4string = fm_CRS(crs))
    df <- data.frame(weight = 1)
    rownames(df) <- "tmp"
    pol <- sp::SpatialPolygonsDataFrame(predpoly, data = df)
  } else {
    po <- sf::st_polygon(
      list(as.matrix(data[c(seq_len(NROW(data)), 1L), cols]))
    )
    po <- sf::st_sfc(po, crs = fm_crs(crs))
    pol <- sf::st_sf(geometry = po, weight = 1, check_ring_dir = FALSE)
  }

  # If requested, change CRS
  if (!is.null(to.crs)) pol <- fm_transform(pol, to.crs)
  pol
}
