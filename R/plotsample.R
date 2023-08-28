# @title Make rectangular \code{SpatialPolygonsDataFrame}.
#
# @description
# Makes a rectangular \code{SpatialPolygonsDataFrame} using its bottom left coordinates,
# its width and its height.
#
# @param start Bottom left coordinates (x,y).
# @param width Width of base of the rectangle.
# @param height Height of the rectangle.
#
# @return A \code{SpatialPolygonsDataFrame}.
#
# @export

makepoly <- function(start, width, height) {
  poly <- matrix(
    c(
      start[1], start[2],
      start[1], start[2] + height,
      start[1] + width, start[2] + height,
      start[1] + width, start[2],
      start[1], start[2]
    ),
    ncol = 2, byrow = TRUE
  )
  return(Polygon(poly))
}

#' @title Create a plot sample.
#'
#' @description
#' Creates a plot sample on a regular grid with a random start location.
#'
#' @param spdf A `SpatialPointsDataFrame` defining the points that are to be sampled
#' by the plot sample.
#' @param boundary A `SpatialPolygonsDataFrame` defining the survey boundary within which
#' the  points occur.
#' @param x.ppn The proportion of the x=axis that is to be included in the plots.
#' @param y.ppn The proportion of the y=axis that is to be included in the plots.
#' @param nx The number of plots in the x-dimension.
#' @param ny The number of plots in the y-dimension.
#'
#' @return A list with three components:
#'  \describe{
#'    \item{`plots`:}{ A `SpatialPolygonsDataFrame` object containing the plots that were
#'    sampled.}
#'    \item{`dets`:}{ A `SpatialPointsDataFrame` object containing the locations of the
#'    points within the plots.}
#'    \item{`counts`:}{ A `dataframe` containing the following columns
#'    \describe{
#'    \item{`x`:}{The x-coordinates of the centres of the plots within the boundary.}
#'    \item{`y`:}{The y-coordinates of the centres of the plots within the boundary.}
#'    \item{`n`:}{The numbers of points in each plot.}
#'    \item{`area`:}{The areas of the plots within the boundary}
#'    }}
#'  }.
#'
#' @examples
#' \donttest{
#' # Some features require the raster package
#' if (bru_safe_sp() &&
#'   require("sp") &&
#'   require("raster", quietly = TRUE) &&
#'   require("ggplot2", quietly = TRUE)) {
#'   data(gorillas, package = "inlabru")
#'   plotpts <- plotsample(gorillas$nests, gorillas$boundary,
#'     x.ppn = 0.4, y.ppn = 0.4, nx = 5, ny = 5
#'   )
#'   ggplot() +
#'     gg(plotpts$plots) +
#'     gg(plotpts$dets, pch = "+", cex = 2) +
#'     gg(gorillas$boundary)
#' }
#' }
#'
#' @export
#'
plotsample <- function(spdf, boundary, x.ppn = 0.25, y.ppn = 0.25, nx = 5, ny = 5) {
  if (x.ppn <= 0 || x.ppn >= 1) stop("'x.ppn' must greater than 0 and less than 1")
  if (y.ppn <= 0 || y.ppn >= 1) stop("'y.ppn' must greater than 0 and less than 1")


  srange <- raster::extent(boundary)
  xrange <- srange[1:2]
  yrange <- srange[3:4]
  nxtot <- round(nx / x.ppn)
  nytot <- round(ny / y.ppn)
  width <- diff(xrange) / nxtot
  height <- diff(yrange) / nytot
  dx <- diff(xrange) / nx
  dy <- diff(yrange) / ny
  startx <- runif(1, xrange[1] - 0.99999 * width, xrange[1] + dx - 0.99999 * width)
  starty <- runif(1, yrange[1] - 0.99999 * height, yrange[1] + dy - 0.99999 * height)
  xs <- startx + (0:nx) * dx
  ys <- starty + (0:ny) * dy

  nxs <- length(xs)
  nys <- length(ys)
  starts <- data.frame(x = rep(xs, nys), y = rep(ys, rep(nxs, nys)))
  nplots <- dim(starts)[1]

  polys <- vector("list", nplots)
  for (i in 1:nplots) {
    polys[[i]] <- Polygons(list(makepoly(as.numeric(starts[i, ]), width, height)), i)
  }
  plots <- SpatialPolygons(polys, proj4string = fm_CRS(spdf))
  plots <- raster::intersect(boundary, plots) # remove bits of plot outside boundary
  dets <- spdf[plots, ] # extract only those nests inside the polygons (neat!)

  return(list(plots = plots, dets = dets))
}



#' @title Convert a plot sample of points into one of counts.
#'
#' @description
#' Converts a plot sample with locations of each point within each plot, into a plot
#' sample with only the count within each plot.
#'
#' @param plots A `SpatialPolygonsDataFrame` object containing the plots that were
#'    sampled.
#' @param dets A `SpatialPointsDataFrame` object containing the locations of the
#'    points within the plots.
#'
#' @return A `SpatialPolygonsDataFrame` with counts in each plot contained in
#'   slot `@data$n`.
#'
#' @examples
#' \donttest{
#' # Some features require the raster package
#' if (bru_safe_sp() &&
#'   require("sp") &&
#'   require("raster", quietly = TRUE) &&
#'   require("ggplot2", quietly = TRUE)) {
#'   data(gorillas, package = "inlabru")
#'   plotpts <- plotsample(gorillas$nests, gorillas$boundary,
#'     x.ppn = 0.4, y.ppn = 0.4, nx = 5, ny = 5
#'   )
#'   p1 <- ggplot() +
#'     gg(plotpts$plots) +
#'     gg(plotpts$dets) +
#'     gg(gorillas$boundary)
#'   countdata <- point2count(plotpts$plots, plotpts$dets)
#'   x <- coordinates(countdata)[, 1]
#'   y <- coordinates(countdata)[, 2]
#'   count <- countdata@data$n
#'   p2 <- ggplot() +
#'     gg(gorillas$boundary) +
#'     gg(plotpts$plots) +
#'     geom_text(aes(label = count, x = x, y = y))
#'   multiplot(p1, p2, cols = 2)
#' }
#' }
#'
#' @export
#'
point2count <- function(plots, dets) {
  np <- length(plots)
  x <- y <- plotarea <- count <- numeric(length = np)
  for (i in 1:np) {
    polylist <- list(Polygons(list(plots@polygons[[i]]@Polygons[[1]]), 1))
    spoly <- SpatialPolygons(polylist, proj4string = CRS(as.character(proj4string(plots))))
    count[i] <- dim(dets[spoly, ])[1]

    plotarea[i] <- plots@polygons[[i]]@area
    xs <- unique(plots@polygons[[i]]@Polygons[[1]]@coords[, 1])
    ys <- unique(plots@polygons[[i]]@Polygons[[1]]@coords[, 2])
    x[i] <- min(xs) + abs(diff(range(xs)) / 2)
    y[i] <- min(ys) + abs(diff(range(ys)) / 2)
  }

  # make a data frame of it
  countdf <- data.frame(n = count, area = plotarea, x = x, y = y)
  # make SpatialPointsDataFrame of it
  plotcounts <- SpatialPointsDataFrame(
    coords = data.frame(x = x, y = y), data = data.frame(n = count, area = plotarea),
    proj4string = CRS(as.character(proj4string(plots)))
  )
  return(plotcounts)
}
