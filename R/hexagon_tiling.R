#' @title Grid partitioning
#'
#' @description
#' `r lifecycle::badge("experimental")` Both the interface and function name
#' may change in future versions.
#'
#' Partitions the region based on the given criteria for calculating residuals
#' in each partition. Parts of this function are taken from concepts in
#' https://rpubs.com/huanfaChen/grid_from_polygon
#'
#' @param samplers A sf object containing region for which
#' partitions to be created
#' @param resolution resolution of the grids that are required
#' @param nrows number of rows of grids that are required
#' @param ncols number of columns of grids that are required
#' @param chess chessboard partitioning
#' @param \dots Unused.
#' @return a partitioned sf object as required or a list of partitioned sf
#'   objects if chess is TRUE
#'
#' @keywords internal
#' @export
#'
#' @author Man Ho Suen
cv_partition <- function(samplers,
                         resolution = NULL,
                         nrows = NULL,
                         ncols = NULL,
                         chess = TRUE,
                         ...) {
  requireNamespace("terra")
  # Create a grid for the given boundary
  if (is.null(resolution)) {
    grid <- terra::rast(
      terra::ext(sf::st_as_sf(
        fmesher::fm_nonconvex_hull(samplers, ...)
      )),
      crs = fmesher::fm_crs(samplers)$input,
      nrows = nrows,
      ncols = ncols
    )
  }

  if (is.null(c(nrows, ncols))) {
    grid <- terra::rast(
      terra::ext(sf::st_as_sf(
        fmesher::fm_nonconvex_hull(samplers, ...)
      )),
      crs = fmesher::fm_crs(samplers)$input,
      resolution = resolution
    )
  }

  gridPolygon <- terra::as.polygons(grid)
  if (isTRUE(chess)) {
    # no idea how it works
    # spatSample(x = gridPolygon, size = 0.5*nrow(gridPolygon),
    #            method="random", strata=NULL, chess="black")
    grid_chess <- terra::init(grid, "chess")
    grid_chess_sf <- sf::st_intersection(
      sf::st_cast(sf::st_as_sf(terra::as.polygons(grid_chess)), "POLYGON"),
      samplers
    )

    # I am not sure which terra version changes how they name for init()
    terra_version <- utils::packageVersion("terra")
    if (terra_version <= "1.8-42") {
      grid_chess_white_sf <- grid_chess_sf[grid_chess_sf$lyr.1 == 1, ]
      grid_chess_black_sf <- grid_chess_sf[grid_chess_sf$lyr.1 == 0, ]
    } else {
      grid_chess_white_sf <- grid_chess_sf[grid_chess_sf$chess == 1, ]
      grid_chess_black_sf <- grid_chess_sf[grid_chess_sf$chess == 0, ]
    }
    return(list(white = grid_chess_white_sf, black = grid_chess_black_sf))
    # ggplot() + gg(data = grid_chess_white_sf, aes(fill = lyr.1)) +
    #   geom_sf(data = nepal_bnd, col = "red", fill = "NA")
    # ggplot() + gg(data = grid_chess_black_sf, aes(fill = lyr.1)) +
    #   geom_sf(data = nepal_bnd, col = "red", fill = "NA")
  } else {
    # Extract the boundary with subpolygons only
    gridPolygon <- sf::st_as_sf(
      terra::intersect(gridPolygon, terra::vect(samplers))
    )
  }
  gridPolygon
}


#' @title Hexagon tiling of region
#'
#' @description
#' `r lifecycle::badge("experimental")` Both the interface and function name
#' may change in future versions.
#'
#' Partitions the region based on the hexagon tiling
#'
#' Input:
#' @param samplers A sf object containing region for which
#' partitions to be created
#' @param cellsize hexagon cellsize, see [sf::st_make_grid()] description
#' @param n_group number of cv folds.
#' @param \dots Passed on to `fm_nonconvex_hull()`, e.g. `resolution`
#' @returns a hexagonal partition covering the `samplers` object, with a `group`
#' integer variable indicating the fold assignment, 1, 2, ..., `n_group`.
#'
#' @keywords internal
#' @export
#'
#' @author Man Ho Suen
#' @examples
#' if (interactive()) {
#'   bnd <- gorillas_sf$boundary
#'   hex_cv <- cv_hex(bnd, cellsize = 0.5, n_group = 3, resolution = 100)
#'   plot(hex_cv)
#'
#'   chess <- cv_partition(bnd, resolution = 0.5, chess = TRUE)
#'   plot(chess$white)
#' }
#'
cv_hex <- function(samplers, cellsize = 0.5, n_group = 3, ...) {
  bnd <- fmesher::fm_nonconvex_hull(samplers, convex = -0.01, ...)
  hex_cell <- sf::st_make_grid(
    bnd,
    cellsize = cellsize,
    square = FALSE
  )
  hex_cell <- sf::st_sf(geometry = hex_cell)

  # get centroids for x positions
  hex_centroids <- sf::st_centroid(hex_cell)
  coords <- sf::st_coordinates(hex_centroids)[, 1]

  # get sorted unique x positions (columns)
  x_unique <- sort(unique(round(coords, 8)))

  # map each hex to a column index
  col_index <- match(round(coords, 8), x_unique)

  # assign groups by cycling 1 → 2 → 3 -> etc -> 1
  col_group <- ((col_index - 1) %% n_group) + 1

  # attach back
  hex_groups <-
    dplyr::mutate(
      hex_cell,
      group = col_group
    )

  sf::st_agr(hex_groups) <- "constant"
  hex_cell_ <- sf::st_intersection(hex_groups, sf::st_geometry(samplers))

  hex_cell_["group"]
}
