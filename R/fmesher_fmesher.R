# Convert loc information to raw matrix coordinates for the mesh
fm_onto_mesh <- function(mesh, loc, crs = NULL) {
  if (!is.matrix(loc) && !fm_crs_is_null(crs)) {
    warning("loc is non-matrix but crs specified; will be ignored")
  }
  if (inherits(loc, c("SpatialPoints", "SpatialPointsDataFrame",
                      "sf", "sfc", "sfg"))) {
    crs <- fm_crs(loc)
  }
  mesh_crs <- fm_crs(mesh)

  loc_needs_normalisation <- FALSE
  if (!fm_crs_is_null(crs) && !fm_crs_is_null(mesh_crs)) {
    if (fm_crs_is_geocent(mesh_crs)) {
      if (!is.matrix(loc)) {
        if (!fm_identical_CRS(crs, mesh_crs)) {
          loc <- fm_transform(loc, crs = mesh_crs, crs0 = crs)
        }
      }
    } else if (!fm_identical_CRS(crs, mesh_crs)) {
      loc <- fm_transform(loc, crs = mesh_crs, crs0 = crs, passthrough = TRUE)
    }
  } else if (identical(mesh$manifold, "S2")) {
    loc_needs_normalisation <- TRUE
  }

  if (inherits(loc, c("SpatialPoints", "SpatialPointsDataFrame"))) {
    loc <- sp::coordinates(loc)
  } else if (inherits(loc, c("sf", "sfc", "sfg"))) {
    loc <- sf::st_coordinates(loc)
    c_names <- colnames(loc)
    c_names <- intersect(c_names, c("X", "Y", "Z"))
    loc <- loc[, c_names, drop = FALSE]
  } else if (!is.matrix(loc)) {
    warning(
      paste0(
        "Unclear if the 'loc' class ('",
        paste0(class(loc), collapse = "', '"),
        "') is of a type we know how to handle."
      ),
      immediate. = TRUE
    )
  }
  if (loc_needs_normalisation) {
    loc <- loc / rowSums(loc^2)^0.5 * mean(rowSums(mesh$loc^2)^0.5)
  }

  loc
}

