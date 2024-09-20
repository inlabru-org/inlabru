#' Shrimp data import
#'
#' Load `shrimp` data stored as file `gamba.Rdata` and construct spatial object
#'
#' @keywords internal
#' @return The [shrimp] data set
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'
import.shrimp <- function() {
  # Load the raw data
  load(file = file.path(
    system.file("extdata", package = "inlabru"),
    "gamba.Rdata"
  ))

  # Use lat/lon which is actually in a utm system
  gamba <- gamba[, c(
    "Peso.Capturado..Kg.",
    "Peso.Retenido..Kg.",
    "Prof.Media",
    "lat",
    "lon"
  )]

  # Rename columns
  colnames(gamba) <- c("catch", "landing", "depth", "northing", "easting")

  original_crs <- fm_crs("+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs")
  new_crs <- fm_crs("+proj=utm +zone=30 +datum=WGS84 +units=km +no_defs")

  # Turn into spatial object
  gamba1 <- sf::st_as_sf(
    gamba,
    coords = c("easting", "northing"),
    crs = original_crs
  )
  gamba1 <- fm_transform(
    gamba1,
    new_crs
  )

  # Make a mesh
  bnd <- fm_nonconvex_hull_inla(gamba1, 20)
  bnd2 <- fm_nonconvex_hull_inla(gamba1, 50)
  mesh <- fm_mesh_2d_inla(
    boundary = list(bnd, bnd2),
    max.edge = c(5, 20),
    crs = fm_crs(gamba1)
  )

  # Final shrimp object
  shrimp <- list(hauls = gamba1, mesh = mesh)

  return(shrimp)
}

# shrimp <- import.shrimp()
# usethis::use_data(shrimp, compress = "xz", overwrite = TRUE)
