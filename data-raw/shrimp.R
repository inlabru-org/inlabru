#' Shrimp data import
#'
#'
#' Load `shrimp` data stored as file `gamba.Rdata` and construct spatial object
#'
#' @aliases import.shrimp
#' @keywords internal
#' @return The [shrimp] data set
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'

import.shrimp <- function() {
  # Load the raw data
  load(file = file.path(system.file("inst/extdata", package = "inlabru"), "gamba.Rdata"))

  # Remove ambiguous coordinates
  gamba <- gamba[, c(1, 2, 3, 6, 7)]

  # Rename columns
  colnames(gamba) <- c("catch", "landing", "depth", "lat", "lon")

  # Turn into spatial object
  coordinates(gamba) <- c("lon", "lat")
  proj4string(gamba) <- CRS("+proj=longlat")

  # Make a mesh
  mesh <- fm_mesh_2d_inla(
    loc.domain = gamba, max.edge = 0.15, offset = 0.15,
    crs = CRS("+proj=longlat")
  )

  # Final shrimp object
  shrimp <- list(hauls = gamba, mesh = mesh)

  return(shrimp)
}

# shrimp <- import.shrimp()
# usethis::use_data(shrimp, overwrite = TRUE)
