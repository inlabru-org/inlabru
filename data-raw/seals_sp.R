#' Seal pup edata import
#'
#' Generate `Spatial` objects from raw seal pup survey data (not inlcuded in
#' [inlabru]). Note that this function
#' will only extract one of the survey transects.
#'
#' @keywords internal
#' @param sealfile Character pointing to the file containing the seal counts and
#'   photo locations
#' @param icefile Character pointing to the `.tif` file containing the ice sheet
#'   covariate
#' @return The [seals_sp] data set
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'

import.seals <- function(
    sealfile = "WestIce2012.csv",
    icefile = "reflectance_0.0025deg_grid_modis_20120328_1310.tif") {
  #' Load seal data

  seals <- utils::read.csv(sealfile)

  #' Turn seal data into spatial object

  sp::coordinates(seals) <- c("longitude", "latitude")
  sp::proj4string(seals) <- sp::CRS("+proj=longlat")

  #' To avoid confusion I remove the x-y coordinates created by Martin

  seals <- seals[, setdiff(names(seals), c("x", "y"))]

  #' Change CRS

  target.p4s <- "+proj=utm +zone=27 +units=km"
  seals <- fm_transform(seals, fm_crs(target.p4s))
  sp::coordnames(seals) <- c("x", "y")

  #' Select a strip

  seals <- seals[(sp::coordinates(seals)[, 2] > 7947) &
                   (sp::coordinates(seals)[, 2] < 7954), ] # strip 9
  # seals = seals[(sp::coordinates(seals)[,2]>7931) &
  #  (sp::coordinates(seals)[,2]<7938.5), ] # strip 12
  # seals = seals[(sp::coordinates(seals)[,2]>7925) &
  #  (sp::coordinates(seals)[,2]<7930), ] # strip 13
  # seals = seals[(sp::coordinates(seals)[,2]>7920) &
  #  (sp::coordinates(seals)[,2]<7924), ] # strip 14
  # seals = seals[(sp::coordinates(seals)[,2]>7910) &
  #  (sp::coordinates(seals)[,2]<7920), ] # strip 15

  #' Add total number of seals

  seals$all <- seals$harps + seals$hooded
  # plot(seals)

  #' Build a mesh. This mesh will be fine at the photo locations but coarse
  #' elsewhere

  bnd <- fm_extensions(sp::coordinates(seals), convex = c(0.5, 0.7))
  mesh <- fm_mesh_2d_inla(
    boundary = bnd,
    max.edge = c(0.2, 3),
    crs = fm_CRS(projargs = sp::CRS(target.p4s))
  )
  # ggplot() + gg(mesh) + gg(seals) + coord_equal()
  # mesh$n

  #' Let's plot the observed counts

  # ggplot() + gg(mesh) + gg(seals, mapping = aes(longitude, latitude, color =
  # log(all/area)), size = 3) + coord_equal()

  #' Ice Covariate

  ice <- sf::read_sf(icefile)
  ice <- fm_transform(ice, fm_crs(target.p4s))
  ii <- fm_is_within(ice, mesh)
  ice <- ice[as.vector(ii), ]

  # ggplot() +  gg(ice, mapping = aes(x, y, color = band1), size = 1) + gg(mesh)
  # + coord_equal() + scale_color_gradientn(colours = topo.colors(100))

  #' Interpolate ice covariate

  stop("Old internal methods needed by import.seals have been removed.")
  # TODO: replace old code in misc/ with modern code
  #  icecv <- covariate(ice, predictor = "band1", mesh = mesh)
  icecv <- NULL
  #  plot(icecv)

  #' Add band1 covariate to seals data frame

  # TODO: replace old code in misc/ with modern code
  # ice.band1 <- evaluator(icecv)
  ice.band1 <- function(x, y) {
    NULL
  }
  seals$ice <- ice.band1(x = sp::coordinates(seals)[, 1],
                         y = sp::coordinates(seals)[, 2])

  #' Plot seal count and ice together
  #
  #   plot.spatial(mesh = mesh,
  #                col = data.frame(mode = icecv$values),
  #                property = "mode", nx = 2000) +
  #     scale_fill_gradientn(colours = topo.colors(100)) +
  #     gg(seals, mapping = aes(x, y, color = log(all/area)), size = 3) +
  #     scale_color_gradientn(colours = heat.colors(100))


  #' Create a data set
  seals <- list(mesh = mesh, points = seals, ice.data = ice, ice.cv = icecv)
}

import.seals_sf <- function() {
  seals <- inlabru::seals_sp

  seals$points <- sf::st_as_sf(seals$points)
  seals$mesh <- fm_as_fm(seals$mesh)
  seals$ice.data <- sf::st_as_sf(seals$ice.data)
  seals$ice.cv$mesh <- fm_as_fm(seals$ice.cv$mesh)

  seals
}


# seals_sf <- import.seals_sf()
# usethis::use_data(seals_sf, overwrite = TRUE, compress = "xz")

# save.seals <- function(...) {
# save.seals("/home/fbachl/devel/r/seals/WestIce2012.csv",
# /home/fbachl/devel/r/seals/reflectance_0.0025deg_grid_modis_20120328_1310.tif
#  seals <- import.seals(...)
#  save("seals",
#       file = paste0(system.file("data", package = "inlabru"),
#       "/seals.RData"))
# }

# seals_sp <- import.seals(...)
#
# seals_sp$mesh <- fm_as_fm(seals_sp$mesh)
# seals_sp$ice.cv$mesh <- fm_as_fm(seals_sp$ice.cv$mesh)
# sp::proj4string(seals_sp$points) <- fm_CRS(seals_sp$mesh)
# sp::proj4string(seals_sp$ice.data) <- fm_CRS(seals_sp$mesh)
#
# use_data(seals_sp, compress = "xz")
