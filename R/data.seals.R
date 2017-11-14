#' @name seals
#' @title Seal pups
#' @docType data
#' @description This is a single transect of an aereal photo seal pup survey in the Greenland Sea
#' 
#' @usage data(seals)
#' 
#' @format The data contain these objects:
#'  \describe{
#'    \item{\code{points}:}{ A \code{SpatialPointsDataFrame} Center locations of the photos}
#'    \item{\code{mesh}:}{ An \code{inla.mesh} enclosing the plane's transect}
#'    \item{\code{ice.data}:}{ An \code{SpatialPointsDataFrame} with MODIS ice concentration estimates}
#'    \item{\code{ice.cv}:}{ An \code{covdata} object with interpolated ice coverage data}
#'  }
#' @source 
#' Martin Jullum <\email{Martin.Jullum@@nr.no}>
#' 
#' @references 
#' Oigard, T. A. (2013) From pup production to quotas: current status of harp seals in the Greenland Sea.
#' ICES Journal of Marine Science, doi.10.1093/icesjms/fst155.
#' 
#' Oigard, T. A. (2014) Current status of hooded seals in the Greenland Sea. Victims of climate change and predation?,
#' Biological Conservation , 2014, 172, 29 - 36.
#' 
#' @examples
#'  data(seals)
#'  ggplot() + gg(seals$mesh) + gg(seals$points)
#'  
NULL

#' Seal pup edata import
#' 
#' Generate \code{Spatial} objects from raw seal pup survey data (not inlcuded in \link{inlabru}). Note that this function
#' will only extract one of the survey transects.
#'
#' @aliases import.seals
#' @keywords internal
#' @importFrom utils read.csv
#' @param sealfile Character pointing to the file containing the seal counts and photo locations
#' @param icefile Character pointing to the .tif file containing the ice sheet covariate
#' @return The \link{seals} data set
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>
#'

import.seals = function(sealfile = "WestIce2012.csv", icefile = "reflectance_0.0025deg_grid_modis_20120328_1310.tif") {
  
  #' Load seal data  
  
  seals = read.csv(sealfile)
  
  #' Turn seal data into spatial object
  
  coordinates(seals) = c("longitude","latitude")
  proj4string(seals) = CRS("+proj=longlat")
  
  #' To avoid confusion I remove the x-y coordinates created by Martin
  
  seals = seals[,setdiff(names(seals), c("x","y"))]
  
  #' Change CRS
  
  target.p4s = "+proj=utm +zone=27 +units=km"
  seals = spTransform(seals, CRS(target.p4s))
  coordnames(seals) = c("x","y")
  
  #' Select a strip
  
  seals = seals[(coordinates(seals)[,2]>7947) & (coordinates(seals)[,2]<7954), ] # strip 9
  # seals = seals[(coordinates(seals)[,2]>7931) & (coordinates(seals)[,2]<7938.5), ] # strip 12
  # seals = seals[(coordinates(seals)[,2]>7925) & (coordinates(seals)[,2]<7930), ] # strip 13
  #seals = seals[(coordinates(seals)[,2]>7920) & (coordinates(seals)[,2]<7924), ] # strip 14
  #seals = seals[(coordinates(seals)[,2]>7910) & (coordinates(seals)[,2]<7920), ] # strip 15
  
  #' Add total number of seals
  
  seals$all = seals$harps + seals$hooded
  # plot(seals)
  
  #' Build a mesh. This mesh will be fine at the photo locations but coarse elsewhere
  
  bnd = INLA::inla.nonconvex.hull(coordinates(seals), resolution = 170, convex = 0.5)
  bnd2 = INLA::inla.nonconvex.hull(coordinates(seals), resolution = 100, convex = 0.7)
  mesh = INLA::inla.mesh.2d(boundary = list(bnd, bnd2), max.edge = c(0.2,3))
  mesh$crs = INLA::inla.CRS(projargs = CRS(target.p4s))
  # ggplot() + gg(mesh) + gg(seals) + coord_equal()
  # mesh$n
  
  #' Let's plot the observed counts
  
  # ggplot() + gg(mesh) + gg(seals, mapping = aes(longitude, latitude, color = log(all/area)), size = 3) + coord_equal()
  
  #' Ice Covariate
  
  ice = rgdal::readGDAL(icefile)
  ice = spTransform(ice, CRS(target.p4s))
  ii = is.inside(mesh, coordinates(ice))
  ice = ice[as.vector(ii),]
  
  # ggplot() +  gg(ice, mapping = aes(x, y, color = band1), size = 1) + gg(mesh) + coord_equal() + scale_color_gradientn(colours = topo.colors(100)) 
  
  #' Interpolate ice covariate
  
  icecv = covariate(ice, predictor = "band1", mesh = mesh)
  plot(icecv)
  
  #' Add band1 covariate to seals data frame
  
  ice.band1 = evaluator(icecv)
  seals$ice = ice.band1(x=coordinates(seals)[,1], y = coordinates(seals)[,2])
  
  #' Plot seal count and ice together
  # 
  #   plot.spatial(mesh = mesh, col = data.frame(mode = icecv$values), property = "mode", nx = 2000) +
  #     scale_fill_gradientn(colours = topo.colors(100)) +
  #     gg(seals, mapping = aes(x, y, color = log(all/area)), size = 3) +
  #     scale_color_gradientn(colours = heat.colors(100))
  
  
  #' Create a data set
  seals = list(mesh = mesh, points = seals, ice.data = ice, ice.cv = icecv)
}

save.seals = function(...) {
  # save.seals("/home/fbachl/devel/r/seals/WestIce2012.csv", "/home/fbachl/devel/r/seals/reflectance_0.0025deg_grid_modis_20120328_1310.tif")
  seals = import.seals(...)
  save("seals", file=paste0(system.file("data",package="inlabru"),"/seals.RData"))
}