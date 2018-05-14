#' @name gorillas
#' @title Gorilla Nesting Sites.
#' @docType data
#' @description This the \code{gorillas} dataset from the package \code{spatstat}, reformatted
#' as point process data for use with \code{inlabru}. 
#' 
#' @usage data(gorillas)
#' 
#' @format The data are a list that contains these elements:
#'  \describe{
#'    \item{\code{nests}:}{ A \code{SpatialPointsDataFrame} object containing the locations of 
#'    the gorilla nests.}
#'    \item{\code{boundary}:}{ An \code{SpatialPolygonsDataFrame} object defining the boundary
#'    of the region that was searched for the nests.}
#'    \item{\code{mesh}:}{ An \code{inla.mesh} object containing a mesh that can be used
#'    with function \code{lgcp} to fit a LGCP to the nest data.}
#'    \item{\code{gcov}:}{ A list of SpatialGridDataFrame objects, one for each of these spatial covariates:
#'     \describe{
#'       \item{\code{aspect}}{ Compass direction of the terrain slope. Categorical, with levels 
#'       N, NE, E, SE, S, SW, W and NW, which are coded as integers 1 to 8.}
#'       \item{\code{elevation}}{ Digital elevation of terrain, in metres.}
#'       \item{\code{heat}}{ Heat Load Index at each point on the surface (Beer's aspect), 
#'       discretised. Categorical with values Warmest (Beer's aspect between 0 and 0.999), 
#'       Moderate (Beer's aspect between 1 and 1.999), Coolest (Beer's aspect equals 2). These are
#'       coded as integers 1, 2 and 3, in that order.}
#'       \item{\code{slopangle}}{ Terrain slope, in degrees.}
#'       \item{\code{slopetype}}{ Type of slope. Categorical, with values Valley, Toe (toe slope), 
#'       Flat, Midslope, Upper and Ridge. These are coded as integers 1 to 6.}
#'       \item{\code{vegetation}}{ Vegetation type: a categorical variable with 6 levels coded as
#'       integers 1 to 6 (in order of increasing expected habitat suitability)}
#'       \item{\code{waterdist}}{ Euclidean distance from nearest water body, in metres.}
#'     }
#'    }
#'    \item{\code{plotsample}}{Plot sample of gorilla nests, sampling 9x9 over the region, with 60\% coverage. Components:
#'    \describe{
#'      \item{\code{counts}}{ A SpatialPointsDataFrame frame with elements \code{x}, \code{y}, \code{count}, 
#'      \code{exposure}, being the x- and y-coordinates of the centre of each plot, the count in each 
#'      plot and the area of each plot.}
#'      \item{\code{plots}}{ A \code{SpatialPolygonsDataFrame} defining the individual plot boundaries.}
#'      \item{\code{nests}}{ A \code{SpatialPointsDataFrame} giving the locations of each detected nest.}
#'    }
#'    }
#'  }
#' @source 
#' Library \code{spatstat}.
#' 
#' 
#' @references 
#' Funwi-Gabga, N. (2008) A pastoralist survey and fire impact assessment in the Kagwene Gorilla 
#' Sanctuary, Cameroon. M.Sc. thesis, Geology and Environmental Science, University of Buea, 
#' Cameroon.
#' 
#' Funwi-Gabga, N. and Mateu, J. (2012) Understanding the nesting spatial behaviour of gorillas 
#' in the Kagwene Sanctuary, Cameroon. Stochastic Environmental Research and Risk Assessment 
#' 26 (6), 793-811.
#' 
#' @examples
#' data(gorillas, package = "inlabru") # get the data
#' # extract all the objects, for convenience:
#'
#' # plot all the nests, mesh and boundary
#' ggplot() + gg(gorillas$mesh) + gg(gorillas$boundary) + gg(gorillas$nests)
#'
#' # Plot the elevation covariate
#' plot(gorillas$gcov$elevation)
#'
#' # Plot the plot sample
#' ggplot() + gg(gorillas$plotsample$plots) + gg(gorillas$plotsample$nests)
NULL

#' Gorilla data import
#'
#' @aliases import.gorillas
#' @keywords internal
#' @importFrom utils data
#' @author Fabian E. Bachl <\email{bachlfab@@gmail.com}>, David L. Borchers <\email{dlb@@st-andrews.ac.uk}> 
#' @return gorilla data

import.gorillas = function() {
  
  # Explicitly load spatstat
  # library(spatstat)
  
  # Load Gorilla data from spatstat
  gorillas = NULL
  gorillas.extra = NULL
  data(gorillas, package="spatstat", envir = environment())
  
  # Create SpatialPoints representing nest locations
  nests = as.data.frame(gorillas)
  coordinates(nests) = c("x","y")
  proj4string(nests) = CRS("+proj=utm +zone=32N +datum=WGS84") # from the Gorillas help file
  
  #' Turn the observation window into spatial polygon
  boundary = spoly(gorillas$window$bdry, crs = CRS("+proj=utm +zone=32N +datum=WGS84"))
  
  #' Build mesh
  bnd = INLA::inla.mesh.segment(loc = data.frame(gorillas$window$bdry[[1]]))
  mesh = INLA::inla.mesh.2d(interior = bnd, max.edge = 222) # ! With higher max.edge we run into various INLA errors/warnings
  mesh$crs = INLA::inla.CRS(proj4string(nests))
  
  #' Turn covariates int SpatialGridDataFrame
  gcov = list()
  for ( nm in names(gorillas.extra) ) { 
    gcov[[nm]] = as(gorillas.extra[[nm]], "SpatialGridDataFrame") 
    proj4string(gcov[[nm]]) = proj4string(nests)
    coordnames(gcov[[nm]]) = c("x","y")
    names(gcov[[nm]]) = nm
  }
  
  #' Hack: change CRS units of the covariates to km
  for ( nm in names(gcov) ) { 
    ga = attributes(gcov[[nm]])$grid
    attributes(ga)$cellcentre.offset = attributes(ga)$cellcentre.offset/1000
    attributes(ga)$cellsize = attributes(ga)$cellsize/1000
    attributes(gcov[[nm]])$grid = ga
    attributes(gcov[[nm]])$proj4string = CRS("+proj=utm +zone=32N +datum=WGS84 +units=km")
  }
  
  
  ### Make final gorilla data set
  gorillas = list(nests = nests,
                  mesh = mesh,
                  boundary = boundary,
                  gcov = gcov)
  
  gorillas = stransform(gorillas, CRS("+proj=utm +zone=32N +datum=WGS84 +units=km"))
  
  # Create a plot sampling data set
  set.seed(121)
  plotpts = plotsample(gorillas$nests, gorillas$boundary, x.ppn=0.6, y.ppn=0.6, nx=5.4, ny=5.4)
  counts = point2count(plotpts$plots,plotpts$dets)
  x = coordinates(counts)[,1]
  y = coordinates(counts)[,2]
  count = counts@data$n
  
  # Make gam data frame
  gnestcount_9x9_60pc = data.frame(x = x, y = y, count = count, exposure = counts$area)
  gnestplots_9x9_60pc = plotpts$plots
  gnestdets_9x9_60pc = plotpts$dets
  sample_9x9_60pc = list( counts = gnestcount_9x9_60pc, 
                          plots = gnestplots_9x9_60pc, 
                          nests = gnestdets_9x9_60pc)
  
  # Attach plotsample to gorilla data
  gorillas$plotsample = sample_9x9_60pc
  
  # Make the count data frame spatial
  coordinates(gorillas$plotsample$counts) = c("x","y")
  proj4string(gorillas$plotsample$counts) = CRS(proj4string(gorillas$nests))
  
  
  # Extrapolate covariate
  pxl = pixels(gorillas$mesh, mask = FALSE, nx = 220, ny = 180)
  for (k in 1:length(gorillas$gcov)) {
    gorillas$gcov[[k]] = sfill(gorillas$gcov[[k]], pxl)
  }
  
  return(gorillas)
}


save.gorillas = function() {
  gorillas = import.gorillas()
  save("gorillas", file=paste0(system.file("data",package="inlabru"),"/gorillas.RData"))
}

