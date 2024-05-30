#' @name mexdolphin
#' @title Pan-tropical spotted dolphins in the Gulf of Mexico
#' @docType data
#' @description This a version of the `mexdolphins` dataset from the package `dsm`, reformatted
#' as point process data for use with `inlabru`. The data are from a combination of several NOAA
#' shipboard surveys conducted on pan-tropical spotted dolphins in the Gulf of Mexico. 47 observations
#' of groups of dolphins were detected. The group size was recorded, as well as the Beaufort sea state at
#' the time of the observation. Transect width is 16 km, i.e. maximal detection
#' distance 8 km (transect half-width 8 km).
#'
#' @format A list of objects:
#'  \describe{
#'    \item{`points`:}{ A `SpatialPointsDataFrame` object containing the locations of
#'    detected dolphin groups, with their size as an attribute.}
#'    \item{`samplers`:}{ A `SpatialLinesDataFrame` object containing the transect lines
#'    that were surveyed.}
#'    \item{`mesh`:}{ An `inla.mesh` object containing a Delaunay triangulation
#'    mesh (a type of discretization of continuous space) covering the survey region.}
#'    \item{`ppoly`:}{ An `SpatialPolygonsDataFrame` object defining the boundary of the
#'    survey region.}
#'    \item{`simulated`:}{ A `SpatialPointsDataFrame` object containing the locations of a
#'    *simulated* population of dolphin groups. The population was simulated from a `inlabru`
#'    model fitted to the actual survey data. Note that the simulated data do not have any associated
#'    size information.}
#'  }
#' @source
#' Library `dsm`.
#'
#' @references
#' Halpin, P.N., A.J. Read, E. Fujioka, B.D. Best, B. Donnelly, L.J. Hazen, C. Kot, K. Urian,
#' E. LaBrecque, A. Dimatteo, J. Cleary, C. Good, L.B. Crowder, and K.D. Hyrenbach. 2009.
#' OBIS-SEAMAP: The world data center for marine mammal, sea bird, and sea turtle distributions.
#' Oceanography 22(2):104-115
#'
#' NOAA Southeast Fisheries Science Center. 1996. Report of a Cetacean Survey of Oceanic and
#' Selected Continental Shelf Waters of the Northern Gulf of Mexico aboard NOAA Ship Oregon II
#' (Cruise 220)
#'
#' @examples
#' \donttest{
#' if (require("ggplot2", quietly = TRUE)) {
#'   data(mexdolphin, package = "inlabru")
#'   ggplot() +
#'     gg(mexdolphin$mesh) +
#'     gg(mexdolphin$ppoly, color = "blue") +
#'     gg(mexdolphin$samplers) +
#'     gg(mexdolphin$points, aes(size = size), color = "red") +
#'     coord_equal()
#'
#'   ggplot() +
#'     gg(mexdolphin$mesh, col = mexdolphin$lambda, mask = mexdolphin$ppoly) +
#'     coord_equal()
#' }
#' }
#' \dontrun{
#' if (requireNamespace("ggmap", quietly = TRUE) &&
#'   require("ggplot2", quietly = TRUE)) {
#'   gmap(mexdolphin$depth) +
#'     gm(mexdolphin$ppoly, color = "blue") +
#'     gm(mexdolphin$samplers) +
#'     gm(mexdolphin$points, aes(size = size), color = "red")
#'
#'   gmap(mexdolphin$depth) +
#'     gm(mexdolphin$depth, aes(col = depth)) +
#'     gm(mexdolphin$ppoly)
#' }
#' }
"mexdolphin"

#' @name mexdolphin_sf
#' @title Pan-tropical spotted dolphins in the Gulf of Mexico
#' @docType data
#' @description This a version of the `mexdolphins` dataset from the package `dsm`, reformatted
#' as point process data for use with `inlabru`, with the parts stored in `sf` format.
#' The data are from a combination of several NOAA
#' shipboard surveys conducted on pan-tropical spotted dolphins in the Gulf of Mexico. 47 observations
#' of groups of dolphins were detected. The group size was recorded, as well as the Beaufort sea state at
#' the time of the observation. Transect width is 16 km, i.e. maximal detection
#' distance 8 km (transect half-width 8 km).
#'
#' @format A list of objects:
#'  \describe{
#'    \item{`points`:}{ An `sf` object containing the locations of
#'    detected dolphin groups, with their size as an attribute.}
#'    \item{`samplers`:}{ An `sf` object containing the transect lines
#'    that were surveyed.}
#'    \item{`mesh`:}{ An `fm_mesh_2d` object containing a Delaunay triangulation
#'    mesh (a type of discretization of continuous space) covering the survey region.}
#'    \item{`ppoly`:}{ An `sf` object defining the boundary of the
#'    survey region.}
#'    \item{`simulated`:}{ A `sf` object containing the locations of a
#'    *simulated* population of dolphin groups. The population was simulated from a `inlabru`
#'    model fitted to the actual survey data. Note that the simulated data do not have any associated
#'    size information.}
#'  }
#' @source
#' Library `dsm`.
#'
#' @references
#' Halpin, P.N., A.J. Read, E. Fujioka, B.D. Best, B. Donnelly, L.J. Hazen, C. Kot, K. Urian,
#' E. LaBrecque, A. Dimatteo, J. Cleary, C. Good, L.B. Crowder, and K.D. Hyrenbach. 2009.
#' OBIS-SEAMAP: The world data center for marine mammal, sea bird, and sea turtle distributions.
#' Oceanography 22(2):104-115
#'
#' NOAA Southeast Fisheries Science Center. 1996. Report of a Cetacean Survey of Oceanic and
#' Selected Continental Shelf Waters of the Northern Gulf of Mexico aboard NOAA Ship Oregon II
#' (Cruise 220)
#'
#' @examples
#' \donttest{
#' if (require("ggplot2", quietly = TRUE)) {
#'   data(mexdolphin_sf, package = "inlabru")
#'   ggplot() +
#'     gg(mexdolphin_sf$mesh) +
#'     gg(mexdolphin_sf$ppoly, color = "blue", alpha = 0, linewidth = 1) +
#'     gg(mexdolphin_sf$samplers) +
#'     gg(mexdolphin_sf$points, aes(size = size), color = "red") +
#'     scale_size_area()
#'
#'   ggplot() +
#'     gg(mexdolphin_sf$mesh, color = mexdolphin_sf$lambda, mask = mexdolphin_sf$ppoly)
#' }
#' }
"mexdolphin_sf"
