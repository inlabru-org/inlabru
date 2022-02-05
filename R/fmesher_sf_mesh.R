#' @title Coercion methods to and from meshes
#' @rdname fm_as
#' @param ... Arguments passed on to other methods
#' @export
fm_as_sfg <- function(x) {
  UseMethod("fm_as_sfg")
}


#' @rdname fm_as
#' @export
fm_as_inla.mesh <- function(...) {
  UseMethod("fm_as_inla.mesh")
}


#' @rdname fm_as
#' @aliases fm_as_sfg fm_as_sfg.inla.mesh
#'
#' @param x An `inla.mesh` mesh object, or an `sfg` `MULTIPOLYGON` object
#' @returns * `fm_as_sfg`: An `sfg` `MULTIPOLYGON` object
#' @exportS3Method fm_as_sfg inla.mesh
#' @export
fm_as_sfg.inla.mesh <- function(x, ...) {
  stopifnot(inherits(x, "inla.mesh"))
  sf::st_multipolygon(
    lapply(seq_len(nrow(x$graph$tv)),
           function(k) {
             sf::st_polygon(
               list(x$loc[x$graph$tv[k, c(1, 2, 3, 1)], , drop = FALSE])
             )
           }
    ),
    dim = "XYZ")
}


#' @rdname fm_as
#' @aliases fm_as_inla.mesh fm_as_inla.mesh.sfg
#'
#' @returns * `fm_as_inla.mesh`: An `inla.mesh` mesh object
#' @exportS3Method fm_as_inla.mesh sfg
#' @export
fm_as_inla.mesh.sfg <- function(x, ...) {
  stopifnot(inherits(x, "MULTIPOLYGON"))
  tv <- matrix(seq_len(3 * length(x)), length(x), 3, byrow = TRUE)
  loc <- do.call(rbind,
                 lapply(as.list(sf_mesh),
                        function(x) as.matrix(x)[1:3,,drop=FALSE]
                        )
                 )
  mesh <- INLA::inla.mesh.create(loc = loc, tv = tv, ...)
  mesh
}
