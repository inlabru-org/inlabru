# General conversion to `bru_obs` and `bru_obs_list` objects ####

#' @title Conversion methods for `bru_obs` and `bru_obs_list` objects
#' @description Methods for converting to `bru_obs` and `bru_obs_list`
#'   objects.
#' @param x An object to convert to [bru_obs] or [bru_obs_list]
#' @param \dots Additional arguments passed to sub-methods.
#' @returns An object of class [bru_obs] or [bru_obs_list].
#' @keywords internal
#' @export
as_bru_obs <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  UseMethod("as_bru_obs")
}

#' @rdname as_bru_obs
#' @export
as_bru_obs_list <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  UseMethod("as_bru_obs_list")
}

#' @rdname as_bru_obs
#' @export
as_bru_obs.bru_obs <- function(x, ...) {
  x
}

#' @rdname as_bru_obs
#' @export
as_bru_obs_list.bru_obs <- function(x, ...) {
  bru_obs_list(list(x), ...)
}

#' @rdname as_bru_obs
#' @export
as_bru_obs_list.list <- function(x, ...) {
  bru_obs_list(x, ...)
}

#' @rdname as_bru_obs
#' @export
as_bru_obs_list.bru_obs_list <- function(x, ...) {
  bru_obs_list(x, ...)
}
