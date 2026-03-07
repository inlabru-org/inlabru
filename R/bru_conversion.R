# General conversion to `bru_obs`, `bru_obs_list` objects ####

#' @title Conversion methods for `bru_obs` and `bru_obs_list` objects
#' @description Methods for converting to `bru_obs` and `bru_obs_list`
#'   objects.
#' @param x An object to convert to [bru_obs] or [bru_obs_list]
#' @param \dots Additional arguments passed to sub-methods.
#' @returns An object of class [bru_obs] or [bru_obs_list].
#' @export
#' @seealso [as_bru_comp_list()]
as_bru_obs <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  UseMethod("as_bru_obs")
}

#' @rdname as_bru_obs
#' @export
as_bru_obs_list <- function(x, .tag = NULL) {
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
#' @param .tag character; optional name for the single observation model in the
#'   returned [bru_obs_list] object. Default is `NULL`, which results in
#'   automatic naming based on the `tag` attribute of `x`, if present.
#' @export
as_bru_obs_list.bru_obs <- function(x, .tag = NULL) {
  bru_obs_list(x, .tag = .tag)
}

#' @rdname as_bru_obs
#' @export
as_bru_obs_list.list <- function(x, .tag = NULL) {
  bru_obs_list(x)
}

#' @rdname as_bru_obs
#' @export
as_bru_obs_list.bru_obs_list <- function(x, .tag = NULL) {
  bru_obs_list(x)
}

#' @rdname as_bru_obs
#' @export
as_bru_obs_list.bru <- function(x, .tag = NULL) {
  as_bru_obs_list(x[["bru_info"]])
}

#' @rdname as_bru_obs
#' @export
as_bru_obs_list.bru_info <- function(x, .tag = NULL) {
  # When/if lhoods is moved to the bru_model class, this code can be simplified.
  if (is.null(x[["lhoods"]])) {
    as_bru_obs_list(x[["model"]])
  } else {
    bru_obs_list(x[["lhoods"]])
  }
}

#' @rdname as_bru_obs
#' @export
as_bru_obs_list.bru_model <- function(x, .tag = NULL) {
  # When/if lhoods is moved to the bru_model class, this code can be simplified.
  if (is.null(x[["lhoods"]])) {
    bru_obs_list(list())
  } else {
    bru_obs_list(x[["lhoods"]])
  }
}

#' @title Conversion methods for `bru_comp` and `bru_comp_list` objects
#' @description Methods for converting to `bru_comp` and `bru_comp_list`
#'   objects.
#' @param x An object to convert to [bru_comp] or [bru_comp_list]
#' @param \dots Additional arguments passed on to [bru_comp_list()].
#' @returns An object of class [bru_comp_list].
#' @export
#' @rdname as_bru_comp
#' @name as_bru_comp
#' @seealso [as_bru_obs()], [as_bru_obs_list()]
as_bru_comp <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  UseMethod("as_bru_comp")
}

#' @rdname as_bru_comp
#' @export
as_bru_comp_list <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }
  UseMethod("as_bru_comp_list")
}

#' @rdname as_bru_comp
#' @export
as_bru_comp.bru_comp <- function(x, ...) {
  x
}

#' @rdname as_bru_comp
#' @export
as_bru_comp_list.bru_comp_list <- function(x, ...) {
  bru_comp_list(x, ...)
}

#' @rdname as_bru_comp
#' @export
as_bru_comp_list.bru_comp <- function(x, ...) {
  bru_comp_list(list(x), ...)
}

#' @describeIn as_bru_comp Extract the component list from a [bru()] object.
#' @export
as_bru_comp_list.bru <- function(x, ...) {
  as_bru_comp_list(x[["bru_info"]], ...)
}

#' @describeIn as_bru_comp Extract the component list from a [bru_info()]
#'   object.
#' @export
as_bru_comp_list.bru_info <- function(x, ...) {
  as_bru_comp_list(x[["model"]], ...)
}

#' @describeIn as_bru_comp Extract the component list from a [bru_model()]
#'   object.
#' @export
as_bru_comp_list.bru_model <- function(x, ...) {
  as_bru_comp_list(x[["effects"]], ...)
}

#' @rdname as_bru_comp
#' @export
as_bru_comp_list.list <- function(x, ...) {
  bru_comp_list(x, ...)
}

#' @rdname as_bru_comp
#' @export
as_bru_comp_list.formula <- function(x, ...) {
  bru_comp_list(x, ...)
}
