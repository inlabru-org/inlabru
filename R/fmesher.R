# Re-export fmesher functions for temporary backwards compatibility.
# These functions are not intended for public use and will be removed in a
# future release.
# For 2.14.1, fm_int and fm_pixels are kept, as they are used
# by intSDM without explicitly importing from fmesher.

# @export
# fmesher::fm_cprod
# @export
# fmesher::fm_crs
#' @export
fmesher::fm_int
#' @export
fmesher::fm_pixels
