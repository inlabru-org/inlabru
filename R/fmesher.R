# Re-export fmesher functions for temporary backwards compatibility.
# These functions are not intended for public use and will be removed in a
# future release.
# For 2.14.1, fm_int and fm_pixels are kept, as they are used
# by intSDM without explicitly importing from fmesher.
# intSDM 2.1.2 still assumes fm_int access, but this is fixed on github so
# the next release should be ok.

# @export
# fmesher::fm_cprod
# @export
# fmesher::fm_crs
#' @export
fmesher::fm_int
# @export
# fmesher::fm_pixels
