#' inlabru
#'
#' Convenient model fitting using (iterated) INLA.
#'
#' @details
#'
#' `inlabru` facilitates Bayesian spatial modeling using integrated nested Laplace approximations. It
#' is heavily based on R-inla (<https://www.r-inla.org>) but adds additional modeling abilities.
#' Tutorials and more information can be found at <http://www.inlabru.org/>.
#'
#' The main function for inference using inlabru is [bru]. For point process inference [lgcp] is
#' a good starting point. The package comes with multiple real world data sets, namely [gorillas],
#' [mexdolphin], [seals]. Plotting these data sets is straight forward using inlabru's extensions
#' to `ggplot2`, e.g. the [gg] function. For educational purposes some simulated data sets are available
#' as well, e.g. [Poisson1_1D], [Poisson2_1D], [Poisson2_1D] and [toygroups].
#'
#' @aliases inlabru
#' @import sp
#' @import stats
#' @import methods
#' @importFrom Matrix diag
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#' @name inlabru
NULL
