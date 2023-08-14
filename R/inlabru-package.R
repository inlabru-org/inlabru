#' inlabru
#'
#' Convenient model fitting using (iterated) INLA.
#'
#' @details
#'
#' `inlabru` facilitates Bayesian spatial modelling using integrated nested
#' Laplace approximations. It is heavily based on R-inla
#' (<https://www.r-inla.org>) but adds additional modelling abilities and simplified
#' syntax for (in particular) spatial models.
#' Tutorials and more information can be found at
#' <https://inlabru-org.github.io/inlabru/> and <http://www.inlabru.org/>.
#' The iterative method used for non-linear predictors is documented in the
#' `method` vignette.
#'
#' The main function for inference using inlabru is [bru()].
#' The general model specification details is documented in [component()] and [like()].
#' Posterior quantities beyond the basic summaries can be calculated with
#' a `predict()` method, documented in [predict.bru()].
#' For point process inference [lgcp()] can be used as a shortcut to `bru(..., like(model="cp", ...))`.
#'
#' The package comes with multiple real world data sets, namely [gorillas],
#' [mexdolphin], [gorillas_sf], [mexdolphin_sf], [seals_sp]. Plotting these data
#'  sets is straight forward using inlabru's extensions
#' to `ggplot2`, e.g. the [gg()] function. For educational purposes some simulated data sets are available
#' as well, e.g. [Poisson1_1D], [Poisson2_1D], [Poisson2_1D] and [toygroups].
#'
#' @aliases inlabru
#' @import sp
#' @import stats
#' @import methods
#' @importFrom Matrix diag
#' @author Fabian E. Bachl \email{bachlfab@@gmail.com}
#'   and Finn Lindgren \email{finn.lindgren@@gmail.com}
"_PACKAGE"

## usethis namespace: start
#' @importFrom lifecycle deprecated
## usethis namespace: end
NULL
