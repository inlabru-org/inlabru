% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inlabru.R
\name{inlabru}
\alias{inlabru}
\title{inlabru}
\description{
Convenient model fitting using (iterated) INLA.
}
\details{
\code{inlabru} facilitates Bayesian spatial modelling using integrated nested
Laplace approximations. It is heavily based on R-inla
(\url{https://www.r-inla.org}) but adds additional modelling abilities and simplified
syntax for (in particular) spatial models.
Tutorials and more information can be found at
\url{https://inlabru-org.github.io/inlabru/} and \url{http://www.inlabru.org/}.
The iterative method used for non-linear predictors is documented in the
\code{method} vignette.

The main function for inference using inlabru is \code{\link[=bru]{bru()}}. For point process inference \code{\link[=lgcp]{lgcp()}} is
a good starting point.
The general model specification details is documented in \code{\link[=component]{component()}} and \code{\link[=like]{like()}}.
Posterior quantities beyond the basic summaries can be calculated with
a \code{predict()} method, documented in \code{\link[=predict.bru]{predict.bru()}}.

The package comes with multiple real world data sets, namely \link{gorillas},
\link{mexdolphin}, \link{seals}. Plotting these data sets is straight forward using inlabru's extensions
to \code{ggplot2}, e.g. the \code{\link[=gg]{gg()}} function. For educational purposes some simulated data sets are available
as well, e.g. \link{Poisson1_1D}, \link{Poisson2_1D}, \link{Poisson2_1D} and \link{toygroups}.
}
\author{
Fabian E. Bachl \email{bachlfab@gmail.com}
and Finn Lindgren \email{finn.lindgren@gmail.com}
}
