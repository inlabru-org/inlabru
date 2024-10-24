% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{bru_safe_inla}
\alias{bru_safe_inla}
\title{Load INLA safely for examples and tests}
\usage{
bru_safe_inla(multicore = NULL, quietly = FALSE, minimum_version = "23.1.31")
}
\arguments{
\item{multicore}{logical; if \code{TRUE}, multiple cores are allowed, and the
INLA \code{num.threads} option is not checked or altered.
If \code{FALSE}, forces \code{num.threads="1:1"}. Default: NULL, checks
if running in testthat or non-interactively, in which case sets
\code{multicore=FALSE}, otherwise \code{TRUE}.}

\item{quietly}{logical; if \code{FALSE} and \code{multicore} is \code{FALSE},
prints a message if the \code{num.threads} option
isn't already "1.1" to alert the user to the change.
Default: FALSE.}

\item{minimum_version}{character; the minimum required INLA version.
Default 23.1.31 (should always match the requirement in the package
DESCRIPTION)}
}
\value{
logical; \code{TRUE} if INLA was loaded safely, otherwise FALSE
}
\description{
Loads the INLA package with \code{requireNamespace("INLA", quietly = TRUE)}, and
optionally checks and sets the multicore \code{num.threads} INLA option.
}
\examples{
\dontrun{
if (bru_safe_inla()) {
  # Run inla dependent calculations
}
}

}
\keyword{internal}
