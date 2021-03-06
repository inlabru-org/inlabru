% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/environment.R
\name{summary.bru_options}
\alias{summary.bru_options}
\alias{print.summary_bru_options}
\title{Print inlabru options}
\usage{
\method{summary}{bru_options}(
  object,
  legend = TRUE,
  include_global = TRUE,
  include_default = TRUE,
  ...
)

\method{print}{summary_bru_options}(x, ...)
}
\arguments{
\item{object}{A \link{bru_options} object to be summarised}

\item{legend}{logical; If \code{TRUE}, include explanatory text, Default: \code{TRUE}}

\item{include_global}{logical; If \code{TRUE}, include global override options}

\item{include_default}{logical; If \code{TRUE}, include default options}

\item{...}{Further parameters, currently ignored}

\item{x}{A \code{summary_bru_options} object to be printed}
}
\description{
Print inlabru options
}
\examples{
if (interactive()) {
  options <- bru_options(verbose = TRUE)

  # Don't print options only set in default:
  print(options, include_default = FALSE)

  # Only include options set in the object:
  print(options, include_default = FALSE, include_global = FALSE)
}
}
