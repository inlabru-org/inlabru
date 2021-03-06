% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{evaluate_predictor}
\alias{evaluate_predictor}
\title{Evaluate component effects or expressions}
\usage{
evaluate_predictor(model, state, data, effects, predictor, format = "auto")
}
\arguments{
\item{state}{A list where each element is a list of named latent state
information, as produced by \code{\link[=evaluate_state]{evaluate_state()}}}

\item{data}{A \code{list}, \code{data.frame}, or \code{Spatial*DataFrame}, with coordinates
and covariates needed to evaluate the model.}

\item{effects}{A list where each element is list of named evaluated effects,
as computed by \code{\link[=evaluate_effect_multi.component_list]{evaluate_effect_multi.component_list()}}}

\item{predictor}{Either a formula or expression}

\item{format}{character; determines the storage format of the output.
Available options:
\itemize{
\item \code{"auto"} If the first evaluated result is a vector or single-column matrix,
the "matrix" format is used, otherwise "list".
\item \code{"matrix"} A matrix where each column contains the evaluated predictor
expression for a state.
\item \code{"list"} A list where each column contains the evaluated predictor
expression for a state.
}

Default: "auto"}
}
\value{
A list or matrix is returned, as specified by \code{format}
}
\description{
Evaluate component effects or expressions, based on a bru model and one or
several states of the latent variables and hyperparameters.
}
\keyword{internal}
