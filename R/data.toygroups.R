#' @name toygroups
#' @title Simulated 1D animal group locations and group sizes
#' @docType data
#' @description
#'
#' This data set serves to teach the concept of modelling species that gather in
#' groups and where the grouping behaviour depends on space.
#'
#' @usage data(toygroups)
#'
#' @format The data are a list that contains these elements:
#'  \describe{
#'    \item{`groups`:}{ A `data.frame` of group locations `x` and size `size`}
#'    \item{`df.size`:}{ IGNORE THIS }
#'    \item{`df.intensity`:}{ A `data.frame` with Poisson process
#'      intensity `d.lambda` at locations `x`}
#'    \item{`df.rate`:}{ A `data.frame` the locations `x` and associated `rate`
#'      which parameterized the exponential distribution from which the group
#'      sizes were drawn.}
#'  }
#'
#'
#' @example inst/examples/data.toygroups.R
"toygroups"
