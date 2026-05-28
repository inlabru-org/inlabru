# 1D LGCP bin count simulation and comparison with data

A common procedure of analyzing the distribution of 1D points is to
chose a binning and plot the data's histogram with respect to this
binning. This function compares the counts that the histogram calculates
to simulations from a 1D log Gaussian Cox process conditioned on the
number of data samples. For each bin this results in a median number of
counts as well as a confidence interval. If the LGCP is a plausible
model for the observed points then most of the histogram counts (number
of points within a bin) should be within the confidence intervals. Note
that a proper comparison is a multiple testing problem which the
function does not solve for you.

## Usage

``` r
bincount(
  result,
  predictor,
  observations,
  breaks,
  nint = 20,
  probs = c(0.025, 0.5, 0.975),
  ...
)
```

## Arguments

- result:

  A result object from a
  [`bru()`](https://inlabru-org.github.io/inlabru/reference/bru.md) or
  [`lgcp()`](https://inlabru-org.github.io/inlabru/reference/lgcp.md)
  call

- predictor:

  A formula describing the prediction of a 1D LGCP via
  [`predict()`](https://rdrr.io/r/stats/predict.html).

- observations:

  A vector of observed values

- breaks:

  A vector of bin boundaries

- nint:

  Number of integration intervals per bin. Increase this if the bins are
  wide and the LGCP is not smooth.

- probs:

  numeric vector of probabilities with values in `[0,1]`

- ...:

  arguments passed on to
  [`predict.bru()`](https://inlabru-org.github.io/inlabru/reference/predict.bru.md)

## Value

An `data.frame` with a ggplot attribute `ggp`

## Examples

``` r
if (FALSE) { # \dontrun{
if (require(ggplot2) && require(fmesher) && bru_safe_inla()) {
  # Load a point pattern
  data(Poisson2_1D)

  # Take a look at the point (and frequency) data

  ggplot(pts2) +
    geom_histogram(
      aes(x = x),
      binwidth = 55 / 20,
      boundary = 0,
      fill = NA,
      color = "black"
    ) +
    geom_point(aes(x), y = 0, pch = "|", cex = 4) +
    coord_fixed(ratio = 1)

  # Fit an LGCP model
  x <- seq(0, 55, length.out = 50)
  mesh1D <- fmesher::fm_mesh_1d(x, boundary = "free")
  matern <- INLA::inla.spde2.pcmatern(mesh1D,
    prior.range = c(1, 0.01),
    prior.sigma = c(1, 0.01),
    constr = TRUE
  )
  mdl <- x ~ spde1D(x, model = matern) + Intercept(1)
  fit.spde <- lgcp(mdl, pts2, domain = list(x = mesh1D))

  # Calculate bin statistics
  bc <- bincount(
    result = fit.spde,
    observations = pts2,
    breaks = seq(0, max(pts2), length.out = 12),
    predictor = x ~ exp(spde1D + Intercept)
  )

  # Plot them!
  attributes(bc)$ggp
}
} # }
```
